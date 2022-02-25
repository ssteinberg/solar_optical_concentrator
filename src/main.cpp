
/*
 *  Copyright Shlomi Steinberg
 */

#include <filesystem>
#include <iostream>
#include <thread>
#include <cassert>
#include <getopt.h>
#include <cstring>
#include <cfenv>
#include <chrono>
#include <future>
#include <functional>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/rotate_vector.hpp>

#include "mirrors.hpp"
#include "raytracer.hpp"

#include "thread_pool.hpp"
#include "random_generator.hpp"

#include "humanreadable_chrono_duration.hpp"
#include "print_progress.hpp"
#include "range.hpp"


enum class target_type {
    unknown,
    cyl,
    flat
};

inline auto create_target(target_type type, point_t offset, f_t size, bool flip) {      // Flip is for PM, where reflected rays arrive from -X and not +X
    return type==target_type::cyl ? 
        std::unique_ptr<primitive_t>(std::make_unique<ball_t>(offset, size)) :
        std::unique_ptr<primitive_t>(std::make_unique<line_t>(offset+point_t{0,-size,0},
                                                              offset+point_t{0, size,0},
                                                              glm::vec3{flip?-1:1,0,0},glm::vec3{flip?-1:1,0,0},
                                                              0,true));
}
inline void write_target_csv(const std::filesystem::path& out, target_type type, point_t offset, f_t alpha, f_t size) {
    std::ofstream f(out);
    f << "type,x,y,z,alpha,size" << std::endl;
    f << (type==target_type::cyl?"cyl":"flat") << "," << serialize<10>(offset) << "," << std::setprecision(10) << std::to_string(alpha) << "," << std::to_string(size) << std::endl;
}

struct optimize_params_t {
    point_t offset{};
    f_t alpha{ 1 };
    rt::statistics_t stats{};
};

struct performance_t {
    f_t epsilon;
    f_t beta_max;
    f_t P;
    f_t Ai, At, Atp;
    point_t offset;
    f_t target_size;
    f_t hitratio;
};

inline auto write_performance_header(const std::filesystem::path &outpath) {
    auto f = std::ofstream(outpath/"performance.csv");
    f << "epsilon,beta_max,P,Ai,At,Atp,target_size,offsetx,offsety,offsetz" << std::endl;
    return f;
}
inline void write_performance(std::ostream& out, const performance_t &p) {
    out << std::setprecision(10) << std::to_string(p.epsilon) << "," 
        << std::setprecision(10) << std::to_string(p.beta_max) << "," 
        << std::setprecision(10) << std::to_string(p.P) << "," 
        << std::setprecision(10) << std::to_string(p.Ai) << "," 
        << std::setprecision(10) << std::to_string(p.At) << "," 
        << std::setprecision(10) << std::to_string(p.Atp) << "," 
        << std::setprecision(10) << std::to_string(p.target_size) << ","
        << serialize<10>(p.offset) << ","
        << std::setprecision(10) << std::to_string(p.hitratio)
        << std::endl;
}
inline void write_performance(const std::filesystem::path &outpath, const performance_t &p) {
    auto o = write_performance_header(outpath);
    write_performance(o, p);
}


performance_t work(
        const std::filesystem::path outpath,
        const bool write_out,
        const f_t f1, 
        const f_t f2,
        const f_t L,
        const f_t dY,
        const f_t beta_max,
        const f_t epsilon,
        const f_t beam_in,
        const point_t target_offset,
        std::optional<f_t> desired_target_size,
        std::size_t p,
        std::size_t total_rays,
        std::size_t mirror_setup,
        bool intersect_back_faces,
        bool only_final_segment,
        bool keep_tracing_after_target,
        bool measure_performance,
        bool optimize,
        bool disable_offset_optimization,
        bool disable_size_optimization,
        const target_type type,
        const target_type src_type
    ) {
    const auto start = std::chrono::high_resolution_clock::now();
    thread_pool tp(p);


    // Numeric integrator -- compute mirrors
    std::cout << "Mirrors..." << std::endl;
    print_progress(.0);

    const auto& I = [=](const f_t a)    { return mirror_setup!=3 || src_type==target_type::cyl ? 1 : glm::cos(a); };
    const auto& S = [=](const f_t beta) { return type==target_type::cyl ? 1 : glm::cos(beta); };
    ::mirrors::mirrors_t mirrors;
    mirrors = mirror_setup == 0 ? ::mirrors::compute_mirror_parabolic(f2,dY,beta_max) :
              mirror_setup != 3 ? ::mirrors::compute_mirrors_design1(f2,L,beam_in,mirror_setup==1,I,S,dY,beta_max) :
                                  ::mirrors::compute_mirrors_design3(f1,f2,L,I,S,dY,beta_max);
    
	std::printf("%c[2K\r", 27);
	std::fflush(stdout);

    
    f_t src_size;
    if (mirror_setup==3) {
        src_size = .5 * f1 * glm::sin(epsilon);
        std::cout << "Using source with radius=" << std::to_string(src_size) << std::endl;
        
        if (!desired_target_size)
            desired_target_size = src_size;
    }

    if (!desired_target_size) {
        desired_target_size = epsilon>0 ? 
                        .9995 * mirrors.D * glm::sin(epsilon/2) / (type==target_type::cyl ? glm::pi<f_t>()*2 : f_t(2)) :
                        1e-3;
    }
    assert(!!desired_target_size && desired_target_size.value()>0);
    std::cout << "Using target with radius=" << std::to_string(desired_target_size.value()) << " offset=<" << serialize<8>(target_offset) << ">" << std::endl;
    
    
    // Populate octree
    std::cout << "Building octree with " << std::to_string(mirrors.m1.size()+mirrors.m2.size()) << " objects..." << std::endl;
    print_progress(.0);
    
    rt::octree_t oct{ mirrors.bb + aabb<f_t>{ point_t{-1},point_t{1} } };
    {
        static constexpr std::size_t chunk = 512;
        const auto sz = mirrors.m1.size()+mirrors.m2.size();
        
        std::vector<std::future<void>> futures;
        futures.reserve(sz/chunk+1);

        std::atomic<std::size_t> done{ 0 };
        for (std::size_t it=0; it<sz; it+=chunk) {
            futures.emplace_back(tp.enqueue([&mirrors,&oct,sz,it,&done]() {
                for (std::size_t i=it;i<it+chunk && i<sz;++i)
                    oct.insert(i<mirrors.m1.size() ? &mirrors.m1[i] : &mirrors.m2[i-mirrors.m1.size()]);
                print_progress(done++*chunk/(float)sz);
            }));
        } 
        for (auto &f : futures) f.get();
    }

	std::printf("%c[2K\r", 27);
	std::fflush(stdout);


    // Trace and optimize
    optimize ?
        std::cout << "Optimizing..." << std::endl :
        std::cout << "Tracing a total of " << std::to_string(2*total_rays) << " paths..." << std::endl;
    const bool printprog = !optimize;

    if (printprog) 
        print_progress(.0);

    const auto target_size = desired_target_size.value();

    
    static constexpr auto opt_levels = 14, opt_size_levels = 18;
    const auto optimize_offset_step_inital = target_size/f_t(64);
    const auto optimize_size_step_max = 1/f_t(2);

    optimize_params_t opt_params_cur{}, opt_params_prev{};
    opt_params_cur.offset = target_offset;
    int optimize_offset_dir   = 0;
    int optimize_offset_level = 0;
    int optimize_size_level   = opt_size_levels;
    auto binsrch_l = f_t(0);
    auto binsrch_h = 2*target_size;
    auto binsrch_thres = target_size*1e-5;

    performance_t perf;

    const bool print_opt_results = optimize;
    bool done = false;
    while (true) {
        // Create target
        const auto cur_target_size = target_size * opt_params_cur.alpha;
        const auto target = create_target(type,opt_params_cur.offset,cur_target_size,mirror_setup==0);
        oct.insert(target.get());
        
        const auto ray_generator = [&](auto rotf, const std::size_t ri) -> std::optional<rt::gen_ray_t> {
            if (ri>=total_rays)
                return std::nullopt;
            
            const auto x = ri/(f_t)total_rays;
            const auto a = rotf();

            if (mirrors.src_type==mirrors::source_type::directional) {
                const auto beam_in = mirrors.src;
                const auto tangent = glm::cross(beam_in, point_t{ 0,0,1 });
                const auto R     = glm::mix(-mirrors.D/2.,mirrors.D/2.,x)*f_t(1.1);
                const auto dst   = (mirrors.centre + f_t(.01456789) * beam_in) + R*tangent;
                const auto ray   = glm::rotateZ(beam_in,a);

                return rt::gen_ray_t{ ray_t<f_t>{ dst - rt::far*ray, ray }, 
                                      glm::abs(R) };
            }
            else {
                const auto alpha = f_t(1.05) * glm::mix(-mirrors.point_src_max_alpha,mirrors.point_src_max_alpha,x);
                const auto ray   = ray_t<f_t>{ mirrors.src, point_t{ -glm::cos(alpha), glm::sin(alpha), 0 } };
                std::optional<f_t> dist;
                oct.trace(ray,[&](const auto *obj) {
                    if (obj->mirror==1) {
                        const auto intr = obj->type==primitive_type::line ?
                            intersect_ray_line(ray, *(line_t*)obj) :
                            intersect_ray_tri(ray, *(tri_t*)obj);
                        if (intr)
                            dist = intr.value().d;
                        return !intr; 
                    }
                    return true;
                });
                assert(!!dist);
                if (!dist) 
                    return {};
                
                const auto mirror_point = ray.p + dist.value() * ray.d;
                point_t src_point = src_type==target_type::cyl ?
                    ray.p + glm::rotateZ(ray.d, 2*a/epsilon * glm::pi<f_t>()/2) * src_size :
                    ray.p + point_t{0,1,0} * f_t(2*a/epsilon * src_size);
                return rt::gen_ray_t{ ray_t<f_t>{ src_point, glm::normalize(mirror_point-src_point) }, 
                                      0 };
            }
        };
        
        const bool montecarlo = measure_performance;
        const auto mc_rotf_gen = [=]() { 
            const auto rand = std::uniform_real_distribution<f_t>{ 0,1 }(random_generator{}());
            return glm::mix(-epsilon/2,+epsilon/2,rand);
        };

        std::unique_ptr<rt::intensity_distribution_t> target_intensity_distribution;
        if (measure_performance)
            target_intensity_distribution = rt::create_intensity_distribution_bins(cur_target_size);
        const auto &raysa = montecarlo ?
            rt::launch_ray_packets(tp, std::bind(ray_generator, mc_rotf_gen, std::placeholders::_1),
                                   printprog ? [](double p) { print_progress(p/2); } : [](double) {},
                                   oct, target_intensity_distribution.get(), intersect_back_faces, only_final_segment, keep_tracing_after_target) :
            rt::launch_ray_packets(tp, std::bind(ray_generator, [=](){ return -epsilon/2; }, std::placeholders::_1),
                                   printprog ? [](double p) { print_progress(p/2); } : [](double) {},
                                   oct, target_intensity_distribution.get(), intersect_back_faces, only_final_segment, keep_tracing_after_target);
        const auto &raysb = montecarlo ?
            rt::launch_ray_packets(tp, std::bind(ray_generator, mc_rotf_gen, std::placeholders::_1),
                                   printprog ? [](double p) { print_progress(p/2+.5); } : [](double) {},
                                   oct, target_intensity_distribution.get(), intersect_back_faces, only_final_segment, keep_tracing_after_target) :
            rt::launch_ray_packets(tp, std::bind(ray_generator, [=](){ return +epsilon/2; }, std::placeholders::_1),
                                   printprog ? [](double p) { print_progress(p/2+.5); } : [](double) {},
                                   oct, target_intensity_distribution.get(), intersect_back_faces, only_final_segment, keep_tracing_after_target);
        opt_params_cur.stats = rt::statistics(std::array<const rt::traced_lines_t*,2>{ &raysa,&raysb });

        if (printprog) {
            std::printf("%c[2K\r", 27);
            std::fflush(stdout);
        }
        

        oct.erase(target.get());
        
        // Optimize
        const auto marginal_improved = opt_params_prev.stats.total_rays == 0 ? true : opt_params_prev.stats.effective_hit<opt_params_cur.stats.effective_hit;
        const auto improved = opt_params_prev.stats.total_rays == 0 ? true : 
            opt_params_prev.stats.effective_hit<1 && 
            (opt_params_cur.stats.effective_hit-opt_params_prev.stats.effective_hit) / glm::max<f_t>(.1,1-opt_params_prev.stats.effective_hit) >= 0.002;
        const auto fullhit = opt_params_cur.stats.total_rays == 0 ? false : opt_params_cur.stats.effective_hit==1;
        auto opt_offset = optimize && !disable_offset_optimization && optimize_offset_level < opt_levels;
        auto opt_size   = optimize && !disable_size_optimization && optimize_size_level>=0;
        
        // Check last iteration
        if (opt_offset) {
            if (!improved) {
                // No improvement? Revert params. 
                opt_params_cur = opt_params_prev;
                if (optimize_offset_dir == 4) {
                    std::cout << 'X';
                    optimize_offset_dir = 0;
                    ++optimize_offset_level;
                    opt_offset = optimize_offset_level < opt_levels;
                    // Continue to size optimization?
                    opt_size   = !fullhit && !disable_size_optimization;
                }
            }
            else if (opt_params_prev.stats.total_rays>0) {
                std::cout << (optimize_offset_dir==2?'v':
                              optimize_offset_dir==3?'^':
                              optimize_offset_dir==1?'<':'>');
                // New params, restart offset search.
                optimize_offset_dir = 0;
            }
        }
        else if (opt_size && !disable_offset_optimization) {
            // If target size is satisfactory, revert and increase level
            // Otherwise, continue as is.
            if (!marginal_improved && !fullhit) {
                std::cout << '=';
                optimize_size_level=glm::max(0,optimize_size_level-1);
            }
            else if (!improved && !fullhit) {
                std::cout << '+';
            }
            else {
                // Go back to offset optimization first
                optimize_offset_level = 0;
                optimize_offset_dir = 0;
                optimize_size_level=opt_size_levels;
                opt_offset = !fullhit && optimize_offset_level < opt_levels;
                opt_size = !fullhit && !disable_size_optimization;
                std::cout << (opt_size ? '*' : '#');
            }
        }
        std::fflush(stdout);

        opt_params_prev = opt_params_cur;
        if (optimize && disable_offset_optimization) {
            // Simple binary search
            const auto mid = binsrch_l+(binsrch_h-binsrch_l)*.5;
            opt_size = !fullhit || (binsrch_h-binsrch_l)>binsrch_thres;
            if (fullhit) { 
                binsrch_h = mid;
                std::cout << '-';
            }
            else {
                binsrch_l = mid;
                std::cout << '+';
            }
            opt_params_cur.alpha = (binsrch_l+(binsrch_h-binsrch_l)*.5)/target_size;
            
            if ((binsrch_h-binsrch_l)<binsrch_thres*1e-1 && binsrch_h==2*target_size && !fullhit) {
                std::cerr << "Initial target too small" << std::endl;
                exit(1);
            }
        }
        else if (opt_offset) {
            // Optimize offset
            ++optimize_offset_dir;
            const point_t dir = 
                optimize_offset_dir == 2 ? point_t{  0,-1, 0 } : 
                optimize_offset_dir == 3 ? point_t{  0, 1, 0 } : 
                optimize_offset_dir == 1 ? point_t{ -1, 0, 0 } : 
                                           point_t{  1, 0, 0 };
            opt_params_cur.offset += dir * f_t(optimize_offset_step_inital / (1<<optimize_offset_level));
        }
        else if (opt_size) {
            // Optimize size
            opt_params_cur.alpha += optimize_size_step_max / (1<<optimize_size_level);
        }

        if (!opt_offset && !opt_size) {
            // Done with tracing/optimizing
            done = true;
            if (optimize) std::cout << "  ";
            if (optimize && measure_performance) {
                optimize = false;
                total_rays *= 10;
                continue;
            }

        }

        // We are done
        if (done) {
            std::cout << "Done." << std::endl;
            
            if (print_opt_results) 
                // Write optimized params
                std::cout << "Optimization results: "
                          << "offset=<" << serialize<10>(opt_params_cur.offset) << ">, "
                          << "alpha=" << std::setprecision(10) << std::to_string(opt_params_cur.alpha) << " "
                          << "(final size=" << std::setprecision(10) << std::to_string(cur_target_size) << ")" << std::endl;

            const auto elapsed = std::chrono::high_resolution_clock::now() - start;
            std::cout << std::endl;
            std::cout << "\t" << std::setprecision(3) << "Total time: " << elapsed << std::endl;
            std::cout << "\t" << std::to_string(opt_params_cur.stats.effective_rays) << " effective segments." << std::endl;
            std::cout << "\tEffective hit: " << std::setprecision(2) << std::to_string(opt_params_cur.stats.effective_hit*100) << "%" << std::endl;
            
            perf.beta_max = beta_max;
            perf.epsilon = epsilon;
            perf.offset = opt_params_cur.offset;
            perf.target_size = cur_target_size;
            perf.P = perf.Ai = perf.At = 0;
            perf.hitratio = opt_params_cur.stats.effective_hit;
            
            if (opt_params_cur.stats.effective_hit==1) {
                const auto effective_beam_area = mirrors.D * opt_params_cur.stats.performance;
                const auto target_area_p = type==target_type::cyl ? 
                    cur_target_size*2*glm::pi<f_t>() :
                    2*cur_target_size;
                const auto target_area = type==target_type::cyl ? 
                    cur_target_size*2*beta_max :
                    2*cur_target_size;
    
                perf.P = effective_beam_area/target_area_p * glm::sin(epsilon/2);
                perf.Ai = effective_beam_area;
                perf.At = target_area;
                perf.Atp = target_area_p;
                
                std::cout << "\tPerformance: " << std::setprecision(6) << std::to_string(perf.P) 
                    << "    (effective beam length entering system: " << std::to_string(effective_beam_area)
                    << ", target length: " << std::to_string(target_area) << ")"
                    << std::endl;
            }
            
            // Write things
            if (write_out) {
                std::cout << std::endl << "Writing data..." <<std::endl;
                std::filesystem::create_directory(outpath);
                
                write_performance(outpath/"performance.csv", { perf });
                if (target_intensity_distribution)
                    rt::write_intdist_data_csv(outpath/"intensity_distribution.csv", *target_intensity_distribution);

                mirrors::write_mirror_data_csv(outpath/"mirror1.csv", mirrors.m1);
                mirrors::write_mirror_data_csv(outpath/"mirror2.csv", mirrors.m2);

                mirrors::write_mirror_R_of_beta_csv(outpath/"R_of_beta.csv", mirrors);

                rt::write_rays_data_csv(outpath/"raysa.csv", raysa);
                rt::write_rays_data_csv(outpath/"raysb.csv", raysb);
                write_target_csv(outpath/"target.csv", type, opt_params_cur.offset, opt_params_cur.alpha, target_size);
            }

            break;
        }
    }
    
    return perf;
}

void print_usage() {
	std::cout << "Usage: tracer <options>" << std::endl;
	std::cout << "-o <output>, --out=<output>\t\t\tSpecify output directory path (required)." << std::endl;
	std::cout << "-f <focal>, --focal=<focal>\t\t\tFocal length of mirror 2, or mirror 1 for setup 3 (required)." << std::endl;
	std::cout << "--focal2=<focal>\t\t\t\tFocal length of mirror 2 (required for setup 3)." << std::endl;
	std::cout << "--setup=<setup>\t\t\t\tMirror setup: 0/1/2/3 (Defaults to 1)." << std::endl;
	std::cout << "--beam_in=<angle>\t\t\tAngle of incident beam, for setup 1/2 only (defaults to 0)." << std::endl;
	std::cout << "-L <dist>, --distance=<dist>\t\t\tDistance of mirror 1 from focus (required)." << std::endl;
	std::cout << "-s <step>, --dY=<step>\t\t\t\tNumeric integration step (defaults to 1e-2)." << std::endl;
	std::cout << "-m <degrees>, --beta_max=<degrees>\t\tNumeric integration end (defaults to 60)." << std::endl;
	std::cout << "-p <threads>, --parallel=<threads>\t\tThread pool threads (defaults to native concurrency)." << std::endl;
	std::cout << "-t <size>, --target=<size>\t\t\tSize of target, radius for cyliderical target, half edge legth for flat target (computed from epsilon and f if unset)." << std::endl;
	std::cout << "-h <shift>, --shift=<shift>\t\t\tOffset of target from origin (Defaults to 0,0,0)." << std::endl;
	std::cout << "-e <angle>, --eosilon=<angle>\t\t\tDiffusivitiy of incident radiation (defaults to 1e-2)." << std::endl;
	std::cout << "-r <count>, --rays=<count>\t\t\tCount of rays to trace (defaults to 512)." << std::endl;
	std::cout << "-y <type>, --type=<type>\t\t\tTarget type (defaults to \"cylinder\")." << std::endl;
	std::cout << "--source_type=<type>\t\t\tSource type, for setup 3 (defaults to \"cylinder\")." << std::endl;
	std::cout << "--performance\t\t\t\t\tPerforms Monte-Carlo ray tracing to measure performance (defaults to FALSE)." << std::endl;
	std::cout << "-b\t\t\t\t\t\tIntersect backfaces (defaults to FALSE)." << std::endl;
	std::cout << "-x\t\t\t\t\t\tOnly record final traced segment (defaults to FALSE)." << std::endl;
	std::cout << "-q\t\t\t\t\t\tNumerically look for target offset and size compensation local minima." << std::endl;
	std::cout << "--no-offset\t\t\t\tDisable shift optimization (defaults to FALSE)." << std::endl;
	std::cout << "--no-size-opt\t\t\t\tDisable size optimization (defaults to FALSE)." << std::endl;
	std::cout << "-k\t\t\t\t\t\tKeep tracing after hitting target." << std::endl;
}

int main(int argc, char* argv[]) {
	feenableexcept(FE_OVERFLOW);
	/* feenableexcept(FE_INVALID | FE_OVERFLOW); */

	// Output opts
	std::filesystem::path outpath{};
    f_t f1{},f2{},L{},
        dY{ 1e-2 },
        beam_in{ 0 };
    range<f_t> beta_maxs{ 60 },
        epsilons{ 1e-2 };
    point_t target_offset{ 0 };
    std::optional<range<f_t>> desired_target_sizes;
    std::size_t p{ std::thread::hardware_concurrency() },
                total_rays{ 512 },
                mirror_setup{ 1 };
    bool intersect_back_faces{ false },
         only_final_segment{ false },
         keep_tracing_after_target{ false },
         optimize{ false },
         measure_performance{ false },
         disable_offset_optimization{ false },
         disable_size_optimization{ false };
    target_type type{ target_type::cyl },
                src_type{ target_type::cyl };

	// Parse intput
	{
		option long_options[] = {
			{"help", no_argument, nullptr, '?'},
			{"performance", no_argument, nullptr, 1},
			{"no-offset", no_argument, nullptr, 10},
			{"no-size-opt", no_argument, nullptr, 11},
			{"out", required_argument, nullptr, 'o'},
			{"focal", required_argument, nullptr, 'f'},
			{"focal2", required_argument, nullptr, 2},
			{"setup", required_argument, nullptr, 3},
			{"distance", required_argument, nullptr, 'L'},
			{"dY", required_argument, nullptr, 's'},
			{"beta_max", required_argument, nullptr, 'm'},
			{"parallel", required_argument, nullptr, 'p'},
			{"target", required_argument, nullptr, 't'},
			{"shift", required_argument, nullptr, 'h'},
			{"epsilon", required_argument, nullptr, 'e'},
			{"rays", required_argument, nullptr, 'r'},
			{"type", required_argument, nullptr, 'y'},
			{"source_type", required_argument, nullptr, 4},
			{"beam_in", required_argument, nullptr, 5},
			{nullptr, no_argument, nullptr, 'b'},
			{nullptr, no_argument, nullptr, 'x'},
			{nullptr, no_argument, nullptr, 'k'},
			{nullptr, no_argument, nullptr, 'q'},
			{nullptr, 0, nullptr, 0}
		};
		char ch;
		while ((ch = getopt_long(argc, argv, "o:f:L:s:h:p:m:t:e:r:y:bxkq?", long_options, nullptr)) != -1) {
			switch (ch) {
			case 'o':
				outpath = optarg; break;
            case 1:
                measure_performance = true; break;
            case 3:
                mirror_setup = atoi(optarg); break;
            case 5:
                beam_in = atof(optarg) * glm::pi<f_t>() / f_t(180); break;
			case 'f':
				f1 = atof(optarg); break;
			case 2:
				f2 = atof(optarg); break;
			case 'L':
				L = atof(optarg); break;
			case 's':
				dY = atof(optarg); break;
			case 'm':
				beta_maxs = range<f_t>(std::string{ optarg }); break;
			case 'e':
				epsilons = range<f_t>(std::string{ optarg }); break;
			case 't':
				desired_target_sizes = range<f_t>(std::string{ optarg }); break;
			case 'h': {
                    auto opt = optarg;
                    f_t x,y,z;
                    x = atof(opt);
                    if (opt=strstr(opt,","), opt)
                        y = atof(opt+1);
                    if (opt=strstr(opt+1,","), opt)
                        z = atof(opt+1);
                    if (!opt) {
                        std::cerr << "Incorrect option to argument shift" << std::endl << std::endl;
                        print_usage();
                        return 1;
                    }
                    target_offset = point_t{ x,y,z };
                } break;
			case 'p':
				p = atoi(optarg); break;
			case 'r':
				total_rays = atoi(optarg); break;
            case 'y': {
                    type = target_type::unknown;
                    if (std::string{ optarg }=="cylinder" || std::string{ optarg }=="cyl")
                        type = target_type::cyl;
                    else if (std::string{ optarg }=="flat")
                        type = target_type::flat;
                } break;
            case 4: {
                    src_type = target_type::unknown;
                    if (std::string{ optarg }=="cylinder" || std::string{ optarg }=="cyl")
                        src_type = target_type::cyl;
                    else if (std::string{ optarg }=="flat")
                        src_type = target_type::flat;
                } break;
            case 'b':
                intersect_back_faces = true; break;
            case 'x':
                only_final_segment = true; break;
            case 'k':
                keep_tracing_after_target = true; break;
            case 'q':
                optimize = true; break;
            case 10:
                disable_offset_optimization = true; break;
            case 11:
                disable_size_optimization = true; break;
			default:
				std::cerr << "Unknown option!" << std::endl << std::endl;
            case '?':
				print_usage();
				return 1;
			}
		}

		if (!outpath.string().length() || type==target_type::unknown || f1<=0 || 
            (mirror_setup==3 && (f2<=0 || src_type==target_type::unknown)) || 
            L<0 || dY<=0 || beta_maxs[0]<=0 || p==0 || epsilons[0]<0) {
			print_usage();
			return 1;
		}
	}
    if (mirror_setup != 3)
        f2=f1;
    if (mirror_setup==0)
        intersect_back_faces = true;
    

    std::filesystem::create_directory(outpath);
    
    const auto has_array_output = epsilons.size()>1 || beta_maxs.size()>1;
    auto perf_out = write_performance_header(outpath);
    for (const auto &desired_target_size : desired_target_sizes.value_or(range<f_t>{.0}))
    for (const auto &epsilon : epsilons)
    for (const auto &beta_max_deg : beta_maxs) {
        const auto beta_max = beta_max_deg * glm::pi<f_t>() / f_t(180);
        std::cout << "epsilon=" << std::to_string(epsilon) << " beta=" + std::to_string(beta_max_deg) << "..." << std::endl;
        const auto r = work(outpath,
                            !has_array_output,
                            f1, 
                            f2,
                            L,
                            dY,
                            beta_max,
                            epsilon,
                            beam_in,
                            target_offset,
                            desired_target_size>0 ? std::optional<f_t>{ desired_target_size } : std::nullopt,
                            p,
                            total_rays,
                            mirror_setup,
                            intersect_back_faces,
                            only_final_segment,
                            keep_tracing_after_target,
                            measure_performance,
                            optimize,
                            disable_offset_optimization,
                            disable_size_optimization,
                            type,
                            src_type);
        write_performance(perf_out, r);
        std::cout << std::endl;
    }
    
    
    return 0;
}
