
#include <fstream>
#include <string>
#include <atomic>
#include <numeric>

#include "common.hpp"
#include "raytracer.hpp"
#include "intersect.hpp"


using namespace rt;

traced_lines_t trace_rays(std::vector<gen_ray_t> &&vin, 
                          const octree_t &oct,
                          intensity_distribution_t *intdist,
                          const bool intersect_back_faces,
                          const bool only_final_segment,
                          const bool trace_after_target) {
    struct traced_segment_t { point_t a,b; };
    
    if (intdist) 
        intdist->bins.reserve(vin.size());

    traced_lines_t traced;
    std::vector<traced_segment_t> path;
    for (const auto &rin : vin) {
        auto r = rin.first;
        const auto R = rin.second;  // Distance from beam center on beam tangent
        if (glm::dot(r.d,r.d)<.1) continue;

        bool finish_tracing{}, hit_ball{}, hit_mirror1{}, hit_mirror2{};
        while (!finish_tracing) {
            std::optional<std::pair<intersection_t<f_t>, const primitive_t*>> shortest;
            f_t m1_hit_r{};
            oct.trace(r,[&](const auto *obj) {
                // Ignore target before reflecting of mirror 2
                if (obj->is_target && !hit_mirror2) return true;
                // Ignore mirrors after reflecting of mirror 2
                if (!obj->is_target && hit_mirror2) return true;
                if (!intersect_back_faces && !hit_mirror1 && obj->mirror!=1) return true;
                if (!intersect_back_faces && hit_mirror1 && !hit_mirror2 && obj->mirror!=2) return true;

                const auto intr = obj->type==primitive_type::ball ? 
                    intersect_ray_sphere(r, *(ball_t*)obj) : obj->type==primitive_type::line ?
                    intersect_ray_line(r, *(line_t*)obj, intersect_back_faces) :
                    intersect_ray_tri(r, *(tri_t*)obj, intersect_back_faces);

                // No intersection?
                if (!intr) return !obj->is_target;

                if (!intersect_back_faces || !shortest || intr.value().d<shortest.value().first.d)
                    shortest = std::make_pair(intr.value(),obj);
                return !intersect_back_faces || !obj->is_target;
            });
            
            const f_t dist = shortest ? shortest.value().first.d : f_t(0);
            const auto p2  = r.p+r.d*dist;
            const auto p2far = r.p+r.d*(hit_mirror2?f_t(2):2*far);
            
            // Handle bookkeeping
            finish_tracing = !shortest;
            if (shortest) {
                if (shortest.value().second->is_target) {
                    hit_ball = true;
                    finish_tracing = true;

                    if (intdist) {
                        // Record intensity distribution point
                        auto theta = glm::atan(p2.y,(hit_mirror1?1:-1)*p2.x);
                        auto gamma = glm::acos(glm::dot(shortest.value().first.n,-r.d));
                        if (theta<0) theta+=2*glm::pi<f_t>();

                        const auto ed = shortest.value().second->type==primitive_type::ball ? 
                                            glm::length(closest_ray_point(r,{})) :
                                            glm::length(p2);
                        
                        std::unique_lock<std::mutex> l(intdist->m);
                        intdist->bins.emplace_back(theta,gamma,ed);
                    }
                }
                else if (shortest.value().second->mirror==1)
                    hit_mirror1 = true;
                else if (shortest.value().second->mirror==2)
                    hit_mirror2 = true;
            }

            path.emplace_back(traced_segment_t{ r.p, trace_after_target ? p2far : p2 });
                
            if (!finish_tracing) {
                r.p = p2;
                r.d = shortest.value().first.refl;
            }
        }

        // Record rays
        const auto result = hit_ball ? traced_line_result::hit : 
                            hit_mirror2 ? traced_line_result::miss_target :
                            hit_mirror1 ? traced_line_result::miss_mirror2 : traced_line_result::miss_mirror1;
        if (!only_final_segment)
            for (const auto &l : path) traced.emplace_back(l.a, l.b, result);
        else if (path.size())          traced.emplace_back(path.back().a, path.back().b, result);
        path.clear();
    }

    return traced;
}

traced_lines_t rt::launch_ray_packets(thread_pool &tp, 
                                      const std::function<std::optional<gen_ray_t>(std::size_t)> &ray_generator, 
                                      const std::function<void(f_t)> &progf,
                                      const octree_t &oct,
                                      intensity_distribution_t *intdist,
                                      const bool intersect_back_faces,
                                      const bool only_final_segment,
                                      const bool trace_after_target) {
    static constexpr auto rays_per_packet = 8;

    // Enqueue packets
    std::vector<std::future<traced_lines_t>> futures;
    std::vector<gen_ray_t> vins{};
    bool end=false;
    std::atomic<std::size_t> done{ 0 }, packets{ 0 };
    std::size_t ri=0;
    for (auto ray=ray_generator(ri);!end;ray=ray_generator(++ri)) {
        end = !ray;
        if (!end)
            vins.emplace_back(ray.value());
        if (end || vins.size() == rays_per_packet) {
            ++packets;
            futures.emplace_back(tp.enqueue(
                [vins=std::move(vins),&done,&packets,&oct,&progf,intdist,intersect_back_faces,only_final_segment,trace_after_target]() mutable { 
                    auto&& ret = trace_rays(std::move(vins),oct,intdist,intersect_back_faces,only_final_segment,trace_after_target);
                    progf(std::min<f_t>(1-1e-10,done++/(f_t)packets.load(std::memory_order_acquire)));
                    return ret;
                }
            ));
            vins.clear();
        }
    }
    
    // Wait for packets
    traced_lines_t traced;
    for (auto &&f : futures) {
        auto &&rs = f.get();
        traced.insert(traced.end(), rs.begin(), rs.end());
    }

    return traced;
}

inline auto& operator<<(std::ostream &o, const traced_line_t &l) {
    const auto status = l.result == traced_line_result::hit ? std::string{ "hit" } :
                        l.result == traced_line_result::miss_target ? std::string{ "misstarget" } :
                        l.result == traced_line_result::miss_mirror1 ? std::string{ "missm1" } : std::string{ "missm2" };
    o << serialize<8>(l.a) << "," << serialize<8>(l.b) << ",\"" << status << "\"";
    return o;
}
void rt::write_rays_data_csv(const std::filesystem::path& out, const traced_lines_t& traced) {
    std::ofstream f(out);
    f << "x1,y1,z1,x2,y2,z2,result" << std::endl;
    if (!traced.size())
        f << traced_line_t{ point_t{-100}, point_t{-100}, traced_line_result::miss_mirror1 } << std::endl;
    for (const auto &l : traced)
        f << l << std::endl;
}

void rt::write_intdist_data_csv(const std::filesystem::path& out, const intensity_distribution_t& intdist) {
    std::ofstream f(out);
    f << "theta,gamma,effective_distance" << std::endl;
    for (const auto& h : intdist.bins)
        f << std::fixed << std::setprecision(10) 
          << std::to_string(h.theta) << "," 
          << std::to_string(h.gamma) << ","
          << std::to_string(h.effective_distance) << std::endl;
}
