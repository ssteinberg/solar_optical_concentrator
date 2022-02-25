
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <filesystem>
#include <vector>
#include <array>
#include <future>
#include <functional>
#include <optional>
#include <atomic>

#include <glm/glm.hpp>

#include "thread_pool.hpp"

#include "primitives.hpp"
#include "ray.hpp"
#include "octree.hpp"

#include "mirrors.hpp"

namespace rt {
    
using octree_t = octree::octree<f_t, primitive_t>;
using gen_ray_t = std::pair<ray_t<f_t>, f_t>;

enum class traced_line_result {
    hit,
    miss_target,
    miss_mirror1,
    miss_mirror2,
};

struct traced_line_t { point_t a,b; traced_line_result result; };
using traced_lines_t = std::vector<traced_line_t>;

struct intensity_distribution_t {
    struct bin_t {
        f_t theta, gamma;
        f_t effective_distance;
    };

    f_t target_size;
    std::mutex m;
    std::vector<bin_t> bins;
};

traced_lines_t launch_ray_packets(thread_pool &tp, 
                                  const std::function<std::optional<gen_ray_t>(std::size_t)> &ray_generator, 
                                  const std::function<void(f_t)> &progf,
                                  const octree_t &octree,
                                  intensity_distribution_t *intdist=nullptr,
                                  const bool intersect_back_faces=false,
                                  const bool only_final_segment=false,
                                  const bool trace_after_target=false);
void write_rays_data_csv(const std::filesystem::path& out, const traced_lines_t& rays);

inline auto create_intensity_distribution_bins(f_t target_size) {
    return std::make_unique<intensity_distribution_t>(target_size);
}
void write_intdist_data_csv(const std::filesystem::path& out, const intensity_distribution_t& intdist);

struct statistics_t {
    std::size_t total_rays{}, rays_entering_system{}, effective_rays{}, effective_hits{};
    f_t performance{}, effective_hit{};
};
template<std::size_t N>
auto statistics(std::array<const traced_lines_t*,N> traces) noexcept {
    statistics_t s{};
    for (const auto *rays : traces) {
        for (const auto &r : *rays) {
            ++s.total_rays;
            if (r.result!=traced_line_result::miss_mirror1)
                ++s.rays_entering_system;
            if (r.result==traced_line_result::hit || r.result==traced_line_result::miss_target)
                ++s.effective_rays;
            if (r.result==traced_line_result::hit)
                ++s.effective_hits;
        }
    }
    s.performance = s.total_rays>0 ? s.effective_rays/(f_t)s.total_rays : f_t(0);
    s.effective_hit = s.effective_rays>0 ? s.effective_hits/(f_t)s.effective_rays : f_t(0);
    return s;
}

static constexpr f_t far = 1000.;

}
