#include <fstream>
#include <omp.h>

#include "design.h"

#include <unordered_set>

void Design::addGeometry(Geometry * const g) {
    geometries.push_back(g);
    if (g->type == Type::SOURCE) source = g;
    else if (g->type == Type::TARGET) target = g;
}

bool Design::intersect(const Ray &ray, HitInfo &minHitInfo, const std::optional<std::reference_wrapper<const std::unordered_set<Type>>> ignoredTypes) const {
    bool hit = false;
    HitInfo tempMinHitInfo;
    minHitInfo.l = std::numeric_limits<float_type>::max();
    for (const auto& geometry : geometries) {
        if (ignoredTypes && ignoredTypes->get().contains(geometry->type)) {
            continue;
        }

        if (geometry->intersect(ray, tempMinHitInfo)) {
            if (tempMinHitInfo.l < minHitInfo.l) {
                hit = true;
                minHitInfo = tempMinHitInfo;
            }
        }
    }
    return hit;
}

Path Design::traceRay(const Ray &ray) const {
    Path path;
    path.addVertex(ray.o);
    Ray r = ray;
    std::unordered_set initialIgnoredCollisionTypes { Type::SOURCE, Type::TARGET };
    std::unordered_set finalIgnoredCollisionTypes { Type::SOURCE };
    for (int i = 0; i < 3; ++i) {
        std::unordered_set<Type> ignoredCollisionTypes = i < 2 ? initialIgnoredCollisionTypes : finalIgnoredCollisionTypes;
        if (HitInfo h; intersect(r, h, ignoredCollisionTypes)) {
            if (i == 0 && h.t != Type::MIRROR_1a && h.t != Type::MIRROR_1b) return Path();
            if (i == 1 && h.t != Type::MIRROR_2a && h.t != Type::MIRROR_2b) return Path();
            path.addVertex(h.p);
            if (h.t == Type::TARGET) break;
            float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
            r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
        }
        else break;
    }
    return path;
}

std::vector<Path> Design::rayTrace(const std::vector<Ray> &rays) const {
    std::vector<std::vector<Path>> threadVectors(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (const auto& ray : rays) {
        int threadID = omp_get_thread_num();
        threadVectors[threadID].push_back(traceRay(ray));
    }
    std::vector<Path> paths;
    for (const auto& v : threadVectors) paths.insert(paths.end(), v.begin(), v.end());
    return paths;
}

void Design::traceMeanRays(const std::string &filePath, const int &numRays) const {
    std::vector<Ray> rays;
    for (int i = 0; i < numRays; ++i) rays.push_back(source->sampleMeanRay());
    traceMeanRays(filePath, rays);
}

void Design::traceMeanRays(const std::string &filePath, const std::vector<Ray> &rays) const {
    std::vector<Path> paths(rayTrace(rays));
    std::ofstream file;
    file.open(filePath + "mean.csv");
    for (const auto& path : paths) path.writePath(file);
    file.close();
}

void Design::traceDiffuseRays(const std::string &filePath, const int &numRays) const {
    std::vector<Ray> rays;
    for (int i = 0; i < numRays; ++i) rays.push_back(source->sampleDiffuseRay().first);
    std::vector<Path> paths(rayTrace(rays));
    std::ofstream file;
    file.open(filePath + "diffuse.csv");
    for (auto path : paths) path.writePath(file);
    file.close();
}

void Design::traceExtremeDiffuseRays(const std::string &filePath) const {
    auto extremeRays = source->generateExtremeDiffuseRays();
    auto paths = std::pair(rayTrace(extremeRays.first), rayTrace(extremeRays.second));
    std::ofstream file;
    file.open(filePath + "extrdiff1.csv");
    for (auto path : paths.first) path.writeFinalSegment(file);
    file.close();
    file.open(filePath + "extrdiff2.csv");
    for (auto path : paths.second) path.writeFinalSegment(file);
    file.close();
}

void Design::traceExtremeInfiniteRays(const std::string &filePath, const int &numRays) {
    auto extremeRays = source->generateExtremeInfiniteRays(numRays);
    auto paths = std::pair(rayTrace(extremeRays.first), rayTrace(extremeRays.second));
    std::ofstream file;
    file.open(filePath + "extrinf1.csv");
    for (auto path : paths.first) path.writeFinalSegment(file);
    file.close();
    file.open(filePath + "extrinf2.csv");
    for (auto path : paths.second) path.writeFinalSegment(file);
    file.close();
}

void Design::traceExtremeDiffuseRaysEllipse(const std::string &filePath) {
    auto extremeRays = source->generateExtremeDiffuseRays();

    std::vector<std::vector<Path>> threadVectors1(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (const auto& ray : extremeRays.first) {
        int threadID = omp_get_thread_num();
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        for (int i = 0; i < 2; ++i) {
            if (HitInfo h; intersect(r, h)) {
                path.addVertex(h.p);
                if (h.t == Type::TARGET) break;
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
        threadVectors1[threadID].push_back(path);
    }
    std::vector<Path> paths1;
    for (const auto& v : threadVectors1) paths1.insert(paths1.end(), v.begin(), v.end());

    std::vector<std::vector<Path>> threadVectors2(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (const auto& ray : extremeRays.second) {
        int threadID = omp_get_thread_num();
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        for (int i = 0; i < 2; ++i) {
            if (HitInfo h; intersect(r, h)) {
                path.addVertex(h.p);
                if (h.t == Type::TARGET) break;
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
        threadVectors2[threadID].push_back(path);
    }
    std::vector<Path> paths2;
    for (const auto& v : threadVectors2) paths2.insert(paths2.end(), v.begin(), v.end());

    std::ofstream file;
    file.open(filePath + "extrdiff1.csv");
    for (auto path : paths1) path.writeFinalSegment(file);
    file.close();
    file.open(filePath + "extrdiff2.csv");
    for (auto path : paths2) path.writeFinalSegment(file);
    file.close();
}

void Design::traceExtremeDiffuseRaysParabola(const std::string &filePath) {
    auto extremeRays = source->generateExtremeDiffuseRays();

    std::vector<std::vector<Path>> threadVectors1(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (const auto& ray : extremeRays.first) {
        int threadID = omp_get_thread_num();
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        for (int i = 0; i < 3; ++i) {
            if (HitInfo h; intersect(r, h)) {
                if (i == 0 && h.t != Type::MIRROR1) {
                    path = Path();
                    break;
                }
                if (i == 1 && h.t != Type::MIRROR2) {
                    path = Path();
                    break;
                }
                path.addVertex(h.p);
                if (h.t == Type::TARGET) break;
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
        threadVectors1[threadID].push_back(path);
    }
    std::vector<Path> paths1;
    for (const auto& v : threadVectors1) paths1.insert(paths1.end(), v.begin(), v.end());

    std::vector<std::vector<Path>> threadVectors2(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (const auto& ray : extremeRays.second) {
        int threadID = omp_get_thread_num();
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        for (int i = 0; i < 3; ++i) {
            if (HitInfo h; intersect(r, h)) {
                if (i == 0 && h.t != Type::MIRROR1) {
                    path = Path();
                    break;
                }
                if (i == 1 && h.t != Type::MIRROR2) {
                    path = Path();
                    break;
                }
                path.addVertex(h.p);
                if (h.t == Type::TARGET) break;
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
        threadVectors2[threadID].push_back(path);
    }
    std::vector<Path> paths2;
    for (const auto& v : threadVectors2) paths2.insert(paths2.end(), v.begin(), v.end());

    std::ofstream file;
    file.open(filePath + "extrdiff1.csv");
    for (auto path : paths1) path.writeFinalSegment(file);
    file.close();
    file.open(filePath + "extrdiff2.csv");
    for (auto path : paths2) path.writeFinalSegment(file);
    file.close();
}

void Design::tracePhaseSpace(const std::string &filePath, const int &numRays) const {

    std::vector<std::vector<std::vector<float_type>>> threadVectorsTarget(omp_get_max_threads());
    std::vector<std::vector<std::vector<float_type>>> threadVectorsMirror1(omp_get_max_threads());
    std::vector<std::vector<std::vector<float_type>>> threadVectorsMirror2(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        auto sample = source->sampleDiffuseRay();
        Ray r = sample.first;
        const float_type source_angle = sample.second;
        const float_vec source_pos = r.o - source->getCentre();

        // for 3 path segments
        for (int j = 0; j < 3; ++j) {

            // check intersection
            if (HitInfo h; intersect(r, h)) {

                // mirror 1
                if (j == 0 && h.t != Type::MIRROR_1a && h.t != Type::MIRROR_1b) break;
                // if (j == 0 && (h.t == Type::MIRROR_1a || h.t == Type::MIRROR_1b)) {
                //     if (source->shape == Shape::FLAT) {
                //         std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, source_pos.y};
                //         threadVectorsMirror1[threadID].emplace_back(std::move(temp));
                //     }
                //     else if (source->shape == Shape::CYLINDRICAL) {
                //         float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                //         if (source_theta < 0) source_theta += 2 * PI;
                //         std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, source_theta * RAD_TO_DEG};
                //         threadVectorsMirror1[threadID].emplace_back(std::move(temp));
                //     }
                // }

                // mirror 2
                if (j == 1 && h.t != Type::MIRROR_2a && h.t != Type::MIRROR_2b) break;
                // if (j == 1 && (h.t == Type::MIRROR_2a || h.t == Type::MIRROR_2b)) {
                //     if (source->shape == Shape::FLAT) {
                //         std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, source_pos.y};
                //         threadVectorsMirror2[threadID].emplace_back(std::move(temp));
                //     }
                //     else if (source->shape == Shape::CYLINDRICAL) {
                //         float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                //         if (source_theta < 0) source_theta += 2 * PI;
                //         std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, source_theta * RAD_TO_DEG};
                //         threadVectorsMirror2[threadID].emplace_back(std::move(temp));
                //     }
                // }

                // target
                if (j == 2 && h.t != Type::TARGET) break;
                if (j == 2 && h.t == Type::TARGET) {
                    if (source->shape == Shape::FLAT && target->shape == Shape::FLAT) {
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.p.y, std::sin(source_angle), source_pos.y};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::CYLINDRICAL && target->shape == Shape::CYLINDRICAL) {
                        float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                        if (source_theta < 0) source_theta += 2 * PI;
                        float_type target_theta = std::atan2(h.p.y, h.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r.d, h.n), target_theta * RAD_TO_DEG, std::sin(source_angle), source_theta * RAD_TO_DEG};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::FLAT && target->shape == Shape::CYLINDRICAL) {
                        float_type target_theta = std::atan2(h.p.y, h.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r.d, h.n), target_theta * RAD_TO_DEG, std::sin(source_angle), source_pos.y};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::CYLINDRICAL && target->shape == Shape::FLAT) {
                        float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                        if (source_theta < 0) source_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.p.y, std::sin(source_angle), source_theta * RAD_TO_DEG};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    break;
                }

                // reflect
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }

    std::ofstream file;

    // target
    std::vector<std::vector<float_type>> dataTarget;
    for (const auto& v : threadVectorsTarget) dataTarget.insert(dataTarget.end(), v.begin(), v.end());
    file.open(filePath + "phasetarget.csv");
    for (const auto& p : dataTarget) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    file.close();

    // mirror 1
    std::vector<std::vector<float_type>> dataMirror1;
    for (const auto& v : threadVectorsMirror1) dataMirror1.insert(dataMirror1.end(), v.begin(), v.end());
    file.open(filePath + "phasemirror1.csv");
    for (const auto& p : dataMirror1) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    file.close();

    // mirror 2
    std::vector<std::vector<float_type>> dataMirror2;
    for (const auto& v : threadVectorsMirror2) dataMirror2.insert(dataMirror2.end(), v.begin(), v.end());
    file.open(filePath + "phasemirror2.csv");
    for (const auto& p : dataMirror2) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    file.close();
}

void Design::tracePhaseSpaceEllipse(const std::string &filePath, const int &numRays) const {
    std::vector<std::vector<std::vector<float_type>>> threadVectorsTarget(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        auto sample = source->sampleDiffuseRay();
        Ray r = sample.first;
        for (int j = 0; j < 2; ++j) {
            if (HitInfo h; intersect(r, h)) {
                if (j == 1 && h.t != Type::TARGET) break;
                if (j == 1 && h.t == Type::TARGET) {
                    if (target->shape == Shape::FLAT) {
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.p.y};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    else if (target->shape == Shape::CYLINDRICAL) {
                        float_type target_theta = std::atan2(h.p.y, h.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r.d, h.n), target_theta * RAD_TO_DEG};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    break;
                }
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }
    std::ofstream file;
    std::vector<std::vector<float_type>> dataTarget;
    for (const auto& v : threadVectorsTarget) dataTarget.insert(dataTarget.end(), v.begin(), v.end());
    file.open(filePath + "phase.csv");
    for (const auto& p : dataTarget) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    file.close();
}

void Design::tracePhaseSpaceParabola(const std::string &filePath, const int &numRays) const {
    std::vector<std::vector<std::vector<float_type>>> threadVectorsTarget(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        auto sample = source->sampleDiffuseRay();
        Ray r = sample.first;
        for (int j = 0; j < 3; ++j) {
            if (HitInfo h; intersect(r, h)) {
                if (j == 0 && h.t != Type::MIRROR1) break;
                if (j == 1 && h.t != Type::MIRROR2) break;
                if (j == 2 && h.t == Type::TARGET) {
                    if (target->shape == Shape::FLAT) {
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.p.y};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    else if (target->shape == Shape::CYLINDRICAL) {
                        float_type target_theta = std::atan2(h.p.y, h.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r.d, h.n), target_theta * RAD_TO_DEG};
                        threadVectorsTarget[threadID].emplace_back(std::move(temp));
                    }
                    break;
                }
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }
    std::ofstream file;
    std::vector<std::vector<float_type>> dataTarget;
    for (const auto& v : threadVectorsTarget) dataTarget.insert(dataTarget.end(), v.begin(), v.end());
    file.open(filePath + "phase.csv");
    for (const auto& p : dataTarget) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    file.close();
}

float_type Design::traceHitData(const int &numRays) const {
    std::vector<int> threadCounts(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        Ray r = source->sampleDiffuseRay().first;
        for (int j = 0; j < 3; ++j) {
            if (HitInfo h; intersect(r, h)) {
                if (h.t == Type::TARGET) {
                    threadCounts[threadID] += 1;
                    break;
                }
                if (h.t == Type::BARRIER) break;
                // if (j == 0 && h.t != Type::MIRROR_1a && h.t != Type::MIRROR_1b) break;
                // if (j == 1 && h.t != Type::MIRROR_2a && h.t != Type::MIRROR_2b) break;
                // if (j == 2 && h.t != Type::TARGET) break;
                // if (j == 2 && h.t == Type::TARGET) {
                //     threadCounts[threadID] += 1;
                //     break;
                // }
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }
    int finalCount = 0;
    for (const auto& count : threadCounts) finalCount += count;
    return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
}

float_type Design::traceHitDataEllipse(const int &numRays) const {
    std::vector<int> threadCounts(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        Ray r = source->sampleDiffuseRay().first;
        for (int j = 0; j < 2; ++j) {
            if (HitInfo h; intersect(r, h)) {
                if (j == 1 && h.t == Type::TARGET) {
                    threadCounts[threadID] += 1;
                    break;
                }
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }
    int finalCount = 0;
    for (const auto& count : threadCounts) finalCount += count;
    return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
}

float_type Design::traceHitDataParabola(const int &numRays) const {
    std::vector<int> threadCounts(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < numRays; ++i) {
        int threadID = omp_get_thread_num();
        Ray r = source->sampleDiffuseRay().first;
        for (int j = 0; j < 3; ++j) {
            if (HitInfo h; intersect(r, h)) {
                if (j == 0 && h.t != Type::MIRROR1) break;
                if (j == 1 && h.t != Type::MIRROR2) break;
                if (j == 2 && h.t == Type::TARGET) {
                    threadCounts[threadID] += 1;
                    break;
                }
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
    }
    int finalCount = 0;
    for (const auto& count : threadCounts) finalCount += count;
    return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
}

void Design::traceFinalPlotRays(const std::string &filePath) const {
    auto rays = source->generateFinalPlotRays();
    auto paths = rayTrace(rays);
    std::ofstream file;
    file.open(filePath + "finalplotrays.csv");
    for (auto path : paths) path.writePath(file);
    file.close();
}
