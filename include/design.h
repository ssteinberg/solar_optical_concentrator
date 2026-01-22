#pragma once

#include "path.h"
#include "geometry/geometry.h"

struct Design {
    Geometry* source;
    Geometry* target;
    std::vector<Geometry*> geometries;

    void addGeometry(Geometry* const g);

    bool intersect(const Ray& ray, HitInfo& minHitInfo) const;


    Path traceRay(const Ray& ray) const;

    std::vector<Path> rayTrace(const std::vector<Ray>& rays) const;


    void traceMeanRays(const std::string& filePath, const int& numRays);

    void traceDiffuseRays(const std::string& filePath, const int& numRays) const;


    void traceExtremeDiffuseRays(const std::string& filePath) const;

    void traceExtremeInfiniteRays(const std::string& filePath, const int& numRays);

    void traceExtremeDiffuseRaysEllipse(const std::string& filePath);

    void traceExtremeDiffuseRaysParabola(const std::string& filePath);


    void tracePhaseSpace(const std::string& filePath, const int& numRays) const;

    void tracePhaseSpaceEllipse(const std::string& filePath, const int& numRays) const;

    void tracePhaseSpaceParabola(const std::string& filePath, const int& numRays) const;


    float_type traceHitData(const int& numRays) const;

    float_type traceHitDataEllipse(const int& numRays) const;

    float_type traceHitDataParabola(const int& numRays) const;


    void traceFinalPlotRays(const std::string& filePath) const;
};
