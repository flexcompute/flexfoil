// mesh2d: mesh a 2D multi-element airfoil from closed contour curves using the
// in-house CavityBasedMesher, driven by an a-priori Spalding boundary-layer metric.
// Part of flexfoil/rans. Input: contours.txt (see rans/contours.py). Output: legacy
// VTK triangle mesh. Build: mesher/build.sh.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Libraries/MeshPrimitives/Vector3.h"
#include "Libraries/MeshPrimitives/Curve.h"
#include "Libraries/GeometryKernel/AnalyticGeometry/PlaneGeometry.h"
#include "Libraries/MeshDataStructures/SurfaceMesh.h"
#include "Libraries/MeshDataStructures/CurveProjection.h"
#include "Libraries/Metric/MetricTensor.h"
#include "Libraries/CavityBasedMesher/CavityMesh/CavityMesh.h"
#include "Libraries/CavityBasedMesher/ConstrainedDelaunayMesher/ConstrainedDelaunayMesher.h"
#include "Libraries/CavityBasedMesher/SnowflakePointGenerator/CavityRemeshPlane.h"

using namespace Flynn360;

// A-priori Spalding-style anisotropic boundary-layer metric.
// Wall-normal spacing grows linearly from h0 with distance (gradation ~ growth),
// capped at hmax; tangential spacing grows from hwall. Anisotropic near walls,
// isotropic far away.
struct SpaldingMetric {
    const std::vector<std::unique_ptr<CurveProjection>>* walls;
    double h0, growth, hwall, hmax;

    SymmetricTensor3 operator()(const Vector3& x) const {
        double bestD2 = 1e300;
        Vector3 proj{0, 0, 0};
        for (const auto& w : *walls) {
            Vector3 p = w->getProjection(x);
            double dx = x[0] - p[0], dy = x[1] - p[1];
            double d2 = dx * dx + dy * dy;
            if (d2 < bestD2) { bestD2 = d2; proj = p; }
        }
        double d = std::sqrt(bestD2);
        double hn = std::min(h0 + (growth - 1.0) * d, hmax);
        double ht = std::min(hwall + (growth - 1.0) * d, hmax);

        double nx = x[0] - proj[0], ny = x[1] - proj[1];
        double nm = std::sqrt(nx * nx + ny * ny);
        if (nm < 1e-12) { nx = 1.0; ny = 0.0; } else { nx /= nm; ny /= nm; }
        double tx = -ny, ty = nx;

        double ln = 1.0 / (hn * hn), lt = 1.0 / (ht * ht);
        double Mxx = ln * nx * nx + lt * tx * tx;
        double Myy = ln * ny * ny + lt * ty * ty;
        double Mxy = ln * nx * ny + lt * tx * ty;
        double Mzz = 1.0 / (ht * ht);  // out-of-plane (irrelevant for planar mesh)
        return SymmetricTensor3{{Mxx}, {Mxy, Myy}, {0.0, 0.0, Mzz}};
    }
};

int main(int argc, char** argv) {
    std::string inPath = argc > 1 ? argv[1] : "/tmp/poc_contours.txt";
    std::string outPath = argc > 2 ? argv[2] : "/tmp/poc_mesh.vtk";

    std::ifstream in(inPath);
    if (!in) { std::cerr << "cannot open " << inPath << "\n"; return 1; }

    double h0 = 0, growth = 0, hwall = 0, hmax = 0;
    std::vector<Vector3> points;
    std::vector<Curve> curves;
    std::vector<int> wallFlag;

    std::string tok;
    while (in >> tok) {
        if (tok == "H0") in >> h0;
        else if (tok == "GROWTH") in >> growth;
        else if (tok == "HWALL") in >> hwall;
        else if (tok == "HMAX") in >> hmax;
        else if (tok == "NPOINTS") {
            size_t n; in >> n;
            points.reserve(n);
            for (size_t i = 0; i < n; i++) { double x, y; in >> x >> y; points.push_back(Vector3{x, y, 0.0}); }
        } else if (tok == "NCURVES") {
            size_t m; in >> m;
            for (size_t c = 0; c < m; c++) {
                int isWall; size_t cnt; in >> isWall >> cnt;
                Curve cv;
                cv.nodeIndices.reserve(cnt);
                for (size_t k = 0; k < cnt; k++) { uint64_t idx; in >> idx; cv.nodeIndices.push_back(idx); }
                curves.push_back(std::move(cv));
                wallFlag.push_back(isWall);
            }
        }
    }
    std::cout << "read " << points.size() << " points, " << curves.size() << " curves; "
              << "h0=" << h0 << " growth=" << growth << " hwall=" << hwall << " hmax=" << hmax << "\n";

    auto pointsPtr = std::make_shared<std::vector<Vector3>>(points);
    SurfaceMesh cdtMesh(pointsPtr);
    SurfaceMesh outMesh(pointsPtr);

    Flow360Geometry::PlaneGeometry planeGeometry{Vector3{0.0, 0.0, 1.0}, 0.0};
    planeGeometry.curves = curves;

    // Wall projections (skip the farfield curve).
    std::vector<std::unique_ptr<CurveProjection>> walls;
    for (size_t i = 0; i < curves.size(); i++)
        if (wallFlag[i]) walls.push_back(std::make_unique<CurveProjection>(curves[i], *pointsPtr));
    SpaldingMetric metric{&walls, h0, growth, hwall, hmax};

    // 1. Initial constrained Delaunay tessellation from the contour curves.
    CavityMesh cavityMesh;
    auto mesher = makeConstrainedDelaunayMesher(cavityMesh, planeGeometry, planeGeometry.curves, *pointsPtr);
    if (!mesher.success()) {
        std::cerr << "CDT failed, status=" << static_cast<int>(mesher.getStatus()) << "\n";
        return 2;
    }
    cavityMesh.extendSurfaceMesh(cdtMesh, mesher.nodeInputMap);
    std::cout << "CDT tessellation: " << cdtMesh.triangles.size() << " triangles\n";

    // 2. Metric-driven interior remesh with the Spalding BL metric.
    cavityRemeshPlane(outMesh, cdtMesh, planeGeometry, metric, "airfoil");
    std::cout << "after Spalding remesh: " << outMesh.triangles.size() << " triangles, "
              << outMesh.points().size() << " points\n";

    // 3. Write legacy VTK (triangles) for ParaView.
    std::ofstream out(outPath);
    out << std::setprecision(17);
    const auto& P = outMesh.points();
    out << "# vtk DataFile Version 3.0\npoc multi-element 2D mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << P.size() << " double\n";
    for (const auto& p : P) out << p[0] << " " << p[1] << " " << p[2] << "\n";
    out << "CELLS " << outMesh.triangles.size() << " " << 4 * outMesh.triangles.size() << "\n";
    for (const auto& t : outMesh.triangles)
        out << "3 " << t.vertices[0] << " " << t.vertices[1] << " " << t.vertices[2] << "\n";
    out << "CELL_TYPES " << outMesh.triangles.size() << "\n";
    for (size_t i = 0; i < outMesh.triangles.size(); i++) out << "5\n";
    std::cout << "wrote " << outPath << "\n";
    return 0;
}
