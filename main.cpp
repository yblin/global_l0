#include <map>
#include <set>

#include "codelibrary/base/time_counter.h"
#include "codelibrary/geometry/kernel/transform_3d.h"
#include "codelibrary/geometry/pca_estimate_normals_3d.h"
#include "codelibrary/point_cloud/global_l0_extractor.h"
#include "codelibrary/point_cloud/io/xyz_io.h"
#include "codelibrary/point_cloud/region_growing.h"
#include "codelibrary/statistics/kernel/mean.h"
#include "codelibrary/visualization/plot.h"

using namespace cl;

void DrawMappedNormals(const Array<RPlane3D>& planes) {
    std::set<RVector3D> normals;
    for (const RPlane3D& plane : planes) {
        if (plane.normal().z < 0.0)
            normals.insert(-plane.normal());
        else
            normals.insert(plane.normal());
    }

    plot::PointPlot point_plot;
    for (const RVector3D& normal : normals) {
        double u = 0.5 + std::atan2(normal.y, normal.x) / (2.0 * M_PI);
        double v = 1.0 - 2.0 * std::asin(normal.z) / M_PI;
        u = Clamp(u, 0.0, 1.0);
        v = Clamp(v, 0.0, 1.0);

        point_plot.Draw(RPoint2D(u, v), RGB32Color(0, 0, 255, 128), 10);
    }
    point_plot.set_axes_ticks_visibility(false);
    point_plot.set_axes_fixed(true);
    point_plot.Save("test.svg", 512, 512);
    system("test.svg");
}

/**
 * Assign each plane a color.
 */
void AssignColor(const Array<RPoint3D>& points,
                 const Array<int>& labels,
                 const Array<RPlane3D>& planes,
                 Array<RGB32Color>* colors) {
    assert(points.size() == labels.size());
    assert(colors);

    Array<RGB32Color> colors1 = cl::ColorMap::Set1().colors();

    std::map<RVector3D, int> hash;
    Array<int> labels2(points.size(), -1);
    int id = 0;
    for (int i = 0; i < points.size(); ++i) {
        if (labels[i] == -1) continue;

        const RPlane3D& plane = planes[labels[i]];
        RVector3D v = plane.normal();
        if (v.z < 0.0) v = -v;

        if (hash.find(v) == hash.end()) {
            hash[v] = id++;
        }
        labels2[i] = hash[v];
    }

    Array<int> cnt(id, 0);
    for (int i = 0; i < labels.size(); ++i) {
        if (labels2[i] == -1) continue;
        ++cnt[labels2[i]];
    }

    Array<int> index;
    IndexSort(cnt.begin(), cnt.end(), &index);
    std::reverse(index.begin(), index.end());
    Array<int> rank(planes.size());
    for (int i = 0; i < index.size(); ++i) {
        rank[index[i]] = i;
    }

    std::mt19937 random(32);
    for (int i = colors1.size(); i < planes.size(); ++i) {
        colors1.push_back(cl::RGB32Color(random()));
    }

    colors->resize(points.size());
    for (int i = 0; i < points.size(); ++i) {
        if (labels2[i] == -1) continue;

        int id = labels2[i];
        (*colors)[i] = colors1[rank[id]];
    }
}

void Metric(const Array<RPoint3D>& points,
            const Array<int>& labels,
            const Array<RPlane3D>& planes) {
    Array<double> dis1;
    for (int i = 0; i < points.size(); ++i) {
        double dis = DBL_MAX;
        for (const RPlane3D& plane : planes) {
            dis = std::min(dis, Distance(points[i], plane));
        }
        dis1.push_back(dis);
    }

    LOG(INFO) << "RMS of distances1:"
              << RootMeanSquare(dis1.begin(), dis1.end());

    Array<double> dis2;
    KDTree<RPoint3D> kd_tree(points.begin(), points.end());
    for (int i = 0; i < points.size(); ++i) {
        int l = labels[i];
        if (l == -1) continue;

        RPoint3D p1 = Project(points[i], planes[l]);
        RPoint3D p;
        kd_tree.FindNearestPoint(p1, &p);
        dis2.push_back(Distance(p, p1));
    }

    LOG(INFO) << "RMS of distances2:"
              << RootMeanSquare(dis2.begin(), dis2.end());

    int cnt = 0;
    for (int i = 0; i < points.size(); ++i) {
        if (labels[i] != -1) ++cnt;
    }
    LOG(INFO) << "Coverage: " << static_cast<double>(cnt) / points.size();
}

int main() {
    LOG_ON(INFO);

    Array<RPoint3D> points;
    point_cloud::ReadXYZPoints("foam_box.xyz", &points);

    int n_points = points.size();
    LOG(INFO) << "Import points: " << n_points;

    LOG(INFO) << "Estimate normals...";
    KDTree<RPoint3D> kd_tree;
    kd_tree.SwapPoints(&points);
    Array<RVector3D> normals;
    geometry::PCAEstimateNormals(kd_tree, 30, &normals);

    LOG(INFO) << "Start plane extracting.";
    TimeCounter timer;
    timer.Start();

    Array<RPlane3D> planes;
    Array<int> labels;
    point_cloud::GlobalL0Extractor extractor(15, 50, 3, 1.0);
    extractor.ExtractPlanes(kd_tree, normals, &planes, &labels);

    timer.Stop();
    LOG(INFO) << "Time: " << timer.elapsed_time();
    LOG(INFO) << "Extracted " << planes.size() << " planes.";

    kd_tree.SwapPoints(&points);
    kd_tree.clear();

    // Print metrics.
    Metric(points, labels, planes);

    DrawMappedNormals(planes);

    // Output points.
    Array<RGB32Color> colors;
    AssignColor(points, labels, planes, &colors);
    cl::point_cloud::WriteXYZPoints("out.xyz", points, colors);

    LOG(INFO) << "Done";

    system("out.xyz");

    return 0;
}
