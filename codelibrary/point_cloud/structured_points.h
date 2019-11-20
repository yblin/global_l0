//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_POINT_CLOUD_STRUCTURED_POINTS_H_
#define GEOMETRY_POINT_CLOUD_STRUCTURED_POINTS_H_

#include <cassert>
#include <cfloat>
#include <cmath>
#include <set>

#include "codelibrary/geometry/kernel/angle.h"
#include "codelibrary/geometry/kernel/center_2d.h"
//#include "codelibrary/geometry/kernel/center_3d.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/intersect_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/transform_3d.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace point_cloud {

/**
 * Point set with structure.
 * This algorithm takes advantage of a set of detected planes: it detects
 * adjacency relationships between planes and resamples the detected planes and
 * edges to produce a structured point set.
 *
 * There are three kind of points are extracted.
 * Planes:  inliers of each detected plane are replaced by sets of noise-free
 *          points sampled at the vertices of a regular grid: this is achieved
 *          by filling an occupancy grid aligned on the principal components of
 *          the inlier sets with a spacing determined by user.
 * Edges:   adjacencies between 2 planes are detected and regularly resampled on
 *          an occupancy array along the edge with a spacing determined by user.
 *
 * The size parameter `epsilon` is used both for detecting adjacencies and for
 * setting the sampling density of the structured point set.
 *
 * This algorithm is well suited to point sets sampled on surfaces with planar
 * sections and sharp edges.
 *
 * The following code is adapted from CGAL(www.cgal.org).
 *
 * Reference:
 *   Lafarge F, Alliez P. Surface Reconstruction through Point Set Structuring.
 *   Computer Graphics Forum, 2013, 32:225-234.
 */
class StructuredPoints {
public:
    struct PlanePoint {
        PlanePoint() = default;

        PlanePoint(int id, const RPoint3D& p)
            : plane_id(id), point(p) {}

        int plane_id = -1;
        RPoint3D point;
    };

    struct EdgePoint {
        EdgePoint() = default;

        EdgePoint(int id, const RPoint3D& p)
            : edge_id(id), point(p) {}

        int edge_id = -1;
        RPoint3D point;
    };

    struct Edge {
        Edge() = default;

        Edge(int id1, int id2, const RSegment3D& l)
            : plane_id1(id1), plane_id2(id2), line(l) {}

        int plane_id1 = -1, plane_id2 = -1;
        RSegment3D line;
    };

    template <typename Point>
    StructuredPoints(const Array<Point>& points,
                     const Array<RPlane3D>& planes,
                     const Array<int>& labels,
                     double epsilon,
                     double attraction_factor = 3.0)
        : planes_(planes),
          epsilon_(epsilon),
          attraction_factor_(attraction_factor) {
        assert(points.size() == labels.size());
        assert(epsilon_ > 0.0);
        assert(attraction_factor_ >= 1.0);

        // Assign each point to its plane.
        assigned_points_.resize(planes_.size());
        for (int i = 0; i < labels.size(); ++i) {
            if (labels[i] == -1) continue;
            assert(labels[i] >= 0 && labels[i] < planes_.size());

            assigned_points_[labels[i]].push_back(i);
        }

        ResamplePoints(points);
        ComputeEdgePoints();
    }

    /**
     * Return planar points.
     */
    const Array<PlanePoint>& plane_points() const {
        return plane_points_;
    }

    /**
     * Return edge points.
     */
    const Array<EdgePoint>& edge_points() const {
        return edge_points_;
    }

    /**
     * Return the planes.
     */
    const Array<RPlane3D>& planes() const {
        return planes_;
    }

    /**
     * Return the edges.
     */
    const Array<Edge>& edges() const {
        return edges_;
    }

private:
    /**
     * Sample noise-free points on each detected plane.
     */
    template <typename Point>
    void ResamplePoints(const Array<Point>& points) {
        const double grid_length = epsilon_;

        for (int k = 0; k < planes_.size(); ++k) {
            const RPlane3D& plane = planes_[k];
            Quaternion rot(plane.normal(), RVector3D(0.0, 0.0, 1.0));
            Quaternion reverse_rot = rot.Reverse();

            RPoint3D q = Rotate(plane.point(), rot);
            double z = q.z;

            Array<RPoint2D> projected_points;
            for (int j : assigned_points_[k]) {
                RPoint3D p = Rotate(points[j], rot);
                projected_points.emplace_back(p.x, p.y);
            }

            RBox2D box(projected_points.begin(), projected_points.end());

            double x_len = box.x_length() / epsilon_;
            double y_len = box.y_length() / epsilon_;
            assert(x_len < INT_MAX && "The value of 'epsilon' is too small.");
            assert(y_len < INT_MAX && "The value of 'epsilon' is too small.");

            int n_x = static_cast<int>(x_len) + 1;
            int n_y = static_cast<int>(y_len) + 1;
            assert(n_x > 0 && n_y > 0 && n_x <= INT_MAX / n_y &&
                   "The value of 'epsilon' is too small.");

            ArrayND<bool> mask(n_x, n_y), mask_border(n_x, n_y);
            mask.fill(false);
            mask_border.fill(false);
            ArrayND<Array<RPoint2D> > point_map(n_x, n_y);

            // Store the projected points in a 2D-grid.
            for (const RPoint2D& p : projected_points) {
                int ind_x = static_cast<int>((p.x - box.x_min()) / grid_length);
                int ind_y = static_cast<int>((p.y - box.y_min()) / grid_length);
                ind_x = Clamp(ind_x, 0, n_x - 1);
                ind_y = Clamp(ind_y, 0, n_y - 1);

                mask(ind_x, ind_y) = true;
                point_map(ind_x, ind_y).push_back(p);
            }

            // Hole filing in mask in 4-connectivity.
            for (int j = 1; j < n_y - 1; ++j)
                for (int i = 1; i < n_x - 1; ++i)
                    if (!mask(i, j) && mask(i - 1, j) && mask(i, j - 1) &&
                        mask(i, j + 1) && mask(i + 1, j))
                        mask(i, j) = true;

            // Finding mask border in 8-connectivity.
            for (int j = 1; j < n_y - 1; ++j)
                for (int i = 1; i < n_x - 1; ++i)
                    if (mask(i, j) && (!mask(i - 1, j - 1) || !mask(i - 1, j) ||
                        !mask(i - 1, j + 1) || !mask(i, j - 1) ||
                        !mask(i, j + 1) || !mask(i + 1, j - 1) ||
                        !mask(i + 1, j) || !mask(i + 1, j + 1)))
                        mask_border(i, j) = true;

            for (int j = 0; j < n_y; ++j) {
                if (mask(0, j)) mask_border(0, j) = true;
                if (mask(n_x - 1, j)) mask_border(n_x - 1, j) = true;
            }

            for (int i = 0; i < n_x; ++i) {
                if (mask(i, 0)) mask_border(i, 0) = true;
                if (mask(i, n_y - 1)) mask_border(i, n_y - 1) = true;
            }

            for (int j = 0; j < n_y; ++j) {
                for (int i = 0; i < n_x; ++i) {
                    if (!point_map(i, j).empty()) {
                        if (!mask_border(i, j)) {
                            // Inside: use the center of the cell.
                            double x = (i + 0.5) * grid_length + box.x_min();
                            double y = (j + 0.4) * grid_length + box.y_min();

                            if (i % 2 == 1) {
                                x = (i + 0.5) * grid_length + box.x_min();
                                y = (j + 0.6) * grid_length + box.y_min();
                            }

                            RPoint3D p(x, y, z);
                            p = Rotate(p, reverse_rot);
                            plane_points_.emplace_back(k, p);
                        } else {
                            // Border: use the barycenter of the points in the
                            // cell.
                            const Array<RPoint2D>& ps = point_map(i, j);
                            RPoint2D p = Centroid(ps);
                            RPoint3D q(p.x, p.y, z);
                            q = Rotate(q, reverse_rot);
                            plane_points_.emplace_back(k, q);
                        }
                    } else if (point_map(i, j).empty() &&
                               !mask_border(i, j) && mask(i, j)) {
                        // Points used to filling 4-connectivity holes.
                        double x = (i + 0.5) * grid_length + box.x_min();
                        double y = (j + 0.49) * grid_length + box.y_min();
                        if (i % 2 == 1) {
                            x = (i + 0.5) * grid_length + box.x_min();
                            y = (j + 0.51) * grid_length + box.y_min();
                        }
                        RPoint3D p(x, y, z);
                        p = Rotate(p, reverse_rot);
                        plane_points_.emplace_back(k, p);
                    }
                }
            }
        }
    }

    /**
     * Find adjacent planes.
     */
    void FindAdjacentPlanes() {
        const double radius = epsilon_ * attraction_factor_;

        Array<RPoint3D> points;
        for (const PlanePoint& p : plane_points_) {
            points.push_back(p.point);
        }
        KDTree<RPoint3D> kd_tree;
        kd_tree.SwapPoints(&points);

        Array<std::set<int> > adjacent(planes_.size());

        Array<int> neighbors;
        for (int i = 0; i < kd_tree.size(); ++i) {
            kd_tree.FindRadiusNeighbors(kd_tree.points()[i], radius * radius,
                                        &neighbors);
            int p1 = plane_points_[i].plane_id;
            for (int j : neighbors) {
                int p2 = plane_points_[j].plane_id;
                if (p1 != p2) {
                    adjacent[p1].insert(p2);
                    adjacent[p2].insert(p1);
                }
            }
        }

        adjacent_planes_.resize(planes_.size());
        for (int i = 0; i < planes_.size(); ++i) {
            for (int j : adjacent[i]) {
                adjacent_planes_[i].push_back(j);
            }
        }
    }

    /**
     * Compute edge points between pairs of adjacent planes.
     */
    void ComputeEdgePoints() {
        FindAdjacentPlanes();

        const double d_delta_edge = epsilon_;
        const double r_edge = d_delta_edge * 0.5;
        const double radius = epsilon_ * attraction_factor_;

        // Indices[i] stores the indices of points of the i-th plane.
        Array<Array<int> > indices(planes_.size());
        for (int i = 0; i < plane_points_.size(); ++i) {
            const PlanePoint& p = plane_points_[i];
            indices[p.plane_id].push_back(i);
        }

        Array<bool> removed(plane_points_.size(), false);

        for (int i = 0; i < adjacent_planes_.size(); ++i) {
            for (int j : adjacent_planes_[i]) {
                if (j <= i) continue;
                const RPlane3D& plane1 = planes_[i];
                const RPlane3D& plane2 = planes_[j];
                int two_planes[2] = {i, j};
                RPlane3D planes[2] = {plane1, plane2};

                RLine3D line;
                if (!geometry::Cross(plane1, plane2, &line)) continue;

                Array<RPoint3D> intersection[2];
                for (int k : indices[i]) {
                    const RPoint3D& p = plane_points_[k].point;
                    if (Distance(p, line) < radius) {
                        intersection[0].push_back(p);
                    }
                }
                if (intersection[0].empty()) continue;

                for (int k : indices[j]) {
                    const RPoint3D& p = plane_points_[k].point;
                    if (Distance(p, line) < radius) {
                        intersection[1].push_back(p);
                    }
                }
                if (intersection[1].empty()) continue;

                PointDirectionCompare3D compare(line.direction());
                RPoint3D p_min1[2], p_max1[2];
                for (int i = 0; i < 2; ++i) {
                    p_min1[i] = *std::min_element(intersection[i].begin(),
                                                  intersection[i].end(),
                                                  compare);
                    p_max1[i] = *std::max_element(intersection[i].begin(),
                                                  intersection[i].end(),
                                                  compare);
                }
                RPoint3D p_min = std::max(p_min1[0], p_min1[1], compare);
                RPoint3D p_max = std::min(p_max1[0], p_max1[1], compare);
                p_min = Project(p_min, line);
                p_max = Project(p_max, line);
                RSegment3D seg(p_min, p_max);
                edges_.emplace_back(i, j, seg);

                int n_division = static_cast<int>(seg.length() / d_delta_edge);
                Array<Array<RPoint3D> > division_tab[2];
                for (int i = 0; i < 2; ++i) {
                    division_tab[i].resize(n_division + 1);
                    for (const RPoint3D& p : intersection[i]) {
                        RPoint3D q = Project(p, line);
                        double dis = Distance(q, seg.lower_point());
                        int tab_index = static_cast<int>(dis / d_delta_edge);
                        tab_index = Clamp(tab_index, 0, n_division);
                        division_tab[i][tab_index].push_back(p);
                    }
                }

                for (int k = 0; k <= n_division; ++k) {
                    if (division_tab[0][k].empty() ||
                        division_tab[1][k].empty()) continue;

                    RPoint3D p = seg.lower_point() + seg.direction() *
                                 (k + 0.5) * (1.0 / (n_division + 1));
                    edge_points_.emplace_back(edges_.size() - 1, p);
                }

                // Remove the plane points that too close to the edge.
                for (int k : indices[i]) {
                    const RPoint3D& p = plane_points_[k].point;
                    if (Distance(p, seg) < d_delta_edge) {
                        removed[k] = true;
                    }
                }
                for (int k : indices[j]) {
                    const RPoint3D& p = plane_points_[k].point;
                    if (Distance(p, seg) < d_delta_edge) {
                        removed[k] = true;
                    }
                }

                // Create the plane points that near the edge.
                for (int k = 0; k < n_division; ++k) {
                    RPoint3D anchor = seg.lower_point() + seg.direction() *
                                      (k + 1) * (1.0 / (n_division + 1));
                    RPlane3D ortho(anchor, seg.direction());

                    for (int t = 0; t < 2; ++t) {
                        if (division_tab[t][k].empty()) continue;

                        RLine3D l;
                        if (!geometry::Cross(planes[t], ortho, &l)) continue;

                        RVector3D v = l.direction();
                        v *= 1.0 / v.norm();
                        RVector3D v1 = Centroid(division_tab[t][k]) - anchor;
                        if (v1 * v < 0.0) v = -v;
                        plane_points_.emplace_back(two_planes[t],
                                                   anchor + v * r_edge);
                    }
                }
            }
        }

        // Filter the plane points.
        int n_plane_points = 0;
        for (int i = 0; i < removed.size(); ++i) {
            if (!removed[i]) plane_points_[n_plane_points++] = plane_points_[i];
        }
        for (int i = removed.size(); i < plane_points_.size(); ++i) {
            plane_points_[n_plane_points++] = plane_points_[i];
        }
        plane_points_.resize(n_plane_points);
    }

    // The input planes.
    Array<RPlane3D> planes_;

    // Used both for detecting adjacencies and for setting the sampling density
    // of the structured point set.
    double epsilon_;

    // Used to attract the points to the edge.
    double attraction_factor_;

    // Assigned points of each plane.
    Array<Array<int> > assigned_points_;

    // Planar points.
    Array<PlanePoint> plane_points_;

    // Edge points.
    Array<EdgePoint> edge_points_;

    // Edges betwee two adjacent planes.
    Array<Edge> edges_;

    // Adjacent graph of planes.
    Array<Array<int> > adjacent_planes_;
};

} // namespace point_cloud
} // namespace cl

#endif // GEOMETRY_POINT_CLOUD_STRUCTURED_POINTS_H_
