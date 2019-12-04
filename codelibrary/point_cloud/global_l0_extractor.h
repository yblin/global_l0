//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_GLOBAL_L0_EXTRACTOR_H_
#define POINT_CLOUD_GLOBAL_L0_EXTRACTOR_H_

#include <cassert>
#include <queue>
#include <unordered_map>

#include "codelibrary/geometry/kernel/distance_2d.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/transform_2d.h"
#include "codelibrary/graph/graph.h"
#include "codelibrary/optimization/discrete/subset_selection/greedy_submodular.h"
#include "codelibrary/point_cloud/k_nearest_neighbors.h"
#include "codelibrary/statistics/kernel/median.h"
#include "codelibrary/util/set/disjoint_set.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace point_cloud {

/**
 * Global-L0 algorithm to extract 2D lines or 3D planes.
 */
class GlobalL0Extractor {
    typedef Graph::Edge Edge;

public:
    /**
     * Constructor.
     *
     * Parameters:
     *  k_neighbors        - number of neighbors for each point.
     *  min_support_points - the minimum number of points to support a plane.
     *  n_constraints      - the number of constraints normals.
     */
    GlobalL0Extractor(int k_neighbors,
                      int min_support_points,
                      int n_constraints,
                      double outlier_penalty)
        : k_neighbors_(k_neighbors),
          min_support_points_(min_support_points),
          n_constraints_(n_constraints),
          outlier_penalty_(outlier_penalty) {
        assert(k_neighbors_ > 0);
        assert(min_support_points_ >= 3);
        assert(n_constraints_ > 0);
        assert(outlier_penalty_ >= 0.0 && outlier_penalty_ <= 10.0);
    }

    /**
     * Extract regular lines from a 2D unstructured point set with unoriented
     * normal vectors.
     *
     * Parameters:
     *  kd_tree - the input point set is stored in KD-tree.
     *  normals - the input unoriented normal vectors.
     *  lines   - the extracted line segments.
     *  labels  - labels[i] denotes the line that own the i-th point, or
     *            labels[i] = -1 if no line of the i-th point is found.
     */
    template <typename T>
    void ExtractLines(const KDTree<Point2D<T>>& kd_tree,
                      const Array<Vector2D<T>>& normals,
                      Array<Segment2D<T>>* lines,
                      Array<int>* labels) const {
        static_assert(std::is_floating_point<T>::value, "");

        assert(lines);
        assert(labels);

        // Segment the points.
        Array<std::pair<Point2D<T>, Vector2D<T>>> models;
        Segment(kd_tree, normals, labels, &models);
        RefineLines(kd_tree.points(), *labels, models, lines);
    }

    /**
     * Extract regular planes from a 3D unstructured point set with unoriented
     * normal vectors.
     *
     * Parameters:
     *  kd_tree - the input point set is stored in KD-tree.
     *  normals - the input unoriented normal vectors.
     *  lines   - the extracted line segments.
     *  labels  - labels[i] denotes the line that own the i-th point, or
     *            labels[i] = -1 if no line of the i-th point is found.
     */
    template <typename T>
    void ExtractPlanes(const KDTree<Point3D<T>>& kd_tree,
                       const Array<Vector3D<T>>& normals,
                       Array<Plane3D<T>>* planes,
                       Array<int>* labels) const {
        static_assert(std::is_floating_point<T>::value, "");

        assert(planes);
        assert(labels);

        // Segment the points.
        Array<std::pair<Point3D<T>, Vector3D<T>>> models;
        Segment(kd_tree, normals, labels, &models);
        RefinePlanes(kd_tree.points(), *labels, models, planes);
    }

    void set_min_support_points(int min_support_points) {
        assert(min_support_points > 0);
        min_support_points_ = min_support_points;
    }

    int min_support_points() const {
        return min_support_points_;
    }

    void set_n_constraints(int n_constraints) {
        assert(n_constraints > 0);
        n_constraints = n_constraints_;
    }

    int n_constraints() const {
        return n_constraints_;
    }

private:
    /**
     * Return distance between two normals.
     */
    template <typename Vector>
    static double Distance(const Vector& v1, const Vector& v2) {
        return 2.0 - 2.0 * std::fabs(v1 * v2);
    }

    /**
     * Estimate the minimal value of lambda.
     */
    template <typename Vector>
    static double MinLambda(const Array<Vector>& normals,
                            const Graph& graph) {
        int n = graph.n_vertices();

        Array<double> dis(n, DBL_MAX);
        for (int i = 0; i < n; ++i) {
            for (const Edge* e : graph.edges_from(i)) {
                int j = e->target();
                if (i != j) {
                    dis[i] = std::min(dis[i], Distance(normals[i], normals[j]));
                }
            }
        }
        return std::max(DBL_EPSILON, Median(dis.begin(), dis.end()));
    }

    /**
     * Segment the points by global-L0 algorithm.
     */
    template <typename Point, typename Vector>
    void Segment(const KDTree<Point>& kd_tree,
                 const Array<Vector>& normals,
                 Array<int>* labels,
                 Array<std::pair<Point, Vector>>* models) const {
        assert(kd_tree.size() == normals.size());
        assert(labels);
        assert(models);

        int n = kd_tree.size();
        const Array<Point>& points = kd_tree.points();
        labels->assign(n, -1);
        if (n <= k_neighbors_ || n < min_support_points_) {
            return;
        }

        Graph graph;
        MutualKNearestNeighbors(kd_tree, k_neighbors_, &graph);

        Array<Vector> normals1 = normals;
        Graph::EdgeProperty<int> weights =
            graph.AddEdgeProperty<int>("weights", 1);
        DisjointSet set(n + 1);
        Array<int> sizes(n, 1);
        Array<double> distances(n, 0.0);

        Array<int> available_vertices(n);
        for (int i = 0; i < n; ++i) {
            available_vertices[i] = i;
        }

        bool terminal = false;
        const double min_lambda = MinLambda(normals1, graph);
        for (double lambda = min_lambda; ; lambda *= 1.2) {
            bool finish = true;
            for (int i : available_vertices) {
                if (i != set.Find(i)) continue;

                std::unordered_map<int, Edge*> hash;
                for (Edge* e : graph.edges_from(i)) {
                    hash[e->target()] = e;
                }

                bool is_outlier = true;
                for (Edge* e : graph.edges_from(i)) {
                    int j = e->target();
                    if (set.Find(j) == n) continue;

                    if (sizes[i] < min_support_points_) finish = false;

                    double c = sizes[i] + sizes[j];
                    double dis = Distance(normals1[i], normals1[j]);
                    double loss = (dis * sizes[i] * sizes[j]) / c;
                    double improvement = std::min(lambda, outlier_penalty_) *
                                         weights[e] - loss;

                    if (normals1[i] != normals1[j] &&
                        lambda > outlier_penalty_) {
                        continue;
                    }
                    is_outlier = false;

                    if (improvement > 0.0) {
                        set.Link(j, i);
                        double a = sizes[i] / c;
                        double b = sizes[j] / c;
                        if (normals1[i] != normals1[j]) {
                            Vector v;
                            if (normals1[i] * normals1[j] > 0.0) {
                                v = a * normals1[i] + b * normals1[j];
                            } else {
                                v = a * normals1[i] - b * normals1[j];
                            }
                            v.Normalize();

                            distances[i] += Distance(normals1[i], v) * sizes[i];
                            distances[i] += Distance(normals1[j], v) * sizes[j];
                            normals1[i] = v;
                        }
                        distances[i] += distances[j];
                        sizes[i] += sizes[j];

                        for (Edge* e1 : graph.edges_from(j)) {
                            int k = e1->target();
                            if (k == i) continue;

                            auto iter = hash.find(k);
                            if (iter == hash.end()) {
                                Edge* new_edge = graph.InsertTwoWayEdge(i, k);
                                hash[k] = new_edge;
                                weights[new_edge] = weights[e1];
                                weights[new_edge->twin()] = weights[e1];
                            } else {
                                const Edge* e2 = iter->second;
                                weights[e2] += weights[e1];
                                weights[e2->twin()] += weights[e1];
                            }
                        }
                        graph.EraseEdgesOfVertex(j);
                    }
                }

                if (is_outlier && sizes[i] < min_support_points_) {
                    set.Link(i, n);
                    graph.EraseEdgesOfVertex(i);
                }
            }
            if (terminal) break;

            int n_available_vertices = 0;
            for (int i : available_vertices) {
                if (set.Find(i) == i) {
                    available_vertices[n_available_vertices++] = i;
                }
            }
            available_vertices.resize(n_available_vertices);

            // Find the representatives via greedy sub-modular.
            Array<int> subset;
            for (int i : available_vertices) {
                if (sizes[i] >= min_support_points_) subset.push_back(i);
            }
            if (subset.size() <= n_constraints_) {
                if (finish) terminal = true;
                continue;
            }

            int m = subset.size();
            Array<int> ns(m);
            for (int i = 0; i < subset.size(); ++i) {
                ns[i] = sizes[subset[i]];
            }

            Array<int> seq;
            IndexSort(ns.begin(), ns.end(), &seq);
            std::reverse(seq.begin(), seq.end());
            RMatrix dis(m, m);
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < m; ++j) {
                    if (i == j) {
                        dis(i, j) = 0.0;
                        continue;
                    }

                    int a = subset[seq[i]];
                    int b = subset[seq[j]];
                    dis(i, j) = Distance(normals1[a], normals1[b]) * sizes[a] +
                                distances[a];
                }
            }

            Array<int> representatives;
            double regularity = lambda;
            optimization::GreedySubmodularSubsetSelection(dis,
                                                          regularity,
                                                          &representatives);

            for (int i = 0; i < m; ++i) {
                double t = DBL_MAX;
                int index = i;
                for (int j : representatives) {
                    double d = dis(i, j);
                    if (d < t) {
                        t = d;
                        index = j;
                    }
                }

                normals1[subset[seq[i]]] = normals1[subset[seq[index]]];
            }
            if (finish && representatives.size() <= n_constraints_)
                terminal = true;
        }

        // Set the label array.
        available_vertices.clear();
        for (int i = 0; i < n; ++i) {
            int a = set.Find(i);
            if (a == n) continue;
            if (sizes[a] >= min_support_points_) {
                (*labels)[i] = a;
                if (a == i) available_vertices.push_back(i);
            }
        }

        // Relabeling, so that the value of labels[i] is in the range of
        // [0, n_planes). (Or -1 if the point is not assigned to a plane).
        Array<int> map(n);
        for (int i = 0; i < available_vertices.size(); ++i) {
            int p = available_vertices[i];
            map[p] = i;
            models->emplace_back(points[p], normals1[p]);
        }
        for (int i = 0; i < n; ++i) {
            if ((*labels)[i] == -1) continue;
            (*labels)[i] = map[(*labels)[i]];
        }
    }

    /**
     * Refine lines.
     */
    template <typename T>
    void RefineLines(const Array<Point2D<T>>& points,
                     const Array<int>& labels,
                     const Array<std::pair<Point2D<T>, Vector2D<T>>>& models,
                     Array<Segment2D<T>>* lines) const {
        assert(lines);

        using Line = Line2D<T>;
        using Point = Point2D<T>;

        Array<Array<int>> clusters(models.size());
        for (int i = 0; i < labels.size(); ++i) {
            if (labels[i] == -1) continue;
            clusters[labels[i]].push_back(i);
        }

        for (int i = 0; i < clusters.size(); ++i) {
            if (clusters[i].empty()) continue;

            Line line(models[i].first, models[i].second);
            line = Line(models[i].first, line.normal());

            Array<double> distances(clusters[i].size());
            Array<int> seq;
            PointDirectionCompare2D<T> compare(line.direction());
            Point p1 = points[clusters[i].front()];
            Point p2 = points[clusters[i].front()];
            for (int j = 0; j < clusters[i].size(); ++j) {
                const Point& p = points[clusters[i][j]];
                distances[j] = SignedDistance(p, line);

                if (compare(p, p1)) p1 = p;
                if (compare(p2, p)) p2 = p;
            }
            IndexSort(distances.begin(), distances.end(), &seq);

            int mid = seq[seq.size() / 2];
            line = Line(points[clusters[i][mid]], line.direction());
            lines->emplace_back(Project(p1, line), Project(p2, line));
        }
    }

    /**
     * Refine planes.
     */
    template <typename T>
    void RefinePlanes(const Array<Point3D<T>>& points,
                      const Array<int>& labels,
                      const Array<std::pair<Point3D<T>, Vector3D<T>>>& models,
                      Array<Plane3D<T>>* planes) const {
        using Point = Point3D<T>;
        using Plane = Plane3D<T>;

        Array<Array<int>> clusters(models.size());
        for (int i = 0; i < labels.size(); ++i) {
            if (labels[i] == -1) continue;
            clusters[labels[i]].push_back(i);
        }

        for (int i = 0; i < clusters.size(); ++i) {
            Plane plane(models[i].first, models[i].second);

            Array<double> distances(clusters[i].size());
            Array<int> seq;
            for (int j = 0; j < clusters[i].size(); ++j) {
                const Point& p = points[clusters[i][j]];
                distances[j] = SignedDistance(p, plane);
            }
            IndexSort(distances.begin(), distances.end(), &seq);

            int mid = seq[seq.size() / 2];
            planes->emplace_back(points[clusters[i][mid]], plane.normal());
        }
    }

    int k_neighbors_;
    int min_support_points_;
    int n_constraints_;
    double outlier_penalty_;
};

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_GLOBAL_L0_EXTRACTOR_H_
