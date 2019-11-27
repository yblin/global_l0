TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    codelibrary/base/algorithm.h \
    codelibrary/base/array.h \
    codelibrary/base/array_nd.h \
    codelibrary/base/bits.h \
    codelibrary/base/equal.h \
    codelibrary/base/format.h \
    codelibrary/base/log.h \
    codelibrary/base/message.h \
    codelibrary/base/object_pool.h \
    codelibrary/base/time_counter.h \
    codelibrary/base/units.h \
    codelibrary/geometry/kernel/angle.h \
    codelibrary/geometry/kernel/box_2d.h \
    codelibrary/geometry/kernel/box_3d.h \
    codelibrary/geometry/kernel/center_2d.h \
    codelibrary/geometry/kernel/center_3d.h \
    codelibrary/geometry/kernel/circle_2d.h \
    codelibrary/geometry/kernel/distance_2d.h \
    codelibrary/geometry/kernel/distance_3d.h \
    codelibrary/geometry/kernel/generator_2d.h \
    codelibrary/geometry/kernel/generator_3d.h \
    codelibrary/geometry/kernel/halfedge_list.h \
    codelibrary/geometry/kernel/intersect_2d.h \
    codelibrary/geometry/kernel/intersect_3d.h \
    codelibrary/geometry/kernel/line_2d.h \
    codelibrary/geometry/kernel/line_3d.h \
    codelibrary/geometry/kernel/multi_polygon_2d.h \
    codelibrary/geometry/kernel/plane_3d.h \
    codelibrary/geometry/kernel/point_2d.h \
    codelibrary/geometry/kernel/point_3d.h \
    codelibrary/geometry/kernel/polygon_2d.h \
    codelibrary/geometry/kernel/predicate_2d.h \
    codelibrary/geometry/kernel/predicate_3d.h \
    codelibrary/geometry/kernel/quaternion.h \
    codelibrary/geometry/kernel/ray_3d.h \
    codelibrary/geometry/kernel/segment_2d.h \
    codelibrary/geometry/kernel/segment_3d.h \
    codelibrary/geometry/kernel/sphere_3d.h \
    codelibrary/geometry/kernel/transform_2d.h \
    codelibrary/geometry/kernel/transform_3d.h \
    codelibrary/geometry/kernel/triangle_2d.h \
    codelibrary/geometry/kernel/triangle_3d.h \
    codelibrary/geometry/kernel/vector_2d.h \
    codelibrary/geometry/kernel/vector_3d.h \
    codelibrary/geometry/pca_estimate_normals_2d.h \
    codelibrary/geometry/pca_estimate_normals_3d.h \
    codelibrary/graph/graph.h \
    codelibrary/math/matrix/eigen/symmetric_eigen.h \
    codelibrary/math/matrix/matrix.h \
    codelibrary/math/number/bigint.h \
    codelibrary/math/number/determinant.h \
    codelibrary/math/number/exact_float.h \
    codelibrary/math/number/interval_float.h \
    codelibrary/math/basic_linear_algebra.h \
    codelibrary/math/vector.h \
    codelibrary/optimization/discrete/subset_selection/greedy_submodular.h \
    codelibrary/point_cloud/io/xyz_io.h \
    codelibrary/point_cloud/sampling/grid_sampling.h \
    codelibrary/point_cloud/sampling/poisson_disk_sample.h \
    codelibrary/point_cloud/sampling/sampling.h \
    codelibrary/point_cloud/bilateral_filter.h \
    codelibrary/point_cloud/estimate_local_planes.h \
    codelibrary/point_cloud/global_l0_extractor.h \
    codelibrary/point_cloud/iterative_closest_point.h \
    codelibrary/point_cloud/k_nearest_neighbors.h \
    codelibrary/point_cloud/mst_orient_normals.h \
    codelibrary/point_cloud/point_octree.h \
    codelibrary/point_cloud/region_growing.h \
    codelibrary/point_cloud/structured_points.h \
    codelibrary/point_cloud/supervoxel_segmentation.h \
    codelibrary/statistics/kernel/covariance.h \
    codelibrary/statistics/kernel/deviation.h \
    codelibrary/statistics/kernel/mean.h \
    codelibrary/statistics/kernel/median.h \
    codelibrary/statistics/kernel/ranking.h \
    codelibrary/statistics/regression/iteratively_reweighted_least_sqaure_fitting.h \
    codelibrary/statistics/regression/least_median_squares_fitting.h \
    codelibrary/statistics/regression/linear_least_squares_fitting.h \
    codelibrary/statistics/regression/linear_regression.h \
    codelibrary/statistics/regression/logistic_regression.h \
    codelibrary/statistics/regression/principal_component_analysis_2d.h \
    codelibrary/statistics/regression/principal_component_analysis_3d.h \
    codelibrary/statistics/regression/ridge_regression.h \
    codelibrary/statistics/regression/softmax_regression.h \
    codelibrary/string/string_split.h \
    codelibrary/util/io/file_map.h \
    codelibrary/util/io/line_reader.h \
    codelibrary/util/list/indexed_list.h \
    codelibrary/util/metric/angular.h \
    codelibrary/util/metric/cosine.h \
    codelibrary/util/metric/euclidean.h \
    codelibrary/util/metric/jaccard.h \
    codelibrary/util/metric/manhattan.h \
    codelibrary/util/metric/squared_euclidean.h \
    codelibrary/util/metric/tanimoto.h \
    codelibrary/util/set/disjoint_set.h \
    codelibrary/util/tree/kd_tree.h \
    codelibrary/visualization/color/cmyk_color.h \
    codelibrary/visualization/color/color_map.h \
    codelibrary/visualization/color/hsi_color.h \
    codelibrary/visualization/color/hsl_color.h \
    codelibrary/visualization/color/hsv_color.h \
    codelibrary/visualization/color/lab_color.h \
    codelibrary/visualization/color/rgb32_color.h \
    codelibrary/visualization/color/rgb_color.h \
    codelibrary/visualization/color/xyz_color.h \
    codelibrary/visualization/plot/axis.h \
    codelibrary/visualization/plot/legend.h \
    codelibrary/visualization/plot/line_plot.h \
    codelibrary/visualization/plot/mix_plot.h \
    codelibrary/visualization/plot/plot.h \
    codelibrary/visualization/plot/point_plot.h \
    codelibrary/visualization/plot/polygon_plot.h \
    codelibrary/visualization/plot/sub_plot.h \
    codelibrary/visualization/terminal/svg_terminal.h \
    codelibrary/visualization/terminal/terminal.h \
    codelibrary/visualization/color.h \
    codelibrary/visualization/font.h \
    codelibrary/visualization/pen.h \
    codelibrary/visualization/plot.h
