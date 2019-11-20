//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_POISSON_DISK_SAMPLING_H_
#define POINT_CLOUD_POISSON_DISK_SAMPLING_H_
#include <set>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/point_cloud/sampling/sampling.h"
#include "codelibrary/util/tree/octree.h"

namespace cl {
namespace point_cloud {

/**
 * Poisson-Disk sampling for point clouds.It ensures that no two points have
 * a distance smaller than the given resolution after sampling.
 */
class PoissonDiskSampling : public Sampling {
public:
    explicit PoissonDiskSampling(double resolution, int random_seed = 0)
        : resolution_(resolution), Sampling(random_seed) {
        assert(resolution_ > 0.0);
    }

    template <typename Iterator>
    void Sample(Iterator first, Iterator last, Array<int>* samples) {
        assert(samples);

        samples->clear();

        Array<int> seq;
        this->Generate(first, last, &seq);

        RBox3D box(first, last);
        double resolution = resolution_ / std::sqrt(3.0);
        assert(box.x_length() / resolution < INT_MAX);
        assert(box.y_length() / resolution < INT_MAX);
        assert(box.z_length() / resolution < INT_MAX);
        int size1 = static_cast<int>(box.x_length() / resolution) + 1;
        int size2 = static_cast<int>(box.y_length() / resolution) + 1;
        int size3 = static_cast<int>(box.z_length() / resolution) + 1;

        const double resolution2 = resolution_ * resolution_;

        Octree<int> octree(size1, size2, size3);
        using LeafNode = typename Octree<int>::LeafNode;

        // Add the voxels into the octree.
        for (int s : seq) {
            const RPoint3D& p = first[s];
            int x = static_cast<int>((p.x - box.x_min()) / resolution);
            int y = static_cast<int>((p.y - box.y_min()) / resolution);
            int z = static_cast<int>((p.z - box.z_min()) / resolution);
            x = Clamp(x, 0, size1 - 1);
            y = Clamp(y, 0, size2 - 1);
            z = Clamp(z, 0, size3 - 1);

            bool found = false;
            for (int x1 = x - 2; !found && x1 <= x + 2; ++x1) {
                if (x1 < 0 || x1 >= size1) continue;
                for (int y1 = y - 2; !found && y1 <= y + 2; ++y1) {
                    if (y1 < 0 || y1 >= size2) continue;
                    for (int z1 = z - 2; !found && z1 <= z + 2; ++z1) {
                        if (z1 < 0 || z1 >= size3) continue;

                        const LeafNode* node = octree.Find(x1, y1, z1);
                        if (!node) continue;

                        if (SquaredDistance(p, first[node->data()]) <
                            resolution2) {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if (!found) {
                octree.Insert(x, y, z, s);
                samples->push_back(s);
            }
        }
    }

private:
    double resolution_;
};

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_POISSON_DISK_SAMPLING_H_
