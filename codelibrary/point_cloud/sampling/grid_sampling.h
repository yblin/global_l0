//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_SAMPLING_GRID_SAMPLE_H_
#define POINT_CLOUD_SAMPLING_GRID_SAMPLE_H_

#include <cassert>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/point_cloud/sampling/sampling.h"
#include "codelibrary/util/tree/octree.h"

namespace cl {
namespace point_cloud {

/**
 * GridSampling considers a regular grid covering the bounding box of the input 
 * point cloud and clusters all points sharing the same cell of the grid by 
 * picking as represent one arbitrarily chosen point.
 */
class GridSampling : public Sampling {
public:
    explicit GridSampling(double resolution, int random_seed = 0)
        : Sampling(random_seed), resolution_(resolution) {
        assert(resolution_ > 0.0);
    }

    template <typename Iterator>
    void Sample(Iterator first, Iterator last, Array<int>* samples) {
        assert(samples);
    
        samples->clear();

        Array<int> seq;
        this->Generate(first, last, &seq);

        RBox3D box(first, last);
        assert(box.x_length() / resolution_ < INT_MAX);
        assert(box.y_length() / resolution_ < INT_MAX);
        assert(box.z_length() / resolution_ < INT_MAX);
        int size1 = static_cast<int>(box.x_length() / resolution_) + 1;
        int size2 = static_cast<int>(box.y_length() / resolution_) + 1;
        int size3 = static_cast<int>(box.z_length() / resolution_) + 1;

        Octree<bool> octree(size1, size2, size3);
        using LeafNode = typename Octree<bool>::LeafNode;
        using Point = typename std::iterator_traits<Iterator>::value_type;

        // Add the voxels into the octree.
        for (int s : seq) {
            const Point& p = first[s];
            int x = static_cast<int>((p.x - box.x_min()) / resolution_);
            int y = static_cast<int>((p.y - box.y_min()) / resolution_);
            int z = static_cast<int>((p.z - box.z_min()) / resolution_);
            x = Clamp(x, 0, size1 - 1);
            y = Clamp(y, 0, size2 - 1);
            z = Clamp(z, 0, size3 - 1);

            std::pair<LeafNode*, bool> pair = octree.Insert(x, y, z, true);
            if (pair.second) {
                samples->push_back(s);
            }
        }
    }

private:
    double resolution_; // The resolution of the girds.
};

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_SAMPLING_GRID_SAMPLE_H_
