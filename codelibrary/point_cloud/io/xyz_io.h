//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_IO_XYZ_IO_H_
#define POINT_CLOUD_IO_XYZ_IO_H_

#include <sstream>
#include <fstream>

#include "codelibrary/base/array.h"
#include "codelibrary/base/format.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/util/io/line_reader.h"
#include "codelibrary/visualization/color/rgb32_color.h"

namespace cl {
namespace point_cloud {

/**
 * Read points from XYZ file.
 * Return false if an error occurs.
 */
template <typename T>
bool ReadXYZPoints(const char* filename, Array<Point3D<T>>* points) {
    assert(points);

    points->clear();

    io::LineReader line_reader;
    if (!line_reader.Open(filename)) return false;
    
    T x, y, z;
    char* line = nullptr;
    while ((line = line_reader.ReadLine())) {
        std::stringstream ss(line);
        if (!(ss >> x) || !(ss >> y) || !(ss >> z)) {
            LOG(INFO) << "Invalid XYZ format at line: " << line_reader.n_line();
            return false;
        }
        points->emplace_back(x, y, z);
    }

    return true;
}

/**
 * Read color points from XYZ file.
 */
template <typename T>
bool ReadXYZPoints(const char* filename, Array<Point3D<T>>* points,
                   Array<RGB32Color>* colors) {
    assert(points);
    assert(colors);

    points->clear();
    colors->clear();

    char* line;
    io::LineReader line_reader(filename);
    T x, y, z;
    int r, g, b;
    while (line = line_reader.ReadLine()) {
        std::stringstream ss(line);
        if (!(ss >> x) || !(ss >> y) || !(ss >> z) ||
            !(ss >> r) || !(ss >> g) || !(ss >> b)) {
            LOG(INFO) << "Invalid XYZ format at line: " << line_reader.n_line();
            return false;
        }
        points->emplace_back(x, y, z);
        colors->emplace_back(r, g, b);
    }

    return true;
}

/**
 * Read oriented points from XYZ file.
 */
template <typename T>
bool ReadXYZPoints(const char* filename, Array<Point3D<T>>* points,
                   Array<Vector3D<T>>* normals) {
    assert(points);
    assert(normals);

    points->clear();
    normals->clear();

    char* line;
    io::LineReader line_reader(filename);
    T x, y, z, nx, ny, nz;
    while (line = line_reader.ReadLine()) {
        std::stringstream ss(line);
        if (!(ss >> x) || !(ss >> y) || !(ss >> z) ||
            !(ss >> nx) || !(ss >> ny) || !(ss >> nz)) {
            LOG(INFO) << "Invalid XYZ format at line: " << line_reader.n_line();
            return false;
        }
        points->emplace_back(x, y, z);
        normals->emplace_back(nx, ny, nz);
    }

    return true;
}

/**
 * Write points into XYZ file.
 */
template <typename T>
bool WriteXYZPoints(const char* filename, const Array<Point3D<T>>& points) {
    FILE* file = std::fopen(filename, "wb");
    if (!file) {
        LOG(INFO) << "Cannot open XYZ file '" << filename << "' for writing.";
        return false;
    }

    for (const Point3D<T>& p : points) {
        fmt::print(file, "{} {} {}\n", p.x, p.y, p.z);
    }

    std::fclose(file);
    return true;
}

/**
 * Write color points into XYZ file.
 */
template <typename T>
bool WriteXYZPoints(const char* filename, const Array<Point3D<T>>& points,
                    const Array<RGB32Color>& colors) {
    assert(points.size() == colors.size());

    FILE* file = std::fopen(filename, "wb");
    if (!file) {
        LOG(INFO) << "Cannot open XYZ file '" << filename << "' for writing.";
        return false;
    }

    for (int i = 0; i < points.size(); ++i) {
        const Point3D<T>& p = points[i];
        const RGB32Color& c = colors[i];
        fmt::print(file, "{} {} {} {} {} {}\n", p.x, p.y, p.z, 
                   c.red(), c.green(), c.blue());
    }

    std::fclose(file);
    return true;
}

/**
 * Write oriented points into XYZ file.
 */
template <typename T>
bool WriteXYZPoints(const char* filename,
                    const Array<Point3D<T>>& points,
                    const Array<Vector3D<T>>& normals) {
    assert(points.size() == normals.size());

    FILE* file = std::fopen(filename, "wb");
    if (!file) {
        LOG(INFO) << "Cannot open XYZ file '" << filename << "' for writing.";
        return false;
    }

    for (int i = 0; i < points.size(); ++i) {
        const Point3D<T>& p = points[i];
        const Vector3D<T>& n = normals[i];
        fmt::print(file, "{} {} {} {} {} {}\n", p.x, p.y, p.z, n.x, n.y, n.z);
    }

    std::fclose(file);
    return true;
}

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_IO_XYZ_IO_H_
