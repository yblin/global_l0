//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_LAB_COLOR_H_
#define VISUALIZATION_COLOR_LAB_COLOR_H_

#include <algorithm>
#include <cmath>

#include "codelibrary/visualization/color/xyz_color.h"

namespace cl {

/**
 * LAB Color.
 *
 * LABColor is a color space designed to have perceptual uniformity; i.e. equal
 * changes in its components will be perceived by a human to have equal effects.
 * LABColor is the basis for computing distance between colors.
 *
 * LABColor is device independent and corresponds to the CIE 1976 L*, a*, b*
 * color space with [l, a, b] = [L*, a*, b*] / 100.0.
 *
 * The parameters have the following interpretation:
 *   l: lightness or approximate luminance;
 *   a: color, green (negative) to magenta (positive);
 *   b: color, blue (negative) to yellow (positive).
 *
 * LABColor allows any real number for l, a, and b.
 */
struct LABColor {
public:
    LABColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    LABColor(double l, double a, double b, double alpha = 1.0)
        : l_(l), a_(a), b_(b), alpha_(alpha) {}

    /**
     * Construct LABColor from XYZColor.
     */
    explicit LABColor(const XYZColor& xyz) {
        // White reference.
        static const double xn = 0.412453 + 0.357580 + 0.180423;
        static const double yn = 0.212671 + 0.715160 + 0.072169;
        static const double zn = 0.019334 + 0.119193 + 0.950227;

        double fx = ConvertLAB(xyz.x() / xn);
        double fy = ConvertLAB(xyz.y() / yn);
        double fz = ConvertLAB(xyz.z() / zn);

        l_ = std::max(0.0, 1.16 * fy - 0.16);
        a_ = 5.0 * (fx - fy);
        b_ = 2.0 * (fy - fz);
    }

    /**
     * Convert LABColor to XYZColor.
     */
    XYZColor ToXYZColor() const {
        // White reference.
        static const double xn = 0.412453 + 0.357580 + 0.180423;
        static const double yn = 0.212671 + 0.715160 + 0.072169;
        static const double zn = 0.019334 + 0.119193 + 0.950227;

        double cy = (l_ + 0.16) / 1.16;
        double y = yn * ConvertXYZ(cy);
        double py = std::pow(y / yn, 1.0 / 3.0);
        double cx = a_ / 5.0 + py;
        double x = xn * cx * cx * cx;
        double cz = py - b_ / 2.0;
        double z = zn * cz * cz * cz;

        return {x, y, z, alpha_};
    }

    void set_l(double l)     { l_     = l; }
    void set_a(double a)     { a_     = a; }
    void set_b(double b)     { b_     = b; }
    void set_alpha(double a) { alpha_ = Clamp(a, 0.0, 1.0); }

    double l()     const { return l_;     }
    double a()     const { return a_;     }
    double b()     const { return b_;     }
    double alpha() const { return alpha_; }

private:
    /**
     * Help function to convert the XYZColor to LABColor.
     */
    static double ConvertLAB(double x) {
        return x >= 0.008856 ? std::pow(x, 1.0 / 3.0) :
                               7.787 * x + 16.0 / 116.0;
    }

    /**
     * Help function to convert the LABColor to XYZColor.
     */
    static double ConvertXYZ(double x) {
        return x >= 0.206893 ? (x * x * x) : (x - 16.0 / 116.0) / 7.787;
    }

    double l_     = 0.0;
    double a_     = 0.0;
    double b_     = 0.0;
    double alpha_ = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_LAB_COLOR_H_
