//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_RGB_COLOR_H_
#define VISUALIZATION_COLOR_RGB_COLOR_H_

#include <algorithm>

#include "codelibrary/visualization/color/rgb32_color.h"

namespace cl {

/**
 * RGB color.
 *
 * ITU-R BT.709 Linear RGB Color space.
 *
 * RGBColor is a standard RGB color space created cooperatively by HP and
 * Microsoft in 1996 for use on monitors, printers and the Internet, and
 * subsequently standardized by the IEC.
 */
struct RGBColor {
public:
    RGBColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    RGBColor(double r, double g, double b, double a = 1.0) {
        red_   = Clamp(r, 0.0, 1.0);
        green_ = Clamp(g, 0.0, 1.0);
        blue_  = Clamp(b, 0.0, 1.0);
        alpha_ = Clamp(a, 0.0, 1.0);
    }

    explicit RGBColor(const RGB32Color& rgb32) {
        red_   = rgb32.red()   / 255.0;
        green_ = rgb32.green() / 255.0;
        blue_  = rgb32.blue()  / 255.0;
        alpha_ = rgb32.alpha() / 255.0;
    }

    /**
     * Convert RGBColor to RGB32Color.
     */
    RGB32Color ToRGB32Color() const {
        return {static_cast<int>(red_   * 255.0),
                static_cast<int>(green_ * 255.0),
                static_cast<int>(blue_  * 255.0),
                static_cast<int>(alpha_ * 255.0)};
    }

    /**
     * Convert the RGB color to the gray scale.
     *
     * Return a value from 0 to 1.
     */
    double ToGrayScale() {
        return (red_ * 11.0 + green_ * 16.0 + blue_ * 5.0) / 32.0;
    }

    void set_red(double r)   { red_   = Clamp(r, 0.0, 1.0); }
    void set_green(double g) { green_ = Clamp(g, 0.0, 1.0); }
    void set_blue(double b)  { blue_  = Clamp(b, 0.0, 1.0); }
    void set_alpha(double a) { alpha_ = Clamp(a, 0.0, 1.0); }

    double red()   const { return red_;   }
    double green() const { return green_; }
    double blue()  const { return blue_;  }
    double alpha() const { return alpha_; }

    double operator[] (int i) const {
        switch(i) {
        case 0: return red_;
        case 1: return green_;
        case 2: return blue_;
        case 3: return alpha_;
        }

        assert(false);
        return 0.0;
    }

private:
    double red_   = 0.0;
    double green_ = 0.0;
    double blue_  = 0.0;
    double alpha_ = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_RGB_COLOR_H_
