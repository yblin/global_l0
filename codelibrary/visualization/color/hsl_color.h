//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_HSL_COLOR_H_
#define VISUALIZATION_COLOR_HSL_COLOR_H_

#include <algorithm>
#include <cassert>

#include "codelibrary/base/equal.h"
#include "codelibrary/visualization/color/rgb_color.h"

namespace cl {

/**
 * HSL color.
 *
 * HSLColor (hue, saturation, and lightness) corresponds a cylindrical
 * transformation of RGBColor, typically used for color picking, allowing for
 * easier interpretation of color parameters.
 */
class HSLColor {
public:
    HSLColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    HSLColor(double h, double s, double l, double a) {
        hue_        = Clamp(h, 0.0, 1.0);
        saturation_ = Clamp(s, 0.0, 1.0);
        lightness_  = Clamp(l, 0.0, 1.0);
        alpha_      = Clamp(a, 0.0, 1.0);
    }

    /**
     * Construct HSLColor from RGBColor.
     */
    explicit HSLColor(const RGBColor& rgb) {
        alpha_ = rgb.alpha();

        double r = rgb.red();
        double g = rgb.green();
        double b = rgb.blue();

        double max = std::max(std::max(r, g), b);
        double min = std::min(std::min(r, g), b);
        double delta1 = max - min;
        double delta2 = max + min;
        lightness_ = delta2 / 2.0;

        if (min < max) {
            if (Equal(r, max)) {
                hue_ = (g - b) / delta1;
            } else if (Equal(g, max)) {
                hue_ = (b - r) / delta1 + 2.0;
            } else {
                hue_ = (r - g) / delta1 + 4.0;
            }
            if (hue_ < 0.0) hue_ += 6.0;
            hue_ /= 6.0;

            saturation_ = (lightness_ <= 0.5) ? (delta1 / delta2) :
                                                (delta1 / (2.0 - delta2));
        }
    }

    /**
     * Convert to RGB color.
     */
    RGBColor ToRGBColor() const {
        if (saturation_ == 0.0) {
            return {lightness_, lightness_, lightness_, alpha_};
        }

        if (lightness_ == 0.0) {
            return {0, 0, 0, alpha_};
        }

        double v = (lightness_ <= 0.5) ? lightness_ * (1.0 + saturation_)
                                       : lightness_ + saturation_ -
                                         lightness_ * saturation_;

        double min = 2.0 * lightness_ - v;
        double sv = (v - min) / v;

        double h = hue_ == 1.0 ? 0.0 : 6.0 * hue_;
        int sextant = static_cast<int>(h);
        double fract = h - sextant;
        double vsf = v * sv * fract;
        double mid1 = min + vsf;
        double mid2 = v - vsf;

        switch (sextant) {
        case 0:
            return {v, mid1, min, alpha_};
        case 1:
            return {mid2, v, min, alpha_};
        case 2:
            return {min, v, mid1, alpha_};
        case 3:
            return {min, mid2, v, alpha_};
        case 4:
            return {mid1, min, v, alpha_};
        case 5:
            return {v, min, mid2, alpha_};
        default:
            assert(false && "Unreachable code.");
        }

        return {0, 0, 0, alpha_};
    }

    void set_hue(double h)        { hue_        = Clamp(h, 0.0, 1.0); }
    void set_saturation(double s) { saturation_ = Clamp(s, 0.0, 1.0); }
    void set_lightness(double l)  { lightness_  = Clamp(l, 0.0, 1.0); }
    void set_alpha(double a)      { alpha_      = Clamp(a, 0.0, 1.0); }

    double hue()        const { return hue_;        }
    double saturation() const { return saturation_; }
    double lightness()  const { return lightness_;  }
    double alpha()      const { return alpha_;      }

private:
    double hue_        = 0.0;
    double saturation_ = 0.0;
    double lightness_  = 0.0;
    double alpha_      = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_HSL_COLOR_H_
