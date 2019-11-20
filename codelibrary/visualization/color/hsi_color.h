//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_HSI_COLOR_H_
#define VISUALIZATION_COLOR_HSI_COLOR_H_

#include <algorithm>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else 
#include <cmath>
#endif // _USE_MATH_DEFINES

#include "codelibrary/visualization/color/rgb_color.h"

namespace cl {

/**
 * HSI color.
 *
 * HSIColor (hue, saturation, and intensity), common in computer vision
 * applications, attempts to balance the advantages and disadvantages of the
 * HSVColor and HSLColor.
 */
class HSIColor {
public:
    HSIColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    HSIColor(double h, double s, double i, double a) {
        hue_        = Clamp(h, 0.0, 1.0);
        saturation_ = Clamp(s, 0.0, 1.0);
        intensity_  = Clamp(i, 0.0, 1.0);
        alpha_      = Clamp(a, 0.0, 1.0);
    }

    /**
     * Construct HSLColor from RGBColor.
     */
    explicit HSIColor(const RGBColor& rgb) {
        double r = rgb.red();
        double g = rgb.green();
        double b = rgb.blue();

        double min = std::min(std::min(r, g), b);
        double t = std::sqrt((r - g) * (r - g) + (r - b) * (g - b));
        double theta = std::acos(0.5 * ((r - g) + (r - b)) / t) * 180.0 / M_PI;
        double sum = r + g + b;

        if (theta > 0.0) hue_ = (b <= g) ? theta / 360.0 : 1.0 - theta / 360.0;
        if (sum > 0.0) saturation_ = 1.0 - 3.0 / sum * min;
        intensity_ = sum / 3.0;
        alpha_ = rgb.alpha();
    }

    /**
     * Convert to RGB color.
     */
    RGBColor ToRGBColor() const {
        double t = intensity_ * (1.0 - saturation_);
        double h = hue_ * 360.0;
        double r = 0.0, g = 0.0, b = 0.0;

        const double c = M_PI / 180.0;

        if (h < 120.0) {
            r = intensity_ * (1.0 + saturation_ * std::cos(h * c) /
                              std::cos((60.0 - h) * c));
            b = t;
            g = 3.0 * intensity_ - (r + b);
        } else if (h < 240.0) {
            h -= 120.0;
            r = t;
            g = intensity_ * (1.0 + saturation_ * std::cos(h * c) /
                              std::cos((60.0 - h) * c));
            b = 3.0 * intensity_ - (r + g);
        } else {
            h -= 240.0;
            g = t;
            b = intensity_ * (1.0 + saturation_ * std::cos(h * c) /
                              std::cos((60.0 - h) * c));
            r = 3.0 * intensity_ - (g + b);
        }

        return {r, g, b, alpha_};
    }

    void set_hue(double h)        { hue_        = Clamp(h, 0.0, 1.0); }
    void set_saturation(double s) { saturation_ = Clamp(s, 0.0, 1.0); }
    void set_intensity(double i)  { intensity_  = Clamp(i, 0.0, 1.0); }
    void set_alpha(double a)      { alpha_      = Clamp(a, 0.0, 1.0); }

    double hue()        const { return hue_;        }
    double saturation() const { return saturation_; }
    double intensity()  const { return intensity_;  }
    double alpha()      const { return alpha_;      }

private:
    double hue_        = 0.0;
    double saturation_ = 0.0;
    double intensity_  = 0.0;
    double alpha_      = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_HSI_COLOR_H_
