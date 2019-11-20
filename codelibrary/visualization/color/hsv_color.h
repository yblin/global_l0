//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_HSV_COLOR_H_
#define VISUALIZATION_COLOR_HSV_COLOR_H_

#include <algorithm>
#include <cassert>

#include "codelibrary/visualization/color/rgb_color.h"

namespace cl {

/**
 * HSV color.
 *
 * HSVColor (hue, saturation, and value), also known as HSB (hue, saturation,
 * and brightness), is similar to the HSLColor.
 */
class HSVColor {
public:
    HSVColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    HSVColor(double h, double s, double v, double a = 1.0) {
        hue_        = Clamp(h, 0.0, 1.0);
        saturation_ = Clamp(s, 0.0, 1.0);
        value_      = Clamp(v, 0.0, 1.0);
        alpha_      = Clamp(a, 0.0, 1.0);
    }

    /**
     * Construct HSVColor from RGBColor.
     */
    explicit HSVColor(const RGBColor& rgb) {
        alpha_ = rgb.alpha();

        double r = rgb.red();
        double g = rgb.green();
        double b = rgb.blue();

        double max = std::max(std::max(r, g), b);
        double min = std::min(std::min(r, g), b);
        double delta = max - min;

        value_ = max;

        if (Equal(min, max)) {
            hue_ = 0.0;
            saturation_ = 0.0;
        } else {
            if (Equal(r, max)) {
                hue_ = (g - b) / delta;
            } else if (Equal(g, max)) {
                hue_ = (b - r) / delta + 2.0;
            } else {
                hue_ = (r - g) / delta + 4.0;
            }
            if (hue_ < 0.0) hue_ += 6.0;
            hue_ /= 6.0;

            saturation_ = delta / max;
        }
    }

    /**
     * Convert to RGB color.
     */
    RGBColor ToRGBColor() const {
        if (saturation_ == 0.0) {
            return {value_, value_, value_, alpha_};
        }

        double h = hue_ == 1.0 ? 0.0 : 6.0 * hue_;
        int sextant = static_cast<int>(h);
        double fract = h - sextant;
        double p = value_ * (1.0 - saturation_);
        double q = value_ * (1.0 - (saturation_ * fract));
        double t = value_ * (1.0 - (saturation_ * (1.0 - fract)));

        switch (sextant) {
        case 0:
            return {value_, t, p, alpha_};
        case 1:
            return {q, value_, p, alpha_};
        case 2:
            return {p, value_, t, alpha_};
        case 3:
            return {p, q, value_, alpha_};
        case 4:
            return {t, p, value_, alpha_};
        case 5:
            return {value_, p, q, alpha_};
        default:
            assert(false && "Unreachable code.");
        }

        return {0, 0, 0, alpha_};
    }

    void set_hue(double h)        { hue_        = Clamp(h, 0.0, 1.0); }
    void set_saturation(double s) { saturation_ = Clamp(s, 0.0, 1.0); }
    void set_value(double v)      { value_      = Clamp(v, 0.0, 1.0); }
    void set_alpha(double a)      { alpha_      = Clamp(a, 0.0, 1.0); }

    double hue()        const { return hue_;        }
    double saturation() const { return saturation_; }
    double value()      const { return value_;      }
    double alpha()      const { return alpha_;      }

private:
    double hue_        = 0.0;
    double saturation_ = 0.0;
    double value_      = 0.0;
    double alpha_      = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_HSV_COLOR_H_
