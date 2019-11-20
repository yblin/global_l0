//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_COLOR_CMYK_COLOR_H_
#define VISUALIZATION_COLOR_CMYK_COLOR_H_

#include <algorithm>

#include "codelibrary/visualization/color/rgb_color.h"

namespace cl {

/**
 * CMYK color.
 *
 * CMYKColor is a subtractive color model, typically used in printing. CMYK
 * refers to the cyan, magenta, yellow, and black inks used in printing.
 */
class CMYKColor {
public:
    CMYKColor() = default;

    /**
     * Color and opacity levels outside the range 0 to 1 will be clamped.
     */
    CMYKColor(double c, double m, double y, double k, double a = 1.0) {
        cyan_    = Clamp(c, 0.0, 1.0);
        magenta_ = Clamp(m, 0.0, 1.0);
        yellow_  = Clamp(y, 0.0, 1.0);
        key_     = Clamp(k, 0.0, 1.0);
        alpha_   = Clamp(a, 0.0, 1.0);
    }

    /**
     * Construct CMYKColor from RGBColor.
     */
    explicit CMYKColor(const RGBColor& rgb) {
        cyan_    = 1.0 - rgb.red();
        magenta_ = 1.0 - rgb.green();
        yellow_  = 1.0 - rgb.blue();
        key_     = std::min(std::min(cyan_, magenta_), yellow_);
        alpha_   = rgb.alpha();

        if (key_ < 1.0) {
            double k = 1.0 - key_;
            cyan_    = (cyan_    - key_) / k;
            magenta_ = (magenta_ - key_) / k;
            yellow_  = (yellow_  - key_) / k;
        }
    }

    /**
     * Convert to RGB color.
     */
    RGBColor ToRGBColor() const {
        double k = 1.0 - key_;
        double c = cyan_    * k + key_;
        double m = magenta_ * k + key_;
        double y = yellow_  * k + key_;

        return {1.0 - c, 1.0 - m, 1.0 - y, alpha_};
    }

    void set_cyan(double c)    { cyan_    = Clamp(c, 0.0, 1.0); }
    void set_magenta(double m) { magenta_ = Clamp(m, 0.0, 1.0); }
    void set_yellow(double y)  { yellow_  = Clamp(y, 0.0, 1.0); }
    void set_key(double k)     { key_     = Clamp(k, 0.0, 1.0); }
    void set_alpha(double a)   { alpha_   = Clamp(a, 0.0, 1.0); }

    double cyan()    const { return cyan_;    }
    double magenta() const { return magenta_; }
    double yellow()  const { return yellow_;  }
    double key()     const { return key_;     }
    double alpha()   const { return alpha_;   }

private:
    double cyan_    = 0.0;
    double magenta_ = 0.0;
    double yellow_  = 0.0;
    double key_     = 0.0;
    double alpha_   = 0.0;
};

} // namespace cl

#endif // VISUALIZATION_COLOR_CMYK_COLOR_H_
