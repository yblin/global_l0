//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PEN_H_
#define VISUALIZATION_PEN_H_

#include <cassert>

#include "codelibrary/visualization/color.h"

namespace cl {

/**
 * The Pen class defines how a Painter should draw lines and outlines of shapes.
 */
struct Pen {
    // Point Style of Pen.
    enum PointStyle {
        CIRCLE    = 0,
        RECTANGLE = 1
    };

    // Line Style of Pen.
    enum LineStyle {
        SOLID_LINE    = 0,  // -----
        DASH_LINE     = 1,  // - - -
        DOT_LINE      = 2,  // . . .
        DASH_DOT_LINE = 3   // -.-.-
    };

    // Fill Style of Pen.
    enum FillStyle {
        NOT_FILL = 0, // Boundary only.
        FILL = 1      // Fill the close shape.
    };

    RGB32Color color       = RGB32Color(0, 0, 0, 255);
    RGB32Color fill_color  = RGB32Color(0, 0, 0, 255);
    double line_width      = 1.0;
    double point_radius    = 3.0;
    PointStyle point_style = CIRCLE;
    LineStyle line_style   = SOLID_LINE;
    FillStyle fill_style   = FILL;

    Pen(const RGB32Color& c = RGB32Color(0, 0, 0, 255))
        : color(c),
          fill_color(color) {}
};

} // namespace cl

#endif // VISUALIZATION_PEN_H_
