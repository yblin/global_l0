//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_TERMINAL_TERMINAL_H_
#define VISUALIZATION_TERMINAL_TERMINAL_H_

#include <cassert>
#include <string>

#include "codelibrary/geometry/kernel/multi_polygon_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/visualization/color.h"
#include "codelibrary/visualization/font.h"
#include "codelibrary/visualization/pen.h"

namespace cl {

/**
 * The terminal class to perform low-level painting.
 */
class Terminal {
public:
    Terminal(int width = 640, int height = 480)
        : width_(width), height_(height), plot_area_(0.0, width, 0.0, height) {
        assert(height_ >= 0 && width_ >= 0);
    }

    Terminal(const Terminal&) = delete;

    Terminal& operator=(const Terminal&) = delete;

    virtual ~Terminal() = default;

    /**
     * Clear the painter.
     */
    virtual void clear() = 0;

    /**
     * Draw circle.
     */
    virtual void DrawCircle(double x, double y, double r) = 0;

    /**
     * Draw line (x1, y1)->(x2, y2).
     */
    virtual void DrawLine(double x1, double y1, double x2, double y2) = 0;

    /**
     * Draw rectangle.
     */
    virtual void DrawRectangle(double x, double y, double w, double h) = 0;

    /**
     * Draw polyline.
     */
    virtual void DrawPolyline(const Array<RPoint2D>& polyline) = 0;

    /**
     * Draw polygon.
     */
    virtual void DrawPolygon(const RMultiPolygon2D& polygon) = 0;

    /**
     * Draw triangle.
     */
    virtual void DrawTriangle(const RPoint2D& p1, const RPoint2D& p2,
                              const RPoint2D& p3) = 0;

    /**
     * Draw text at position(x, y).
     */
    virtual void DrawText(double x, double y, const std::string& text) = 0;

    /**
     * Draw vertical text at position(x, y).
     */
    virtual void DrawVerticalText(double x, double y,
                                  const std::string& text) = 0;

    /**
     * Save to the file.
     */
    virtual void SaveToFile(const std::string& file) const = 0;

    /**
     * Resize the terminal.
     */
    virtual void Resize(int height, int width) = 0;

    /**
     * Return the height of terminal.
     */
    int height() const {
        return height_;
    }

    /**
     * Return the width of terminal.
     */
    int width() const {
        return width_;
    }

    /**
     * Return the current pen.
     */
    const Pen& pen() const {
        return pen_;
    }

    /**
     * Set the current pen.
     */
    void set_pen(const Pen& pen) {
        pen_ = pen;
    }

    /**
     * Return the current font.
     */
    const Font& font() const {
        return font_;
    }

    /**
     * Set font.
     */
    void set_font(const Font& font) {
        font_ = font;
    }

    /**
     * Set plot area.
     */
    void set_plot_area(const RBox2D& plot_area) {
        assert(plot_area.x_min() >= 0.0);
        assert(plot_area.y_min() >= 0.0);
        assert(plot_area.x_max() <= width_);
        assert(plot_area.y_max() <= height_);

        plot_area_ = plot_area;
    }

    /**
     * Return the plot area.
     */
    const RBox2D& plot_area() const {
        return plot_area_;
    }

protected:
    int width_;        // The width of terminal.
    int height_;       // The height of terminal.
    Font font_;        // The font used for drawing text.
    Pen pen_;          // The current pen.
    RBox2D plot_area_; // The area to plot.
};

} // namespace cl

#endif // VISUALIZATION_TERMINAL_TERMINAL_H_
