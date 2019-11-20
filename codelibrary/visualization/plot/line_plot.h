//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_LINE_PLOT_H_
#define VISUALIZATION_PLOT_LINE_PLOT_H_

#include <cassert>
#include <string>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/polygon_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"
#include "codelibrary/visualization/plot/plot.h"
#include "codelibrary/visualization/terminal/terminal.h"

namespace cl {
namespace plot {

/**
 * 2D line plot.
 *
 * Sample usage:
 *
 *   // Input data.
 *   Array<RPoint2D> data1(100), data2(100);
 *   for (int i = 0; i < 100; ++i) {
 *       double x = 2.0 * M_PI * i / 100.0;
 *       data1[i].x = data2[i].x = x;
 *       data1[i].y = std::sin(x);
 *       data2[i].y = std::cos(x);
 *   }
 *
 *   // Plotting.
 *   plot::LinePlot line_plot;
 *   line_plot.Draw(data1);
 *   line_plot.Draw(data2);
 *
 *   // Save the plot to the file.
 *   line_plot.Save("test.svg");
 */
class LinePlot : public Plot {
    // Object for LinePlot.
    class Object {
        friend class LinePlot;

    public:
        Object() = default;

        template <typename T>
        Object(const Array<Point2D<T> >& points, const std::string& name) {
            SetData(points, name);
        }

        template <typename T>
        Object(const Polygon2D<T>& polygon, const std::string& name) {
            if (!polygon.empty()) {
                Array<Point2D<T> > points = polygon.vertices();
                points.push_back(points.front());
                SetData(points, name);
            }
        }

        template <typename T>
        Object(const Box2D<T>& box, const std::string& name) {
            Array<Point2D<T> > points(5);
            points[0] = Point2D<T>(box.x_min(), box.y_min());
            points[1] = Point2D<T>(box.x_min(), box.y_max());
            points[2] = Point2D<T>(box.x_max(), box.y_max());
            points[3] = Point2D<T>(box.x_max(), box.y_min());
            points[4] = Point2D<T>(box.x_min(), box.y_min());
            SetData(points, name);
        }

        template <typename T>
        Object(const Segment2D<T>& line, const std::string& name) {
            Array<Point2D<T> > points(2);
            points[0] = line.lower_point();
            points[1] = line.upper_point();
            SetData(points, name);
        }

        bool empty() const {
            return data_.empty();
        }

    private:
        template <typename T>
        void SetData(const Array<Point2D<T> >& data, const std::string& name) {
            data_.reserve(data.size());
            for (const Point2D<T>& p : data) {
                data_.emplace_back(static_cast<double>(p.x),
                                   static_cast<double>(p.y));
            }
            name_ = name;
        }
        void SetData(const Array<RPoint2D>& data, const std::string& name) {
            data_ = data;
            name_ = name;
        }

        Array<RPoint2D> data_;
        std::string name_;
    };

public:
    explicit LinePlot(const ColorMap& color_map = ColorMap::Lines())
        : color_map_(color_map) {}

    /**
     * Return true if this plot is empty.
     */
    virtual bool empty() const override {
        return objects_.empty();
    }

    /**
     * Clear this plot.
     */
    virtual void clear() override {
        objects_.clear();
        pens_.clear();
    }

    /**
     * Set the default line width.
     */
    void set_line_width(double line_width) {
        assert(line_width > 0.0);

        line_width_ = line_width;
    }

    double line_width() const {
        return line_width_;
    }

    /**
     * Draw lines with the given stroke color and line style.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object,
              const RGB32Color& color,
              const Pen::LineStyle& line_style = Pen::SOLID_LINE,
              const std::string& name = "") {
        Object o(object, name);
        Pen pen(color);
        pen.line_width = line_width_;
        pen.line_style = line_style;
        DrawObject(o, pen);
    }

    /**
     * Draw lines with the given line style.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object,
              const Pen::LineStyle& line_style,
              const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], line_style, name);
    }

    /**
     * Draw lines with the given stroke color.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object,
              const RGB32Color& color,
              const std::string& name) {
        Draw(object, color, Pen::SOLID_LINE, name);
    }

    /**
     * Draw lines with the given pen.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const Pen& pen,
              const std::string& name = "") {
        assert(pen.line_width > 0.0);

        Object o(object, name);
        DrawObject(o, pen);
    }

    /**
     * Draw lines with the default pen.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], name);
    }

    /**
     * Draw plottable data on the terminal.
     */
    virtual void DrawData(Terminal* terminal) override {
        // Draw ploylines.
        for (int i = 0; i < objects_.size(); ++i) {
            terminal->set_pen(pens_[i]);

            // Convert the coordinate for plotting.
            Array<RPoint2D> polyline;
            for (const RPoint2D& p : objects_[i].data_) {
                const RPoint2D q = ToPlotPosition(p);
                polyline.emplace_back(q.x, q.y);
            }
            terminal->DrawPolyline(polyline);
        }
    }

protected:
    /**
     * Draw polyline determined by the given points with specific pen.
     */
    void DrawObject(const Object& object, const Pen& pen) {
        if (object.empty()) return;

        // Update the data range.
        RBox2D box(object.data_.begin(), object.data_.end());
        if (objects_.empty())
            data_range_ = box;
        else
            data_range_.Join(box);

        objects_.push_back(object);
        pens_.push_back(pen);

        if (!object.name_.empty())
            legend_.InsertItem(Legend::LINE, object.name_, pen);
    }

    // Color map to draw the lines.
    const ColorMap& color_map_;

    // The input data for plotting.
    Array<Object> objects_;

    // The pens for each line.
    Array<Pen> pens_;

    // Line width.
    double line_width_ = 2.0;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_LINE_PLOT_H_
