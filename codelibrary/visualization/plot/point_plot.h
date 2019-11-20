//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_POINT_PLOT_H_
#define VISUALIZATION_PLOT_POINT_PLOT_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/box_2d.h"
#include "codelibrary/geometry/kernel/circle_2d.h"
#include "codelibrary/visualization/plot/plot.h"
#include "codelibrary/visualization/terminal/terminal.h"

namespace cl {
namespace plot {

/**
 * 2D point plot.
 *
 * Sample usage:
 *
 *   // Generate data.
 *   std::mt19937 random;
 *   RBox2D box(0.0, 1.0, 0.0, 1.0);
 *   geometry::RandomPointInBox2D random_point(box);
 *   Array<RPoint2D> points;
 *   for (int i = 0; i < 1000; ++i) {
 *       points.push_back(random_point(&random));
 *   }
 *
 *   // Draw and show the scatter.
 *   plot::PointPlot point_plot;
 *   point_plot.Draw(points);
 *
 *   // Save the plot to the file.
 *   point_plot.Save("test.svg");
 */
class PointPlot : public Plot {    
    // Object for PointPlot.
    class Object {
        friend class PointPlot;
    public:
        Object() = default;

        template <typename T>
        explicit Object(const Array<Point2D<T> >& points,
                        const std::string& name)
            : name_(name) {
            SetData(points);
        }

        template <typename T>
        explicit Object(const Point2D<T>& p, const std::string& name)
            : name_(name) {
            data_.emplace_back(p.x, p.y);
        }

        bool empty() const {
            return data_.empty();
        }

    private:
        template <typename T>
        void SetData(const Array<Point2D<T> >& data) {
            data_.reserve(data.size());
            for (const Point2D<T>& p : data) {
                data_.emplace_back(p.x, p.y);
            }
        }
        void SetData(const Array<RPoint2D>& data) {
            data_ = data;
        }

        Array<RPoint2D> data_;
        std::string name_;
    };

public:
    explicit PointPlot(const ColorMap& color_map = ColorMap::Set2())
        : color_map_(color_map) {
        data_range_ = RBox2D(0.0, 1.0, 0.0, 1.0);
    }

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

    double point_radius() const {
        return point_radius_;
    }

    /**
     * Set default point radius.
     */
    void set_point_radius(double point_radius) {
        assert(point_radius > 0.0);

        point_radius_ = point_radius;
    }

    /**
     * Draw point.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object,
              const RGB32Color& color,
              const Pen::PointStyle& style,
              double radius,
              const std::string& name = "") {
        Object o(object, name);
        Pen pen(color);
        pen.point_style = style;
        pen.point_radius = radius;
        pen.fill_color.set_alpha(128);
        DrawObject(o, pen);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const RGB32Color& color,
              const std::string& name = "") {
        Draw(object, color, Pen::CIRCLE, point_radius_, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const RGB32Color& color,
              const Pen::PointStyle& style,
              const std::string& name = "") {
        Draw(object, color, style, point_radius_, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const RGB32Color& color,
              double radius, const std::string& name = "") {
        Draw(object, color, Pen::CIRCLE, radius, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const Pen::PointStyle& style,
              const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], style, point_radius_, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const Pen::PointStyle& style,
              double radius, const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], style, radius, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, double radius,
              const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], Pen::CIRCLE, radius, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const std::string& name = "") {
        Draw(object, color_map_[pens_.size()], Pen::PointStyle::CIRCLE, name);
    }
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const Pen& pen,
              const std::string& name = "") {
        Object o(object, name);
        DrawObject(o, pen);
    }

    /**
     * Draw plottable data on the terminal.
     */
    virtual void DrawData(Terminal* terminal) override {
        // Draw points.
        for (int i = 0; i < objects_.size(); ++i) {
            const Pen& pen = pens_[i];
            terminal->set_pen(pen);

            // Convert the coordinate for plotting.
            for (const RPoint2D& p : objects_[i].data_) {
                const RPoint2D q = ToPlotPosition(p);

                switch(pen.point_style) {
                case Pen::CIRCLE:
                    terminal->DrawCircle(q.x, q.y, pen.point_radius);
                    break;
                case Pen::RECTANGLE:
                    terminal->DrawRectangle(q.x - pen.point_radius,
                                            q.y + pen.point_radius,
                                            pen.point_radius * 2.0,
                                            pen.point_radius * 2.0);
                    break;
                }
            }
        }
    }

protected:
    /**
     * Draw points with specific pen.
     */
    void DrawObject(const Object& object, const Pen& pen) {
        assert(pen.point_radius > 0.0);

        if (object.empty()) return;

        // Update the bounding box.
        RBox2D box(object.data_.begin(), object.data_.end());
        if (objects_.empty())
            data_range_ = box;
        else
            data_range_.Join(box);

        objects_.push_back(object);
        pens_.push_back(pen);

        if (!object.name_.empty())
            legend_.InsertItem(Legend::POINT, object.name_, pen);
    }

    // Color map to draw the lines.
    const ColorMap& color_map_;

    // The input object for plotting.
    Array<Object> objects_;

    // The pens for each point data.
    Array<Pen> pens_;

    // The default radius for each point.
    double point_radius_ = 3.0;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_POINT_PLOT_H_
