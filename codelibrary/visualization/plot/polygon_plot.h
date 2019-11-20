//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_POLYGON_PLOT_H_
#define VISUALIZATION_PLOT_POLYGON_PLOT_H_

#include <cassert>

#include "codelibrary/visualization/plot/plot.h"

namespace cl {
namespace plot {

/**
 * Polygon plot.
 */
class PolygonPlot : public Plot {
    // Object for PolygonPlotter.
    class Object {
        friend class PolygonPlot;
    public:
        Object() = default;

        template <typename T>
        explicit Object(const Polygon2D<T>& polygon) {
            SetData(MultiPolygon2D<T>(polygon));
        }

        template <typename T>
        explicit Object(const Box2D<T>& box) {
            SetData(MultiPolygon2D<T>(Polygon2D<T>(box)));
        }

        template <typename T>
        explicit Object(const MultiPolygon2D<T>& polygon) {
            SetData(polygon);
        }

        bool empty() const {
            return polygon_.empty();
        }

    private:
        template <typename T>
        void SetData(const MultiPolygon2D<T>& polygon) {
            using Boundary = typename MultiPolygon2D<T>::Boundary;
            for (const Boundary& b : polygon.boundaries()) {
                Array<RPoint2D> points;
                points.reserve(b.polygon.size());
                for (const Point2D<T>& p : b.polygon) {
                    points.emplace_back(p.x, p.y);
                }

                RPolygon2D polygon(points);
                if (!polygon.empty()) polygon_.Insert(polygon, b.is_outer);
            }
        }
        void SetData(const RMultiPolygon2D& polygon) {
            polygon_ = polygon;
        }

        RMultiPolygon2D polygon_; // The input polygon.
    };

public:
    explicit PolygonPlot(const ColorMap& color_map = ColorMap::Lines())
        : color_map_(color_map) {
        data_range_ = RBox2D(0.0, 1.0, 0.0, 1.0);
    }

    /**
     * Return true if this plotter is empty.
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
     * Draw polygon with fill color.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const RGB32Color& fill_color) {
        Object o(object);
        Pen pen(fill_color);
        pen.line_width = 0.0;
        DrawObject(o, pen);
    }

    /**
     * Draw polygon with specific pen.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object, const Pen& pen) {
        Object o(object);
        DrawObject(o, pen);
    }

    /**
     * Draw polygon with color_map's color.
     */
    template <typename PlottableObject>
    void Draw(const PlottableObject& object) {
        Draw(object, color_map_[pens_.size()]);
    }

    /**
     * Draw plottable data on the terminal.
     */
    virtual void DrawData(Terminal* terminal) override {
        // Draw polygons.
        for (int i = 0; i < objects_.size(); ++i) {
            Object obj = objects_[i];
            terminal->set_pen(pens_[i]);

            RMultiPolygon2D poly;
            for (int i = 0; i < obj.polygon_.n_boundaries(); ++i) {
                Array<RPoint2D> points;
                for (const RPoint2D& p : obj.polygon_.boundaries()[i].polygon) {
                    RPoint2D q = ToPlotPosition(p);
                    points.push_back(q);
                }
                poly.Insert(RPolygon2D(points),
                            obj.polygon_.boundaries()[i].is_outer);
            }

            terminal->DrawPolygon(poly);
        }
    }

protected:
    /**
     * Draw polygon determined with specific pen.
     */
    void DrawObject(const Object& object, const Pen& pen) {
        if (object.empty()) return;

        // Update the bounding box.
        RBox2D box = object.polygon_.bounding_box();
        if (objects_.empty())
            data_range_ = box;
        else
            data_range_.Join(box);

        objects_.push_back(object);
        pens_.push_back(pen);
    }

    // Color map to draw the lines.
    const ColorMap& color_map_;

    // The objects for plotting.
    Array<Object> objects_;

    // The pens for each line.
    Array<Pen> pens_;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_POLYGON_PLOT_H_
