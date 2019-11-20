//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_LEGEND_H_
#define VISUALIZATION_PLOT_LEGEND_H_

#include <string>

#include "codelibrary/visualization/pen.h"
#include "codelibrary/visualization/terminal/terminal.h"

namespace cl {
namespace plot {

/**
 * Legend for plot.
 */ 
class Legend {
public:
    // Type for item.
    enum ItemType {
        LINE, POINT
    };

    // Where to draw the legend in the plot.
    enum Position {
        RIGHT_TOP, LEFT_TOP, LEFT_BOTTOM, RIGHT_BOTTOM
    };

    // Item of Plot.
    struct Item {
        ItemType type;    // Type of this item.
        std::string name; // Item name.
        Pen pen;          // Pen to draw this item.

        Item() {}

        Item(const ItemType& t, const std::string& n, const Pen& p)
            : type(t), name(n), pen(p) {}

        /**
         * Draw item legend in the given box.
         */
        void DrawLegend(const RBox2D& box, Terminal* terminal) const {
            switch (type) {
            case LINE:
            {
                Pen pen1 = pen;
                pen1.line_width = 1;
                terminal->set_pen(pen1);
                terminal->DrawLine(box.x_min(), (box.y_min() +
                                                 box.y_max()) * 0.5,
                                   box.x_max(), (box.y_min() +
                                                 box.y_max()) * 0.5);
                break;
            }
            case POINT:
            {
                Pen pen1 = pen;
                pen1.line_width = 1;
                terminal->set_pen(pen1);
                double radius = 0.2 * box.y_length();

                cl::RPoint2D c1(box.x_min() + radius * 2.0,
                                box.y_min() + radius * 3.0);
                cl::RPoint2D c2((box.x_min() + box.x_max()) * 0.5,
                                box.y_max() - radius);
                cl::RPoint2D c3(box.x_max() - radius * 2.0,
                                box.y_min() + radius * 2.5);

                switch (pen.point_style) {
                case Pen::CIRCLE:
                    terminal->DrawCircle(c1.x, c1.y, radius);
                    terminal->DrawCircle(c2.x, c2.y, radius);
                    terminal->DrawCircle(c3.x, c3.y, radius);
                    break;
                case Pen::RECTANGLE:
                    terminal->DrawRectangle(c1.x - radius, c1.y + radius,
                                            2.0 * radius, 2.0 * radius);
                    terminal->DrawRectangle(c2.x - radius, c2.y + radius,
                                            2.0 * radius, 2.0 * radius);
                    terminal->DrawRectangle(c3.x - radius, c3.y + radius,
                                            2.0 * radius, 2.0 * radius);
                    break;
                default:
                    break;
                }
                break;
            }
            default:
                break;
            }
        }
    };

    /**
     * Add an item into the legend.
     */
    void InsertItem(const ItemType& type, const std::string& name,
                 const Pen& pen) {
        items_.emplace_back(type, name, pen);
    }

    /**
     * Insert one legend into this legend.
     */
    void InsertLegend(const Legend& legend) {
        for (const Item& item : legend.items_) {
            items_.push_back(item);
        }
    }

    void clear() {
        items_.clear();
    }

    const Array<Item>& items() const {
        return items_;
    }

    void set_position(const Position& position) {
        position_ = position;
    }

    const Position& position() const {
        return position_;
    }

private:
    Position position_ = RIGHT_TOP;
    Array<Item> items_;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_LEGEND_H_
