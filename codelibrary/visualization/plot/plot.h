//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_PLOT_H_
#define VISUALIZATION_PLOT_PLOT_H_

#include <cassert>
#include <string>

#include "codelibrary/string/string_split.h"
#include "codelibrary/visualization/plot/axis.h"
#include "codelibrary/visualization/plot/legend.h"
#include "codelibrary/visualization/terminal/svg_terminal.h"

namespace cl {
namespace plot {

/**
 * Plot geometric objects on XY axes.
 *
 * The following is the framework of the plot.
 *
 *    -----------------------------------
 *    |              Title              |
 *    -----------------------------------
 *    |   |   |                         |
 *    | V |   |                         |
 *    | e |   |                         |
 *    | r | Y |                         |
 *    | t |   |                         |
 *    | i | A |                         |
 *    | c | x |       Data area         |
 *    | a | i |                         |
 *    | l | x |                         |
 *    |   |   |                         |
 *    | t |   |                         |
 *    | i |   |                         |
 *    | t |   |                         |
 *    | l |   |                         |
 *    | e |   --------------------------|
 *    |   |         X Axis              |
 *     --  -----------------------------|
 *    |   |      Horizental title       |
 *     ---------------------------------
 */
class Plot {
public:
    Plot() {
        Initialize();
    }

    virtual ~Plot() = default;

    /**
     * Return true if this plot is empty.
     */
    virtual bool empty() const = 0;

    /**
     * Clear this plot.
     */
    virtual void clear() = 0;

    /**
     * Show the plot on the given terminal.
     */
    virtual void Show(Terminal* terminal) {
        assert(terminal);

        terminal->clear();

        UpdateAxes();
        EvaluatePlotAreas(terminal);
        DrawData(terminal);
        DrawAxes(terminal);
        DrawTitle(terminal);
        DrawLegend(terminal);
    }

    /**
     * Save the plot into the given file.
     */
    virtual void Save(const std::string& filename,
                      int width = 640,
                      int height = 480) {
        Array<std::string> splits;
        StringSplit(filename, '.', &splits);
        assert(!splits.empty());

        const std::string& suffix = splits.back();
        if (suffix == "svg" || suffix == "SVG") {
            SVGTerminal terminal(width, height);
            Show(&terminal);
            terminal.SaveToFile(filename);
        } else {
            // Currently we only support svg file.
            assert(false);
        }
    }

    /**
     * Return the bounding box of the input data.
     */
    const RBox2D& data_range() const {
        return data_range_;
    }

    /**
     * Return the title of the plot.
     */
    const std::string& title() const {
        return title_;
    }

    /**
     * Return the vertical title of the plot.
     */
    const std::string& vertical_title() const {
        return vertical_title_;
    }

    /**
     * Return the horizental title of the plot.
     */
    const std::string& horizental_title() const {
        return horizental_title_;
    }

    /**
     * Return the position of legend.
     */
    const Legend::Position& legend_position() const {
        return legend_.position();
    }

    /**
     * Return the legend of this plot.
     */
    const Legend& legend() const {
        return legend_;
    }

    /**
     * Set the font of plot.
     */
    void set_font(const Font& font) {
        font_ = font;
    }

    /**
     * Set plot title.
     */
    void set_title(const std::string& title) {
        title_ = title;
    }

    /**
     * Set horizental title.
     */
    void set_horizental_title(const std::string& title) {
        horizental_title_ = title;
    }

    /**
     * Set vertical title.
     */
    void set_vertical_title(const std::string& title) {
        vertical_title_ = title;
    }

    /**
     * Show or hide the axes.
     * Note that, we will refresh the visibility states of ticks and mirror.
     */
    void set_axes_visibility(bool visibility) {
        axes_visibility_ = visibility;
        axes_ticks_visibility_ = visibility;
        axes_mirror_visibility_ = visibility;
    }

    /**
     * Show or hide the axes ticks and labels.
     */
    void set_axes_ticks_visibility(bool visibility) {
        axes_ticks_visibility_ = visibility;
    }

    /**
     * Show or hide the mirror axes lines.
     */
    void set_mirror_axes_visibility(bool visibility) {
        axes_mirror_visibility_ = visibility;
    }

    /**
     * Set axes fixed.
     */
    void set_axes_fixed(bool axes_fixed) {
        axes_fixed_ = axes_fixed;
    }

    /**
     * Set data range of axes.
     */
    void set_data_range(const RBox2D& data_range) {
        data_range_ = data_range;
        UpdateAxes();
    }

    /**
     * Set data area for plotting data.
     */
    void set_data_area(const RBox2D& data_area) {
        data_area_ = data_area;
    }

    /**
     * Set the position of legend.
     */
    void set_legend_position(const Legend::Position& position) {
        legend_.set_position(position);
    }

    /**
     * Reverse the direction of X axis.
     */
    void ReverseXAxis() {
        x_axis_.Reverse();
    }

    /**
     * Reverse the direction of Y axis.
     */
    void ReverseYAxis() {
        y_axis_.Reverse();
    }

    /**
     * Draw plottable data on the terminal.
     */
    virtual void DrawData(Terminal* terminal) = 0;

    /**
     * Draw the title on the terminal.
     */
    void DrawTitle(Terminal* terminal) const {
        assert(terminal);

        double title_font_size = 2.0 * font_.size();
        if (!title_.empty()) {
            title_font_size = std::min(title_area_.y_length(), title_font_size);
        }

        double subtitle_font_size = title_font_size * 0.7;
        if (!vertical_title_.empty()) {
            subtitle_font_size = std::min(vertical_title_area_.x_length(),
                                          subtitle_font_size);
        }
        if (!horizental_title_.empty()) {
            subtitle_font_size = std::min(horizental_title_area_.y_length(),
                                          subtitle_font_size);
        }

        if (!title_.empty()) {
            // Show title.
            Font font = font_;
            font.set_size(title_font_size);
            font.set_aligment(Font::MIDDLE);
            font.set_bold(true);
            terminal->set_font(font);

            double x = (title_area_.x_min() + title_area_.x_max()) * 0.5;
            double y = title_area_.y_max();
            terminal->DrawText(x, y, title_);
        }

        if (!vertical_title_.empty()) {
            // Show vertical title.
            Font font = font_;
            font.set_size(subtitle_font_size);
            font.set_aligment(Font::MIDDLE);
            terminal->set_font(font);

            double x = vertical_title_area_.x_min() + subtitle_font_size;
            double y = 0.5 * (vertical_title_area_.y_min() +
                              vertical_title_area_.y_max()) +
                       subtitle_font_size;
            terminal->DrawVerticalText(x, y, vertical_title_);
        }

        if (!horizental_title_.empty()) {
            // Show horizental title.
            Font font = font_;
            font.set_size(subtitle_font_size);
            font.set_aligment(Font::MIDDLE);
            terminal->set_font(font);

            double x = 0.5 * (horizental_title_area_.x_min() +
                              horizental_title_area_.x_max());
            double y = horizental_title_area_.y_min() + subtitle_font_size;
            terminal->DrawText(x, y, horizental_title_);
        }
    }

    /**
     * Draw axes on the terminal.
     */
    void DrawAxes(Terminal* terminal) {
        assert(terminal);

        if (!axes_visibility_) return;

        // The tick length: default is 5.0 (5 pixels).
        const double tick_length = 5.0;

        // Set terminal pen and font.
        Pen pen(RGB32Color::Black());
        terminal->set_pen(pen);

        Font font = axes_font_;

        // Draw X axis.
        font.set_aligment(Font::MIDDLE);
        terminal->set_font(font);

        // Draw ticks and labels.
        if (axes_ticks_visibility_) {
            for (int i = 0; i < x_axis_.ticks().size(); ++i) {
                const std::string& label = x_axis_.tick_labels()[i];
                double x = ToXPlotPosition(x_axis_.ticks()[i]);

                // The Y position of tick on the axis.
                double tick = data_area_.y_min();
                terminal->DrawLine(x, tick, x, tick + tick_length);
                terminal->DrawText(x, tick, label);

                if (axes_mirror_visibility_) {
                    // The Y position of mirror tick on the axis.
                    double mirror_tick = data_area_.y_max();
                    terminal->DrawLine(x, mirror_tick,
                                       x, mirror_tick - tick_length);
                }
            }
        }

        // Draw X axis lines.
        terminal->DrawLine(data_area_.x_min(), data_area_.y_min(),
                           data_area_.x_max(), data_area_.y_min());
        if (axes_mirror_visibility_) {
            terminal->DrawLine(data_area_.x_min(), data_area_.y_max(),
                               data_area_.x_max(), data_area_.y_max());
        }

        // Draw horizonal axis.
        font.set_aligment(Font::END);
        terminal->set_font(font);

        // Adjust the label position.
        const double label_y_shift = font.size() * 0.75;
        const double label_x_shift = font.size() * 0.75;

        // Draw Y ticks and labels.
        if (axes_ticks_visibility_) {
            for (int i = 0; i < y_axis_.ticks().size(); ++i) {
                const std::string& label = y_axis_.tick_labels()[i];
                double y = ToYPlotPosition(y_axis_.ticks()[i]);

                // The X position of tick on the axis.
                double tick = data_area_.x_min();
                terminal->DrawLine(tick, y, tick + tick_length, y);
                terminal->DrawText(tick - label_x_shift, y + label_y_shift,
                                   label);

                if (axes_mirror_visibility_) {
                    // The X position of mirror tick on the axis.
                    double mirror_tick = data_area_.x_max();
                    terminal->DrawLine(mirror_tick, y,
                                       mirror_tick - tick_length, y);
                }
            }
        }

        // Draw Y axis line.
        terminal->DrawLine(data_area_.x_min(), data_area_.y_min(),
                           data_area_.x_min(), data_area_.y_max());
        if (axes_mirror_visibility_) {
            terminal->DrawLine(data_area_.x_max(), data_area_.y_min(),
                               data_area_.x_max(), data_area_.y_max());
        }
    }

    /**
     * Convert the axis point to the plot position.
     */
    RPoint2D ToPlotPosition(const RPoint2D& p) const {
        return {ToXPlotPosition(p.x), ToYPlotPosition(p.y)};
    }

    /**
     * Convert the data value to the X plot position.
     */
    double ToXPlotPosition(double value) const {
        if (data_area_.x_min() >= data_area_.x_max()) return data_area_.x_min();
        if (x_axis_.min() >= x_axis_.max()) return x_axis_.min();

        double plot_length = data_area_.x_length();
        if (x_axis_.is_increasing()) {
            return (value - x_axis_.min()) / x_axis_.length() * plot_length +
                   data_area_.x_min();
        }

        return (x_axis_.max() - value) / x_axis_.length() * plot_length +
               data_area_.x_min();
    }

    /**
     * Convert the value to the Y plot position.
     */
    double ToYPlotPosition(double value) const {
        if (data_area_.y_min() >= data_area_.y_max()) return data_area_.y_min();
        if (y_axis_.min() >= y_axis_.max()) return y_axis_.min();

        double plot_length = data_area_.y_length();
        if (y_axis_.is_increasing()) {
            return (value - y_axis_.min()) / y_axis_.length() * plot_length +
                   data_area_.y_min();
        }

        return (y_axis_.max() - value) / y_axis_.length() * plot_length +
               data_area_.y_min();
    }

    /**
     * Evaluate plot areas.
     */
    void EvaluatePlotAreas(Terminal* terminal) {
        const RBox2D& plot_area = terminal->plot_area();
        EvaluateTitleArea(plot_area);
        EvaluateAxesArea(plot_area);
        EvaluateAxesFontSize();
        EvaluateDataArea();
        horizental_title_area_ = RBox2D(data_area_.x_min(),
                                        data_area_.x_max(),
                                        horizental_title_area_.y_min(),
                                        horizental_title_area_.y_max());
    }

    /**
     * Update XY axes.
     */
    void UpdateAxes() {
        const int max_ticks = 10;
        x_axis_.Reset(data_range_.x_min(), data_range_.x_max(), max_ticks);
        y_axis_.Reset(data_range_.y_min(), data_range_.y_max(), max_ticks);
    }

protected:
    /**
     * Set default option values.
     */
    void Initialize() {
        font_ = Font("times", 12.0);
        data_range_ = RBox2D(0.0, 1.0, 0.0, 1.0);
        axes_visibility_ = true;
        axes_ticks_visibility_ = true;
        axes_mirror_visibility_ = true;
        axes_fixed_ = false;
    }

    /**
     * Evaluate plot area for title.
     */
    void EvaluateTitleArea(const RBox2D& plot_area) {
        double font_size = font_.size() * 2.0;
        double title_height = std::min(0.15 * plot_area.y_length(), font_size);
        if (title_.empty()) title_height = 0.0;

        double font_size1 = font_.size() * 1.5;
        double vertical_title_width = std::min(0.1 * plot_area.x_length(),
                                               font_size1);
        if (vertical_title_.empty()) vertical_title_width = 0.0;

        double horizental_title_height = std::min(0.1 * plot_area.y_length(),
                                                  font_size1);
        if (horizental_title_.empty()) horizental_title_height = 0.0;

        title_area_ = RBox2D(plot_area.x_min() + vertical_title_width * 2,
                             plot_area.x_max(),
                             plot_area.y_max() - title_height,
                             plot_area.y_max());

        double x = plot_area.x_min() + vertical_title_width;
        double y = plot_area.y_min() + horizental_title_height;
        vertical_title_area_ = RBox2D(plot_area.x_min(), x,
                                      y, plot_area.y_max());
        horizental_title_area_ = RBox2D(x, plot_area.x_max(),
                                        plot_area.y_min(), y);
    }

    /**
     * Evaluate axes area.
     */
    void EvaluateAxesArea(const RBox2D& plot_area) {
        double x_spacing = 0.2 * vertical_title_area_.x_length();
        double y_spacing = 0.5 * title_area_.y_length();
        if (vertical_title_.empty()) x_spacing = 0.0;
        if (horizental_title_.empty()) y_spacing = 0.0;

        axes_area_ = RBox2D(vertical_title_area_.x_max() + x_spacing,
                            plot_area.x_max(),
                            horizental_title_area_.y_max() + y_spacing * 0.5,
                            title_area_.y_min() - 0.5 * title_area_.y_length());
    }

    /**
     * Evaluate data area.
     */
    void EvaluateDataArea() {
        double x_axis_height = 0.0;
        double y_axis_width = 0.0;
        if (axes_ticks_visibility_) {
            size_t s = 0;
            for (const std::string& label : y_axis_.tick_labels()) {
                s = std::max(label.size(), s);
            }
            ++s;
            double font_size = axes_font_.size();
            x_axis_height = std::min(0.2 * axes_area_.y_length(), font_size);
            y_axis_width = std::min(0.2 * axes_area_.x_length(), s * font_size);
        }

        data_area_ = RBox2D(axes_area_.x_min() + y_axis_width,
                            axes_area_.x_max() - y_axis_width * 0.5,
                            axes_area_.y_min() + x_axis_height,
                            axes_area_.y_max() - x_axis_height);
    }

    /**
     * Evaluate the font size of axes.
     */
    void EvaluateAxesFontSize() {
        double s1 = 0.0;
        for (const std::string& label : x_axis_.tick_labels()) {
            s1 += label.size();
        }

        double s2 = y_axis_.tick_labels().size();

        double size = std::min(font_.size(), 1.5 * axes_area_.x_length() / s1);
        size = std::min(size, axes_area_.y_length() / s2);
        axes_font_.set_size(size);
    }

    /**
     * Draw legend.
     */
    void DrawLegend(Terminal* terminal) {
        double height = font_.size() * 2.0;
        double width  = 2.0 * height;

        const Legend::Position& position = legend_.position();

        double x, y;
        if (position == Legend::LEFT_BOTTOM || position == Legend::LEFT_TOP)
            x = this->data_area_.x_min() + width * 0.5;
        else
            x = this->data_area_.x_max() - width * 1.5;
        if (position == Legend::LEFT_TOP || position == Legend::RIGHT_TOP)
            y = this->data_area_.y_max() - height * 1.5;
        else
            y = this->data_area_.y_min() + height * 1.5;

        Font font = font_;
        font.set_size(font_.size() * 1.2);
        if (position == Legend::RIGHT_TOP || position == Legend::RIGHT_BOTTOM) {
            font.set_aligment(Font::END);
        } else {
            font.set_aligment(Font::START);
        }
        terminal->set_font(font);

        int n_items = legend_.items().size();

        if (legend_.position() == Legend::RIGHT_TOP) {
            for (int i = 0; i < n_items; ++i) {
                const Legend::Item& item = legend_.items()[i];
                double y1 = y - height * i;
                double y2 = y1 + height * 0.5;
                double y3 = y1 - height * 0.5;

                item.DrawLegend(cl::RBox2D(x, x + width, y3, y2), terminal);
                terminal->DrawText(x - 0.5 * height, y2, item.name);
            }
        } else if (legend_.position() == Legend::RIGHT_BOTTOM) {
            for (int i = n_items - 1; i >= 0; --i) {
                const Legend::Item& item = legend_.items()[i];
                double y1 = y + height * (n_items - 1 - i);
                double y2 = y1 + height * 0.5;
                double y3 = y1 - height * 0.5;

                item.DrawLegend(cl::RBox2D(x, x + width, y3, y2), terminal);
                terminal->DrawText(x - 0.5 * height, y2, item.name);
            }
        } else if (legend_.position() == Legend::LEFT_BOTTOM) {
            for (int i = n_items - 1; i >= 0; --i) {
                const Legend::Item& item = legend_.items()[i];
                double y1 = y + height * (n_items - 1 - i);
                double y2 = y1 + height * 0.5;
                double y3 = y1 - height * 0.5;

                item.DrawLegend(cl::RBox2D(x, x + width, y3, y2), terminal);
                terminal->DrawText(x + width + 0.5 * height, y2, item.name);
            }
        } else {
            for (int i = 0; i < n_items; ++i) {
                const Legend::Item& item = legend_.items()[i];
                double y1 = y - height * i;
                double y2 = y1 + height * 0.5;
                double y3 = y1 - height * 0.5;

                item.DrawLegend(cl::RBox2D(x, x + width, y3, y2), terminal);
                terminal->DrawText(x + width + 0.5 * height, y2, item.name);
            }
        }
    }

    // Font for plot title.
    Font font_;

    // Font for plot axes.
    Font axes_font_;

    // Bounding box of the input data.
    RBox2D data_range_;

    // The title of plotter.
    std::string title_;

    // The horizental title.
    std::string horizental_title_;

    // The vertical title.
    std::string vertical_title_;

    // XY axes.
    Axis x_axis_, y_axis_;

    // Area to plot vertical title.
    RBox2D vertical_title_area_;

    // Area to plot horizental title.
    RBox2D horizental_title_area_;

    // Area to plot title.
    RBox2D title_area_;

    // Area to plot XY axes and data.
    RBox2D axes_area_;

    // Area to plot data.
    RBox2D data_area_;

    // Set true to show the major axes.
    bool axes_visibility_ = true;

    // Set true to show the axes ticks and labels.
    bool axes_ticks_visibility_ = true;

    // Set true to show the mirror lines of axes.
    bool axes_mirror_visibility_ = true;

    // Set true to fix the axes.
    bool axes_fixed_ = false;

    // Plot legend.
    Legend legend_;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_PLOT_H_
