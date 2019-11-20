//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_SUB_PLOT_H_
#define VISUALIZATION_PLOT_SUB_PLOT_H_

#include <cassert>
#include <memory>

#include "codelibrary/base/array_nd.h"
#include "codelibrary/visualization/plot/plot.h"

namespace cl {
namespace plot {

/**
 * Sub Plot.
 *
 * It organizes multiple plots through a 2D grid layout.
 */ 
class SubPlot : public Plot {
public:
    /**
     * Subplot divides the current figure into an m-by-n grid.
     */ 
    SubPlot(int n_rows, int n_columns) {
        assert(n_rows > 0);
        assert(n_columns > 0);

        plots_.reshape(n_rows, n_columns);
        this->axes_visibility_ = false;
        this->axes_ticks_visibility_ = false;
        this->axes_mirror_visibility_ = false;
        sub_font_size_ = font_.size();
        x_spacing_ = 10.0;
        y_spacing_ = 10.0;
    }

    /**
     * Return true if the plotter is empty.
     */ 
    virtual bool empty() const override {
        return plots_.empty();
    }

    // Clear this plot.
    virtual void clear() override {
        plots_.clear();
    }

    // Assign a plot to the cell (m, n).
    template <typename PlotType>
    void SetPlot(int m, int n, const PlotType& plot) {
        static_assert(std::is_base_of<Plot, PlotType>::value,
                      "The input parameter must be derived form Plot");
        assert(m >= 0 && m < plots_.shape(0));
        assert(n >= 0 && n < plots_.shape(1));

        if (plot.empty()) return;
        plots_(m, n) = std::make_shared<PlotType>(plot);
    }

    // Set the X spacing between plots.
    void set_x_spacing(double x_spacing) {
        assert(x_spacing_ >= 0.0);

        x_spacing_ = x_spacing;
    }

    // Set the Y spacing between plots.
    void set_y_spacing(double y_spacing) {
        assert(y_spacing_ >= 0.0);

        y_spacing_ = y_spacing;
    }

    // Set font size for sub plots.
    void set_sub_font_size(double sub_font_size) {
        assert(sub_font_size_ >= 0.0);

        sub_font_size_ = sub_font_size;
    }

protected:
    // Draw all subplots on terminal.
    virtual void DrawData(Terminal* terminal) override {
        assert(terminal);

        for (int i = 0; i < plots_.shape(0); ++i) {
            for (int j = 0; j < plots_.shape(1); ++j) {
                if (plots_(i, j) != nullptr) {
                    terminal->set_plot_area(ComputeSubplotArea(i, j, terminal));
                    Font font = font_;
                    font.set_size(sub_font_size_);
                    plots_(i, j)->set_font(font);
                    plots_(i, j)->UpdateAxes();
                    plots_(i, j)->EvaluatePlotAreas(terminal);
                    plots_(i, j)->DrawData(terminal);
                    plots_(i, j)->DrawAxes(terminal);
                    plots_(i, j)->DrawTitle(terminal);
                }
            }
        }
        terminal->set_plot_area(RBox2D(0.0, terminal->width(),
                                       0.0, terminal->height()));
    }

    // Compute subplot area of cell (m, n)
    RBox2D ComputeSubplotArea(int m, int n, Terminal* terminal) const {
        m = plots_.shape(0) - m - 1;
        double x = axes_area_.x_length() / plots_.shape(1);
        double y = axes_area_.y_length() / plots_.shape(0);
        double x_min = axes_area_.x_min() + x * n;
        double y_min = axes_area_.y_min() + y * m;

        double x_min1 = x_min + x_spacing_ * 0.5;
        double x_max1 = std::max(x_min + x - x_spacing_ * 0.5, x_min1);
        double y_min1 = y_min + y_spacing_ * 0.5;
        double y_max1 = std::max(y_min + y - y_spacing_ * 0.5, y_min1);

        x_max1 = std::min(x_max1, static_cast<double>(terminal->width()));
        y_max1 = std::min(y_max1, static_cast<double>(terminal->height()));
        x_min1 = std::min(x_min1, x_max1);
        y_min1 = std::min(y_min1, y_max1);

        return {x_min1, x_max1, y_min1, y_max1};
    }

    double x_spacing_;
    double y_spacing_;
    double sub_font_size_;
    ArrayND<std::shared_ptr<Plot> > plots_;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_SUB_PLOT_H_
