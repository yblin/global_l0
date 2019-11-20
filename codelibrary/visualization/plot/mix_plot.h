//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_MIX_PLOT_H_
#define VISUALIZATION_PLOT_MIX_PLOT_H_

#include <cassert>
#include <memory>
#include <string>
#include <type_traits>

#include "codelibrary/visualization/plot/plot.h"

namespace cl {
namespace plot {

/**
 * Mix plot.
 *
 * This class is used to plot several different plots on the same axes.
 *
 * Sample usage:
 *
 *   plot::MixPlot mix_plot;
 *
 *   // Add plots.
 *   mix_plot.AddPlot(line_plot);
 *   mix_plot.AddPlot(scatter_plot);
 *
 *   // Save to SVG.
 *   mix_plot.Save("test.svg");
 */
class MixPlot : public Plot {
public:
    MixPlot() = default;

    /**
     * Return true if this plot is empty.
     */
    virtual bool empty() const override {
        return plots_.empty();
    }

    /**
     * Clear this plot.
     */
    virtual void clear() override {
        plots_.clear();
    }

    /**
     * Add a plot.
     */
    template <typename PlotType>
    void AddPlot(const PlotType& plot) {
        static_assert(std::is_base_of<Plot, PlotType>::value,
                      "The input parameter must be derived form Plot");

        if (plot.empty()) return;

        if (!axes_fixed_) {
            if (plots_.empty()) {
                data_range_ = plot.data_range();
            } else {
                if (!axes_fixed_) data_range_.Join(plot.data_range());
            }
        }
        plots_.push_back(std::make_shared<PlotType>(plot));
    }

    /**
     * Draw plottable data on the terminal.
     */
    virtual void DrawData(Terminal* terminal) override {
        assert(terminal);

        for (const auto& plot : plots_) {
            plot->set_data_range(data_range_);
            plot->set_data_area(data_area_);
            plot->DrawData(terminal);
            legend_.InsertLegend(plot->legend());
        }
    }

private:
    Array<std::shared_ptr<Plot> > plots_;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_MIX_PLOT_H_
