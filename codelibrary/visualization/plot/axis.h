//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_PLOT_AXIS_H_
#define VISUALIZATION_PLOT_AXIS_H_

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <string>

#include "codelibrary/base/format.h"
#include "codelibrary/visualization/terminal/terminal.h"

namespace cl {
namespace plot {

/**
 * Axis for Plot.
 */
class Axis {
public:
    /**
     * Construct a axis by given range and maximal number of ticks.
     * Note that, the data range not necessarily the same as the axis range.
     *
     * Parameters:
     *  min       - the expected minimal value of axis.
     *  max       - the expected maximal value of axis.
     *  max_ticks - the expected maximal number of ticks on the axis.
     */
    Axis(double min = 0.0, double max = 1.0, int max_ticks = 10) {
        Reset(min, max, max_ticks);
    }

    /**
     * Reset this axis.
     */
    void Reset(double min, double max, int max_ticks = 10) {
        assert(min <= max);
        assert(max_ticks > 0);

        ChooseInterval(min, max, max_ticks);
        GenerateTicks(min, max);
    }

    /**
     * Return the minimum value of tick on the axis.
     */
    double min() const {
        return ticks_.front();
    }

    /**
     * Return the maximal value of tick on the axis.
     */
    double max() const {
        return ticks_.back();
    }

    /**
     * Return the length of axis.
     */
    double length() const {
        return max() - min();
    }

    /**
     * Return the interval between two neighbor ticks.
     */
    double interval() const {
        return interval_;
    }

    /**
     * Manually set an interval.
     */
    void SetInterval(double interval) {
        assert(interval > 0.0);
        assert(length() < INT_MAX * interval &&
               "The given interval can not be too small");

        interval_ = interval;
        GenerateTicks(min(), max());
    }

    /**
     * Return the tick values.
     */
    const Array<double>& ticks() const {
        return ticks_;
    }

    /**
     * Return the tick labels.
     */
    const Array<std::string>& tick_labels() const {
        return tick_labels_;
    }

    /**
     * Reverse the direction of axis.
     */
    void Reverse() {
        is_increasing_ ^= 1;
    }

    /**
     * Return true if the axis value is increasing.
     */
    bool is_increasing() const {
        return is_increasing_;
    }

private:
    /**
     * Automatically choose a usable ticking interval.
     */
    void ChooseInterval(double min, double max, int max_ticks) {
        // Order of magnitude of argument.
        double length = max - min;
        if (length == 0.0) length = 1.0;
        double power = std::pow(10.0, std::floor(std::log10(length)));

        // Approximate number of decades, we expect that 1 <= xnorm <= 10.
        double xnorm = length / power;
        xnorm = Clamp(xnorm, 1.0, 10.0);

        // Approximate number of ticks per decade.
        double n_ticks = max_ticks / xnorm;
        if (n_ticks > 20.0) {
            interval_ = 0.05;         // e.g. 0, 0.05, 0.10, ...
        } else if (n_ticks > 10.0) {
            interval_ = 0.1;          // e.g. 0, 0.1, 0.2, ...
        } else if (n_ticks > 5.0) {
            interval_ = 0.2;          // e.g. 0, 0.2, 0.4, ...
        } else if (n_ticks > 2.0) {
            interval_ = 0.5;          // e.g. 0, 0.5, 1, ...
        } else if (n_ticks > 1.0) {
            interval_ = 1.0;          // e.g. 0, 1, 2, ...
        } else if (n_ticks > 0.5) {
            interval_ = 2.0;          // e.g. 0, 2, 4, ...
        } else {
            interval_ = std::ceil(xnorm);
        }

        interval_ *= power;
    }

    /**
     * Generate tick values and labels.
     */
    void GenerateTicks(double min, double max) {
        int start = static_cast<int>(std::floor(min / interval_));
        int end   = static_cast<int>(std::ceil(max / interval_));
        assert(start <= end);

        ticks_.clear();
        for (int i = start; i <= end; ++i) {
            double x = i * interval_;
            ticks_.push_back(x);
        }
        if (ticks_.size() == 1) {
            ticks_.push_back((end + 1) * interval_);
        }

        // Generate tick labels.
        tick_labels_.clear();
        for (double x : ticks_) {
            tick_labels_.push_back(fmt::format("{:g}", x));
        }
    }

    // The interval between two neighbor ticks.
    double interval_ = 0.0;

    // Position of ticks.
    Array<double> ticks_;

    // Labels for every tick.
    Array<std::string> tick_labels_;

    // True if the values are increasing.
    bool is_increasing_ = true;
};

} // namespace plot
} // namespace cl

#endif // VISUALIZATION_PLOT_AXIS_H_
