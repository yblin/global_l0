//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_TIME_COUNTER_H_
#define BASE_TIME_COUNTER_H_

#include <cstdint>
#include <ctime>
#include <string>

#include "codelibrary/base/units.h"

namespace cl {

/**
 * A timer for calculating user's process time.
 */
class TimeCounter {
public:
    TimeCounter() = default;

    /**
     * Reset the timer.
     */
    void Reset() {
        elapsed_ = 0;
        if (running_) {
            started_ = clock(); // clock() is faster than std::chrono.
        }
    }

    /**
     * Timer start.
     */
    void Start() {
        running_ = true;
        started_ = clock();
    }

    /**
     * Timer stop.
     */
    void Stop() {
        if (running_) {
            running_ = false;
            elapsed_ += clock() - started_;
        }
    }

    /**
     * Return the elapsed time in seconds.
     */
    double elapsed_seconds() const {
        return static_cast<double>(elapsed_) / CLOCKS_PER_SEC;
    }

    /**
     * Return the human readable elapsed time.
     */
    std::string elapsed_time() const {
        return unit::time::ToString(elapsed_seconds());
    }

private:
    // True if timer is running.
    bool running_ = false;

    // The elapsed time.
    int64_t elapsed_ = 0;

    // The time of started.
    int64_t started_ = 0;
};

} // namespace cl

#endif // BASE_TIME_COUNTER_H_
