//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_IO_LINE_READER_H_
#define UTIL_IO_LINE_READER_H_

#include <cassert>
#include <cstdio>
#include <cstring>
#include <memory>

#include "codelibrary/base/log.h"

namespace cl {
namespace io {

/**
 * A class to efficiently read large files line by line.
 *
 * For VS complier, one should add _CRT_SECURE_NO_WARNINGS option.
 */
class LineReader {
    // The size of the buffer.
    // We assume that the line length is no longer than BUFFER_SIZE.
    static const size_t BUFFER_SIZE;

public:
    LineReader() {}

    explicit LineReader(const char* filename) {
        if (OpenFile(filename)) 
            Initialize();
    }

    virtual ~LineReader() {
        if (file_) fclose(file_);
    }
     
    /**
     * Open file for reading.
     */
    bool Open(const char* filename) {
        if (!OpenFile(filename)) return false;

        Initialize();
        return true;
    }

    char* ReadLine() {
        if (data_begin_ == data_end_) return nullptr;

        ++n_line_;

        if (data_begin_ >= BUFFER_SIZE) {
            data_begin_ -= BUFFER_SIZE;
            data_end_   -= BUFFER_SIZE;
            std::memcpy(buffer_.get(), buffer_.get() + BUFFER_SIZE, 
                        BUFFER_SIZE);
            data_end_ += std::fread(buffer_.get() + BUFFER_SIZE, 1, BUFFER_SIZE,
                                    file_);
        }

        size_t i = data_begin_;
        while (i < data_end_ && buffer_[i] != '\n') ++i;

        if (i - data_begin_ >= BUFFER_SIZE) {
            LOG(INFO) << "Reading error in line " << n_line_ << std::endl
                      << "Line length limit exceeded: " 
                      << i - data_begin_ << " vs " << BUFFER_SIZE << ".";
            assert(false && "Line length limit exceeded.");
        }

        if (i != data_end_ && buffer_[i] == '\n') {
            buffer_[i] = '\0';
        } else {
            // We did not found '\n'.
            ++data_end_;
            buffer_[i] = '\0';
        }

        // Handle \r\n-line breaks.
        if (i != data_begin_ && buffer_[i - 1] == '\r') {
            buffer_[i - 1] = '\0';
        }

        char* result = buffer_.get() + data_begin_;
        data_begin_ = i + 1;
        return result;
    }

    int n_line() const {
        return n_line_;
    }

private:
    /**
     * Open file for reading.
     */
    bool OpenFile(const char* filename) {
        file_ = std::fopen(filename, "rb");
        
        if (!file_) {
            LOG(INFO) << "Cannot open file '" << filename << "' for reading.";
            return false;
        }

        return true;
    }

    /**
     * Initialize.
     */
    void Initialize() {
        n_line_ = 0;
        data_begin_ = 0;
        buffer_.reset(new char[2 * BUFFER_SIZE + 1]);
        data_end_ = std::fread(buffer_.get(), 1, 2 * BUFFER_SIZE, file_);
    }

    // The handle of the file.
    FILE* file_ = nullptr;

    // Data buffer.
    std::unique_ptr<char[]> buffer_;

    // The current number of lines.
    int n_line_ = 0;

    size_t data_begin_ = 0;
    size_t data_end_ = 0;
};

const size_t LineReader::BUFFER_SIZE = 1 << 24;

} // namespace io
} // namespace cl

#endif // UTIL_IO_LINE_READER_H_
