//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_IO_FILE_MAP_H_
#define UTIL_IO_FILE_MAP_H_

#include <cassert>
#include <cstdio>
#include <cinttypes>
#include <string>

#include "codelibrary/base/log.h"

namespace cl {
namespace io {

/**
 * Map file to string.
 *
 * Return false if an error occurs.
 */
inline bool MapFile(const std::string& filename, Array<char>* data) {
    assert(data);

    FILE* file = std::fopen(filename.c_str(), "rb");
    if (!file) {
        LOG(INFO) << "Can't open input file: " << filename;
        return false;
    }

    off_t buffer_size = 8192;
    if (std::fseek(file, 0, SEEK_END) == 0) {
        buffer_size = std::max<off_t>(std::ftell(file), 1);
        if (std::fseek(file, 0, SEEK_SET) != 0) {
            LOG(INFO) << "fseek error, when reading file: " << filename;
            return false;
        }
    } else if (std::ferror(file)) {
        LOG(INFO) << "fseek error, when reading file: " << filename;
        return false;
    }

    Array<char> buf(buffer_size);
    while (!std::feof(file)) {
        std::size_t read_bytes = std::fread(buf.data(), 1, buffer_size, file);
        if (std::ferror(file)) {
            LOG(INFO) << "fread error, when reading file: " << filename;
            return false;
        }
        data->insert(buf.begin(), buf.begin() + read_bytes);
    }

    std::fclose(file);
    return true;
}

} // namespace io
} // namespace cl

#endif // UTIL_IO_FILE_MAP_H_
