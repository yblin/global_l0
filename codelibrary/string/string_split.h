//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STRING_STRING_SPLIT_H_
#define STRING_STRING_SPLIT_H_

#include <cassert>
#include <cstring>
#include <string>

#include "codelibrary/base/array.h"

namespace cl {

/**
 * Split 'str' into a vector of strings delimited by 'c', placing the results
 * in 'result'. If several instances of 'c' are contiguous, or if 'str' begins
 * with or ends with 'c', then an empty string is inserted.
 */
template <typename String>
void StringSplit(const String& str, const typename String::value_type& c,
                 Array<String>* result) {
    assert(result);

    result->clear();
    size_t last = 0;
    size_t size = str.size();
    for (size_t i = 0; i <= size; ++i) {
        if (i == size || str[i] == c) {
            String tmp(str, last, i - last);

            // Avoid converting an empty or all-whitespace source string into a
            // vector of one empty string.
            if (i != size || !result->empty() || !tmp.empty()) {
                result->push_back(tmp);
            }
            last = i + 1;
        }
    }
}

/**
 * Split 'str' into a vector of strings by whitespace.
 * Whitespace is defined as: space, tab, LF, line tab, FF, or CR.
 */
template <typename String>
void StringSplit(const String& str, Array<String>* result) {
    assert(result);

    result->clear();
    const size_t length = str.size();
    if (!length) return;

    bool last_was_ws = false;
    size_t last_non_ws_start = 0;
    for (size_t i = 0; i < length; ++i) {
        switch (str[i]) {
        case L' ':
        case L'\t':
        case L'\xA':
        case L'\xB':
        case L'\xC':
        case L'\xD':
            if (!last_was_ws) {
                if (i > 0) {
                    result->push_back(str.substr(last_non_ws_start,
                                                 i - last_non_ws_start));
                }
                last_was_ws = true;
            }
            break;
        default:  // Not a space character.
            if (last_was_ws) {
                last_was_ws = false;
                last_non_ws_start = i;
            }
            break;
        }
    }
    if (!last_was_ws) {
        result->push_back(str.substr(last_non_ws_start,
                                     length - last_non_ws_start));
    }
}

/**
 * The same as SplitString, but use a substring delimiter instead of a char.
 */
template <typename String>
void StringSplit(const String& str, const String& s, Array<String>* r) {
    r->clear();

    size_t begin_index = 0;
    while (true) {
        const size_t end_index = str.find(s, begin_index);
        if (end_index == String::npos) {
            r->push_back(str.substr(begin_index));
            return;
        }
        r->push_back(str.substr(begin_index, end_index - begin_index));
        begin_index = end_index + s.size();
    }
}

/**
 * Split 'str' with 'c'. Special treat for wchar_t*.
 */
inline void StringSplit(const wchar_t* str, wchar_t c,
                        Array<std::wstring>* res) {
    StringSplit(std::wstring(str, str + wcslen(str)), c, res);
}

/**
 * Split 'str' with whitespace. Special treat for wchar_t*.
 */
inline void StringSplit(const wchar_t* str, Array<std::wstring>* res) {
    StringSplit(std::wstring(str, str + wcslen(str)), res);
}

/**
 * Split 'str' with a substring. Special treat for wchar_t*.
 */
inline void StringSplit(const wchar_t* str, const wchar_t* cs,
                        Array<std::wstring>* res) {
    StringSplit(std::wstring(str, str + wcslen(str)),
                std::wstring(cs, cs + wcslen(cs)),
                res);
}

/**
 * Split 'str' with 'c'. Special treat for char*.
 */
inline void StringSplit(const char* str, char c,
                        Array<std::string>* res) {
    StringSplit(std::string(str, str + strlen(str)), c, res);
}

/**
 * Split 'str' with whitespace. Special treat for char*.
 */
inline void StringSplit(const char* str, Array<std::string>* res) {
    StringSplit(std::string(str, str + strlen(str)), res);
}

/**
 * Split 'str' with a substring. Special treat for char*.
 */
inline void StringSplit(const char* str, const char* cs,
                        Array<std::string>* res) {
    StringSplit(std::string(str, str + strlen(str)),
                std::string(cs, cs + strlen(cs)),
                res);
}

/**
 * Split 'str' with a substring. Special treat for char*.
 */
inline void StringSplit(const std::string& s, const char* cs,
                        Array<std::string>* res) {
    StringSplit(s, std::string(cs), res);
}

} // namespace cl

#endif // STRING_STRING_SPLIT_H_
