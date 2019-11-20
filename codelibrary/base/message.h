//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_MESSAGE_H_
#define BASE_MESSAGE_H_

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include "codelibrary/base/array.h"
#include "codelibrary/base/format.h"

namespace cl {

/**
 * The Message class works like an ostream.
 *
 * Typical usage:
 *
 *   1. First, stream a bunch of values to a Message object.
 *      It will store the text in a stringstream.
 *   2. Then stream the Message object to an ostream.
 *      This causes the text in the Message to be streamed to the ostream.
 *
 * For example;
 *
 *   Message foo;
 *   foo << 1 << " != " << 2;
 *   std::cout << foo;
 *
 * will print "1 != 2".
 */
class Message {
    // The type of basic IO manipulators (endl, ends, and flush) for narrow
    // streams.
    using BasicNarrowIoManip = std::ostream &(*)(std::ostream &);

public:
    /**
     * Constructs an empty Message.
     */
    Message() = default;

    /**
     * Copy constructor.
     */
    Message(const Message& message) {
        stream_ << message.ToString();
    }

    /**
     * Construct a Message from a C-string.
     */
    template <typename T>
    explicit Message(const T& n) {
        this->operator<<(n);
    }

    /**
     * Stream a non-pointer value to this object.
     */
    template <typename T>
    Message& operator <<(const T& value) {
        stream_ << fmt::format("{}", value);
        return *this;
    }

    /**
     * Stream a pointer value to this object.
     * This function is an overload of the previous one.
     */
    template <typename T>
    Message& operator <<(T* pointer) {
        if (pointer == nullptr) {
            stream_ << "(nullptr)";
        } else { 
            stream_ << pointer;
        }
        return *this;
    }

    /**
     * Since the basic IO manipulators are overloaded for both narrow and wide
     * streams, we have to provide this specialized definition of operator <<,
     * even though its body is the same as the templated version above.
     *
     * Without this definition, streaming endl or other basic IO manipulators to
     * Message will confuse the compiler.
     */
    Message& operator <<(BasicNarrowIoManip value) {
        stream_ << value;
        return *this;
    }

    /**
     * Instead of 1/0, we want to see true/false for bool values.
     */
    Message& operator <<(bool b) {
        return *this << (b ? "true" : "false");
    }

    /**
     * These two overloads allow streaming a wide C string to a Message using
     * the UTF-8 encoding.
     */
    Message& operator <<(const wchar_t* wide_c_str) {
        return *this << WStringToUtf8(wide_c_str);
    }
    Message& operator <<(wchar_t* wide_c_str) {
        return *this << WStringToUtf8(wide_c_str);
    }

    /**
     * Convert the given wide string to a narrow string using the UTF-8
     * encoding and stream the result to this Message object.
     */
    Message& operator <<(const std::wstring& wstring) {
        return *this << WStringToUtf8(wstring);
    }

    /**
     * Message of STL pair.
     */
    template <class T1, class T2>
    Message& operator <<(const std::pair<T1, T2>& p) {
        return *this << fmt::format("({}, {})", p.first, p.second);
    }

    /**
     * Get the text streamed to this object so far as an std::string.
     * When the input value is a std::string or std::wstring object, each NUL
     * character in it is replaced with "\\0".
     */
    std::string ToString() const {
        const std::string& str = stream_.str();
        std::string result;
        result.reserve(str.size());
        for (char i : str) {
            if (i == '\0') {
                result += "\\0";
            } else {
                result += i;
            }
        }

        return result;
    }

    /**
     * Streams a Message to an ostream.
     */
    friend std::ostream& operator <<(std::ostream& os, const Message& message) {
        return os << message.ToString();
    }

private:
    /** Convert a Unicode code point to a narrow string in UTF-8 encoding.
     *
     * code_point parameter is of unsigned 32-bit integer type because wchar_t
     * may not be wide enough to contain a code point.
     *
     * A Unicode code-point can have up to 21 bits and is encoded in UTF-8 like
     * this:
     *
     * Code-point length   Encoding
     *   0 -  7 bits       0xxxxxxx
     *   8 - 11 bits       110xxxxx 10xxxxxx
     *  12 - 16 bits       1110xxxx 10xxxxxx 10xxxxxx
     *  17 - 21 bits       11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
     *
     * If the code_point is not a valid Unicode code point
     * (i.e. outside of Unicode range U+0 to U+10FFFF) it will be converted to
     * "(Invalid Unicode 0xXXXXXXXX)".
     */
    static std::string UnicodeCodePointToUtf8(uint32_t c) {
        // The maximum code-point a one-byte UTF-8 sequence can represent.
        const uint32_t max_code_point1 = (static_cast<uint32_t>(1) << 7) - 1;

        // The maximum code-point a two-byte UTF-8 sequence can represent.
        const uint32_t max_code_point2 = (static_cast<uint32_t>(1) << 11) - 1;

        // The maximum code-point a three-byte UTF-8 sequence can represent.
        const uint32_t max_code_point3 = (static_cast<uint32_t>(1) << 16) - 1;

        // The maximum code-point a four-byte UTF-8 sequence can represent.
        const uint32_t max_code_point4 = (static_cast<uint32_t>(1) << 21) - 1;

        if (c > max_code_point4) {
            std::stringstream ss;
            ss << std::hex << std::uppercase << c;
            return "(Invalid Unicode 0x" + ss.str() + ")";
        }

        char str[5]; // Big enough for the largest valid code point.
        if (c <= max_code_point1) {
            str[1] = '\0';
            str[0] = static_cast<char>(c);                         // 0xxxxxxx
        } else if (c <= max_code_point2) {
            str[2] = '\0';
            str[1] = static_cast<char>(0x80 | (c & 0x3F));         // 10xxxxxx
            str[0] = static_cast<char>(0xC0 | (c >> 6));           // 110xxxxx
        } else if (c <= max_code_point3) {
            str[3] = '\0';
            str[2] = static_cast<char>(0x80 | (c & 0x3F));         // 10xxxxxx
            str[1] = static_cast<char>(0x80 | ((c >> 6) & 0x3F));  // 10xxxxxx
            str[0] = static_cast<char>(0xE0 | (c >> 12));          // 1110xxxx
        } else {
            str[4] = '\0';
            str[3] = static_cast<char>(0x80 | (c & 0x3F));         // 10xxxxxx
            str[2] = static_cast<char>(0x80 | ((c >> 6) & 0x3F));  // 10xxxxxx
            str[1] = static_cast<char>(0x80 | ((c >> 12) & 0x3F)); // 10xxxxxx
            str[0] = static_cast<char>(0xF0 | (c >> 18));          // 11110xxx
        }

        return str;
    }

    /**
     * Convert a wide string to a narrow string in UTF-8 encoding.
     * The wide string is assumed to have the following encoding:
     *   UTF-16 if sizeof(wchar_t) == 2 (on Windows, Cygwin, Symbian OS)
     *   UTF-32 if sizeof(wchar_t) == 4 (on Linux)
     *
     * Note If the string contains code points that are not valid Unicode code
     * points (i.e. outside of Unicode range U+0 to U+10FFFF) they will be
     * output as '(Invalid Unicode 0xXXXXXXXX)'.
     * If the string is in UTF16 encoding and contains invalid UTF-16 surrogate
     * pairs, values in those pairs will be encoded as individual Unicode
     * characters from Basic Normal Plane.
     */
    static std::string WStringToUtf8(const std::wstring& str) {
        size_t len = str.size();

        std::string result;
        for (size_t i = 0; i < len; ++i) {
            uint32_t code_point;

            if (str[i] == L'\0') break;

            if (i + 1 < len && sizeof(wchar_t) == 2 &&
                (str[i] & 0xFC00) == 0xD800 &&
                (str[i + 1] & 0xFC00) == 0xDC00) {
                const uint32_t mask = (1 << 10) - 1;
                code_point = (((str[i] & mask) << 10) |
                    (str[i + 1] & mask)) + 0x10000;
                ++i;
            } else {
                code_point = static_cast<uint32_t>(str[i]);
            }

            result += UnicodeCodePointToUtf8(code_point);
        }

        return result;
    }

    // We'll hold the text streamed to this object here.
    std::stringstream stream_;

    // The precision for floating type.
    int n_precision_ = 0;
};

} // namespace cl

#endif // BASE_MESSAGE_H_
