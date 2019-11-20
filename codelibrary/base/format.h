//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_FORMAT_H_
#define BASE_FORMAT_H_

#include <algorithm>
#include <cassert>
#include <vector>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

// Formating library. It can be replaced by c++20<format> in future.
#include "third_party/fmt/format.h"
#include "third_party/fmt/ostream.h"

namespace fmt {

/**
 * Support to format every thing.
 */
template <typename Char, typename CharTraits, typename T>
std::basic_ostream<Char, CharTraits>& operator<<(
    std::basic_ostream<Char, CharTraits>& os, const T& object) {
    auto bytes = reinterpret_cast<const unsigned char*>(&object);

    int count = sizeof(object);

    // Tell the user how big the object is.
    os << count << "-byte object <";

    int i = 0;
    for (; i < std::max(128, count); ++i) {
        if (i != 0) {
            // Organizes the bytes into groups of 2 for easy parsing by human.
            if ((i % 2) == 0)
                os << ' ';
            else
                os << '-';
        }
        auto t = static_cast<int>(bytes[i]);
        if (t < 16) {
            os << std::hex << 0 << t;
        } else {
            os << std::hex << t;
        }
    }
    if (i < count) os << " ... ";

    os << ">";

    return os;
}

} // namespace fmt

namespace cl {

/**
 * Format printing.
 */
class Printer {
public:
    static Printer* GetInstance() {
        static Printer printer;
        return &printer;
    }

    /**
     * Python like print function.
     */
    template <typename T>
    void Print(const std::string& var, const T& a) {
        std::string str = fmt::format("{}", a);
        std::vector<std::string> lines;

        SplitLine(str, &lines);

        size_t width = 0;
        for (auto line : lines) {
            width = std::max(width, line.size());
        }

        if (width + var.size() + 3 > scrren_width_) {
            printf("%s = \n", var.c_str());
            printf("%s\n", str.c_str());
        } else {
            printf("%s = ", var.c_str());
            if (lines.empty()) return;

            std::string spaces(var.size() + 3, ' ');

            printf("%s", lines[0].c_str());
            for (size_t i = 1; i < lines.size(); ++i) {
                printf("%s%s", spaces.c_str(), lines[i].c_str());
            }
            printf("\n");
        }
    }
       
    /**
     * Format a sequence in one line, do not consider alignment.
     */
    template <typename Iterator>
    std::string Format(Iterator first, Iterator last) {
        std::size_t n = std::distance(first, last);
        if (n == 0) return "[]";

        std::string line = "[";
        if (n <= n_max_print_elements_ || n_max_print_elements_ == 0) {
            Iterator p = first;
            line += fmt::format("{}", *p++);
            for (; p != last; ++p) {
                line += fmt::format(", {}", *p);
            }
        } else {
            Iterator p = first;
            line += fmt::format("{}", *p++);
            for (int k = 1; k < (n_max_print_elements_ + 1) / 2; ++k) {
                line += fmt::format(", {}", *p++);
            }
            line += ", ...";
            std::advance(p, n - n_max_print_elements_);
            for (; p != last; ++p) {
                line += fmt::format(", {}", *p);
            }
        }
        return line + "]";
    }

    /**
     * Format a sequence in multi-lines with alignment.
     * 
     * This function is especially useful for printing matrix or vector.
     * It can greatly improve the human readability.
     *
     * For example,
     *   [ 1,  2,  3,
     *     4,  5,  6,
     *     7,  8,  9,
     *    10, 11, 12]
     * is much better than
     *   [1, 2, 3,
     *   4, 5, 6,
     *   7, 8, 9,
     *   10, 11, 12]
     * 
     * Parameters:
     *  [first, last)    - the input sequence.
     *  width            - the width for each element.
     *                     The default value is 0 which means that we will set
     *                     width automatically.
     *  n_leading_spaces - the number of leading spaces.
     *  n_tailing_spaces - the number of tailing spaces.
     */
    template <typename Iterator>
    std::string AlignFormat(Iterator first, Iterator last,
                            int width = 0,
                            int n_leading_spaces = 0,
                            int n_tailing_spaces = 0) {
        assert(width >= 0);
        assert(n_leading_spaces >= 0);
        assert(n_tailing_spaces >= 0);

        if (first == last) return "[]";

        if (width == 0) {
            // If width = 0 and the string can be printed in a single line.
            // We do not compute the width.
            std::string line = Format(first, last);
            if (line.size() <= scrren_width_) {
                return line;
            }
        }

        size_t n = std::distance(first, last);

        std::vector<std::string> token;
        int split_pos = -1;
        if (n <= n_max_print_elements_ || n_max_print_elements_ == 0) {
            for (Iterator p = first; p != last; ++p) {
                token.push_back(fmt::format("{}", *p));
            }
        } else {
            Iterator p = first;
            for (split_pos = 0; 
                 p != last && split_pos < (n_max_print_elements_ + 1) / 2;
                 ++p, ++split_pos) {
                token.push_back(fmt::format("{}", *p));
            }
            std::advance(p, n - n_max_print_elements_);
            for (; p != last; ++p) {
                token.push_back(fmt::format("{}", *p));
            }
        }

        if (width == 0) {
            for (const std::string& t : token) {
                assert(t.size() <= INT_MAX);
                width = std::max(static_cast<int>(t.size()), width);
            }
        }

        const int max_len_in_line = scrren_width_ - 7 - n_leading_spaces -
                                    n_tailing_spaces;
        int n_elements_per_line = std::max(1, max_len_in_line / (width + 2));

        std::string buf;
        std::string space;
        for (int i = 0; i < n_leading_spaces; ++i) 
            space.push_back(' ');
        
        // We do not add the leading space for the first line. Because we want
        // to get the following result:
        // 
        //  a = [1, 2, 3
        //       4, 5, 6]
        // 
        // but not
        //  
        //  a =     [1, 2, 3
        //       4, 5, 6]
        std::string line = "[";
        bool multi_line = false;
        int n_elements_in_line = 1;
        line += fmt::format("{:>{}}", token[0], width);
        for (int i = 1; i < token.size(); ++i) {
            line += ",";

            if (i == split_pos) {
                if (multi_line || line.size() + 5 > max_len_in_line) {
                    buf += fmt::format("{}\n{} {:>{}},\n", line, space, "...",
                                       width);
                    line = space;
                    n_elements_in_line = 0;
                } else {
                    line += " ...,";
                }
            }

            // If current token is not the last token.
            if (n_elements_in_line && 
                n_elements_in_line % n_elements_per_line == 0) {
                buf += line + "\n";
                multi_line = true;
                line = space;
            }

            line += fmt::format(" {:>{}}", token[i], width);
            ++n_elements_in_line;
        }

        // Specific for n_max_print = 1.
        if (n_max_print_elements_ == 1 && split_pos != -1) {
            if (line.size() + 5 > max_len_in_line) {
                line += ",\n  ...";
            } else {
                line += ", ...";
            }
        }

        buf += line + "]";

        return buf;
    }
    
    void set_max_print_elements(int max_print_elements) {
        assert(max_print_elements >= 0);

        n_max_print_elements_ = max_print_elements;
    }

private:
    /**
     * Split a string into multi lines.
     */
    void SplitLine(const std::string& str, std::vector<std::string>* lines) {
        assert(lines);

        lines->clear();
        std::string line;
        for (size_t i = 0; i < str.size(); ++i) {
            line += str[i];
            if (str[i] == '\n') {
                lines->push_back(line);
                line = "";
            }
        }
        if (!line.empty()) lines->push_back(line);
    }

    // The width of the screen (the maximal characters in one line).
    int scrren_width_ = 100;

    // The maximal number of elements of the array that need to be printed.
    // n_max_print_elements_ = 0 means that all elements will be printed.
    int n_max_print_elements_ = 10;
};

} // namespace cl

// Python like print function.
#define println(a) cl::Printer::GetInstance()->Print(#a, a)

#endif // BASE_FORMAT_H_
