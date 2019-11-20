//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_FONT_H_
#define VISUALIZATION_FONT_H_

#include <string>

namespace cl {

/** 
 * Generic font class.
 */
class Font {
public:
    /**
     * The Alignment property is used to align (start-, middle- or
     * end-alignment) a string of text relative to a given point.
     */
    enum Alignment {
        START  = 0,
        MIDDLE = 1,
        END    = 2
    };

    /**
     * A default empty font.
     */
    Font() = default;

    /**
     * Construct a font object with the specified name, size, weight and italic
     * settings.
     */
    Font(const std::string& name, double size = 12.0, bool bold = false,
         bool italic = false, Alignment alignment = START)
        : name_(name),
          size_(size),
          bold_(bold),
          italic_(italic),
          alignment_(alignment) {}

    const std::string& name() const {
        return name_;
    }

    void set_name(const std::string& name) {
        name_ = name;
    }

    double size() const {
        return size_;
    }

    void set_size(double size) {
        size_ = size;
    }

    bool bold() const {
        return bold_;
    }

    void set_bold(bool bold) {
        bold_ = bold;
    }

    bool italic() const {
        return italic_;
    }

    void set_italic(bool italic) {
        italic_ = italic;
    }

    Alignment aligment() const {
        return alignment_;
    }

    void set_aligment(Alignment alignment) {
        alignment_ = alignment;
    }

private:
    std::string name_;            // Name of font.
    double size_ = 0.0;           // Size of font.
    bool bold_ = false;           // Check if the font is bold.
    bool italic_ = false;         // Check if the font is italic.
    Alignment alignment_ = START; // Alignment of text.
};

} // namespace cl

#endif // VISUALIZATION_FONT_H_
