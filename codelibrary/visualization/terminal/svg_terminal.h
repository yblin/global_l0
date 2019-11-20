//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef VISUALIZATION_TERMINAL_SVG_TERMINAL_H_
#define VISUALIZATION_TERMINAL_SVG_TERMINAL_H_

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include "codelibrary/base/format.h"
#include "codelibrary/geometry/kernel/multi_polygon_2d.h"
#include "codelibrary/string/string_split.h"
#include "codelibrary/visualization/terminal/terminal.h"

namespace cl {

/**
 * SVG terminal
 *
 * Visualizing the figures on the Scalable Array Graphics File.
 */
class SVGTerminal : public Terminal {
public:
    SVGTerminal(int width = 640, int height = 480)
        : Terminal(width, height) {
        font_.set_name("Times");
        Initialize();
    }

    SVGTerminal(const SVGTerminal&) = delete;

    SVGTerminal& operator=(const SVGTerminal&) = delete;

    /**
     * Clear the content of SVG.
     */
    virtual void clear() override {
        content_.clear();
    }

    /**
     * Draw circle.
     */
    virtual void DrawCircle(double x, double y, double r) override {
        content_ += ElementStart("circle") +
                    Attribute("cx", x) +
                    Attribute("cy", height_ - y) +
                    Attribute("r", r) +
                    PenAttribute() +
                    EmptyElementEnd();
    }

    /**
     * Draw rectangle.
     */
    virtual void DrawRectangle(double x, double y, double w, double h) {
        content_ += ElementStart("rect") +
                    PenAttribute() +
                    Attribute("x", x) +
                    Attribute("y", height_ - y) +
                    Attribute("width", w) +
                    Attribute("height", h) +
                    EmptyElementEnd();
    }

    /**
     * Draw line (x1, y1)->(x2, y2).
     */
    virtual void DrawLine(double x1, double y1, double x2, double y2) override {
        content_ += ElementStart("line") +
                    Attribute("x1", x1) +
                    Attribute("y1", height_ - y1) +
                    Attribute("x2", x2) +
                    Attribute("y2", height_ - y2) +
                    LineAttribute() +
                    EmptyElementEnd();
    }

    /**
     * Draw polyline.
     */
    virtual void DrawPolyline(const Array<RPoint2D>& polyline) override {
        std::string points;
        for (const RPoint2D& p : polyline) {
            points += fmt::format("{:g}", p.x) + "," +
                      fmt::format("{:g}", height_ - p.y) + " ";
        }

        content_ += ElementStart("polyline") +
                    Attribute("fill", "none") +
                    LineAttribute() +
                    Attribute("points", points) +
                    EmptyElementEnd();
    }

    /**
     * Draw multiple-polygon.
     */
    virtual void DrawPolygon(const RMultiPolygon2D& polygon) override {
        std::string path;
        for (int i = 0; i < polygon.n_boundaries(); ++i) {
            path += PathAttribute(polygon.boundaries()[i].polygon.vertices());
            path += "Z ";
        }
        path = ElementStart("path") + Attribute("d", path) + EmptyElementEnd();

        content_ += ElementStart("g") +
                    PenAttribute() +
                    Attribute("fill-rule", "evenodd") +
                    + ">\n" +
                    path +
                    ElementEnd("g");
    }

    /**
     * Draw triangle.
     */
    virtual void DrawTriangle(const RPoint2D& p1, const RPoint2D& p2,
                              const RPoint2D& p3) override {
        std::string points;
        points += fmt::format("{:g}", p1.x) + "," +
                  fmt::format("{:g}", height_ - p1.y) + " ";
        points += fmt::format("{:g}", p2.x) + "," +
                  fmt::format("{:g}", height_ - p2.y) + " ";
        points += fmt::format("{:g}", p3.x) + "," +
                  fmt::format("{:g}", height_ - p3.y) + " ";

        content_ += ElementStart("polygon") +
                    PenAttribute() +
                    Attribute("fill-rule", "evenodd") +
                    Attribute("points", points) +
                    EmptyElementEnd();
    }

    /**
     * Draw text at position (x, y).
     */
    virtual void DrawText(double x, double y,
                          const std::string& text) override {
        DrawText(x, y, false, text);
    }

    /**
     * Draw vertical text at position (x, y).
     */
    virtual void DrawVerticalText(double x, double y,
                                  const std::string& text) override {
        DrawText(x, y, true, text);
    }

    /**
     * Save to the SVG file.
     */
    virtual void SaveToFile(const std::string& file) const override {
        std::ofstream fout(file);
        fout << head_ << content_ << ElementEnd("svg");
    }

    /**
     * Resize the terminal.
     */
    virtual void Resize(int height, int width) override {
        assert(height > 0 && width > 0);

        if (height_ == height && width_ == width) return;

        height_ = height;
        width_  = width;
        Initialize();
    }

private:
    /**
     * Convert a value into XML attribute.
     */
    template <typename T>
    static const std::string Attribute(const std::string& attribute_name,
                                       const T& value) {
        std::stringstream ss;
        ss << attribute_name << "=\"" << value << "\" ";
        return ss.str();
    }
    static const std::string Attribute(const std::string& attribute_name,
                                       const RGB32Color& color) {
        std::stringstream ss;
        ss << attribute_name << "=\"" << ColorToString(color) << "\" ";
        return ss.str();
    }

    /**
     * Return a string that represents starting of a SVG element.
     */
    static const std::string ElementStart(const std::string& element_name) {
        return "\t<" + element_name + " ";
    }

    /**
     * Return a string that represents ending of a SVG element.
     */
    static const std::string ElementEnd(const std::string& element_name) {
        return "\t</" + element_name + ">\n";
    }

    /**
     * Return a string that represents empty ending of a SVG element.
     */
    static const std::string EmptyElementEnd() {
        return "/>\n";
    }

    /**
     * Convert the RGB Color into string.
     */
    static const std::string ColorToString(const RGB32Color& c) {
        if (c.alpha() == 255) {
            return "rgb(" + std::to_string(c.red())   + "," +
                            std::to_string(c.green()) + "," +
                            std::to_string(c.blue())  + ")";
        }

        return "rgba(" + std::to_string(c.red())   + "," +
                         std::to_string(c.green()) + "," +
                         std::to_string(c.blue())  + "," +
                         std::to_string(c.alpha() / 255.0) + ")";
    }

    /**
     * Initialize the SVG terminal.
     */
    void Initialize() {
        // Initialize the SVG header.
        head_ = "<?xml " + Attribute("version", "1.0") + "?>\n" +
                "<svg " +
                Attribute("xmlns", "http://www.w3.org/2000/svg") +
                Attribute("xmlns:xlink", "http://www.w3.org/1999/xlink") +
                Attribute("width",  width_) +
                Attribute("height", height_) +
                ">\n\n";
    }

    /**
     * Return pen attribute.
     */
    std::string PenAttribute() const {
        if (pen_.line_width == 0.0) {
            if (pen_.fill_style == Pen::NOT_FILL)
                return Attribute("fill", "none");

            return Attribute("fill", pen_.fill_color);
        }

        if (pen_.fill_style == Pen::NOT_FILL) {
            return LineAttribute() + Attribute("fill", "none");
        }

        return LineAttribute() + Attribute("fill", pen_.fill_color);
    }

    /**
     * Get the pen attribute for line drawing.
     */
    std::string LineAttribute() const {
        std::string dash;
        switch (pen_.line_style) {
        case Pen::SOLID_LINE:
            break;
        case Pen::DASH_LINE:
            dash = Attribute("stroke-dasharray", std::string("10,10"));
            break;
        case Pen::DOT_LINE:
            dash = Attribute("stroke-dasharray", std::string("2,2"));
            break;
        case Pen::DASH_DOT_LINE:
            dash = Attribute("stroke-dasharray", std::string("2,10,2"));
            break;
        }

        return Attribute("stroke", pen_.color) +
               Attribute("stroke-width", pen_.line_width) +
               dash;
    }

    /**
     * DrawText at position(x, y) with the vertical mode.
     */
    void DrawText(double x, double y, bool is_vertical,
                  const std::string& text) {
        y -= font_.size();

        std::string aligment;
        switch (font_.aligment()) {
        case Font::START:
            aligment = "start";
            break;
        case Font::END:
            aligment = "end";
            break;
        case Font::MIDDLE:
            aligment = "middle";
            break;
        }

        std::string font_weight = font_.bold() ? "bold" : "normal";

        std::string rotate = fmt::format("rotate(-90, {:g} {:g})",
                                         x, height_ - y);

        content_ += ElementStart("text") +
                    Attribute("x", x) +
                    Attribute("y", height_ - y) +
                    (is_vertical ? Attribute("transform", rotate)
                                 : "") +
                    Attribute("font-size", font_.size()) +
                    Attribute("text-anchor", aligment) +
                    Attribute("font-weight", font_weight) +
                    Attribute("font-family", font_.name()) + ">";
        content_ += ElementStart("tspan") +
                    Attribute("x", x) +
                    Attribute("y", height_ - y) +
                    ">" +
                    text +
                    ElementEnd("tspan");
        content_ += ElementEnd("text");
    }

    /**
     * Compute a close SVG path element from the given vertices list.
     */
    std::string ClosePath(const Array<RPoint2D>& points) const {
        std::string attribute = PathAttribute(points) + "Z";
        return ElementStart("path") + Attribute("d", attribute) +
               EmptyElementEnd();
    }

    /**
     * Compute a SVG path attribute from the given vertices list.
     */
    std::string PathAttribute(const Array<RPoint2D>& points) const {
        if (points.empty()) return "";

        std::string attribute = "M ";
        attribute += std::to_string(points.front().x) + " ";
        attribute += std::to_string(height_ - points.front().y) + " ";
        attribute += "L ";
        for (int i = 1; i < points.size(); ++i) {
            attribute += std::to_string(points[i].x) + " " +
                         std::to_string(height_ - points[i].y) + " ";
        }

        return attribute;
    }

    std::string head_;    // SVG head.
    std::string content_; // SVG content.
};

} // namespace cl

#endif // VISUALIZATION_TERMINAL_SVG_TERMINAL_H_
