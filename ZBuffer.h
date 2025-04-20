#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include "easy_image.h"
#include "Figure.h"
#include "vector3d.h"
#include "Line2D.h"
#include <vector>

// void triangulate(Figures3D &figs);

class ZBuffer: public std::vector<std::vector<double>> {
public:
    /**
     * @brief Constructor: Creates a Z-Buffer of the correct size and initializes all fields.
     * Initializes buffer values to the lowest possible double value (~negative infinity)
     * to correctly compare 1/z values (larger 1/z means closer).
     * @param width The width of the Z-buffer (image width).
     * @param height The height of the Z-buffer (image height).
     */
    ZBuffer(const int width, const int height);

    /**
     * @brief Draws a 2D line segment onto the image using Z-buffering.
     * Performs Z-interpolation using 1/z and updates the buffer and image
     * if the current line segment is closer than what's already buffered.
     * This version aims for pixel rasterization consistency with EasyImage::draw_line.
     * @param image The EasyImage object to draw on.
     * @param line The Line2D object containing endpoints (p1, p2),
     * depth values (z1, z2), and color. Assumes coordinates are
     * already scaled and translated, and the line is clipped.
     */
    void draw_zbuf_line(img::EasyImage &image, Line2D line);

    // static Triangles doProjectTriangle(const Figures3D &figures, const double &d);
    // void draw_zbuf_triangle(img::EasyImage &image, const Triangle& triangle, const double &d);
};


#endif //ENGINE_ZBUFFER_H