#include "ZBuffer.h"
#include "Figure.h"
#include "Projection.h"
#include "vector3d.h"
#include "easy_image.h"
#include "Line2D.h"
#include <limits>       // Required for std::numeric_limits
#include <cmath>        // Required for std::round, std::sqrt, std::isinf
#include <sstream>      // Required for std::stringstream
#include <vector>
#include <algorithm>    // Required for std::min, std::max, std::swap

using namespace std;

// Constructor: Initializes the Z-buffer
ZBuffer::ZBuffer(const int width, const int height) {
    // Initialize with the smallest possible double value (~negative infinity)
    double negInf = std::numeric_limits<double>::lowest();

    if (width <= 0 || height <= 0) {
        int valid_width = std::max(0, width);
        int valid_height = std::max(0, height);
        this->resize(valid_width);
        if (valid_width > 0) {
             for (int i = 0; i < valid_width; ++i) {
                this->at(i).resize(valid_height, negInf);
            }
        }
         std::cerr << "Warning: ZBuffer created with non-positive dimensions (" << width << "x" << height << "). Size set to " << valid_width << "x" << valid_height << "." << std::endl;
        return;
    }

    this->resize(width);
    for (int i = 0; i < width; ++i) {
        this->at(i).resize(height, negInf);
    }
}

void triangulate(Figures3D &figs) {

    for (auto &figure: figs) {
        vector<Face> faces;

        for (auto &face: figure.faces) {
            auto vertexCount = face.point_indexes.size();

            if (vertexCount < 3) {
                // Invalid or a line
                faces.push_back(face);
            } else if (vertexCount == 3) {
                // The face is already a triangle
                faces.push_back(face);
            }/* else if (vertexCount == 4) {
                // The face is a rectangle/square => we will split it into 2 triangles
                faces.push_back({{
                                         face.point_indexes[0],
                                         face.point_indexes[1],
                                         face.point_indexes[2]
                                 }});
                faces.push_back({{
                                         face.point_indexes[0],
                                         face.point_indexes[2],
                                         face.point_indexes[3]
                                 }});
            }*/ else {
                // The face is an n-gon where n > 4. We'll split it into (n-2) triangles.
                for (int i = 0; i < vertexCount - 2; i++) {
                    faces.push_back({{
                                             face.point_indexes[0],
                                             face.point_indexes[i + 1],
                                             face.point_indexes[i + 2]
                                     }});
                }
            }
        }
        figure.faces = faces;
    }
}

/*
Triangles ZBuffer::doProjectTriangle(const Figures3D &figures, const double &d) {
    Triangles triangles;

    for (const auto &figure : figures) {
        for (const auto &face : figure.faces) {
            // Points A, B, C in 3D
            const Vector3D &A = figure.points[face.point_indexes[0]];
            const Vector3D &B = figure.points[face.point_indexes[1]];
            const Vector3D &C = figure.points[face.point_indexes[2]];

            // Project points A, B, C to 2D using the projection formula
            Point2D projectA;
            projectA.x = 1 * A.x / -A.z;
            projectA.y = 1 * A.y / -A.z;
            double z1 = -A.z;

            Point2D projectB;
            projectB.x = 1 * B.x / -B.z;
            projectB.y = 1 * B.y / -B.z;
            double z2 = -B.z;

            Point2D projectC;
            projectC.x = 1 * C.x / -C.z;
            projectC.y = 1 * C.y / -C.z;
            double z3 = -C.z;

            // Add the projected 2D triangle to the collection
            triangles.push_back({
                projectA, projectB, projectC,
                z1, z2, z3, figure.color
            });
        }
    }

    return triangles;
}


void ZBuffer::draw_zbuf_triangle(img::EasyImage &image, const Triangle& triangle, const double &d) {
    Vector3D A = Vector3D::point(triangle.p1.x, triangle.p1.y, triangle.z1);
    Vector3D B = Vector3D::point(triangle.p2.x, triangle.p2.y, triangle.z2);
    Vector3D C = Vector3D::point(triangle.p3.x, triangle.p3.y,  triangle.z3);

    Vector3D U = B - A;
    Vector3D V = C - A;
    Vector3D W = Vector3D::cross(U, V);

    Vector3D middel_punt = Vector3D::point(
            A.x/3 + B.x/3 + C.x/3,
            A.y/3 + B.y/3 + C.y/3,
            A.z/3 + B.z/3 + C.z/3
    );

    double k = W.x + A.x + W.y + A.y + W.z + A.z;
    double dzdx = W.x / (-d * k);
    double dzdy = W.y / (-d * k);

    int max_y = lround(std::max(std::max(A.y, B.y), C.y) - 0.5);
    int min_y = lround(std::min(std::min(A.y, B.y), C.y) + 0.5);

    for (int y = min_y; y <= max_y; y++) {
        double AB_l = numeric_limits<double>::infinity();
        double AB_r = -numeric_limits<double>::infinity();

        double AC_l = numeric_limits<double>::infinity();
        double AC_r = -numeric_limits<double>::infinity();

        double BC_l = numeric_limits<double>::infinity();
        double BC_r = -numeric_limits<double>::infinity();

        if ((y - A.y) * (y - B.y) <= 0 and A.y != B.y) {
            double x = B.x + (A.x - B.x) * (y - B.y) / (A.y - B.y);
            AB_l = AB_r = x;
        }

        if ((y - A.y) * (y - C.y) <= 0 and A.y != C.y) {
            double x = C.x + (A.x - C.x) * (y - C.y) / (A.y - C.y);
            AC_l = AC_r = x;
        }

        if ((y - B.y) * (y - C.y) <= 0 and B.y != C.y) {
            double x = C.x + (B.x - C.x) * (y - C.y) / (B.y - C.y);
            BC_l = BC_r = x;
        }

        int max_x = lround(std::max(std::max(AB_r, AC_r), BC_r) - 0.5);
        int min_x = lround(std::min(std::min(AB_l, AC_l), BC_l) + 0.5);

        for (int x = min_x; x <= max_x; x++) {
            double z = 1.0001 * middel_punt.z + (x - middel_punt.x) * dzdx + (y - middel_punt.y) * dzdy;

            if (z <= this->at(x)[y]) {
                image(x, y) = {
                    static_cast<uint8_t>(triangle.color.red * 255),
                    static_cast<uint8_t>(triangle.color.green * 255),
                    static_cast<uint8_t>(triangle.color.blue * 255)
                };
                this->at(x)[y] = z;
            }
        }
    }
}
*/

// Draws a 2D line segment with Z-buffering.
// Assumes line coordinates are already transformed and clipped.
void ZBuffer::draw_zbuf_line(img::EasyImage &image, Line2D line) {
    // Use doubles for precision during calculations
    double x0_d = line.p1.x;
    double y0_d = line.p1.y;
    double x1_d = line.p2.x;
    double y1_d = line.p2.y;

    // Handle potential infinite Z values
    double inv_z0 = (line.z1 > 0 && !std::isinf(line.z1)) ? 1.0 / line.z1 : std::numeric_limits<double>::lowest();
    double inv_z1 = (line.z2 > 0 && !std::isinf(line.z2)) ? 1.0 / line.z2 : std::numeric_limits<double>::lowest();

    // Convert color
    img::Color color = img::Color(static_cast<uint8_t>(round(line.color.red * 255)),
                                  static_cast<uint8_t>(round(line.color.green * 255)),
                                  static_cast<uint8_t>(round(line.color.blue * 255)));

    // Get integer coordinates by rounding for pixel access and range calculations
    // *** CHANGED lround to round ***
    int x0_i = round(x0_d);
    int y0_i = round(y0_d);
    int x1_i = round(x1_d);
    int y1_i = round(y1_d);

    // Ensure line goes left-to-right or bottom-to-top
    if ((x0_i > x1_i) || ((x0_i == x1_i) && (y0_i > y1_i))) {
        std::swap(x0_d, x1_d);
        std::swap(y0_d, y1_d);
        std::swap(inv_z0, inv_z1);
        // Update integer versions after swap
        // *** CHANGED lround to round ***
        x0_i = round(x0_d);
        y0_i = round(y0_d);
        x1_i = round(x1_d);
        y1_i = round(y1_d);
    }

    int imgWidth = image.get_width();
    int imgHeight = image.get_height();

    // Calculate integer differences
    int dx_i = x1_i - x0_i;
    int dy_i = y1_i - y0_i;

    // --- Rasterization Logic ---
    if (dx_i == 0) { // Vertical line
        if (dy_i == 0) { // Single point
             if (x0_i >= 0 && x0_i < imgWidth && y0_i >= 0 && y0_i < imgHeight) {
                 if (inv_z0 > this->at(x0_i)[y0_i]) {
                     image(x0_i, y0_i) = color;
                     this->at(x0_i)[y0_i] = inv_z0;
                 }
             }
        } else { // dy_i > 0
            for (int y = y0_i; y <= y1_i; ++y) {
                 if (x0_i >= 0 && x0_i < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dy_i == 0) ? 0.0 : (double)(y - y0_i) / dy_i;
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                    if (current_inv_z > this->at(x0_i)[y]) {
                         image(x0_i, y) = color;
                         this->at(x0_i)[y] = current_inv_z;
                    }
                }
            }
        }
    } else if (dy_i == 0) { // Horizontal line (dx_i > 0)
        for (int x = x0_i; x <= x1_i; ++x) {
            if (x >= 0 && x < imgWidth && y0_i >= 0 && y0_i < imgHeight) {
                double t = (dx_i == 0) ? 0.0 : (double)(x - x0_i) / dx_i;
                double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                if (current_inv_z > this->at(x)[y0_i]) {
                     image(x, y0_i) = color;
                     this->at(x)[y0_i] = current_inv_z;
                }
            }
        }
    } else { // Diagonal line
        double m = (y1_d - y0_d) / (x1_d - x0_d); // Use double slope

        if (-1.0 <= m && m <= 1.0) { // More horizontal
            for (int x = x0_i; x <= x1_i; ++x) {
                // *** CHANGED lround to round ***
                int y = round(y0_d + m * (x - x0_i));
                 if (x >= 0 && x < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dx_i == 0) ? 0.0 : (double)(x - x0_i) / dx_i;
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                     if (current_inv_z > this->at(x)[y]) {
                         image(x, y) = color;
                         this->at(x)[y] = current_inv_z;
                     }
                }
            }
        } else if (m > 1.0) { // Steep positive slope
            for (int y = y0_i; y <= y1_i; ++y) {
                // *** CHANGED lround to round ***
                int x = round(x0_d + (y - y0_i) / m);
                 if (x >= 0 && x < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dy_i == 0) ? 0.0 : (double)(y - y0_i) / dy_i;
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                     if (current_inv_z > this->at(x)[y]) {
                         image(x, y) = color;
                         this->at(x)[y] = current_inv_z;
                     }
                }
            }
        } else { // m < -1.0, Steep negative slope
             for (int y = y0_i; y >= y1_i; --y) {
                 // *** CHANGED lround to round ***
                 int x = round(x0_d + (y - y0_i) / m);
                 if (x >= 0 && x < imgWidth && y >= 0 && y < imgHeight) {
                     double t = (dy_i == 0) ? 0.0 : (double)(y - y0_i) / dy_i;
                     double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                     if (current_inv_z > this->at(x)[y]) {
                          image(x, y) = color;
                          this->at(x)[y] = current_inv_z;
                     }
                 }
             }
        }
    }
}
