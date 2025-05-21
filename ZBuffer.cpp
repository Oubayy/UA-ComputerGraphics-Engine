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

/*
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

            if (z > this->at(x)[y]) {
            //if (z <= this->at(x)[y]) {
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

void ZBuffer::draw_zbuf_triangle(img::EasyImage &image, const Triangle& triangle_details, const double &d_projection_plane_distance_unused /* This parameter is not used */) {
    Point2D pA_screen = triangle_details.p1; // Already scaled and translated screen coordinates
    Point2D pB_screen = triangle_details.p2;
    Point2D pC_screen = triangle_details.p3;

    // Original positive Z-depths from the eye
    double depth_A = triangle_details.z1;
    double depth_B = triangle_details.z2;
    double depth_C = triangle_details.z3;

    // Calculate 1/depth for each vertex.
    double inv_depth_A = (depth_A > 1e-9) ? 1.0 / depth_A : std::numeric_limits<double>::lowest();
    double inv_depth_B = (depth_B > 1e-9) ? 1.0 / depth_B : std::numeric_limits<double>::lowest();
    double inv_depth_C = (depth_C > 1e-9) ? 1.0 / depth_C : std::numeric_limits<double>::lowest();

    // If all points are effectively at/behind camera (invalid depths), skip.
    if (inv_depth_A == std::numeric_limits<double>::lowest() &&
        inv_depth_B == std::numeric_limits<double>::lowest() &&
        inv_depth_C == std::numeric_limits<double>::lowest()) {
        return;
    }

    // Calculate the screen-space bounding box of the triangle
    int min_x_bb = lround(std::min({pA_screen.x, pB_screen.x, pC_screen.x}));
    int max_x_bb = lround(std::max({pA_screen.x, pB_screen.x, pC_screen.x}));
    int min_y_bb = lround(std::min({pA_screen.y, pB_screen.y, pC_screen.y}));
    int max_y_bb = lround(std::max({pA_screen.y, pB_screen.y, pC_screen.y}));

    // Clip bounding box to image dimensions
    min_x_bb = std::max(0, min_x_bb);
    max_x_bb = std::min((int)image.get_width() - 1, max_x_bb);
    min_y_bb = std::max(0, min_y_bb);
    max_y_bb = std::min((int)image.get_height() - 1, max_y_bb);

    // Denominator for barycentric coordinates (2 * area of triangle)
    // Using pA, pB, pC for screen coordinates (triangle_details.p1 etc.)
    double bary_denom = (pB_screen.y - pC_screen.y) * (pA_screen.x - pC_screen.x) +
                        (pC_screen.x - pB_screen.x) * (pA_screen.y - pC_screen.y);

    if (std::abs(bary_denom) < 1e-9) { // Degenerate triangle on screen
        return;
    }

    for (int y_pixel = min_y_bb; y_pixel <= max_y_bb; ++y_pixel) {
        for (int x_pixel = min_x_bb; x_pixel <= max_x_bb; ++x_pixel) {
            // Calculate barycentric coordinates for the center of the current pixel
            // It's often fine to use (x_pixel, y_pixel) directly if vertices were rounded,
            // but using pixel centers (x_pixel + 0.5, y_pixel + 0.5) is more accurate.
            // Let's use integer pixel coordinates for simplicity matching the example.
            double px_test = static_cast<double>(x_pixel); // or x_pixel + 0.5;
            double py_test = static_cast<double>(y_pixel); // or y_pixel + 0.5;

            double lambda1 = ((pB_screen.y - pC_screen.y) * (px_test - pC_screen.x) +
                              (pC_screen.x - pB_screen.x) * (py_test - pC_screen.y)) / bary_denom;
            double lambda2 = ((pC_screen.y - pA_screen.y) * (px_test - pC_screen.x) +
                              (pA_screen.x - pC_screen.x) * (py_test - pC_screen.y)) / bary_denom;
            double lambda3 = 1.0 - lambda1 - lambda2;

            // Check if pixel is inside or on the edge of the triangle
            // A small epsilon can help with floating point inaccuracies at edges
            // double epsilon = 1e-7; // Adjust if needed
            // if (lambda1 >= -epsilon && lambda1 <= 1.0 + epsilon &&
            //     lambda2 >= -epsilon && lambda2 <= 1.0 + epsilon &&
            //     lambda3 >= -epsilon && lambda3 <= 1.0 + epsilon)
            // A simpler check (pixel center strictly inside or on edge):
            if (lambda1 >= 0.0 && lambda1 <= 1.0 &&
                lambda2 >= 0.0 && lambda2 <= 1.0 &&
                lambda3 >= 0.0 && lambda3 <= 1.0)
            {
                // Interpolate 1/depth using barycentric coordinates
                double current_inv_depth = inv_depth_A * lambda1 +
                                           inv_depth_B * lambda2 +
                                           inv_depth_C * lambda3;

                // Optional: Small bias for Z-fighting (consistent with previous attempt)
                // current_inv_depth *= 1.0001; // Re-evaluate if this is needed/helpful

                // Z-Buffer check: larger 1/depth value means closer.
                // this->at(x_pixel)[y_pixel] is ZBuffer's internal storage.
                if (x_pixel >= 0 && x_pixel < image.get_width() && y_pixel >=0 && y_pixel < image.get_height()){ // Redundant due to BB clip, but safe
                    if (current_inv_depth > this->at(x_pixel)[y_pixel]) {
                        image(x_pixel, y_pixel) = img::Color(
                            static_cast<uint8_t>(round(triangle_details.color.red * 255)),
                            static_cast<uint8_t>(round(triangle_details.color.green * 255)),
                            static_cast<uint8_t>(round(triangle_details.color.blue * 255))
                        );
                        this->at(x_pixel)[y_pixel] = current_inv_depth;
                    }
                }
            }
        }
    }
}


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
