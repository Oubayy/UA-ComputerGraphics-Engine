#include "ZBuffer.h"
#include "Figure.h"       // Figures3D, Figure, Face
#include "Projection.h"
#include "vector3d.h"
#include "easy_image.h"
#include "Line2D.h"       // Triangle, Point2D, Color
#include <limits>         // std::numeric_limits
#include <cmath>          // std::round, std::abs, std::min, std::max, lround
#include <sstream>        // std::stringstream
#include <vector>
#include <algorithm>      // std::min, std::max, std::swap
#include <iostream>       // std::cerr

using namespace std;

// Constructor: Initializes the Z-buffer
ZBuffer::ZBuffer(const int width, const int height) {
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
        this->at(i).resize(height, negInf); // Initialize with a very small value (far away for 1/z)
    }
}

// Projects 3D figures to a list of 2D triangles with depth info
Triangles ZBuffer::doProjectTriangle(const Figures3D &figures, const double &d_projection_plane) {
    Triangles triangles;

    for (const auto &figure : figures) {
        if (figure.faces.empty() || figure.points.empty()) continue;

        for (const auto &face : figure.faces) {
            if (face.point_indexes.size() < 3) continue; // Skip lines or points, only triangles

            // Assuming triangulation already happened, so all faces are triangles.
            // If not, this loop needs to handle triangulation or expect only triangles.
            // For simplicity, we assume face.point_indexes already refers to a triangle.
            const Vector3D &P1_3D_eye = figure.points[face.point_indexes[0]];
            const Vector3D &P2_3D_eye = figure.points[face.point_indexes[1]];
            const Vector3D &P3_3D_eye = figure.points[face.point_indexes[2]];

            // Project points (d_projection_plane is 'd' from projection formula x' = d*x_eye / -z_eye)
            Point2D p1_2D, p2_2D, p3_2D;
            double z1_eye, z2_eye, z3_eye; // Positive depth values

            // Projection for P1
            if (P1_3D_eye.z >= 0) { // Point at or behind eye
                p1_2D.x = (P1_3D_eye.x > 0) ? 1e6 : -1e6; p1_2D.y = (P1_3D_eye.y > 0) ? 1e6 : -1e6;
                z1_eye = std::numeric_limits<double>::infinity();
            } else {
                p1_2D.x = d_projection_plane * P1_3D_eye.x / -P1_3D_eye.z;
                p1_2D.y = d_projection_plane * P1_3D_eye.y / -P1_3D_eye.z;
                z1_eye = -P1_3D_eye.z;
            }

            // Projection for P2
            if (P2_3D_eye.z >= 0) {
                p2_2D.x = (P2_3D_eye.x > 0) ? 1e6 : -1e6; p2_2D.y = (P2_3D_eye.y > 0) ? 1e6 : -1e6;
                z2_eye = std::numeric_limits<double>::infinity();
            } else {
                p2_2D.x = d_projection_plane * P2_3D_eye.x / -P2_3D_eye.z;
                p2_2D.y = d_projection_plane * P2_3D_eye.y / -P2_3D_eye.z;
                z2_eye = -P2_3D_eye.z;
            }

            // Projection for P3
            if (P3_3D_eye.z >= 0) {
                p3_2D.x = (P3_3D_eye.x > 0) ? 1e6 : -1e6; p3_2D.y = (P3_3D_eye.y > 0) ? 1e6 : -1e6;
                z3_eye = std::numeric_limits<double>::infinity();
            } else {
                p3_2D.x = d_projection_plane * P3_3D_eye.x / -P3_3D_eye.z;
                p3_2D.y = d_projection_plane * P3_3D_eye.y / -P3_3D_eye.z;
                z3_eye = -P3_3D_eye.z;
            }

            // Only add triangle if all points are in front of the eye (or implement clipping)
            if (z1_eye != std::numeric_limits<double>::infinity() &&
                z2_eye != std::numeric_limits<double>::infinity() &&
                z3_eye != std::numeric_limits<double>::infinity()) {
                triangles.push_back({
                    p1_2D, p2_2D, p3_2D,
                    z1_eye, z2_eye, z3_eye, figure.color
                });
            }
        }
    }
    return triangles;
}


// d_proj is the projection plane distance used in doProjectTriangle (typically 1.0)
// It is NOT the image scaling factor 'd' from the calculate() function.
void ZBuffer::draw_zbuf_triangle(img::EasyImage &image, const Triangle& triangle, const double &d_proj_unused) {
    // triangle.p1.x, .p1.y etc. are ALREADY scaled and translated screen coordinates.
    // triangle.z1, .z2, .z3 are original positive eye-space depths.

    Point2D p1_screen = triangle.p1;
    Point2D p2_screen = triangle.p2;
    Point2D p3_screen = triangle.p3;

    double inv_z1_eye = 1.0 / triangle.z1;
    double inv_z2_eye = 1.0 / triangle.z2;
    double inv_z3_eye = 1.0 / triangle.z3;

    // Centroid of the screen-projected triangle
    double G_x_screen = (p1_screen.x + p2_screen.x + p3_screen.x) / 3.0;
    double G_y_screen = (p1_screen.y + p2_screen.y + p3_screen.y) / 3.0;
    // 1/z value at the centroid (interpolated in a perspective-correct manner)
    double G_inv_z_eye = (inv_z1_eye + inv_z2_eye + inv_z3_eye) / 3.0;

    // Vectors for calculating screen-space gradients of 1/z
    // Using (screen_x, screen_y, 1/z_eye)
    Vector3D vA_scr_invZ = Vector3D::point(p1_screen.x, p1_screen.y, inv_z1_eye);
    Vector3D vB_scr_invZ = Vector3D::point(p2_screen.x, p2_screen.y, inv_z2_eye);
    Vector3D vC_scr_invZ = Vector3D::point(p3_screen.x, p3_screen.y, inv_z3_eye);

    Vector3D U_scr_invZ = vB_scr_invZ - vA_scr_invZ;
    Vector3D V_scr_invZ = vC_scr_invZ - vA_scr_invZ;
    Vector3D W_scr_invZ = Vector3D::cross(U_scr_invZ, V_scr_invZ); // Normal to the plane in (screen_x, screen_y, 1/z_eye) space

    double dzdx_screen, dzdy_screen;
    // W_scr_invZ.z is proportional to the signed area of the triangle on screen.
    // If W_scr_invZ.z is zero, the triangle is degenerate (projects to a line or point).
    if (std::abs(W_scr_invZ.z) < 1e-9) { // Threshold for degenerate triangle
        // Handle as a line or skip. For now, skip to prevent division by zero.
        // A more robust solution might draw the dominant edges as lines using draw_zbuf_line.
        return;
    }

    // Gradients of 1/z_eye with respect to screen_x and screen_y
    // Plane equation: W_scr_invZ.x * x_s + W_scr_invZ.y * y_s + W_scr_invZ.z * (1/z_eye) + D_const = 0
    // So, d(1/z_eye)/dx_s = -W_scr_invZ.x / W_scr_invZ.z
    // So, d(1/z_eye)/dy_s = -W_scr_invZ.y / W_scr_invZ.z
    dzdx_screen = -W_scr_invZ.x / W_scr_invZ.z;
    dzdy_screen = -W_scr_invZ.y / W_scr_invZ.z;

    // Determine bounding box of the triangle in screen coordinates
    // Using lround and +/- 0.5 for conservative pixel bounds (similar to your original approach)
    long y_min_tri = lround(std::min({p1_screen.y, p2_screen.y, p3_screen.y}) + 0.5);
    long y_max_tri = lround(std::max({p1_screen.y, p2_screen.y, p3_screen.y}) - 0.5);
    long x_min_tri_overall = lround(std::min({p1_screen.x, p2_screen.x, p3_screen.x}) + 0.5); // For early exit
    long x_max_tri_overall = lround(std::max({p1_screen.x, p2_screen.x, p3_screen.x}) - 0.5); // For early exit

    int img_width = image.get_width();
    int img_height = image.get_height();

    // Clip bounding box to image dimensions
    y_min_tri = std::max(0L, y_min_tri);
    y_max_tri = std::min((long)img_height - 1, y_max_tri);
    x_min_tri_overall = std::max(0L, x_min_tri_overall);
    x_max_tri_overall = std::min((long)img_width - 1, x_max_tri_overall);


    // Scanline rasterization
    for (long y_current = y_min_tri; y_current <= y_max_tri; ++y_current) {
        // Calculate x_L and x_R for the current scanline y_current
        // Intersect scanline y_current with triangle edges (p1p2, p2p3, p3p1)
        double x_intersect[3];
        int intersect_count = 0;

        // Edge p1-p2
        if ((p1_screen.y <= y_current && p2_screen.y > y_current) || (p2_screen.y <= y_current && p1_screen.y > y_current)) {
            if (std::abs(p2_screen.y - p1_screen.y) > 1e-9) { // Avoid division by zero for horizontal edge
                x_intersect[intersect_count++] = p1_screen.x + (p2_screen.x - p1_screen.x) * (y_current - p1_screen.y) / (p2_screen.y - p1_screen.y);
            }
        }
        // Edge p2-p3
        if ((p2_screen.y <= y_current && p3_screen.y > y_current) || (p3_screen.y <= y_current && p2_screen.y > y_current)) {
             if (std::abs(p3_screen.y - p2_screen.y) > 1e-9) {
                x_intersect[intersect_count++] = p2_screen.x + (p3_screen.x - p2_screen.x) * (y_current - p2_screen.y) / (p3_screen.y - p2_screen.y);
            }
        }
        // Edge p3-p1
        if ((p3_screen.y <= y_current && p1_screen.y > y_current) || (p1_screen.y <= y_current && p3_screen.y > y_current)) {
            if (std::abs(p1_screen.y - p3_screen.y) > 1e-9) {
                x_intersect[intersect_count++] = p3_screen.x + (p1_screen.x - p3_screen.x) * (y_current - p3_screen.y) / (p1_screen.y - p3_screen.y);
            }
        }

        if (intersect_count < 2) continue; // Not a valid span for this scanline

        long x_L = lround(std::min(x_intersect[0], x_intersect[1]) + 0.5);
        long x_R = lround(std::max(x_intersect[0], x_intersect[1]) - 0.5);
        if (intersect_count == 3) { // Should not happen with convex triangles unless an edge is horizontal
            x_L = lround(std::min({x_intersect[0], x_intersect[1], x_intersect[2]}) + 0.5);
            x_R = lround(std::max({x_intersect[0], x_intersect[1], x_intersect[2]}) - 0.5);
        }


        // Clip x_L and x_R to image bounds
        x_L = std::max(x_min_tri_overall, x_L); // Use overall min/max for x to avoid going too far out
        x_R = std::min(x_max_tri_overall, x_R);


        for (long x_current = x_L; x_current <= x_R; ++x_current) {
            // Interpolate 1/z for the current pixel (x_current, y_current)
            double current_inv_z_eye = G_inv_z_eye + (x_current - G_x_screen) * dzdx_screen + (y_current - G_y_screen) * dzdy_screen;

            // Z-buffer stores 1/z values; larger 1/z means closer
            if (x_current >=0 && x_current < img_width && y_current >=0 && y_current < img_height) { // Double check bounds before access
                if (current_inv_z_eye > this->at(x_current)[y_current]) {
                    this->at(x_current)[y_current] = current_inv_z_eye;
                    image(x_current, y_current) = img::Color(
                        static_cast<uint8_t>(round(triangle.color.red * 255)),
                        static_cast<uint8_t>(round(triangle.color.green * 255)),
                        static_cast<uint8_t>(round(triangle.color.blue * 255))
                    );
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

    // Handle potential infinite Z values and ensure z > 0 for 1/z
    double inv_z0 = (line.z1 > 1e-9 && !std::isinf(line.z1)) ? 1.0 / line.z1 : std::numeric_limits<double>::lowest();
    double inv_z1 = (line.z2 > 1e-9 && !std::isinf(line.z2)) ? 1.0 / line.z2 : std::numeric_limits<double>::lowest();


    // Convert color
    img::Color color = img::Color(static_cast<uint8_t>(round(line.color.red * 255)),
                                  static_cast<uint8_t>(round(line.color.green * 255)),
                                  static_cast<uint8_t>(round(line.color.blue * 255)));

    int x0_i = round(x0_d);
    int y0_i = round(y0_d);
    int x1_i = round(x1_d);
    int y1_i = round(y1_d);

    if ((x0_i > x1_i) || ((x0_i == x1_i) && (y0_i > y1_i))) {
        std::swap(x0_d, x1_d);
        std::swap(y0_d, y1_d);
        std::swap(inv_z0, inv_z1);
        x0_i = round(x0_d);
        y0_i = round(y0_d);
        x1_i = round(x1_d);
        y1_i = round(y1_d);
    }

    int imgWidth = image.get_width();
    int imgHeight = image.get_height();

    int dx_i = x1_i - x0_i;
    int dy_i = y1_i - y0_i;

    if (dx_i == 0) {
        if (dy_i == 0) {
             if (x0_i >= 0 && x0_i < imgWidth && y0_i >= 0 && y0_i < imgHeight) {
                 if (inv_z0 > this->at(x0_i)[y0_i]) {
                     image(x0_i, y0_i) = color;
                     this->at(x0_i)[y0_i] = inv_z0;
                 }
             }
        } else {
            for (int y = y0_i; y <= y1_i; ++y) {
                 if (x0_i >= 0 && x0_i < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dy_i == 0) ? 0.0 : static_cast<double>(y - y0_i) / dy_i;
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                    if (current_inv_z > this->at(x0_i)[y]) {
                         image(x0_i, y) = color;
                         this->at(x0_i)[y] = current_inv_z;
                    }
                }
            }
        }
    } else if (dy_i == 0) {
        for (int x = x0_i; x <= x1_i; ++x) {
            if (x >= 0 && x < imgWidth && y0_i >= 0 && y0_i < imgHeight) {
                double t = (dx_i == 0) ? 0.0 : static_cast<double>(x - x0_i) / dx_i;
                double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                if (current_inv_z > this->at(x)[y0_i]) {
                     image(x, y0_i) = color;
                     this->at(x)[y0_i] = current_inv_z;
                }
            }
        }
    } else {
        double m = static_cast<double>(dy_i) / dx_i; // Use integer slope for this rasterization part

        if (-1.0 <= m && m <= 1.0) {
            for (int x = x0_i; x <= x1_i; ++x) {
                int y = round(y0_d + m * (x - x0_i)); // Use original doubles for y calculation for precision
                 if (x >= 0 && x < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dx_i == 0) ? 0.5 : static_cast<double>(x - x0_i) / dx_i; // t based on dominant axis
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                     if (current_inv_z > this->at(x)[y]) {
                         image(x, y) = color;
                         this->at(x)[y] = current_inv_z;
                     }
                }
            }
        } else { // More vertical
            // Ensure y iterates in the correct direction
            int start_y = y0_i, end_y = y1_i;
            double start_x_d = x0_d;
            if (y0_i > y1_i) { // if steep and going "up" visually (y decreases)
                 std::swap(start_y, end_y);
                 // We already swapped x0_d,y0_d,inv_z0 with x1_d,y1_d,inv_z1 if original x0_i > x1_i or (x0_i==x1_i && y0_i > y1_i)
                 // For this branch (m > 1 or m < -1), the primary iteration is over y.
                 // The initial swap logic prioritizes x. If steep, y changes more.
                 // Let's re-ensure y0_d is less than y1_d for this loop structure if m is positive.
                 // Or iterate from y0_d towards y1_d regardless.
            }

            for (int y = y0_i; ; (y0_i <= y1_i ? y++ : y--)) {
                int x = round(x0_d + (y - y0_d) / m); // Use original doubles for x calculation
                 if (x >= 0 && x < imgWidth && y >= 0 && y < imgHeight) {
                    double t = (dy_i == 0) ? 0.5 : static_cast<double>(y - y0_i) / dy_i; // t based on y-axis
                    double current_inv_z = inv_z0 + t * (inv_z1 - inv_z0);
                     if (current_inv_z > this->at(x)[y]) {
                         image(x, y) = color;
                         this->at(x)[y] = current_inv_z;
                     }
                }
                if (y == y1_i) break;
            }
        }
    }
}