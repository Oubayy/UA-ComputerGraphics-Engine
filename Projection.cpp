/**
 * @file Projection.cpp
 * @brief Implementation of 3D to 2D projection functions.
 */

#include "Projection.h"
#include <iostream>
#include <limits>

/**
 * @brief Projects a single 3D point onto the 2D projection plane.
 *
 * Uses perspective projection with the projection plane at distance d.
 * Assumes the point is in eye coordinates (eye at origin, looking down -Z axis).
 *
 * @param point The 3D point (in eye coordinates) to project.
 * @param d The distance from the eye (origin) to the projection plane.
 * @return The projected 2D point.
 */
Point2D doProjection(const Vector3D &point, const double d) {
    Point2D projectedPoint;
    if (point.z >= 0) {
        projectedPoint.x = (point.x > 0) ? 1e6 : -1e6;
        projectedPoint.y = (point.y > 0) ? 1e6 : -1e6;
    } else {
        projectedPoint.x = d * point.x / -point.z;
        projectedPoint.y = d * point.y / -point.z;
    }
    return projectedPoint;
}

/**
 * @brief Projects a list of 3D figures into a list of 2D lines.
 *
 * Iterates through each figure and its faces, projecting the edges (lines)
 * onto the 2D plane. Assumes figures are in eye coordinates.
 * Stores the projected 2D coordinates and the original depth (-z) for Z-buffering.
 *
 * @param figs The list of 3D figures (in eye coordinates) to project.
 * @return A list of Line2D objects representing the projected lines with depth info.
 */
Lines2D doProjection(const Figures3D &figs) {
    Lines2D lines;
    const double d = 1.0;

    for (const auto &fig : figs) {
        // The 'fig.color' member is the old single color, used for non-lit wireframes.
        // For lit rendering, material properties are used by ZBuffer::doProjectTriangle.
        Color current_fig_line_color = fig.color; // Use the fallback/old color member

        for (const auto &face : fig.faces) {
            if (face.point_indexes.size() == 2) { // Line segment
                 int p1_idx = face.point_indexes[0];
                 int p2_idx = face.point_indexes[1];

                 if (p1_idx < 0 || static_cast<size_t>(p1_idx) >= fig.points.size() ||
                     p2_idx < 0 || static_cast<size_t>(p2_idx) >= fig.points.size()) {
                     std::cerr << "Warning: Invalid point index in line face. Skipping line." << std::endl;
                     continue;
                 }

                 const Vector3D& P1_3D = fig.points[p1_idx];
                 const Vector3D& P2_3D = fig.points[p2_idx];

                 Point2D p1_2D = doProjection(P1_3D, d);
                 Point2D p2_2D = doProjection(P2_3D, d);

                 double z1 = (P1_3D.z < 0) ? -P1_3D.z : std::numeric_limits<double>::infinity();
                 double z2 = (P2_3D.z < 0) ? -P2_3D.z : std::numeric_limits<double>::infinity();

                 if (z1 != std::numeric_limits<double>::infinity() || z2 != std::numeric_limits<double>::infinity()) {
                     lines.push_back(Line2D(p1_2D, p2_2D, z1, z2, current_fig_line_color));
                 }
            }
            else if (face.point_indexes.size() > 2) { // Polygon edges
                for (size_t i = 0; i < face.point_indexes.size(); ++i) {
                    int p1_idx = face.point_indexes[i];
                    int p2_idx = face.point_indexes[(i + 1) % face.point_indexes.size()];

                    if (p1_idx < 0 || static_cast<size_t>(p1_idx) >= fig.points.size() ||
                        p2_idx < 0 || static_cast<size_t>(p2_idx) >= fig.points.size()) {
                        std::cerr << "Warning: Invalid point index in polygon face. Skipping line segment." << std::endl;
                        continue;
                    }

                    const Vector3D& P1_3D = fig.points[p1_idx];
                    const Vector3D& P2_3D = fig.points[p2_idx];

                    Point2D p1_2D = doProjection(P1_3D, d);
                    Point2D p2_2D = doProjection(P2_3D, d);

                    double z1 = (P1_3D.z < 0) ? -P1_3D.z : std::numeric_limits<double>::infinity();
                    double z2 = (P2_3D.z < 0) ? -P2_3D.z : std::numeric_limits<double>::infinity();

                     if (z1 != std::numeric_limits<double>::infinity() || z2 != std::numeric_limits<double>::infinity()) {
                         lines.push_back(Line2D(p1_2D, p2_2D, z1, z2, current_fig_line_color));
                     }
                }
            }
        }
    }
    return lines;
}
