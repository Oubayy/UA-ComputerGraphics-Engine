/**
 * @file Projection.cpp
 * @brief Implementation of 3D to 2D projection functions.
 */

#include "Projection.h"
#include <iostream> // For debug

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

    // Basic check for points at or behind the projection plane (z >= 0 in eye coords)
    // A more robust implementation would handle clipping.
    if (point.z >= 0) {
         // Handle point at or behind the eye - return a point far away or handle as error
         // std::cerr << "Warning: Point at or behind the eye (z=" << point.z << "). Projection is undefined or infinite." << std::endl;
         // Returning a point at "infinity" or a default value might be options.
         // For simplicity, let's project to a large coordinate, but proper clipping is needed.
         projectedPoint.x = (point.x > 0) ? 1e6 : -1e6;
         projectedPoint.y = (point.y > 0) ? 1e6 : -1e6;

    } else {
        // Standard perspective projection formula
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
    const double d = 1.0; // Distance to the projection plane (often set to 1)

    for (const auto &fig : figs) {
        for (const auto &face : fig.faces) {
            // Handle faces defined by lines (e.g., from LineDrawing or 3DLSystem)
            if (face.point_indexes.size() == 2) {
                 int p1_idx = face.point_indexes[0];
                 int p2_idx = face.point_indexes[1];

                 // Check index bounds
                 if (p1_idx < 0 || p1_idx >= fig.points.size() || p2_idx < 0 || p2_idx >= fig.points.size()) {
                     std::cerr << "Warning: Invalid point index in face. Skipping line." << std::endl;
                     continue;
                 }

                 const Vector3D& P1_3D = fig.points[p1_idx];
                 const Vector3D& P2_3D = fig.points[p2_idx];

                 // Project points
                 Point2D p1_2D = doProjection(P1_3D, d);
                 Point2D p2_2D = doProjection(P2_3D, d);

                 // Store depth (distance from eye, which is -z after eye transformation)
                 // Ensure z is not zero or positive before calculating depth
                 double z1 = (P1_3D.z < 0) ? -P1_3D.z : std::numeric_limits<double>::infinity(); // Use positive distance
                 double z2 = (P2_3D.z < 0) ? -P2_3D.z : std::numeric_limits<double>::infinity();

                 // Add line if depths are valid (or handle clipping properly)
                 if (z1 != std::numeric_limits<double>::infinity() || z2 != std::numeric_limits<double>::infinity()) {
                     lines.push_back(Line2D(p1_2D, p2_2D, z1, z2, fig.color));
                 } else {
                     // Both points might be behind the camera, skip line (or clip)
                     //std::cerr << "Skipping line with both points behind the eye." << std::endl;
                 }

            }
            // Handle faces defined by polygons (triangles, quads, etc.)
            else if (face.point_indexes.size() > 2) {
                for (size_t i = 0; i < face.point_indexes.size(); ++i) {
                    int p1_idx = face.point_indexes[i];
                    int p2_idx = face.point_indexes[(i + 1) % face.point_indexes.size()]; // Connect last to first

                    // Check index bounds
                    if (p1_idx < 0 || p1_idx >= fig.points.size() || p2_idx < 0 || p2_idx >= fig.points.size()) {
                        std::cerr << "Warning: Invalid point index in face. Skipping line segment." << std::endl;
                        continue;
                    }

                    const Vector3D& P1_3D = fig.points[p1_idx];
                    const Vector3D& P2_3D = fig.points[p2_idx];

                    // Project points
                    Point2D p1_2D = doProjection(P1_3D, d);
                    Point2D p2_2D = doProjection(P2_3D, d);

                    // Store depth (distance from eye, which is -z after eye transformation)
                    double z1 = (P1_3D.z < 0) ? -P1_3D.z : std::numeric_limits<double>::infinity();
                    double z2 = (P2_3D.z < 0) ? -P2_3D.z : std::numeric_limits<double>::infinity();

                     // Add line if depths are valid (or handle clipping properly)
                     if (z1 != std::numeric_limits<double>::infinity() || z2 != std::numeric_limits<double>::infinity()) {
                         lines.push_back(Line2D(p1_2D, p2_2D, z1, z2, fig.color));
                     } else {
                         // Both points might be behind the camera, skip line (or clip)
                         //std::cerr << "Skipping line segment with both points behind the eye." << std::endl;
                     }
                }
            }
             else {
                 // Handle cases with less than 2 points if necessary
                 //std::cerr << "Warning: Face with less than 2 points encountered." << std::endl;
             }
        }
    }
    return lines;
}
