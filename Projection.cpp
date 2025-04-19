#include "Projection.h"

Point2D doProjection(const Vector3D &point, const double d) {
    //std::cout << "Upper doProjection called" << std::endl;

    Point2D projectedPoint;

    //std::cout << "3D Point: (" << point.x << ", " << point.y << ", " << point.z << ")\n";

    projectedPoint.x = d * point.x / -point.z;
    projectedPoint.y = d * point.y / -point.z;

    //std::cout << "Projected 2D Point: (" << projectedPoint.x << ", " << projectedPoint.y << ")\n";

    return projectedPoint;
}

Lines2D doProjection(const Figures3D &figs) {
    //std::cout << "Lower doProjection called" << std::endl;

    Lines2D lines;
    //std::cout << "Starting projection of figures\n";

    for (const auto &fig : figs) {
        for (const auto &face : fig.faces) {
            for (size_t i = 0; i < face.point_indexes.size() - 1; ++i) {
                int begin = face.point_indexes[i];
                int end = face.point_indexes[i + 1];

                // Punt 1 & 2
                Point2D p1;
                p1.x = 1 * fig.points[begin].x / -fig.points[begin].z;
                p1.y = 1 * fig.points[begin].y / -fig.points[begin].z;

                Point2D p2;
                p2.x = 1 * fig.points[end].x / -fig.points[end].z;
                p2.y = 1 * fig.points[end].y / -fig.points[end].z;

                double z1 = fig.points[begin].z;
                double z2 = fig.points[end].z;

                Line2D line(p1, p2, z1, z2, fig.color);
                lines.push_back(line);

                //std::cout << "Projected line from (" << p1.x << ", " << p1.y << ") to ("
                //          << p2.x << ", " << p2.y << ")\n";
            }

            int begin = face.point_indexes.back();
            int end = face.point_indexes.front();

            // Punt 1 & 2
            Point2D p1;
            p1.x = 1 * fig.points[begin].x / -fig.points[begin].z;
            p1.y = 1 * fig.points[begin].y / -fig.points[begin].z;

            Point2D p2;
            p2.x = 1 * fig.points[end].x / -fig.points[end].z;
            p2.y = 1 * fig.points[end].y / -fig.points[end].z;

            double z1 = fig.points[begin].z;
            double z2 = fig.points[end].z;

            Line2D line(p1, p2, z1, z2, fig.color);
            lines.push_back(line);

            //std::cout << "Projected line from (" << p1.x << ", " << p1.y << ") to ("
            //          << p2.x << ", " << p2.y << ")\n";
        }
    }
    //std::cout << "Finished projection, number of lines: " << lines.size() << std::endl;
    return lines;
}
