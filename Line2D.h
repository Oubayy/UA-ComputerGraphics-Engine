#ifndef LINE2D_H
#define LINE2D_H

#include <list>

class Color {
public:
    double red;
    double green;
    double blue;

    Color() : red(0), green(0), blue(0) {}
    Color(double r = 0.0, double g = 0.0, double b = 0.0) : red(r), green(g), blue(b) {}
};

class Point2D {
public:
    double x;
    double y;

    Point2D(double xCoord = 0.0, double yCoord = 0.0) : x(xCoord), y(yCoord) {}
};

class Line2D {
public:
    Point2D p1;
    Point2D p2;

    double z1; // z waarde = depth
    double z2;

    Color color;

    Line2D(Point2D start, Point2D end, double z1, double z2, Color lineColor) : p1(start), p2(end), z1(z1), z2(z2), color(lineColor) {}
};
using Lines2D = std::list<Line2D>;


class Triangle {
public:
    Point2D p1;
    Point2D p2;
    Point2D p3;

    double z1;
    double z2;
    double z3;

    Color color;
};
using Triangles = std::list<Triangle>;


#endif // LINE2D_H