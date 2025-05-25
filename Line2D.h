// Line2D.h
#ifndef LINE2D_H
#define LINE2D_H

#include "vector3d.h" // For Vector3D in Triangle struct
#include <list>

// Color class defined here
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

    double z1;
    double z2;

    Color color; // For simple colored lines (e.g., 2D L-systems, wireframes)

    Line2D(Point2D start, Point2D end, double _z1, double _z2, Color lineColor) :
        p1(start), p2(end), z1(_z1), z2(_z2), color(lineColor) {}
};
using Lines2D = std::list<Line2D>;

class FigurePropertiesForTriangle {
public:
    Color ambientReflection;
    Color diffuseReflection;
    Color specularReflection;
    double reflectionCoefficient;

    FigurePropertiesForTriangle() :
        ambientReflection(0.0,0.0,0.0), // Initialize to prevent uninitialized reads
        diffuseReflection(0.0,0.0,0.0),
        specularReflection(0.0,0.0,0.0),
        reflectionCoefficient(1.0) {}
};

class Triangle {
public:
    Point2D p1;
    Point2D p2;
    Point2D p3;

    double z1_eye;
    double z2_eye;
    double z3_eye;

    Vector3D p1_eye_space;
    Vector3D p2_eye_space;
    Vector3D p3_eye_space;

    FigurePropertiesForTriangle material;

    Triangle(const Point2D& _p1, const Point2D& _p2, const Point2D& _p3,
             double _z1_eye, double _z2_eye, double _z3_eye,
             const Vector3D& _p1_eye_space, const Vector3D& _p2_eye_space, const Vector3D& _p3_eye_space,
             const FigurePropertiesForTriangle& _mat) :
        p1(_p1), p2(_p2), p3(_p3),
        z1_eye(_z1_eye), z2_eye(_z2_eye), z3_eye(_z3_eye),
        p1_eye_space(_p1_eye_space), p2_eye_space(_p2_eye_space), p3_eye_space(_p3_eye_space),
        material(_mat)
    {}
};
using Triangles = std::list<Triangle>;

#endif // LINE2D_H