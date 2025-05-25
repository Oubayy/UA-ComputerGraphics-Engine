// Figure.h
#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>
#include <string>
#include "vector3d.h" // Vector3D must be defined before use
#include "Line2D.h"   // This is where Color is defined

// Face class
class Face {
public:
    std::vector<int> point_indexes;
};

// Figure class
class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    // Lighting properties
    Color ambientReflection;    // Ka
    Color diffuseReflection;    // Kd
    Color specularReflection;   // Ks
    double reflectionCoefficient; // m_s (shininess/specular exponent)

    // Old color member for backward compatibility (e.g., simple wireframes, L-systems, fractal base)
    Color color;

    Figure() :
        ambientReflection(0.1, 0.1, 0.1),
        diffuseReflection(0.7, 0.7, 0.7),
        specularReflection(0.0, 0.0, 0.0),
        reflectionCoefficient(1.0),
        color(1.0, 1.0, 1.0) // Default for the old 'color' member
    {}

    Vector3D center() const;

    static Figure createCube();
    static Figure createTetrahedron();
    static Figure createOctahedron();
    static Figure createIcosahedron();
    static Figure createDodecahedron();
    static Figure createCylinder(int n, double height);
    static Figure createCone(int n, double height);
    static Figure createSphere(int n);
    static Figure createTorus(double R, double r, int n, int m);
    // line_color_fallback can be used to set default ambient for L-system lines if no other material defined
    static Figure generate3DLSystem(const std::string &inputFile, const Color &line_color_fallback);
    static Figure createBuckyBall();
};

typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H