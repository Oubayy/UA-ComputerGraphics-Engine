#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>
#include "vector3d.h"
// Removed l_parser.h as it's not directly needed here.
// #include "l_parser.h"
#include "Line2D.h" // Contains Color definition
#include <cmath>
#include <string> // Added for std::string

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
    Color color; // Uses Color from Line2D.h

    Figure() : color(1.0, 1.0, 1.0) {} // Default color to white or another default

    // Method to calculate the geometric center
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
    static Figure generate3DLSystem(const std::string &inputFile, const Color &color);

    static Figure createBuckyBall();
};

// Define a list of Figures3D
// Ensure this typedef doesn't clash if defined elsewhere. Usually put with Figure.h
typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H