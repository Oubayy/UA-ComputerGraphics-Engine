#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>
//#include "../vector3d.h"
#include "l_parser.h"
#include "Line2D.h"
#include <cmath>

// Face class
class Face {
public:
    std::vector<int> point_indexes;
};

// Figure class
class Figure {
public:
    //std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;

    Figure() : color(1.0, 0.0, 0.0) {} // Default constructor with explicit initialization

    static Figure createCube();
    static Figure createTetrahedron();
    static Figure createOctahedron();
    static Figure createIcosahedron();
    static Figure createDodecahedron();
    static Figure createCylinder(int n, double height);
    static Figure createCone(int n, double height);

    static Figure createSphere(int n);
    //static Vector3D normalize(const Vector3D& v);
    //void subdivideTriangle(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3, int depth, std::vector<Vector3D>& vertices, std::vector<Face>& faces);

    static Figure createTorus(double R, double r, int n, int m);

    static Figure generate3DLSystem(const std::string &inputFile, const Color &color);
};

// Define a list of Figures3D
typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
