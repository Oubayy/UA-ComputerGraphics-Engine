#include "Figure.h"
#include "Transformations.h"
#include "Projection.h"
#include "vector3d.h"
#include "l_parser.h"
#include <vector>
#include <list>
#include <map>
#include "Line2D.h"
#include <functional>
#include <fstream>
#include <stack>

#ifndef Pi
#define Pi 3.14159265358979323846
#endif

using namespace std;

Figure Figure::createCube() {
    // Figuur aanmaken
    Figure cube;
    // Voeg punten toe
    cube.points = {
            Vector3D::point(-1, -1, -1), Vector3D::point(1, -1, -1), Vector3D::point(1, 1, -1), Vector3D::point(-1, 1, -1),
            Vector3D::point(-1, -1, 1), Vector3D::point(1, -1, 1), Vector3D::point(1, 1, 1), Vector3D::point(-1, 1, 1)
    };
    // Voeg vlakken toe
    cube.faces = {
            Face{{0, 1, 2, 3}}, // front face
            Face{{4, 5, 6, 7}}, // back face
            Face{{0, 1, 5, 4}}, // bottom face
            Face{{2, 3, 7, 6}}, // top face
            Face{{0, 3, 7, 4}}, // left face
            Face{{1, 2, 6, 5}}  // right face
    };
    return cube;
}

Figure Figure::createTetrahedron() {
    // Figuur aanmaken
    Figure tetrahedron;

    // Voeg punten toe (=vertices)
    /*tetrahedron.points = {
            Vector3D::point(1, 1, 1),     // Vertex 0
            Vector3D::point(-1, -1, 1),   // Vertex 1
            Vector3D::point(-1, 1, -1),   // Vertex 2
            Vector3D::point(1, -1, -1)    // Vertex 3
    };*/
    tetrahedron.points = {
            Vector3D::point(1, -1, -1),     // Vertex 0
            Vector3D::point(-1, 1, -1),   // Vertex 1
            Vector3D::point(1, 1, 1),   // Vertex 2
            Vector3D::point(-1, -1, 1)    // Vertex 3
    };

    // Voeg vlakken toe (=faces)
    /*tetrahedron.faces = {
            Face{{0, 1, 2}}, // Face connecting vertices 0, 1, and 2
            Face{{0, 1, 3}}, // Face connecting vertices 0, 1, and 3
            Face{{0, 2, 3}}, // Face connecting vertices 0, 2, and 3
            Face{{1, 2, 3}}  // Face connecting vertices 1, 2, and 3
    };*/
    tetrahedron.faces = {
            Face{{0, 1, 2}}, // Face connecting vertices 0, 1, and 2
            Face{{1, 3, 2}}, // Face connecting vertices 0, 1, and 3
            Face{{0, 3, 1}}, // Face connecting vertices 0, 2, and 3
            Face{{0, 2, 3}}  // Face connecting vertices 1, 2, and 3
    };

    return tetrahedron;
}

Figure Figure::createOctahedron() {
    // Figuur aanmaken
    Figure octahedron;

    // Voeg punten toe (=vertices)
    octahedron.points = {
            Vector3D::point(1, 0, 0),    // Vertex 0
            Vector3D::point(-1, 0, 0),   // Vertex 1
            Vector3D::point(0, 1, 0),    // Vertex 2
            Vector3D::point(0, -1, 0),   // Vertex 3
            Vector3D::point(0, 0, 1),    // Vertex 4
            Vector3D::point(0, 0, -1)    // Vertex 5
    };

    // Voeg vlakken toe (=faces)
    octahedron.faces = {
            Face{{0, 2, 4}}, // Face connecting vertices 0, 2, and 4
            Face{{0, 2, 5}}, // Face connecting vertices 0, 2, and 5
            Face{{0, 3, 4}}, // Face connecting vertices 0, 3, and 4
            Face{{0, 3, 5}}, // Face connecting vertices 0, 3, and 5
            Face{{1, 2, 4}}, // Face connecting vertices 1, 2, and 4
            Face{{1, 2, 5}}, // Face connecting vertices 1, 2, and 5
            Face{{1, 3, 4}}, // Face connecting vertices 1, 3, and 4
            Face{{1, 3, 5}}  // Face connecting vertices 1, 3, and 5
    };

    return octahedron;
}

/*
Figure Figure::createIcosahedron() {
    // Figuur aanmaken
    Figure icosahedron;

    const double phi = (1.0 + sqrt(5.0)) / 2.0; // The golden ratio
    const double a = 1.0;
    const double b = 1.0 / phi;

    // Voeg punten toe (vertices of the icosahedron)
    icosahedron.points = {
            Vector3D::point(-a,  b,  0), Vector3D::point( a,  b,  0),
            Vector3D::point(-a, -b,  0), Vector3D::point( a, -b,  0),

            Vector3D::point( 0, -a,  b), Vector3D::point( 0,  a,  b),
            Vector3D::point( 0, -a, -b), Vector3D::point( 0,  a, -b),

            Vector3D::point( b,  0, -a), Vector3D::point( b,  0,  a),
            Vector3D::point(-b,  0, -a), Vector3D::point(-b,  0,  a)
    };

    // Voeg vlakken toe (faces of the icosahedron)
    icosahedron.faces = {
            Face{{0, 11, 5}}, Face{{0, 5, 1}}, Face{{0, 1, 7}}, Face{{0, 7, 10}}, Face{{0, 10, 11}},
            Face{{1, 5, 9}},  Face{{5, 11, 4}}, Face{{11, 10, 2}}, Face{{10, 7, 6}}, Face{{7, 1, 8}},
            Face{{3, 9, 4}},  Face{{3, 4, 2}},  Face{{3, 2, 6}},  Face{{3, 6, 8}},  Face{{3, 8, 9}},
            Face{{4, 9, 5}},  Face{{2, 4, 11}}, Face{{6, 2, 10}}, Face{{8, 6, 7}},  Face{{9, 8, 1}}
    };

    return icosahedron;
} */

Figure Figure::createIcosahedron() {
    Figure icosahedron;

    // Define vertices
    Vector3D point1 = Vector3D::point(0, 0, sqrt(5) / 2);
    icosahedron.points.push_back(point1);

    for (int i = 2; i < 7; ++i) {
        Vector3D point = Vector3D::point(cos((i - 2) * 2 * Pi / 5), sin((i - 2) * 2 * Pi / 5), 0.5);
        icosahedron.points.push_back(point);
    }

    for (int i = 7; i < 12; ++i) {
        Vector3D point = Vector3D::point(cos(Pi / 5 + (i - 7) * 2 * Pi / 5), sin(Pi / 5 + (i - 7) * 2 * Pi / 5), -0.5);
        icosahedron.points.push_back(point);
    }

    Vector3D point12 = Vector3D::point(0, 0, -sqrt(5) / 2);
    icosahedron.points.push_back(point12);

    // Define faces (correctly connecting vertices)
    icosahedron.faces = {
            Face{{0, 1, 2}}, Face{{0, 2, 3}}, Face{{0, 3, 4}}, Face{{0, 4, 5}}, Face{{0, 5, 1}},
            Face{{1, 6, 2}}, Face{{2, 6, 7}}, Face{{2, 7, 3}}, Face{{3, 7, 8}}, Face{{3, 8, 4}},
            Face{{4, 8, 9}}, Face{{4, 9, 5}}, Face{{5, 9, 10}}, Face{{5, 10, 1}}, Face{{1, 10, 6}},
            Face{{11, 6, 7}}, Face{{11, 7, 8}}, Face{{11, 8, 9}}, Face{{11, 9, 10}}, Face{{11, 10, 6}}
    };

    return icosahedron;
}

/*
Figure Figure::createDodecahedron() {
    Figure dodecahedron;

    // Golden ratio
    const double phi = (1 + sqrt(5)) / 2;

    // Vertices of a dodecahedron
    dodecahedron.points = {
            Vector3D::point( 1,  1,  1),   // 0
            Vector3D::point( 1,  1, -1),   // 1
            Vector3D::point( 1, -1,  1),   // 2
            Vector3D::point( 1, -1, -1),   // 3
            Vector3D::point(-1,  1,  1),   // 4
            Vector3D::point(-1,  1, -1),   // 5
            Vector3D::point(-1, -1,  1),   // 6
            Vector3D::point(-1, -1, -1),   // 7
            Vector3D::point( 0,  phi,  1/phi),  // 8
            Vector3D::point( 0,  phi, -1/phi),  // 9
            Vector3D::point( 0, -phi,  1/phi),  // 10
            Vector3D::point( 0, -phi, -1/phi),  // 11
            Vector3D::point( 1/phi, 0,  phi),   // 12
            Vector3D::point( 1/phi, 0, -phi),   // 13
            Vector3D::point(-1/phi, 0,  phi),   // 14
            Vector3D::point(-1/phi, 0, -phi),   // 15
            Vector3D::point( phi,  1/phi, 0),   // 16
            Vector3D::point( phi, -1/phi, 0),   // 17
            Vector3D::point(-phi,  1/phi, 0),   // 18
            Vector3D::point(-phi, -1/phi, 0)    // 19
    };

    // Faces of a dodecahedron (elk vlak is een vijfhoek/pentagon)
    dodecahedron.faces = {
            Face{{0, 8, 10, 2, 12}},  // First face
            Face{{0, 12, 17, 1, 16}},
            Face{{1, 16, 9, 4, 8}},
            Face{{1, 9, 5, 15, 13}},
            Face{{2, 10, 11, 3, 17}},
            Face{{3, 17, 16, 1, 13}},
            Face{{3, 13, 15, 7, 19}},
            Face{{5, 9, 8, 0, 4}},
            Face{{5, 4, 18, 6, 15}},
            Face{{7, 19, 14, 6, 15}},
            Face{{7, 14, 12, 2, 10}},
            Face{{6, 18, 4, 9, 8}}
    };

    return dodecahedron;
}
*/

Figure Figure::createDodecahedron() {
    Figure dodecahedron;
    Figure icosahedron = createIcosahedron();

    // Calculate centroids of each face of the icosahedron to use as vertices of the dodecahedron
    for (const auto& face : icosahedron.faces) {
        Vector3D p1 = icosahedron.points[face.point_indexes[0]];
        Vector3D p2 = icosahedron.points[face.point_indexes[1]];
        Vector3D p3 = icosahedron.points[face.point_indexes[2]];

        // Calculate the centroid of the face
        Vector3D centroid = Vector3D::point(
                (p1.x + p2.x + p3.x) / 3.0,
                (p1.y + p2.y + p3.y) / 3.0,
                (p1.z + p2.z + p3.z) / 3.0
        );


        /* Normalize the centroids
        double magnitude = sqrt(centroid.x * centroid.x + centroid.y * centroid.y + centroid.z * centroid.z);
        if (magnitude != 0) {
            centroid.x /= magnitude;
            centroid.y /= magnitude;
            centroid.z /= magnitude;
        }*/

        // Add the centroids to the dodecahedron points
        dodecahedron.points.push_back(centroid);
    }

    // Define the faces of the dodecahedron based on the centroids calculated
    dodecahedron.faces = {
            Face{{0, 1, 2, 3, 4}},    // Face connecting vertices 0, 1, 2, 3, and 4
            Face{{0, 5, 6, 7, 1}},    // enz..
            Face{{1, 7, 8, 9, 2}},
            Face{{2, 9, 10, 11, 3}},
            Face{{3, 11, 12, 13, 4}},
            Face{{4, 13, 14, 5, 0}},
            Face{{19, 18, 17, 16, 15}},
            Face{{19, 14, 13, 12, 18}},
            Face{{18, 12, 11, 10, 17}},
            Face{{17, 10, 9, 8, 16}},
            Face{{16, 8, 7, 6, 15}},
            Face{{15, 6, 5, 14, 19}}
    };

    return dodecahedron;
}

Figure Figure::createCylinder(int n, double height) {
    Figure cylinder;

    vector<Vector3D> baseVertices;
    vector<Vector3D> topVertices;

    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * Pi * i / n;
        double x = cos(angle);
        double y = sin(angle);
        baseVertices.push_back(Vector3D::point(x, y, 0));
        topVertices.push_back(Vector3D::point(x, y, height));
    }

    // Combine base and top vertices
    cylinder.points.insert(cylinder.points.end(), baseVertices.begin(), baseVertices.end());
    cylinder.points.insert(cylinder.points.end(), topVertices.begin(), topVertices.end());

    // Create base face
    vector<int> baseFaceIndices;
    for (int i = 0; i < n; ++i) {
        baseFaceIndices.push_back(i);
    }
    cylinder.faces.push_back(Face{baseFaceIndices});

    // Create top face
    vector<int> topFaceIndices;
    for (int i = 0; i < n; ++i) {
        topFaceIndices.push_back(n + i); // top vertices start from index n
    }
    cylinder.faces.push_back(Face{topFaceIndices});

    // Create side faces
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        cylinder.faces.push_back(Face{{i, next, n + next, n + i}});
    }

    // extra uitleg: baseVertices and topVertices vectors slaan de points op de basis en op de top (op deze cirkels dus)

    return cylinder;
}

Figure Figure::createCone(int n, double height) {
    Figure cone;
    std::vector<Vector3D> vertices;

    // Create the vertices for the base
    for (int i = 0; i < n; ++i) {
        double angle = 2 * Pi * i / n;
        vertices.push_back(Vector3D::point(cos(angle), sin(angle), 0));
    }

    // Add the apex of the cone
    Vector3D apex = Vector3D::point(0, 0, height);
    vertices.push_back(apex);

    // Store the vertices in the cone object
    cone.points = vertices;

    // Create the base face (n-sided polygon)
    std::vector<int> baseFace;
    for (int i = 0; i < n; ++i) {
        baseFace.push_back(i);
    }
    cone.faces.push_back(Face{baseFace});

    // Create the triangular side faces
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        cone.faces.push_back(Face{{i, next, n}});
    }

    return cone;
}

/*
Figure Figure::createSphere(int n) {
    Figure sphere;

    // Constants for the initial icosahedron vertices
    const double X = 0.525731112119133606;
    const double Z = 0.850650808352039932;

    // Create initial icosahedron vertices
    std::vector<Vector3D> icosahedronVertices = {
            Vector3D::point(-X, 0.0, Z), Vector3D::point(X, 0.0, Z), Vector3D::point(-X, 0.0, -Z), Vector3D::point(X, 0.0, -Z),
            Vector3D::point(0.0, Z, X), Vector3D::point(0.0, Z, -X), Vector3D::point(0.0, -Z, X), Vector3D::point(0.0, -Z, -X),
            Vector3D::point(Z, X, 0.0), Vector3D::point(-Z, X, 0.0), Vector3D::point(Z, -X, 0.0), Vector3D::point(-Z, -X, 0.0)
    };

    // Define the initial icosahedron faces
    std::vector<Face> icosahedronFaces = {
            Face{{0, 1, 4}}, Face{{0, 4, 9}}, Face{{9, 4, 5}}, Face{{4, 8, 5}}, Face{{4, 1, 8}},
            Face{{8, 1, 10}}, Face{{8, 10, 3}}, Face{{5, 8, 3}}, Face{{5, 3, 2}}, Face{{2, 3, 7}},
            Face{{7, 3, 10}}, Face{{7, 10, 6}}, Face{{7, 6, 11}}, Face{{11, 6, 0}}, Face{{0, 6, 1}},
            Face{{6, 10, 1}}, Face{{9, 5, 2}}, Face{{9, 2, 11}}, Face{{11, 2, 7}}, Face{{0, 9, 11}}
    };

    std::vector<Vector3D> vertices;
    std::vector<Face> faces;

    // Define helper functions
    auto midpoint = [](const Vector3D& a, const Vector3D& b) {
        return Vector3D::point(
                (a.x + b.x) / 2,
                (a.y + b.y) / 2,
                (a.z + b.z) / 2
        );
    };

    auto normalize = [](const Vector3D& v) {
        double length = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        return Vector3D::point(v.x / length, v.y / length, v.z / length);
    };

    std::function<void(const Vector3D&, const Vector3D&, const Vector3D&, int)> subdivideAndNormalize;

    subdivideAndNormalize = [&](const Vector3D& v1, const Vector3D& v2, const Vector3D& v3, int depth) {
        if (depth == 0) {
            int index = vertices.size();
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            faces.push_back(Face{{index, index + 1, index + 2}});
            return;
        }

        // Calculate midpoints and normalize them to lie on the sphere surface
        Vector3D v12 = normalize(midpoint(v1, v2));
        Vector3D v23 = normalize(midpoint(v2, v3));
        Vector3D v31 = normalize(midpoint(v3, v1));

        // Recursively subdivide further
        subdivideAndNormalize(v1, v12, v31, depth - 1);
        subdivideAndNormalize(v2, v23, v12, depth - 1);
        subdivideAndNormalize(v3, v31, v23, depth - 1);
        subdivideAndNormalize(v12, v23, v31, depth - 1);
    };

    // Subdivide each triangle of the icosahedron
    for (const auto& face : icosahedronFaces) {
        subdivideAndNormalize(
                icosahedronVertices[face.point_indexes[0]],
                icosahedronVertices[face.point_indexes[1]],
                icosahedronVertices[face.point_indexes[2]],
                n
        );
    }

    // Rescale all points to lie exactly on the sphere's surface
    for (auto& vertex : vertices) {
        vertex = normalize(vertex);
    }

    sphere.points = vertices;
    sphere.faces = faces;

    return sphere;
}
*/

// Helper function to normalize a vector to lie on the sphere surface
Vector3D normalize(const Vector3D& v) {
    double length = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return Vector3D::point(v.x / length, v.y / length, v.z / length);
}

Figure Figure::createSphere(int n) {
    Figure figure;

    // #1: Generate an initial icosahedron
    figure = createIcosahedron();

    // #2: Subdivide each face of the icosahedron n times
    for (int counter = 0; counter < n; ++counter) {
        Figure tempFigure;

        for (const auto& face : figure.faces) {
            Vector3D pointA = figure.points[face.point_indexes[0]];
            Vector3D pointB = figure.points[face.point_indexes[1]];
            Vector3D pointC = figure.points[face.point_indexes[2]];

            Vector3D pointD = normalize((pointA + pointB) / 2);
            Vector3D pointE = normalize((pointA + pointC) / 2);
            Vector3D pointF = normalize((pointB + pointC) / 2);

            tempFigure.points.push_back(pointA);
            tempFigure.points.push_back(pointB);
            tempFigure.points.push_back(pointC);
            tempFigure.points.push_back(pointD);
            tempFigure.points.push_back(pointE);
            tempFigure.points.push_back(pointF);

            tempFigure.faces.push_back(Face{{(int)tempFigure.points.size() - 6, (int)tempFigure.points.size() - 3, (int)tempFigure.points.size() - 2}});
            tempFigure.faces.push_back(Face{{(int)tempFigure.points.size() - 5, (int)tempFigure.points.size() - 1, (int)tempFigure.points.size() - 3}});
            tempFigure.faces.push_back(Face{{(int)tempFigure.points.size() - 4, (int)tempFigure.points.size() - 2, (int)tempFigure.points.size() - 1}});
            tempFigure.faces.push_back(Face{{(int)tempFigure.points.size() - 3, (int)tempFigure.points.size() - 1, (int)tempFigure.points.size() - 2}});
        }

        figure = tempFigure;
    }

    // #4: Normalize all points to lie on the sphere's surface
    for (auto& point : figure.points) {
        point = normalize(point);
    }

    return figure;
}

Figure Figure::createTorus(double R, double r, int n, int m) {
    Figure torus;

    std::vector<Vector3D> vertices;
    std::vector<Face> faces;

    // Generate vertices
    for (int i = 0; i < n; ++i) {
        double theta = 2 * Pi * i / n;
        for (int j = 0; j < m; ++j) {
            double phi = 2 * Pi * j / m;

            double x = (R + r * cos(phi)) * cos(theta);
            double y = (R + r * cos(phi)) * sin(theta);
            double z = r * sin(phi);

            vertices.push_back(Vector3D::point(x, y, z));
        }
    }

    // Generate faces
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int current = i * m + j;
            int nextJ = (j + 1) % m;
            int nextI = (i + 1) % n;

            int currentNextJ = i * m + nextJ;
            int nextINextJ = nextI * m + nextJ;
            int nextICurrentJ = nextI * m + j;

            faces.push_back(Face{{current, currentNextJ, nextINextJ, nextICurrentJ}});
        }
    }

    torus.points = vertices;
    torus.faces = faces;

    return torus;
}

/*
std::string generateLSystemString(const std::string &axiom, const std::map<char, std::string> &rules, int n) {
    std::string currentString = axiom;
    for (int i = 0; i < n; ++i) {
        std::string nextString;
        for (char c : currentString) {
            if (rules.find(c) != rules.end()) {
                nextString += rules.at(c);
            } else {
                nextString += c;
            }
        }
        currentString = nextString;
    }
    return currentString;
}


std::string generateLSystemString(const LParser::LSystem3D &lSystem) {
    std::string currentString = lSystem.get_initiator();
    unsigned int iterations = lSystem.get_nr_iterations();

    // Perform the iterative replacement process
    for (unsigned int i = 0; i < iterations; ++i) {
        std::string nextString;

        for (char symbol : currentString) {
            // Check if the symbol is part of the L-system's alphabet and has a replacement rule
            if (lSystem.get_alphabet().find(symbol) != lSystem.get_alphabet().end()) {
                nextString += lSystem.get_replacement(symbol);
            } else {
                // If there's no replacement rule for the symbol, it remains unchanged
                nextString += symbol;
            }
        }

        currentString = nextString;
    }

    return currentString;
}

// Function to generate the 3D L-System
Figure Figure::generate3DLSystem(const string &inputFile, const Color &color) {
    LParser::LSystem3D l_system;

    //std::ifstream inputStream(inputFile);

    // Adjust the path to include the correct directory
    std::string filePath = "configs/3d_lichamen/" + inputFile;

    std::ifstream inputStream(filePath);
    if (!inputStream.is_open()) {
        throw std::runtime_error("Could not open L3D file: " + filePath);
    }

    inputStream >> l_system;
    inputStream.close();

    Figure figure;
    figure.color = color;

    double angle = l_system.get_angle() * Pi / 180.0;
    //double bangle = l_system.get_replacement()

    Vector3D currentPosition = Vector3D::point(0, 0, 0);
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);

    std::stack<Vector3D> positionStack;
    std::stack<std::vector<Vector3D>> angleStack;

    // Generate the L-System string
    std::string lsystemString = generateLSystemString(l_system);

    for (char command : lsystemString) {
        if (command == '+') {
            H = H * cos(angle) + L * sin(angle);
            L = -H * sin(angle) + L * cos(angle);
        } else if (command == '-') {
            H = H * cos(-angle) + L * sin(-angle);
            L = -H * sin(-angle) + L * cos(-angle);
        } else if (command == '&') {
            H = H * cos(angle) + U * sin(angle);
            U = -H * sin(angle) + U * cos(angle);
        } else if (command == '^') {
            H = H * cos(-angle) + U * sin(angle);
            U = -H * sin(-angle) + U * cos(-angle);
        } else if (command == '\\') {
            L = L * cos(angle) - U * sin(angle);
            U = L * sin(angle) + U * cos(angle);
        } else if (command == '/') {
            L = L * cos(-angle) - U * sin(angle);
            U = L * sin(-angle) + U * cos(-angle);
        } else if (command == '|') {
            H = -H;
            L = -L;
        } else if (command == '(') {
            positionStack.push(currentPosition);
            angleStack.push({H, L, U});
        } else if (command == ')') {
            currentPosition = positionStack.top();
            positionStack.pop();
            std::vector<Vector3D> angles = angleStack.top();
            angleStack.pop();
            H = angles[0];
            L = angles[1];
            U = angles[2];
        } else if (l_system.draw(command)) {
            Vector3D newPosition = currentPosition + H;
            figure.points.push_back(currentPosition);
            figure.points.push_back(newPosition);
            figure.faces.push_back(Face{{(int)figure.points.size() - 2, (int)figure.points.size() - 1}});
            currentPosition = newPosition;
        } else {
            currentPosition += H;
        }
    }

    return figure;
}
*/


Figure Figure::generate3DLSystem(const std::string &inputFile, const Color &color) {
    LParser::LSystem3D l_system;

    // Adjust the path to include the correct directory
    std::string filePath = inputFile;

    std::ifstream inputStream(filePath);
    if (!inputStream.is_open()) {
        throw std::runtime_error("Could not open L3D file: " + filePath);
    }

    inputStream >> l_system;
    inputStream.close();

    Figure figure;
    figure.color = color;

    double angle = l_system.get_angle() * Pi / 180.0;
    Vector3D currentPosition = Vector3D::point(0, 0, 0);
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);

    std::stack<Vector3D> positionStack;
    std::stack<std::vector<Vector3D>> angleStack;

    // Generate the L-System string
    std::set<char> alphabet = l_system.get_alphabet();
    std::string lsystemString = l_system.get_initiator();
    std::string nextString;

    for (unsigned int i = 0; i < l_system.get_nr_iterations(); ++i) {
        nextString.clear();
        for (char symbol : lsystemString) {
            if (alphabet.find(symbol) != alphabet.end()) {
                nextString += l_system.get_replacement(symbol);
            } else {
                nextString += symbol;
            }
        }
        lsystemString = nextString;
    }

    Vector3D _H = H;
    Vector3D _L = L;
    Vector3D _U = U;

    for (char command : lsystemString) {
        switch (command) {
            case '+':
                _H = H;
                _L = L;

                H = _H * cos(angle) + _L * sin(angle);
                L = -_H * sin(angle) + _L * cos(angle);
                break;
            case '-':
                _H = H;
                _L = L;

                H = _H * cos(-angle) + _L * sin(-angle);
                L = -_H * sin(-angle) + _L * cos(-angle);
                break;
            case '^':
                _H = H;
                _U = U;

                H = _H * cos(angle) + _U * sin(angle);
                U = -_H * sin(angle) + _U * cos(angle);
                break;
            case '&':
                _H = H;
                _U = U;

                H = _H * cos(-angle) + _U * sin(-angle);
                U = -_H * sin(-angle) + _U * cos(-angle);
                break;
            case '\\':
                _U = U;
                _L = L;

                U = _U * cos(angle) + _L * sin(angle);
                L = -_U * sin(angle) + _L * cos(angle);
                break;
            case '/':
                _U = U;
                _L = L;

                U = _U * cos(-angle) + _L * sin(-angle);
                L = -_U * sin(-angle) + _L * cos(-angle);
                break;
            case '|':
                H = -H;
                L = -L;
                break;
            case '(':
                positionStack.push(currentPosition);
                angleStack.push({H, L, U});
                break;
            case ')': {
                currentPosition = positionStack.top();
                positionStack.pop();
                std::vector<Vector3D> angles = angleStack.top();
                angleStack.pop();
                H = angles[0];
                L = angles[1];
                U = angles[2];
                break;
            }
            default:
                if (l_system.draw(command)) {
                    Vector3D newPosition = currentPosition + H;
                    figure.points.push_back(currentPosition);
                    figure.points.push_back(newPosition);
                    figure.faces.push_back(Face{{(int)figure.points.size() - 2, (int)figure.points.size() - 1}});
                    currentPosition = newPosition;
                } else {
                    currentPosition += H;
                }
                break;
        }
    }

    return figure;
}
