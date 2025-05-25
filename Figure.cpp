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
#include <cmath> // For sqrt, cos, sin
#include <stdexcept> // For runtime_error
#include <numeric>   // For center calculation sum (optional)

#ifndef Pi
#define Pi 3.14159265358979323846
#endif

using namespace std;

Figure Figure::createBuckyBall() {
    // Placeholder:
    Figure bucky = createIcosahedron();
    bucky.color = Color(0.9, 0.9, 0.9);
    // TODO: Implement actual Buckyball generation (Truncated Icosahedron)
    return bucky;
}

// Definition for the center() method
Vector3D Figure::center() const {
    if (points.empty()) {
        return Vector3D::point(0, 0, 0); // Return origin if no points
    }

    double sumX = 0, sumY = 0, sumZ = 0;
    for (const Vector3D& p : points) {
        sumX += p.x;
        sumY += p.y;
        sumZ += p.z;
    }
    size_t numPoints = points.size();
    return Vector3D::point(sumX / numPoints, sumY / numPoints, sumZ / numPoints);
}


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
            Face{{0, 1, 2, 3}}, // front face (-Z)
            Face{{7, 6, 5, 4}}, // back face (+Z) - Ensure consistent winding order if needed (e.g., CCW from outside)
            Face{{4, 5, 1, 0}}, // bottom face (-Y)
            Face{{3, 2, 6, 7}}, // top face (+Y)
            Face{{4, 0, 3, 7}}, // left face (-X)
            Face{{1, 5, 6, 2}}  // right face (+X)
    };
    return cube;
}

Figure Figure::createTetrahedron() {
    // Figuur aanmaken
    Figure tetrahedron;

    tetrahedron.points = {
            Vector3D::point(1, -1, -1),     // Vertex 0
            Vector3D::point(-1, 1, -1),   // Vertex 1
            Vector3D::point(-1, -1, 1),   // Vertex 2 - Changed to form a standard tetrahedron base
            Vector3D::point(1, 1, 1)      // Vertex 3 - Apex
    };

    // Faces (ensure consistent winding order, e.g., CCW from outside)
    tetrahedron.faces = {
            Face{{0, 1, 2}}, // Base triangle
            Face{{0, 3, 1}}, // Side face
            Face{{1, 3, 2}}, // Side face
            Face{{2, 3, 0}}  // Side face
    };

    return tetrahedron;
}

Figure Figure::createOctahedron() {
    Figure octahedron;
    octahedron.points = {
        Vector3D::point( 1,  0,  0), Vector3D::point(-1,  0,  0), // 0, 1 (+x, -x)
        Vector3D::point( 0,  1,  0), Vector3D::point( 0, -1,  0), // 2, 3 (+y, -y)
        Vector3D::point( 0,  0,  1), Vector3D::point( 0,  0, -1)  // 4, 5 (+z, -z) - Top, Bottom
    };

    // Faces (ensure consistent winding order, e.g., CCW from outside)
    octahedron.faces = {
        // Top faces (connected to vertex 4)
        Face{{0, 2, 4}}, Face{{2, 1, 4}}, Face{{1, 3, 4}}, Face{{3, 0, 4}},
        // Bottom faces (connected to vertex 5)
        Face{{0, 5, 2}}, Face{{2, 5, 1}}, Face{{1, 5, 3}}, Face{{3, 5, 0}}
    };
    return octahedron;
}


Figure Figure::createIcosahedron() {
    Figure icosahedron;

    // default color, can be overridden later
    icosahedron.color = Color(0.5, 0.5, 1.0);

    icosahedron.points.push_back(Vector3D::point(0, 0, std::sqrt(5.0) / 2.0));

    for (int i = 0; i < 5; ++i) {
        icosahedron.points.push_back(Vector3D::point(
            std::cos(i * 2.0 * Pi / 5.0),
            std::sin(i * 2.0 * Pi / 5.0),
            0.5
        ));
    }

    for (int i = 0; i < 5; ++i) {
        icosahedron.points.push_back(Vector3D::point(
            std::cos(Pi / 5.0 + i * 2.0 * Pi / 5.0),
            std::sin(Pi / 5.0 + i * 2.0 * Pi / 5.0),
            -0.5
        ));
    }

    icosahedron.points.push_back(Vector3D::point(0, 0, -std::sqrt(5.0) / 2.0));

    icosahedron.faces.push_back(Face({{0, 1, 2}}));
    icosahedron.faces.push_back(Face({{0, 2, 3}}));
    icosahedron.faces.push_back(Face({{0, 3, 4}}));
    icosahedron.faces.push_back(Face({{0, 4, 5}}));
    icosahedron.faces.push_back(Face({{0, 5, 1}}));

    icosahedron.faces.push_back(Face({{1, 6, 2}}));
    icosahedron.faces.push_back(Face({{2, 6, 7}}));
    icosahedron.faces.push_back(Face({{2, 7, 3}}));
    icosahedron.faces.push_back(Face({{3, 7, 8}}));
    icosahedron.faces.push_back(Face({{3, 8, 4}}));
    icosahedron.faces.push_back(Face({{4, 8, 9}}));
    icosahedron.faces.push_back(Face({{4, 9, 5}}));
    icosahedron.faces.push_back(Face({{5, 9, 10}}));
    icosahedron.faces.push_back(Face({{5, 10, 1}}));
    icosahedron.faces.push_back(Face({{1, 10, 6}}));

    icosahedron.faces.push_back(Face({{11, 7, 6}}));
    icosahedron.faces.push_back(Face({{11, 8, 7}}));
    icosahedron.faces.push_back(Face({{11, 9, 8}}));
    icosahedron.faces.push_back(Face({{11, 10, 9}}));
    icosahedron.faces.push_back(Face({{11, 6, 10}}));

    for(auto& p : icosahedron.points) {
        double len = p.length();
        if (len > 1e-9) {
             p.x /= len; p.y /= len; p.z /= len;
        }
    }

    return icosahedron;
}


Figure Figure::createDodecahedron() {
    Figure dodecahedron;

    // default color, can be overridden later
    dodecahedron.color = Color(0.9, 0.7, 0.3);

    Figure icosahedron = createIcosahedron();

    if (icosahedron.points.empty() || icosahedron.faces.empty()) {
        std::cerr << "Error: Base icosahedron for dodecahedron generation is empty." << std::endl;
        return dodecahedron;
    }
    if (icosahedron.faces.size() != 20) {
        std::cerr << "Error: Base icosahedron does not have 20 faces for dodecahedron generation." << std::endl;
        return dodecahedron;
    }

    dodecahedron.points.clear();
    for (const auto& ico_face : icosahedron.faces) {
        if (ico_face.point_indexes.size() == 3) {
            const Vector3D& p1_ico = icosahedron.points[ico_face.point_indexes[0]];
            const Vector3D& p2_ico = icosahedron.points[ico_face.point_indexes[1]];
            const Vector3D& p3_ico = icosahedron.points[ico_face.point_indexes[2]];

            Vector3D centroid = (p1_ico + p2_ico + p3_ico) / 3.0;

            double len = centroid.length();
            if (len > 1e-9) {
                centroid.x /= len;
                centroid.y /= len;
                centroid.z /= len;
            }
            dodecahedron.points.push_back(centroid);
        } else {
            std::cerr << "Warning: Icosahedron face encountered that is not a triangle during dodecahedron generation." << std::endl;
        }
    }

    if (dodecahedron.points.size() != 20) {
        std::cerr << "Error: Dodecahedron vertex generation resulted in " << dodecahedron.points.size() << " points, expected 20." << std::endl;
        return Figure();
    }

    dodecahedron.faces = {
        Face({{0, 1, 2, 3, 4}}),
        Face({{0, 5, 6, 7, 1}}),
        Face({{1, 7, 8, 9, 2}}),
        Face({{2, 9, 10, 11, 3}}),
        Face({{3, 11, 12, 13, 4}}),
        Face({{4, 13, 14, 5, 0}}),

        Face({{19, 18, 17, 16, 15}}),
        Face({{19, 14, 13, 12, 18}}),
        Face({{18, 12, 11, 10, 17}}),
        Face({{17, 10, 9, 8, 16}}),
        Face({{16, 8, 7, 6, 15}}),
        Face({{15, 6, 5, 14, 19}})
    };

    return dodecahedron;
}


Figure Figure::createCylinder(int n, double height) {
    Figure cylinder;
    if (n < 3) n = 3; // Need at least 3 points for a base

    // Bottom circle points (z=0)
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * Pi * i / n;
        cylinder.points.push_back(Vector3D::point(cos(angle), sin(angle), 0));
    }
    // Top circle points (z=height)
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * Pi * i / n;
        cylinder.points.push_back(Vector3D::point(cos(angle), sin(angle), height));
    }

    // Create bottom face (indices 0 to n-1) - Reverse order for CCW from outside? Check convention.
    vector<int> baseFaceIndices;
    // If viewed from below, CCW is 0, 1, 2... If viewed from above, CW is 0, 1, 2...
    // Let's assume CCW from outside means CCW when looking towards origin. So base face (z=0) needs CW winding.
    for (int i = n - 1; i >= 0; --i) {
         baseFaceIndices.push_back(i);
    }
    // Or simply CCW if that's the convention:
    // for (int i = 0; i < n; ++i) { baseFaceIndices.push_back(i); }
    cylinder.faces.push_back(Face{baseFaceIndices});


    // Create top face (indices n to 2n-1) - Needs CCW winding when viewed from above.
    vector<int> topFaceIndices;
    for (int i = 0; i < n; ++i) {
        topFaceIndices.push_back(n + i);
    }
    cylinder.faces.push_back(Face{topFaceIndices});

    // Create side faces (quads: bottom_i, bottom_next, top_next, top_i)
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        // Indices: i, next, n + next, n + i (CCW winding from outside)
        cylinder.faces.push_back(Face{{i, next, n + next, n + i}});
    }

    return cylinder;
}

Figure Figure::createCone(int n, double height) {
    Figure cone;
    if (n < 3) n = 3;

    // Base circle points (z=0)
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * Pi * i / n;
        cone.points.push_back(Vector3D::point(cos(angle), sin(angle), 0));
    }
    // Apex point (index n)
    cone.points.push_back(Vector3D::point(0, 0, height));
    int apexIndex = n;

    // Create base face (indices 0 to n-1) - Needs CW winding if viewed from above for CCW outside.
    vector<int> baseFaceIndices;
    for (int i = n - 1; i >= 0; --i) { // CW order
        baseFaceIndices.push_back(i);
    }
    // Or CCW if convention demands:
    // for (int i = 0; i < n; ++i) { baseFaceIndices.push_back(i); }
    cone.faces.push_back(Face{baseFaceIndices});


    // Create side faces (triangles: base_i, base_next, apex)
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        // Indices: i, next, apexIndex (CCW winding from outside)
        cone.faces.push_back(Face{{i, next, apexIndex}});
    }

    return cone;
}

// Helper function to normalize a vector - place inside Figure.cpp or globally if needed elsewhere
static Vector3D normalize(const Vector3D& v) {
    double length = v.length();
    if (length < 1e-9) return v; // Avoid division by zero, return original vector
    return v / length; // Assuming Vector3D supports division by scalar
}

// Helper for sphere: recursively subdivide triangle faces
static void subdivideAndNormalize(Figure& figure, int n_iterations) {
    for (int iter = 0; iter < n_iterations; ++iter) {
        Figure subdividedFigure;
        subdividedFigure.color = figure.color; // Preserve color
        std::map<std::pair<int, int>, int> midpoint_cache; // Cache midpoints to avoid duplicates

        auto get_midpoint_index = [&](int p1_idx, int p2_idx) -> int {
            // Ensure order for cache key
            if (p1_idx > p2_idx) std::swap(p1_idx, p2_idx);
            auto key = std::make_pair(p1_idx, p2_idx);

            auto it = midpoint_cache.find(key);
            if (it != midpoint_cache.end()) {
                return it->second; // Return cached index
            }

            // Calculate, normalize, add midpoint, and cache index
            const Vector3D& p1 = figure.points[p1_idx];
            const Vector3D& p2 = figure.points[p2_idx];
            Vector3D midpoint = normalize((p1 + p2) * 0.5); // Use overloaded ops if available
            subdividedFigure.points.push_back(midpoint);
            int new_index = subdividedFigure.points.size() - 1;
            midpoint_cache[key] = new_index;
            return new_index;
        };

        // Add original points to the new figure first
        subdividedFigure.points = figure.points;

        // Process each face of the old figure
        for (const auto& face : figure.faces) {
             if (face.point_indexes.size() != 3) continue; // Only subdivide triangles

             int a_idx = face.point_indexes[0];
             int b_idx = face.point_indexes[1];
             int c_idx = face.point_indexes[2];

             // Get indices of midpoints (calculates and adds them if not cached)
             int ab_mid_idx = get_midpoint_index(a_idx, b_idx);
             int bc_mid_idx = get_midpoint_index(b_idx, c_idx);
             int ca_mid_idx = get_midpoint_index(c_idx, a_idx);

             // Create 4 new faces
             subdividedFigure.faces.push_back(Face{{a_idx, ab_mid_idx, ca_mid_idx}});
             subdividedFigure.faces.push_back(Face{{b_idx, bc_mid_idx, ab_mid_idx}});
             subdividedFigure.faces.push_back(Face{{c_idx, ca_mid_idx, bc_mid_idx}});
             subdividedFigure.faces.push_back(Face{{ab_mid_idx, bc_mid_idx, ca_mid_idx}}); // Central triangle
        }
        figure = std::move(subdividedFigure); // Replace the old figure with the subdivided one
    }

     // Final normalization of all points after all subdivisions
     for (auto& point : figure.points) {
         point = normalize(point);
     }
}


Figure Figure::createSphere(int n) {
    // Start with an Icosahedron as the base shape
    Figure sphere = createIcosahedron();
    sphere.color = Color(1.0, 1.0, 1.0); // Default sphere color

    // Subdivide and normalize
    subdivideAndNormalize(sphere, n);

    return sphere;
}


Figure Figure::createTorus(double R, double r, int n, int m) {
    Figure torus;
     if (n < 3) n = 3;
     if (m < 3) m = 3;

    // Generate vertices
    for (int i = 0; i < n; ++i) { // Main ring segments
        double u = 2.0 * Pi * i / n; // Angle around Y-axis (or Z)
        for (int j = 0; j < m; ++j) { // Tube segments
            double v = 2.0 * Pi * j / m; // Angle around tube center

            double x = (R + r * cos(v)) * cos(u);
            double y = (R + r * cos(v)) * sin(u);
            double z = r * sin(v);

            torus.points.push_back(Vector3D::point(x, y, z));
        }
    }

    // Generate faces (quads)
    for (int i = 0; i < n; ++i) {
        int i_next = (i + 1) % n;
        for (int j = 0; j < m; ++j) {
            int j_next = (j + 1) % m;

            // Indices of the 4 corners of the quad
            int idx00 = i * m + j;
            int idx10 = i_next * m + j;
            int idx11 = i_next * m + j_next;
            int idx01 = i * m + j_next;

            // Ensure CCW winding from outside
            torus.faces.push_back(Face{{idx00, idx10, idx11, idx01}});
        }
    }

    return torus;
}

Figure Figure::generate3DLSystem(const std::string &inputFile, const Color &color) {
    LParser::LSystem3D l_system;

    // Adjust the path to include the correct directory if needed
    std::string filePath = inputFile; // Assume inputFile is the full or relative path

    std::ifstream inputStream(filePath);
    if (!inputStream.is_open()) {
        // It's often better to throw an exception than print to cerr and continue
        throw std::runtime_error("Could not open L-System file: " + filePath);
    }

    try {
        inputStream >> l_system; // Use the LParser library to read the L-system
    } catch (const std::exception& e) {
         inputStream.close();
         throw std::runtime_error("Error parsing L-System file '" + filePath + "': " + e.what());
    }
    inputStream.close();

    Figure figure;
    figure.color = color; // Assign the provided color

    // Initial state
    double angle_rad = l_system.get_angle() * Pi / 180.0; // Angle in radians
    Vector3D currentPosition = Vector3D::point(0, 0, 0);
    Vector3D H = Vector3D::vector(1, 0, 0); // Heading
    Vector3D L = Vector3D::vector(0, 1, 0); // Left
    Vector3D U = Vector3D::vector(0, 0, 1); // Up

    std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D>> stateStack; // Stack for position, H, L, U

    // Generate the full instruction string
    std::string lsystemString = l_system.get_initiator();
    std::set<char> alphabet = l_system.get_alphabet();
    bool ignoreUnknown = true; // Decide how to handle symbols not in alphabet or commands

    for (unsigned int i = 0; i < l_system.get_nr_iterations(); ++i) {
        std::string nextString = "";
        for (char symbol : lsystemString) {
            if (alphabet.count(symbol)) { // Check if symbol is in the alphabet for replacement
                nextString += l_system.get_replacement(symbol);
            } else {
                // Keep commands (+, -, ^, &, \, /, |, (, )) and potentially other known symbols
                 if (string("+-^&\\/|()").find(symbol) != string::npos || !ignoreUnknown) {
                     nextString += symbol;
                 } else if (!alphabet.count(symbol)) {
                     // If it's not a command and not in alphabet, keep it only if it should be drawn
                     if (l_system.draw(symbol)) {
                         nextString += symbol;
                     }
                     // else: ignore unknown non-drawable symbols silently
                 }
            }
        }
        lsystemString = nextString;
    }

    // Process the final instruction string
    int pointIndexOffset = 0; // To keep track of indices for faces

    for (char command : lsystemString) {
        Vector3D _H, _L, _U; // Temporary holders for rotations

        switch (command) {
            case '+': // Turn Left (+)
                _H = H; _L = L;
                H = _H * cos(angle_rad) + _L * sin(angle_rad);
                L = -_H * sin(angle_rad) + _L * cos(angle_rad);
                break;
            case '-': // Turn Right (-)
                 _H = H; _L = L;
                 H = _H * cos(-angle_rad) + _L * sin(-angle_rad);
                 L = -_H * sin(-angle_rad) + _L * cos(-angle_rad);
                 break;
            case '^': // Pitch Up (^)
                 _H = H; _U = U;
                 H = _H * cos(angle_rad) + _U * sin(angle_rad);
                 U = -_H * sin(angle_rad) + _U * cos(angle_rad);
                 break;
            case '&': // Pitch Down (&)
                 _H = H; _U = U;
                 H = _H * cos(-angle_rad) + _U * sin(-angle_rad);
                 U = -_H * sin(-angle_rad) + _U * cos(-angle_rad);
                 break;
            case '\\': // Roll Left (\) - Note: \ needs escaping in string literals if used directly
                 _L = L; _U = U;
                 L = _L * cos(angle_rad) + _U * sin(angle_rad);
                 U = -_L * sin(angle_rad) + _U * cos(angle_rad);
                 break;
            case '/': // Roll Right (/)
                 _L = L; _U = U;
                 L = _L * cos(-angle_rad) + _U * sin(-angle_rad);
                 U = -_L * sin(-angle_rad) + _U * cos(-angle_rad);
                 break;
            case '|': // Turn Around (|)
                 H = -H;
                 L = -L; // Roll 180 degrees
                 break;
            case '(': // Push state
                stateStack.push({currentPosition, H, L, U});
                break;
            case ')': // Pop state
                if (!stateStack.empty()) {
                     auto state = stateStack.top();
                     stateStack.pop();
                     currentPosition = std::get<0>(state);
                     H = std::get<1>(state);
                     L = std::get<2>(state);
                     U = std::get<3>(state);
                } // else: stack underflow, ignore or warn
                 break;
            default: // Check if it's a symbol to be drawn
                if (l_system.draw(command)) {
                     Vector3D nextPosition = currentPosition + H; // Move forward by H vector
                     // Add start and end points
                     figure.points.push_back(currentPosition);
                     figure.points.push_back(nextPosition);
                     // Create a face representing this line segment
                     figure.faces.push_back(Face{{(int)figure.points.size() - 2, (int)figure.points.size() - 1}});
                     // Update current position
                     currentPosition = nextPosition;
                } else {
                    // If symbol is not drawn but is in alphabet (like a variable used for growth), just move
                    // Or if it's an unknown symbol we decided not to ignore, also just move? Or do nothing?
                    // Standard L-System interpretation: If not draw, just move.
                    // If it's truly unknown, perhaps do nothing.
                     // Let's assume move if not drawing.
                     currentPosition = currentPosition + H;
                }
                break;
        }
    }

    return figure;
}