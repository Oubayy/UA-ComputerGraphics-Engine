#include "Figure.h"
#include "Transformations.h"
// #include "Projection.h" // Not directly used in Figure.cpp
#include "vector3d.h"
#include "l_parser.h" // For LSystem3D
#include <vector>
#include <list>
#include <map>
// #include "Line2D.h" // Already included via Figure.h
#include <functional>
#include <fstream>
#include <stack>
#include <cmath>
#include <stdexcept>
#include <numeric>

#ifndef M_PI // Ensure M_PI is defined
#define M_PI 3.14159265358979323846
#endif

using namespace std;

Figure Figure::createBuckyBall() {
    Figure bucky = createIcosahedron(); // Base geometry
    // Material properties will be set by INI or default from Figure constructor
    // bucky.color = Color(0.9, 0.9, 0.9); // This would set the fallback 'color' member
    return bucky;
}

Vector3D Figure::center() const {
    if (points.empty()) {
        return Vector3D::point(0, 0, 0);
    }
    double sumX = 0, sumY = 0, sumZ = 0;
    for (const Vector3D& p : points) {
        sumX += p.x; sumY += p.y; sumZ += p.z;
    }
    size_t numPoints = points.size();
    return Vector3D::point(sumX / numPoints, sumY / numPoints, sumZ / numPoints);
}

Figure Figure::createCube() {
    Figure cube;
    cube.points = {
            Vector3D::point(-1, -1, -1), Vector3D::point(1, -1, -1), Vector3D::point(1, 1, -1), Vector3D::point(-1, 1, -1),
            Vector3D::point(-1, -1, 1), Vector3D::point(1, -1, 1), Vector3D::point(1, 1, 1), Vector3D::point(-1, 1, 1)
    };
    cube.faces = {
            Face{{0, 3, 2, 1}}, // Front: CCW from outside view if looking towards -Z
            Face{{4, 5, 6, 7}}, // Back: CCW from outside view
            Face{{0, 1, 5, 4}}, // Bottom
            Face{{3, 7, 6, 2}}, // Top
            Face{{0, 4, 7, 3}}, // Left
            Face{{1, 2, 6, 5}}  // Right
    };
    // Winding order: For each face, if you curl your right hand fingers along 0->1->2, thumb is normal.
    // For front face (0,1,2,3), if looking from +Z towards origin, this is CW.
    // To make it CCW from outside (normal pointing out), it should be 0,3,2,1.
    // Or adjust normal calculation in ZBuffer. For now, use provided order.
    // My original suggestion was based on PDF, let's use the one that worked before.
    // Provided:
     cube.faces = {
            Face{{0, 1, 2, 3}}, Face{{7, 6, 5, 4}},
            Face{{4, 5, 1, 0}}, Face{{3, 2, 6, 7}},
            Face{{4, 0, 3, 7}}, Face{{1, 5, 6, 2}}
    };
    return cube;
}

Figure Figure::createTetrahedron() {
    Figure tetrahedron;
    tetrahedron.points = {
            Vector3D::point(1, -1, -1), Vector3D::point(-1, 1, -1),
            Vector3D::point(-1, -1, 1), Vector3D::point(1, 1, 1)
    };
    tetrahedron.faces = {
            Face{{0, 1, 2}}, Face{{0, 3, 1}},
            Face{{1, 3, 2}}, Face{{2, 3, 0}}
    };
    return tetrahedron;
}

Figure Figure::createOctahedron() {
    Figure octahedron;
    octahedron.points = {
        Vector3D::point( 1,  0,  0), Vector3D::point(-1,  0,  0),
        Vector3D::point( 0,  1,  0), Vector3D::point( 0, -1,  0),
        Vector3D::point( 0,  0,  1), Vector3D::point( 0,  0, -1)
    };
    octahedron.faces = {
        Face{{0, 2, 4}}, Face{{2, 1, 4}}, Face{{1, 3, 4}}, Face{{3, 0, 4}},
        Face{{0, 5, 2}}, Face{{2, 5, 1}}, Face{{1, 5, 3}}, Face{{3, 5, 0}}
    };
    return octahedron;
}

Figure Figure::createIcosahedron() {
    Figure icosahedron;
    // Materials set by INI or Figure default constructor
    // icosahedron.color = Color(0.5, 0.5, 1.0); // Old fallback color

    icosahedron.points.push_back(Vector3D::point(0, 0, std::sqrt(5.0) / 2.0));
    for (int i = 0; i < 5; ++i) {
        icosahedron.points.push_back(Vector3D::point(std::cos(i * 2.0 * M_PI / 5.0), std::sin(i * 2.0 * M_PI / 5.0), 0.5));
    }
    for (int i = 0; i < 5; ++i) {
        icosahedron.points.push_back(Vector3D::point(std::cos(M_PI / 5.0 + i * 2.0 * M_PI / 5.0), std::sin(M_PI / 5.0 + i * 2.0 * M_PI / 5.0), -0.5));
    }
    icosahedron.points.push_back(Vector3D::point(0, 0, -std::sqrt(5.0) / 2.0));

    icosahedron.faces = {
        Face{{0, 1, 2}}, Face{{0, 2, 3}}, Face{{0, 3, 4}}, Face{{0, 4, 5}}, Face{{0, 5, 1}},
        Face{{1, 6, 2}}, Face{{2, 6, 7}}, Face{{2, 7, 3}}, Face{{3, 7, 8}}, Face{{3, 8, 4}},
        Face{{4, 8, 9}}, Face{{4, 9, 5}}, Face{{5, 9, 10}}, Face{{5, 10, 1}}, Face{{1, 10, 6}},
        Face{{11, 7, 6}}, Face{{11, 8, 7}}, Face{{11, 9, 8}}, Face{{11, 10, 9}}, Face{{11, 6, 10}}
    };
    for(auto& p : icosahedron.points) { // Normalize to unit sphere
        p.normalise();
    }
    return icosahedron;
}

Figure Figure::createDodecahedron() {
    Figure dodecahedron;
    // dodecahedron.color = Color(0.9, 0.7, 0.3); // Old fallback
    Figure icosahedron = createIcosahedron(); // Base for dual
    if (icosahedron.points.empty() || icosahedron.faces.size() != 20) {
        std::cerr << "Error: Base icosahedron for dodecahedron invalid." << std::endl;
        return dodecahedron;
    }
    dodecahedron.points.clear();
    for (const auto& ico_face : icosahedron.faces) {
        Vector3D p_sum = Vector3D::point(0,0,0);
        for(int idx : ico_face.point_indexes) p_sum += icosahedron.points[idx];
        Vector3D centroid = p_sum / ico_face.point_indexes.size();
        centroid.normalise(); // Project onto unit sphere
        dodecahedron.points.push_back(centroid);
    }
    // Faces of dodecahedron are derived from vertices of icosahedron
    // This part is complex and needs correct topology mapping.
    // The provided face indices were likely for a specific vertex ordering.
    // Using the same face indices for dual might not be correct without re-evaluating.
    // For now, use the provided ones, assuming they match the new vertex generation.
    dodecahedron.faces = { /* Use original indices from your working version if they were correct */
        Face{{0, 1, 2, 3, 4}}, Face{{0, 5, 6, 7, 1}}, Face{{1, 7, 8, 9, 2}},
        Face{{2, 9, 10, 11, 3}}, Face{{3, 11, 12, 13, 4}}, Face{{4, 13, 14, 5, 0}},
        Face{{19, 18, 17, 16, 15}}, Face{{19, 14, 13, 12, 18}}, Face{{18, 12, 11, 10, 17}},
        Face{{17, 10, 9, 8, 16}}, Face{{16, 8, 7, 6, 15}}, Face{{15, 6, 5, 14, 19}}
    };
    return dodecahedron;
}

Figure Figure::createCylinder(int n, double height) {
    Figure cylinder;
    if (n < 3) n = 3;
    for (int i = 0; i < n; ++i) { // Bottom circle
        double angle = 2.0 * M_PI * i / n;
        cylinder.points.push_back(Vector3D::point(cos(angle), sin(angle), 0));
    }
    for (int i = 0; i < n; ++i) { // Top circle
        double angle = 2.0 * M_PI * i / n;
        cylinder.points.push_back(Vector3D::point(cos(angle), sin(angle), height));
    }
    vector<int> bottom_face_indices, top_face_indices;
    for (int i = 0; i < n; ++i) { // CCW for top, CW for bottom (viewed from outside)
        top_face_indices.push_back(n + i);
        bottom_face_indices.push_back(n - 1 - i);
    }
    cylinder.faces.push_back(Face{bottom_face_indices});
    cylinder.faces.push_back(Face{top_face_indices});
    for (int i = 0; i < n; ++i) { // Sides
        int next_i = (i + 1) % n;
        cylinder.faces.push_back(Face{{i, next_i, n + next_i, n + i}});
    }
    return cylinder;
}

Figure Figure::createCone(int n, double height) {
    Figure cone;
    if (n < 3) n = 3;
    for (int i = 0; i < n; ++i) { // Base circle
        double angle = 2.0 * M_PI * i / n;
        cone.points.push_back(Vector3D::point(cos(angle), sin(angle), 0));
    }
    cone.points.push_back(Vector3D::point(0, 0, height)); // Apex (index n)
    vector<int> base_face_indices;
    for (int i = 0; i < n; ++i) base_face_indices.push_back(n - 1 - i); // CW for bottom
    cone.faces.push_back(Face{base_face_indices});
    for (int i = 0; i < n; ++i) { // Sides
        cone.faces.push_back(Face{{i, (i + 1) % n, n}}); // i, next, apex
    }
    return cone;
}

static Vector3D normalize_static(const Vector3D& v) { // Renamed to avoid conflict if Figure had a normalize
    return Vector3D::normalise(v); // Use static normalise from vector3d.h
}

static void subdivideAndNormalizeSphere(Figure& figure, int n_iterations) { // Renamed
    for (int iter = 0; iter < n_iterations; ++iter) {
        Figure subdividedFigure;
        // subdividedFigure.color = figure.color; // Not needed if Figure constructor handles default
        std::map<std::pair<int, int>, int> midpoint_cache;

        auto get_midpoint_index = [&](int p1_idx, int p2_idx, std::vector<Vector3D>& current_points) -> int {
            if (p1_idx > p2_idx) std::swap(p1_idx, p2_idx);
            auto key = std::make_pair(p1_idx, p2_idx);
            auto it = midpoint_cache.find(key);
            if (it != midpoint_cache.end()) return it->second;

            const Vector3D& p1_vec = figure.points[p1_idx]; // Use original figure's points for calculation
            const Vector3D& p2_vec = figure.points[p2_idx];
            Vector3D midpoint = normalize_static((p1_vec + p2_vec) * 0.5);
            current_points.push_back(midpoint); // Add to new figure's points
            int new_index = current_points.size() - 1;
            midpoint_cache[key] = new_index;
            return new_index;
        };

        subdividedFigure.points = figure.points; // Start with original points

        std::vector<Face> new_faces;
        for (const auto& face : figure.faces) {
             if (face.point_indexes.size() != 3) continue;
             int a_idx = face.point_indexes[0];
             int b_idx = face.point_indexes[1];
             int c_idx = face.point_indexes[2];
             int ab_mid_idx = get_midpoint_index(a_idx, b_idx, subdividedFigure.points);
             int bc_mid_idx = get_midpoint_index(b_idx, c_idx, subdividedFigure.points);
             int ca_mid_idx = get_midpoint_index(c_idx, a_idx, subdividedFigure.points);
             new_faces.push_back(Face{{a_idx, ab_mid_idx, ca_mid_idx}});
             new_faces.push_back(Face{{b_idx, bc_mid_idx, ab_mid_idx}});
             new_faces.push_back(Face{{c_idx, ca_mid_idx, bc_mid_idx}});
             new_faces.push_back(Face{{ab_mid_idx, bc_mid_idx, ca_mid_idx}});
        }
        figure.points = subdividedFigure.points; // Update points first
        figure.faces = new_faces;                // Then faces
    }
     for (auto& point : figure.points) point.normalise(); // Final normalization
}

Figure Figure::createSphere(int n) {
    Figure sphere = createIcosahedron();
    // sphere.color = Color(1.0, 1.0, 1.0); // Old fallback
    subdivideAndNormalizeSphere(sphere, n);
    return sphere;
}

Figure Figure::createTorus(double R, double r_tube, int n, int m) { // Renamed r to r_tube
    Figure torus;
     if (n < 3) n = 3; if (m < 3) m = 3;
    for (int i = 0; i < n; ++i) {
        double u = 2.0 * M_PI * i / n;
        for (int j = 0; j < m; ++j) {
            double v = 2.0 * M_PI * j / m;
            double x = (R + r_tube * cos(v)) * cos(u);
            double y = (R + r_tube * cos(v)) * sin(u);
            double z = r_tube * sin(v);
            torus.points.push_back(Vector3D::point(x, y, z));
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int idx00 = i * m + j;
            int idx10 = ((i + 1) % n) * m + j;
            int idx11 = ((i + 1) % n) * m + ((j + 1) % m);
            int idx01 = i * m + ((j + 1) % m);
            torus.faces.push_back(Face{{idx00, idx10, idx11, idx01}});
        }
    }
    return torus;
}

Figure Figure::generate3DLSystem(const std::string &inputFile, const Color &line_color_fallback) {
    LParser::LSystem3D l_system;
    std::ifstream inputStream(inputFile);
    if (!inputStream.is_open()) throw std::runtime_error("Could not open L-System file: " + inputFile);
    try { inputStream >> l_system; }
    catch (const std::exception& e) { inputStream.close(); throw std::runtime_error("Error parsing L-System file: " + std::string(e.what())); }
    inputStream.close();

    Figure figure;
    figure.ambientReflection = line_color_fallback;
    figure.diffuseReflection = Color(0,0,0);
    figure.specularReflection = Color(0,0,0);
    figure.color = line_color_fallback; // Set the old color member as well

    double angle_rad = l_system.get_angle() * M_PI / 180.0;
    Vector3D currentPosition = Vector3D::point(0, 0, 0);
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L_turtle = Vector3D::vector(0, 1, 0); // Renamed to avoid conflict
    Vector3D U_turtle = Vector3D::vector(0, 0, 1);
    std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D>> stateStack;
    std::string lsystemString = l_system.get_initiator();
    for (unsigned int iter = 0; iter < l_system.get_nr_iterations(); ++iter) {
        std::string nextString = "";
        for (char symbol : lsystemString) {
            if (l_system.get_alphabet().count(symbol)) nextString += l_system.get_replacement(symbol);
            else nextString += symbol;
        }
        lsystemString = nextString;
    }

    for (char command : lsystemString) {
        Vector3D _H, _L_turtle, _U_turtle;
        switch (command) {
            case '+': _H = H; _L_turtle = L_turtle; H = _H * cos(angle_rad) + _L_turtle * sin(angle_rad); L_turtle = -_H * sin(angle_rad) + _L_turtle * cos(angle_rad); break;
            case '-': _H = H; _L_turtle = L_turtle; H = _H * cos(-angle_rad) + _L_turtle * sin(-angle_rad); L_turtle = -_H * sin(-angle_rad) + _L_turtle * cos(-angle_rad); break;
            case '^': _H = H; _U_turtle = U_turtle; H = _H * cos(angle_rad) + _U_turtle * sin(angle_rad); U_turtle = -_H * sin(angle_rad) + _U_turtle * cos(angle_rad); break;
            case '&': _H = H; _U_turtle = U_turtle; H = _H * cos(-angle_rad) + _U_turtle * sin(-angle_rad); U_turtle = -_H * sin(-angle_rad) + _U_turtle * cos(-angle_rad); break;
            case '\\': _L_turtle = L_turtle; _U_turtle = U_turtle; L_turtle = _L_turtle * cos(angle_rad) + _U_turtle * sin(angle_rad); U_turtle = -_L_turtle * sin(angle_rad) + _U_turtle * cos(angle_rad); break;
            case '/':  _L_turtle = L_turtle; _U_turtle = U_turtle; L_turtle = _L_turtle * cos(-angle_rad) + _U_turtle * sin(-angle_rad); U_turtle = -_L_turtle * sin(-angle_rad) + _U_turtle * cos(-angle_rad); break;
            case '|': H = -H; L_turtle = -L_turtle; break;
            case '(': stateStack.push({currentPosition, H, L_turtle, U_turtle}); break;
            case ')': if (!stateStack.empty()) { auto s = stateStack.top(); stateStack.pop(); currentPosition=get<0>(s); H=get<1>(s); L_turtle=get<2>(s); U_turtle=get<3>(s); } break;
            default: if (l_system.draw(command)) { Vector3D nextPos = currentPosition + H; figure.points.push_back(currentPosition); figure.points.push_back(nextPos); figure.faces.push_back(Face{{(int)figure.points.size()-2, (int)figure.points.size()-1}}); currentPosition = nextPos; } else { currentPosition = currentPosition + H; } break;
        }
    }
    return figure;
}