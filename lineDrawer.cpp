#include "lineDrawer.h"
#include "Figure.h"
#include "Transformations.h"
#include "Projection.h"
#include "fractals3D.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include "Light.h"           // Added for light generation
#include <cmath>
#include <limits>            // Added for numeric_limits
#include <stack>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#ifndef M_PI // Define M_PI if not already defined
#define M_PI 3.14159265358979323846
#endif

using namespace std;

map<string, double> calculateMinMax(const Triangles &triangles) {
    map<string, double> resultaat;
    if (triangles.empty()) {
        resultaat["xmin"] = -1.0; resultaat["xmax"] = 1.0;
        resultaat["ymin"] = -1.0; resultaat["ymax"] = 1.0;
        return resultaat;
    }
    double xmin = numeric_limits<double>::max();
    double xmax = numeric_limits<double>::lowest();
    double ymin = numeric_limits<double>::max();
    double ymax = numeric_limits<double>::lowest();
    for (const Triangle &triangle : triangles) {
        xmin = std::min({xmin, triangle.p1.x, triangle.p2.x, triangle.p3.x});
        xmax = std::max({xmax, triangle.p1.x, triangle.p2.x, triangle.p3.x});
        ymin = std::min({ymin, triangle.p1.y, triangle.p2.y, triangle.p3.y});
        ymax = std::max({ymax, triangle.p1.y, triangle.p2.y, triangle.p3.y});
    }
    resultaat["xmin"] = xmin; resultaat["xmax"] = xmax;
    resultaat["ymin"] = ymin; resultaat["ymax"] = ymax;
    return resultaat;
}

map<string, double> calculateMinMax(const Lines2D &lines) {
     map<string, double> resultaat;
     if (lines.empty()) {
         resultaat["xmin"] = -1.0; resultaat["xmax"] = 1.0;
         resultaat["ymin"] = -1.0; resultaat["ymax"] = 1.0;
         return resultaat;
     }
     double xmin = numeric_limits<double>::max();
     double xmax = numeric_limits<double>::lowest();
     double ymin = numeric_limits<double>::max();
     double ymax = numeric_limits<double>::lowest();
     for (const Line2D &line : lines) {
         xmin = std::min({xmin, line.p1.x, line.p2.x});
         xmax = std::max({xmax, line.p1.x, line.p2.x});
         ymin = std::min({ymin, line.p1.y, line.p2.y});
         ymax = std::max({ymax, line.p1.y, line.p2.y});
     }
     resultaat["xmin"] = xmin; resultaat["xmax"] = xmax;
     resultaat["ymin"] = ymin; resultaat["ymax"] = ymax;
     return resultaat;
}

std::map<std::string, double> calculateParametersStrict(const std::map<std::string, double>& minMax, const int size) {
    std::map<std::string, double> resultaten = minMax;
    double xmin = minMax.at("xmin");
    double xmax = minMax.at("xmax");
    double ymin = minMax.at("ymin");
    double ymax = minMax.at("ymax");
    double Xrange = std::abs(xmax - xmin);
    double Yrange = std::abs(ymax - ymin);
    const double minRange = 1e-9;
    if (Xrange < minRange) Xrange = 1.0;
    if (Yrange < minRange) Yrange = 1.0;
    resultaten["Xrange"] = Xrange;
    resultaten["Yrange"] = Yrange;
    int safeSize = std::max(1, size);
    int imageX_calc = static_cast<int>(static_cast<double>(safeSize) * (Xrange / std::max(Xrange, Yrange)));
    int imageY_calc = static_cast<int>(static_cast<double>(safeSize) * (Yrange / std::max(Xrange, Yrange)));
    imageX_calc = std::max(1, imageX_calc);
    imageY_calc = std::max(1, imageY_calc);
    resultaten["imageX"] = static_cast<double>(imageX_calc);
    resultaten["imageY"] = static_cast<double>(imageY_calc);
    double d_2D_scaling = 0.95 * (static_cast<double>(imageX_calc) / Xrange);
    resultaten["d"] = d_2D_scaling;
    double DCx = d_2D_scaling * (xmin + xmax) / 2.0;
    double DCy = d_2D_scaling * (ymin + ymax) / 2.0;
    resultaten["DCx"] = DCx;
    resultaten["DCy"] = DCy;
    double dx_trans = (static_cast<double>(imageX_calc) / 2.0) - DCx;
    double dy_trans = (static_cast<double>(imageY_calc) / 2.0) - DCy;
    resultaten["dx"] = dx_trans;
    resultaten["dy"] = dy_trans;
    return resultaten;
}

std::map<std::string, double> calculate(const Lines2D &lines, const int size) {
    auto minMax = calculateMinMax(lines);
    return calculateParametersStrict(minMax, size);
}

std::map<std::string, double> calculate(const Triangles &triangles, const int size) {
    auto minMax = calculateMinMax(triangles);
    return calculateParametersStrict(minMax, size);
}

img::EasyImage draw2DLines(Lines2D &lines, int size, const Color &backgroundColor) {
    auto resultaten = calculate(lines, size);
    int finalImageWidth = static_cast<int>(resultaten.at("imageX"));
    int finalImageHeight = static_cast<int>(resultaten.at("imageY"));
    if (finalImageWidth <= 0 || finalImageHeight <= 0) {
         finalImageWidth = std::max(1, size);
         finalImageHeight = std::max(1, size);
    }
    img::EasyImage image(finalImageWidth, finalImageHeight, img::Color(static_cast<uint8_t>(backgroundColor.red * 255), static_cast<uint8_t>(backgroundColor.green * 255), static_cast<uint8_t>(backgroundColor.blue * 255)));
    if (lines.empty()) return image;

    double d = resultaten.at("d");
    double dx = resultaten.at("dx");
    double dy = resultaten.at("dy");
    for (auto &line : lines) {
        line.p1.x = line.p1.x * d + dx; line.p1.y = line.p1.y * d + dy;
        line.p2.x = line.p2.x * d + dx; line.p2.y = line.p2.y * d + dy;
    }
    for (const auto &line : lines) {
        image.draw_line(
            std::lround(line.p1.x), std::lround(line.p1.y),
            std::lround(line.p2.x), std::lround(line.p2.y),
            img::Color(static_cast<uint8_t>(std::round(line.color.red * 255)), static_cast<uint8_t>(std::round(line.color.green * 255)), static_cast<uint8_t>(std::round(line.color.blue * 255)))
        );
    }
    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system, const Color &lineColor) {
    Lines2D lines;
    std::stack<std::pair<Point2D, double>> stateStack;
    Point2D currentPosition(0.0, 0.0);
    double currentAngle_rad = l_system.get_starting_angle() * M_PI / 180.0;
    double angle_change_rad = l_system.get_angle() * M_PI / 180.0;
    std::string instructions = l_system.get_initiator();
    for (int i = 0; i < l_system.get_nr_iterations(); ++i) {
        std::string nextInstructions;
        for (char c : instructions) {
            if (l_system.get_alphabet().count(c)) nextInstructions += l_system.get_replacement(c);
            else nextInstructions += c;
        }
        instructions = nextInstructions;
    }
    for (char c : instructions) {
        if (c == '+') currentAngle_rad += angle_change_rad;
        else if (c == '-') currentAngle_rad -= angle_change_rad;
        else if (c == '(') stateStack.push({currentPosition, currentAngle_rad});
        else if (c == ')') { if (!stateStack.empty()) { currentPosition = stateStack.top().first; currentAngle_rad = stateStack.top().second; stateStack.pop(); }}
        else {
            Point2D nextPosition = currentPosition;
            nextPosition.x += std::cos(currentAngle_rad);
            nextPosition.y += std::sin(currentAngle_rad);
            if (l_system.draw(c)) lines.push_back(Line2D(currentPosition, nextPosition, 0,0, lineColor));
            currentPosition = nextPosition;
        }
    }
    return lines;
}

Lights3D generateLights(const ini::Configuration &configuration, const Matrix& viewMatrix) {
    Lights3D lights_list;
    int nrLights = 0;
    if (configuration["General"]["nrLights"].exists()) {
        nrLights = configuration["General"]["nrLights"].as_int_or_die();
    }

    for (int i = 0; i < nrLights; ++i) {
        std::string light_key = "Light" + std::to_string(i);
        const auto& lightConfig = configuration[light_key];

        // Basic check if section has any of the core light color keys
        if (!lightConfig["ambientLight"].exists() &&
            !lightConfig["diffuseLight"].exists() &&
            !lightConfig["specularLight"].exists()) {
            std::cerr << "Warning: Light section '" << light_key << "' missing core color fields (ambient/diffuse/specular)." << std::endl;
            continue;
        }

        Color ambientL(0,0,0), diffuseL(0,0,0), specularL(0,0,0);
        if (lightConfig["ambientLight"].exists()) {
            ini::DoubleTuple ambTuple = lightConfig["ambientLight"].as_double_tuple_or_die();
            ambientL = Color(ambTuple[0], ambTuple[1], ambTuple[2]);
        }
        if (lightConfig["diffuseLight"].exists()) {
            ini::DoubleTuple diffTuple = lightConfig["diffuseLight"].as_double_tuple_or_die();
            diffuseL = Color(diffTuple[0], diffTuple[1], diffTuple[2]);
        }
        if (lightConfig["specularLight"].exists()) {
            ini::DoubleTuple specTuple = lightConfig["specularLight"].as_double_tuple_or_die();
            specularL = Color(specTuple[0], specTuple[1], specTuple[2]);
        }

        bool isInfinite = false;
        if (lightConfig["infinity"].exists()) {
            isInfinite = lightConfig["infinity"].as_bool_or_die();
        }

        if (isInfinite) {
            DirectionalLight* dl = new DirectionalLight();
            dl->ambientLight = ambientL; dl->diffuseLight = diffuseL; dl->specularLight = specularL;
            if (lightConfig["direction"].exists()) {
                ini::DoubleTuple dirTuple = lightConfig["direction"].as_double_tuple_or_die();
                // dl->direction is Ld (from light to origin) in world space. Transform to eye space.
                // A direction vector is transformed by the main 3x3 part of the matrix.
                // Assuming viewMatrix has translation in 4th row/col, which is fine for points.
                // For direction vectors (w=0), translation part of matrix mult doesn't apply.
                dl->direction = Vector3D::vector(dirTuple[0], dirTuple[1], dirTuple[2]) * viewMatrix;
            } else {
                 dl->direction = Vector3D::vector(0,0,-1.0) * viewMatrix; // Default world direction transformed
            }
            dl->direction.normalise(); // making sure its normalized after transformation
            lights_list.push_back(dl);
        } else { // Point Light
            PointLight* pl = new PointLight();
            pl->ambientLight = ambientL; pl->diffuseLight = diffuseL; pl->specularLight = specularL;
            if (lightConfig["location"].exists()) {
                ini::DoubleTuple locTuple = lightConfig["location"].as_double_tuple_or_die();
                // Transform point light location from world to eye space
                pl->location = Vector3D::point(locTuple[0], locTuple[1], locTuple[2]) * viewMatrix;
            } else {
                 pl->location = Vector3D::point(0,0,0) * viewMatrix; // Default world location transformed
            }

            if (lightConfig["spotAngle"].exists()) {
                pl->spotAngleDegrees = lightConfig["spotAngle"].as_double_or_die();
                // Spot direction also needs to be in eye space.
                // Default: from light location (now in eye space) towards eye-space origin (0,0,0).
                if (pl->location.length() > 1e-9) {
                    pl->spotDirection = Vector3D::normalise(Vector3D::point(0,0,0) - pl->location);
                } else {
                    // If light is at eye-space origin, default spot direction (e.g., along -Z eye)
                    pl->spotDirection = Vector3D::vector(0,0,-1.0);
                    // This default spot direction is already in eye-space, no further transform needed
                }
            } else {
                pl->spotAngleDegrees = 181.0; // Mark as omni, spotDirection won't be used
            }
            lights_list.push_back(pl);
        }
    }
    return lights_list;
}

Figures3D generateFigures(const ini::Configuration &configuration) {
    Figures3D figures_list;
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    for (int i = 0; i < nrFigures; ++i) {
        std::string fig_key = "Figure" + std::to_string(i);
        if (!configuration[fig_key]["type"].exists()) {
            std::cerr << "Warning: Figure section '" << fig_key << "' or its 'type' not found." << std::endl;
            continue;
        }
        const auto& figureConfig = configuration[fig_key];
        std::string type = figureConfig["type"].as_string_or_die();

        Figure currentFigure; // Uses default constructor for materials

        bool newMaterialFormat = false;
        if (figureConfig["ambientReflection"].exists()) {
            newMaterialFormat = true;
            ini::DoubleTuple ar = figureConfig["ambientReflection"].as_double_tuple_or_die();
            currentFigure.ambientReflection = Color(ar[0], ar[1], ar[2]);
        }
        if (figureConfig["diffuseReflection"].exists()) {
            newMaterialFormat = true;
            ini::DoubleTuple dr = figureConfig["diffuseReflection"].as_double_tuple_or_die();
            currentFigure.diffuseReflection = Color(dr[0], dr[1], dr[2]);
        }
        if (figureConfig["specularReflection"].exists()) {
            newMaterialFormat = true;
            ini::DoubleTuple sr = figureConfig["specularReflection"].as_double_tuple_or_die();
            currentFigure.specularReflection = Color(sr[0], sr[1], sr[2]);
        }
        if (figureConfig["reflectionCoefficient"].exists()) {
            newMaterialFormat = true;
            currentFigure.reflectionCoefficient = figureConfig["reflectionCoefficient"].as_double_or_die();
        }

        if (!newMaterialFormat && figureConfig["color"].exists()) {
            ini::DoubleTuple c = figureConfig["color"].as_double_tuple_or_die();
            currentFigure.ambientReflection = Color(c[0], c[1], c[2]);
            currentFigure.color = Color(c[0], c[1], c[2]); // Set old color member
            currentFigure.diffuseReflection = Color(0,0,0);
            currentFigure.specularReflection = Color(0,0,0);
            currentFigure.reflectionCoefficient = 1.0;
        } else if (newMaterialFormat) {
            // If new format, set the old 'color' member to something reasonable, e.g., diffuse
            currentFigure.color = currentFigure.diffuseReflection;
        }
        // If no color info at all, Figure constructor defaults are used for both new and old color.


        bool figureGenerated = true;
        Figure geomOnlyFigure; // To hold geometry from static creators

        if (type == "LineDrawing") {
            int nrPoints = figureConfig["nrPoints"].as_int_or_die();
            int nrLines = figureConfig["nrLines"].as_int_or_die();
            for (int j = 0; j < nrPoints; ++j) currentFigure.points.push_back(Vector3D::point(figureConfig["point" + std::to_string(j)].as_double_tuple_or_die()[0], figureConfig["point" + std::to_string(j)].as_double_tuple_or_die()[1],figureConfig["point" + std::to_string(j)].as_double_tuple_or_die()[2]));
            for (int j = 0; j < nrLines; ++j) { ini::IntTuple ld = figureConfig["line"+std::to_string(j)].as_int_tuple_or_die(); if(ld.size()>=2) currentFigure.faces.push_back(Face{{ld[0], ld[1]}}); }
        } else if (type == "3DLSystem") {
            std::string inputFile = figureConfig["inputfile"].as_string_or_die();
            try { geomOnlyFigure = Figure::generate3DLSystem(inputFile, currentFigure.ambientReflection); } // Pass a fallback
            catch (const std::runtime_error& e) { std::cerr << "Error 3DLSystem: " << e.what() << std::endl; figureGenerated = false; }
             if(figureGenerated) { currentFigure.points = geomOnlyFigure.points; currentFigure.faces = geomOnlyFigure.faces; }
        } else if (type == "Cube") { geomOnlyFigure = Figure::createCube(); }
          else if (type == "Tetrahedron") { geomOnlyFigure = Figure::createTetrahedron(); }
          else if (type == "Octahedron") { geomOnlyFigure = Figure::createOctahedron(); }
          else if (type == "Icosahedron") { geomOnlyFigure = Figure::createIcosahedron(); }
          else if (type == "Dodecahedron") { geomOnlyFigure = Figure::createDodecahedron(); }
          else if (type == "Cylinder") { geomOnlyFigure = Figure::createCylinder(figureConfig["n"].as_int_or_die(), figureConfig["height"].as_double_or_die()); }
          else if (type == "Cone") { geomOnlyFigure = Figure::createCone(figureConfig["n"].as_int_or_die(), figureConfig["height"].as_double_or_die()); }
          else if (type == "Sphere") { geomOnlyFigure = Figure::createSphere(figureConfig["n"].as_int_or_die()); }
          else if (type == "Torus") { geomOnlyFigure = Figure::createTorus(figureConfig["R"].as_double_or_die(), figureConfig["r"].as_double_or_die(), figureConfig["n"].as_int_or_die(), figureConfig["m"].as_int_or_die()); }
          else if (type == "MengerSponge") { geomOnlyFigure = generateMengerSponge(figureConfig["nrIterations"].as_int_or_die());}
          else if (type == "FractalCube") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalCube(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else if (type == "FractalTetrahedron") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalTetrahedron(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else if (type == "FractalOctahedron") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalOctahedron(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else if (type == "FractalIcosahedron") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalIcosahedron(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else if (type == "FractalDodecahedron") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalDodecahedron(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else if (type == "BuckyBall") { geomOnlyFigure = Figure::createBuckyBall(); }
          else if (type == "FractalBuckyBall") { double scale = figureConfig["fractalScale"].as_double_or_die(); geomOnlyFigure = generateFractalBuckyBall(figureConfig["nrIterations"].as_int_or_die(), (scale > 1e-9 ? 1.0/scale : 1.0));}
          else { std::cerr << "Unknown figure type: " << type << std::endl; figureGenerated = false; }

        if (!figureGenerated) continue;

        // For types other than LineDrawing and 3DLSystem (which populate currentFigure directly or via geomOnlyFigure already)
        if (type != "LineDrawing" && type != "3DLSystem") {
            currentFigure.points = geomOnlyFigure.points;
            currentFigure.faces = geomOnlyFigure.faces;
            // currentFigure.color = geomOnlyFigure.color; // If geomOnlyFigure has a meaningful fallback color
        }

        Matrix scaleMatrix = scaleFigure(figureConfig["scale"].as_double_or_die());
        Matrix rotateXMatrix = rotateX(figureConfig["rotateX"].as_double_or_die());
        Matrix rotateYMatrix = rotateY(figureConfig["rotateY"].as_double_or_die());
        Matrix rotateZMatrix = rotateZ(figureConfig["rotateZ"].as_double_or_die());
        ini::DoubleTuple cd = figureConfig["center"].as_double_tuple_or_die();
        Matrix translateMatrix = translate(Vector3D::point(cd[0], cd[1], cd[2]));
        Matrix transformationMatrix = scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * translateMatrix;
        applyTransformation(currentFigure, transformationMatrix);

        figures_list.push_back(currentFigure);
    }
    return figures_list;
}