#include "lineDrawer.h"
#include "Figure.h"
#include "Transformations.h" // Includes Matrix, Vector3D ops
#include "Projection.h"      // For doProjection
#include "fractals3D.h"      // For fractal generators
#include "ini_configuration.h"
#include "Line2D.h"
#include <cmath>
#include <limits>
#include <stack>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm> // for std::max, std::min
#include <map>       // for std::map

#ifndef Pi
#define Pi 3.14159265358979323846
#endif

using namespace std;

// Calculate Min/Max for Triangles
map<string, double> calculateMinMax(const Triangles &triangles) {
    map<string, double> resultaat;
    if (triangles.empty()) {
        // Provide default bounds if empty to avoid issues in calculateParameters
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

    resultaat["xmin"] = xmin;
    resultaat["xmax"] = xmax;
    resultaat["ymin"] = ymin;
    resultaat["ymax"] = ymax;
    return resultaat;
}

// Calculate Min/Max for Lines
map<string, double> calculateMinMax(const Lines2D &lines) {
     map<string, double> resultaat;
     if (lines.empty()) {
         // Provide default bounds if empty
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

     resultaat["xmin"] = xmin;
     resultaat["xmax"] = xmax;
     resultaat["ymin"] = ymin;
     resultaat["ymax"] = ymax;
     return resultaat;
}

// Calculate Scaling/Translation Parameters (Unified Logic)
// Takes a map containing min/max values and calculates parameters
std::map<std::string, double> calculateParameters(const std::map<std::string, double>& minMax, const int size) {
    std::map<std::string, double> resultaten = minMax; // Copy min/max values

    double xmin = minMax.at("xmin");
    double xmax = minMax.at("xmax");
    double ymin = minMax.at("ymin");
    double ymax = minMax.at("ymax");

    double Xrange = xmax - xmin;
    double Yrange = ymax - ymin;

    // Handle zero range cases
    const double minRange = 1e-9; // Use a small epsilon for comparison
    if (Xrange < minRange) Xrange = 1.0;
    if (Yrange < minRange) Yrange = 1.0;


    resultaten["Xrange"] = Xrange;
    resultaten["Yrange"] = Yrange;

    // Calculate image dimensions based on aspect ratio, round to integer
    // Ensure size is positive before using it
    int safeSize = std::max(1, size);
    int imageX = std::lround(safeSize * (Xrange / std::max(Xrange, Yrange)));
    int imageY = std::lround(safeSize * (Yrange / std::max(Xrange, Yrange)));

    // Ensure minimum size of 1x1
    imageX = std::max(1, imageX);
    imageY = std::max(1, imageY);

    resultaten["imageX"] = static_cast<double>(imageX); // Store as double, but derived from int
    resultaten["imageY"] = static_cast<double>(imageY);

    // Calculate scale factor 'd' using the integer image width
    // Use Xrange here as per friend's code logic `d = 0.95*(imagex/Xrange)`
    double d = 0.95 * (static_cast<double>(imageX) / Xrange);
    resultaten["d"] = d;

    // Calculate center of drawing in original coordinates, scaled
    double DCx = d * (xmin + xmax) / 2.0;
    double DCy = d * (ymin + ymax) / 2.0;
    resultaten["DCx"] = DCx;
    resultaten["DCy"] = DCy;

    // Calculate translation needed to center the drawing in the image (using integer dimensions)
    double dx = (static_cast<double>(imageX) / 2.0) - DCx;
    double dy = (static_cast<double>(imageY) / 2.0) - DCy;
    resultaten["dx"] = dx;
    resultaten["dy"] = dy;

    return resultaten;
}

// Overload for Lines2D
std::map<std::string, double> calculate(const Lines2D &lines, const int size) {
    auto minMax = calculateMinMax(lines); // Will return defaults if lines is empty
    return calculateParameters(minMax, size);
}

// Overload for Triangles
std::map<std::string, double> calculate(const Triangles &triangles, const int size) {
    auto minMax = calculateMinMax(triangles); // Will return defaults if triangles is empty
    return calculateParameters(minMax, size);
}


// Draw 2D Lines (Wireframe mode)
img::EasyImage draw2DLines(Lines2D &lines, int size, const Color &backgroundColor) {
    // Calculate parameters *before* checking if lines is empty,
    // as calculate handles empty inputs now.
    auto resultaten = calculate(lines, size);

    int imageWidth = static_cast<int>(resultaten["imageX"]);
    int imageHeight = static_cast<int>(resultaten["imageY"]);

    // Create the image (dimensions are guaranteed to be >= 1)
    img::EasyImage image(imageWidth, imageHeight, img::Color(backgroundColor.red * 255, backgroundColor.green * 255, backgroundColor.blue * 255));

    // Check if there are any lines to draw *after* creating the image
    if (lines.empty()) {
        std::cerr << "Warning: draw2DLines called with no lines. Returning blank image." << std::endl;
        return image; // Return the blank image
    }

    double d = resultaten["d"];
    double dx = resultaten["dx"];
    double dy = resultaten["dy"];

    // Apply scaling and translation *to the line coordinates*
    for (auto &line : lines) {
        // Scale points
        line.p1.x *= d; line.p1.y *= d;
        line.p2.x *= d; line.p2.y *= d;
        // Translate points
        line.p1.x += dx; line.p1.y += dy;
        line.p2.x += dx; line.p2.y += dy;
    }

    // Draw the transformed lines onto the image
    for (const auto &line : lines) {
        img::Color drawColor(line.color.red * 255, line.color.green * 255, line.color.blue * 255);
        image.draw_line(
            std::lround(line.p1.x), std::lround(line.p1.y),
            std::lround(line.p2.x), std::lround(line.p2.y),
            drawColor
        );
    }

    return image;
}

// Draw L-System
Lines2D drawLSystem(const LParser::LSystem2D &l_system, const Color &lineColor) {
    Lines2D lines;
    std::stack<std::pair<Point2D, double>> stateStack; // Pair: {position, angle_radians}

    Point2D currentPosition(0.0, 0.0);
    double currentAngle_rad = l_system.get_starting_angle() * Pi / 180.0; // Start angle in radians
    double angle_change_rad = l_system.get_angle() * Pi / 180.0; // Angle change per +/- command

    // Generate the full instruction string
    std::string instructions = l_system.get_initiator();
    std::set<char> alphabet = l_system.get_alphabet();
    bool ignoreUnknown = true; // Standard behavior

    for (int i = 0; i < l_system.get_nr_iterations(); ++i) {
        std::string nextInstructions;
        for (char c : instructions) {
            if (alphabet.count(c)) { // If char is in alphabet, replace it
                nextInstructions += l_system.get_replacement(c);
            } else { // Keep commands and potentially other symbols
                 if (string("+-()").find(c) != string::npos || !ignoreUnknown || l_system.draw(c)) {
                    nextInstructions += c;
                 } // Else: ignore silently
            }
        }
        instructions = nextInstructions;
    }

    // Process the final instruction string
    for (char c : instructions) {
        if (c == '+') {
            currentAngle_rad += angle_change_rad;
        } else if (c == '-') {
            currentAngle_rad -= angle_change_rad;
        } else if (c == '(') {
            stateStack.push({currentPosition, currentAngle_rad});
        } else if (c == ')') {
            if (!stateStack.empty()) {
                currentPosition = stateStack.top().first;
                currentAngle_rad = stateStack.top().second;
                stateStack.pop();
            } // Else: stack underflow, ignore
        } else { // Could be a symbol to draw or just move
            // Check if the L-System definition says to draw this symbol
            bool shouldDraw = l_system.draw(c);

            // Calculate next position (always move forward by step size 1)
            Point2D nextPosition = currentPosition;
            nextPosition.x += std::cos(currentAngle_rad); // Assume step size is 1
            nextPosition.y += std::sin(currentAngle_rad);

            if (shouldDraw) {
                // Add line segment if drawing is required
                lines.push_back(Line2D(currentPosition, nextPosition, 0, 0, lineColor));
            }
            // Update position regardless of drawing (turtle moves)
            currentPosition = nextPosition;
        }
    }
    return lines;
}


// Generate 3D figures based on the configuration
Figures3D generateFigures(const ini::Configuration &configuration) {
    Figures3D figures;
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    for (int i = 0; i < nrFigures; ++i) {
        std::string fig_key = "Figure" + std::to_string(i);
        // Check if figure section exists before accessing it
        // Check if section exists with try-catch
        bool figureExists = true;
        try {
            configuration[fig_key];
        } catch (const std::exception &e) {
            figureExists = false;
        }

        if (!figureExists) {
            std::cerr << "Warning: Figure section '" << fig_key << "' not found in configuration." << std::endl;
            continue;
        }

        std::string type = configuration[fig_key]["type"].as_string_or_die();
        //std::cout << "Generating figure " << i << ": " << type << std::endl; // Debug
        Figure figure; // Create figure for this index

        // Read color first, as it might be needed by generators or applied later
        ini::DoubleTuple colorData = configuration[fig_key]["color"].as_double_tuple_or_die();
        Color color(colorData[0], colorData[1], colorData[2]);

        bool figureGenerated = true; // Flag to track if figure was successfully generated

        // Handle different figure types
        if (type == "LineDrawing") {
            int nrPoints = configuration[fig_key]["nrPoints"].as_int_or_die();
            int nrLines = configuration[fig_key]["nrLines"].as_int_or_die();
            // Read points
            for (int j = 0; j < nrPoints; ++j) {
                ini::DoubleTuple pointData = configuration[fig_key]["point" + std::to_string(j)].as_double_tuple_or_die();
                figure.points.push_back(Vector3D::point(pointData[0], pointData[1], pointData[2]));
            }
            // Read lines (as faces with 2 points)
            for (int j = 0; j < nrLines; ++j) {
                ini::IntTuple lineData = configuration[fig_key]["line" + std::to_string(j)].as_int_tuple_or_die();
                // Basic validation for line indices
                if (lineData.size() >= 2 && lineData[0] >= 0 && lineData[0] < nrPoints && lineData[1] >= 0 && lineData[1] < nrPoints) {
                     figure.faces.push_back(Face{{lineData[0], lineData[1]}});
                } else {
                     std::cerr << "Warning: Invalid line indices for line" << j << " in " << fig_key << std::endl;
                }
            }
        } else if (type == "3DLSystem") {
            std::string inputFile = configuration[fig_key]["inputfile"].as_string_or_die();
            try {
                 figure = Figure::generate3DLSystem(inputFile, color); // Color is passed here
            } catch (const std::runtime_error& e) {
                 std::cerr << "Error generating 3DLSystem for " << fig_key << ": " << e.what() << std::endl;
                 figureGenerated = false;
            }
        } else if (type == "Cube") { figure = Figure::createCube(); }
          else if (type == "Tetrahedron") { figure = Figure::createTetrahedron(); }
          else if (type == "Octahedron") { figure = Figure::createOctahedron(); }
          else if (type == "Icosahedron") { figure = Figure::createIcosahedron(); }
          else if (type == "Dodecahedron") { figure = Figure::createDodecahedron(); }
          else if (type == "Cylinder") {
              int n = configuration[fig_key]["n"].as_int_or_die();
              double h = configuration[fig_key]["height"].as_double_or_die();
              figure = Figure::createCylinder(n, h);
          } else if (type == "Cone") {
              int n = configuration[fig_key]["n"].as_int_or_die();
              double h = configuration[fig_key]["height"].as_double_or_die();
              figure = Figure::createCone(n, h);
          } else if (type == "Sphere") {
              int n = configuration[fig_key]["n"].as_int_or_die();
              figure = Figure::createSphere(n);
          } else if (type == "Torus") {
              double R = configuration[fig_key]["R"].as_double_or_die();
              double r = configuration[fig_key]["r"].as_double_or_die();
              int n = configuration[fig_key]["n"].as_int_or_die();
              int m = configuration[fig_key]["m"].as_int_or_die();
              figure = Figure::createTorus(R, r, n, m);
          }
        // --- Basic Shapes & L-Systems above ---
        // --- Fractals and Buckyballs below ---
          else if (type == "MengerSponge") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              figure = generateMengerSponge(nrIterations); // Uses fixed 1/3 scale internally
          } else if (type == "FractalCube") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
              double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; // Avoid division by zero
              figure = generateFractalCube(nrIterations, actualScaleFactor);
          } else if (type == "FractalTetrahedron") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
              double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0;
              figure = generateFractalTetrahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalOctahedron") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
              double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0;
              figure = generateFractalOctahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalIcosahedron") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
              double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0;
              figure = generateFractalIcosahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalDodecahedron") {
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
              double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
              double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0;
              figure = generateFractalDodecahedron(nrIterations, actualScaleFactor);
          } else if (type == "BuckyBall") {
                figure = Figure::createBuckyBall(); // Using placeholder for now
          } else if (type == "FractalBuckyBall") {
                int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die();
                double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die();
                double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0;
                figure = generateFractalBuckyBall(nrIterations, actualScaleFactor); // Using placeholder
          }
          // --- Unknown Type ---
          else {
               std::cerr << "Warning: Unknown figure type '" << type << "' specified for " << fig_key << std::endl;
               figureGenerated = false; // Mark as not generated
          }

        // If figure generation failed or type was unknown, skip transformations and adding it
        if (!figureGenerated) {
            continue;
        }

        // Assign the color to the generated figure (overwrites default or LSystem color if applicable)
        figure.color = color;

        // Apply geometric transformations (scale, rotate, translate) *after* generation
        double scale = configuration[fig_key]["scale"].as_double_or_die();
        double rotateXAngle = configuration[fig_key]["rotateX"].as_double_or_die();
        double rotateYAngle = configuration[fig_key]["rotateY"].as_double_or_die();
        double rotateZAngle = configuration[fig_key]["rotateZ"].as_double_or_die();
        ini::DoubleTuple centerData = configuration[fig_key]["center"].as_double_tuple_or_die();
        Vector3D centerVec = Vector3D::point(centerData[0], centerData[1], centerData[2]);

        // Create transformation matrices
        Matrix scaleMatrix = scaleFigure(scale);
        Matrix rotateXMatrix = rotateX(rotateXAngle);
        Matrix rotateYMatrix = rotateY(rotateYAngle);
        Matrix rotateZMatrix = rotateZ(rotateZAngle);
        Matrix translateMatrix = translate(centerVec);

        // Combine transformations: Scale -> RotateX -> RotateY -> RotateZ -> Translate
        // Order matters! Typically scale first, then rotate, then translate.
        Matrix transformationMatrix = scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * translateMatrix;

        // Apply the combined transformation to the figure
        applyTransformation(figure, transformationMatrix);

        figures.push_back(figure); // Add the fully generated and transformed figure
        //std::cout << "Figure " << i << " (" << type << ") added." << std::endl; // Debug
    }
    //std::cout << "Total figures generated: " << figures.size() << std::endl; // Debug
    return figures;
}