#include "lineDrawer.h"
#include "Figure.h"
#include "Transformations.h" // Includes Matrix, Vector3D ops
#include "Projection.h"      // For doProjection
#include "fractals3D.h"      // For fractal generators
#include "ini_configuration.h"
#include "Line2D.h"
#include <cmath>     // For std::round, std::lround, std::abs, std::max
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

// Calculate Min/Max for Triangles (Keep as is)
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

// Calculate Min/Max for Lines (Keep as is)
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
    //std::cout << "--- calculateParametersStrict ---" << std::endl; // Marker
    std::map<std::string, double> resultaten = minMax; // Copy min/max values

    double xmin = minMax.at("xmin");
    double xmax = minMax.at("xmax");
    double ymin = minMax.at("ymin");
    double ymax = minMax.at("ymax");

    //std::cout << "Input minMax: xmin=" << xmin << ", xmax=" << xmax << ", ymin=" << ymin << ", ymax=" << ymax << std::endl;

    double Xrange = std::abs(xmax - xmin);
    double Yrange = std::abs(ymax - ymin);

    // Handle zero range cases
    const double minRange = 1e-9;
    if (Xrange < minRange) Xrange = 1.0;
    if (Yrange < minRange) Yrange = 1.0;

    resultaten["Xrange"] = Xrange;
    resultaten["Yrange"] = Yrange;
    //std::cout << "Calculated Xrange=" << Xrange << ", Yrange=" << Yrange << std::endl;


    // Calculate image dimensions as INTEGERS first, using truncation (implicit cast)
    int safeSize = std::max(1, size);
    int imageX_calc = static_cast<int>(static_cast<double>(safeSize) * (Xrange / std::max(Xrange, Yrange)));
    int imageY_calc = static_cast<int>(static_cast<double>(safeSize) * (Yrange / std::max(Xrange, Yrange)));

    // Ensure minimum size of 1x1
    imageX_calc = std::max(1, imageX_calc);
    imageY_calc = std::max(1, imageY_calc);

    //std::cout << "Calculated canvas (before storing in map): imageX_calc=" << imageX_calc << ", imageY_calc=" << imageY_calc << std::endl;

    // Store the INTEGER dimensions (cast to double for map storage if needed, but track the int source)
    resultaten["imageX"] = static_cast<double>(imageX_calc);
    resultaten["imageY"] = static_cast<double>(imageY_calc);

    // Calculate scale factor 'd' using the INTEGER image width and Xrange
    // Cast imageX_calc to double for the division
    double d_2D_scaling = 0.95 * (static_cast<double>(imageX_calc) / Xrange);
    resultaten["d"] = d_2D_scaling;
    //std::cout << "Calculated 2D scaling factor d=" << d_2D_scaling << std::endl;


    // Calculate center of drawing in original coordinates, scaled by 'd'
    double DCx = d_2D_scaling * (xmin + xmax) / 2.0;
    double DCy = d_2D_scaling * (ymin + ymax) / 2.0;
    resultaten["DCx"] = DCx;
    resultaten["DCy"] = DCy;
    //std::cout << "Scaled drawing center DCx=" << DCx << ", DCy=" << DCy << std::endl;

    // Calculate translation needed using INTEGER image dimensions (cast to double)
    double dx_trans = (static_cast<double>(imageX_calc) / 2.0) - DCx;
    double dy_trans = (static_cast<double>(imageY_calc) / 2.0) - DCy;
    resultaten["dx"] = dx_trans;
    resultaten["dy"] = dy_trans;
    //std::cout << "Calculated 2D translation dx=" << dx_trans << ", dy=" << dy_trans << std::endl;
    //std::cout << "--- End calculateParametersStrict ---" << std::endl;

    return resultaten;
}

// Overload for Lines2D - Calls new strict parameter function
std::map<std::string, double> calculate(const Lines2D &lines, const int size) {
    auto minMax = calculateMinMax(lines);
    return calculateParametersStrict(minMax, size);
}

// Overload for Triangles - Calls new strict parameter function
std::map<std::string, double> calculate(const Triangles &triangles, const int size) {
    auto minMax = calculateMinMax(triangles);
    return calculateParametersStrict(minMax, size);
}


img::EasyImage draw2DLines(Lines2D &lines, int size, const Color &backgroundColor) {
    //std::cout << "--- draw2DLines ---" << std::endl; // Marker
    auto resultaten = calculate(lines, size); // Uses calculateParametersStrict internally

    int finalImageWidth = static_cast<int>(resultaten.at("imageX"));
    int finalImageHeight = static_cast<int>(resultaten.at("imageY"));
    //std::cout << "Final canvas size for EasyImage: width=" << finalImageWidth << ", height=" << finalImageHeight << std::endl;


    // Double check dimensions are valid (should be >= 1 due to logic in calculateParametersStrict)
    if (finalImageWidth <= 0 || finalImageHeight <= 0) {
         std::cerr << "Warning: Calculated image dimensions are invalid (" << finalImageWidth << "x" << finalImageHeight << "). Using default size." << std::endl;
         finalImageWidth = std::max(1, size);
         finalImageHeight = std::max(1, size);
    }

    img::EasyImage image(finalImageWidth, finalImageHeight, img::Color(backgroundColor.red * 255, backgroundColor.green * 255, backgroundColor.blue * 255));

    if (lines.empty()) {
        //std::cout << "--- End draw2DLines (empty lines) ---" << std::endl;
        return image;
    }

    double d = resultaten.at("d");
    double dx = resultaten.at("dx");
    double dy = resultaten.at("dy");

    // Apply scaling and translation (using doubles)
    for (auto &line : lines) {
        line.p1.x *= d; line.p1.y *= d;
        line.p2.x *= d; line.p2.y *= d;
        line.p1.x += dx; line.p1.y += dy;
        line.p2.x += dx; line.p2.y += dy;
    }

    // Draw the transformed lines onto the image using std::round()
    //std::cout << "--- Drawing Scaled Lines ---" << std::endl;
    for (const auto &line : lines) {

        //std::cout << "Drawing line from (" << std::round(line.p1.x) << "," << std::round(line.p1.y)
        //          << ") to (" << std::round(line.p2.x) << "," << std::round(line.p2.y) << ")" << std::endl;

        img::Color drawColor(
            static_cast<uint8_t>(std::round(line.color.red * 255)),
            static_cast<uint8_t>(std::round(line.color.green * 255)),
            static_cast<uint8_t>(std::round(line.color.blue * 255))
        );
        image.draw_line(
            std::round(line.p1.x), std::round(line.p1.y),
            std::round(line.p2.x), std::round(line.p2.y),
            drawColor
        );
    }
    //std::cout << "--- Finished Drawing Scaled Lines ---" << std::endl;
    //std::cout << "--- End draw2DLines ---" << std::endl;
    return image;
}

// Draw L-System (Keep previous version - likely okay)
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
            bool shouldDraw = l_system.draw(c);
            Point2D nextPosition = currentPosition;
            nextPosition.x += std::cos(currentAngle_rad); // Assume step size is 1
            nextPosition.y += std::sin(currentAngle_rad);
            if (shouldDraw) {
                lines.push_back(Line2D(currentPosition, nextPosition, 0, 0, lineColor));
            }
            currentPosition = nextPosition;
        }
    }
    return lines;
}


// Generate 3D figures (Keep previous version - was okay)
Figures3D generateFigures(const ini::Configuration &configuration) {
    Figures3D figures;
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    for (int i = 0; i < nrFigures; ++i) {
        std::string fig_key = "Figure" + std::to_string(i);

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
        Figure figure;
        ini::DoubleTuple colorData = configuration[fig_key]["color"].as_double_tuple_or_die();
        Color color(colorData[0], colorData[1], colorData[2]);
        bool figureGenerated = true;

        // Handle different figure types (simplified)
        if (type == "LineDrawing") { /* ... keep implementation ... */
            int nrPoints = configuration[fig_key]["nrPoints"].as_int_or_die();
            int nrLines = configuration[fig_key]["nrLines"].as_int_or_die();
            for (int j = 0; j < nrPoints; ++j) {
                ini::DoubleTuple pointData = configuration[fig_key]["point" + std::to_string(j)].as_double_tuple_or_die();
                figure.points.push_back(Vector3D::point(pointData[0], pointData[1], pointData[2]));
            }
            for (int j = 0; j < nrLines; ++j) {
                ini::IntTuple lineData = configuration[fig_key]["line" + std::to_string(j)].as_int_tuple_or_die();
                if (lineData.size() >= 2 && lineData[0] >= 0 && lineData[0] < nrPoints && lineData[1] >= 0 && lineData[1] < nrPoints) {
                     figure.faces.push_back(Face{{lineData[0], lineData[1]}});
                } else { std::cerr << "Warning: Invalid line indices for line" << j << " in " << fig_key << std::endl; }
            }
        } else if (type == "3DLSystem") { /* ... keep implementation ... */
            std::string inputFile = configuration[fig_key]["inputfile"].as_string_or_die();
            try { figure = Figure::generate3DLSystem(inputFile, color); }
            catch (const std::runtime_error& e) { std::cerr << "Error generating 3DLSystem for " << fig_key << ": " << e.what() << std::endl; figureGenerated = false; }
        } else if (type == "Cube") { figure = Figure::createCube(); }
          else if (type == "Tetrahedron") { figure = Figure::createTetrahedron(); }
          else if (type == "Octahedron") { figure = Figure::createOctahedron(); }
          else if (type == "Icosahedron") { figure = Figure::createIcosahedron(); }
          else if (type == "Dodecahedron") { figure = Figure::createDodecahedron(); }
          else if (type == "Cylinder") { /* ... keep implementation ... */
              int n = configuration[fig_key]["n"].as_int_or_die(); double h = configuration[fig_key]["height"].as_double_or_die(); figure = Figure::createCylinder(n, h);
          } else if (type == "Cone") { /* ... keep implementation ... */
              int n = configuration[fig_key]["n"].as_int_or_die(); double h = configuration[fig_key]["height"].as_double_or_die(); figure = Figure::createCone(n, h);
          } else if (type == "Sphere") { /* ... keep implementation ... */
              int n = configuration[fig_key]["n"].as_int_or_die(); figure = Figure::createSphere(n);
          } else if (type == "Torus") { /* ... keep implementation ... */
              double R = configuration[fig_key]["R"].as_double_or_die(); double r = configuration[fig_key]["r"].as_double_or_die(); int n = configuration[fig_key]["n"].as_int_or_die(); int m = configuration[fig_key]["m"].as_int_or_die(); figure = Figure::createTorus(R, r, n, m);
          } else if (type == "MengerSponge") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); figure = generateMengerSponge(nrIterations);
          } else if (type == "FractalCube") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalCube(nrIterations, actualScaleFactor);
          } else if (type == "FractalTetrahedron") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalTetrahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalOctahedron") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalOctahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalIcosahedron") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalIcosahedron(nrIterations, actualScaleFactor);
          } else if (type == "FractalDodecahedron") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalDodecahedron(nrIterations, actualScaleFactor);
          } else if (type == "BuckyBall") { /* ... keep implementation ... */
              figure = Figure::createBuckyBall();
          } else if (type == "FractalBuckyBall") { /* ... keep implementation ... */
              int nrIterations = configuration[fig_key]["nrIterations"].as_int_or_die(); double fractalScaleInput = configuration[fig_key]["fractalScale"].as_double_or_die(); double actualScaleFactor = (fractalScaleInput > 1e-9) ? (1.0 / fractalScaleInput) : 1.0; figure = generateFractalBuckyBall(nrIterations, actualScaleFactor);
          } else {
               std::cerr << "Warning: Unknown figure type '" << type << "' specified for " << fig_key << std::endl;
               figureGenerated = false;
          }

        if (!figureGenerated) continue; // Skip applying color/transforms if generation failed

        figure.color = color; // Apply the color

        // Apply geometric transformations (scale, rotate, translate)
        double scale = configuration[fig_key]["scale"].as_double_or_die();
        double rotateXAngle = configuration[fig_key]["rotateX"].as_double_or_die();
        double rotateYAngle = configuration[fig_key]["rotateY"].as_double_or_die();
        double rotateZAngle = configuration[fig_key]["rotateZ"].as_double_or_die();
        ini::DoubleTuple centerData = configuration[fig_key]["center"].as_double_tuple_or_die();
        Vector3D centerVec = Vector3D::point(centerData[0], centerData[1], centerData[2]);

        Matrix scaleMatrix = scaleFigure(scale);
        Matrix rotateXMatrix = rotateX(rotateXAngle);
        Matrix rotateYMatrix = rotateY(rotateYAngle);
        Matrix rotateZMatrix = rotateZ(rotateZAngle);
        Matrix translateMatrix = translate(centerVec);
        Matrix transformationMatrix = scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * translateMatrix;

        applyTransformation(figure, transformationMatrix);
        figures.push_back(figure);
    }
    return figures;
}