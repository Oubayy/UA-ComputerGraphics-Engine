#include "lineDrawer.h"
#include "Figure.h"
#include "Transformations.h"
#include "Projection.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include <cmath>
#include <limits>
#include <stack>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef Pi
#define Pi 3.14159265358979323846
#endif

using namespace std;

map<string, double> calculateMinMax(const Triangles &triangles) {
    map<string, double> resultaat;

    // initialiseren
    double xmin = numeric_limits<double>::max();
    double xmax = numeric_limits<double>::lowest();
    double ymin = numeric_limits<double>::max();
    double ymax = numeric_limits<double>::lowest();

    // Doorloop alle lijnen om de minimale en maximale waarden te vinden
    for (const Triangle &triangle : triangles) {
        if (triangle.p1.x < xmin) xmin = triangle.p1.x;
        if (triangle.p1.x > xmax) xmax = triangle.p1.x;
        if (triangle.p2.x < xmin) xmin = triangle.p2.x;
        if (triangle.p2.x > xmax) xmax = triangle.p2.x;
        if (triangle.p1.y < ymin) ymin = triangle.p1.y;
        if (triangle.p1.y > ymax) ymax = triangle.p1.y;
        if (triangle.p2.y < ymin) ymin = triangle.p2.y;
        if (triangle.p2.y > ymax) ymax = triangle.p2.y; //
        if (triangle.p3.x < xmin) xmin = triangle.p3.x;
        if (triangle.p3.x > xmax) xmax = triangle.p3.x;
        if (triangle.p3.y < ymin) ymin = triangle.p3.y;
        if (triangle.p3.y > ymax) ymax = triangle.p3.y;
    }

    // resultaten opslaan in de map
    resultaat["xmin"] = xmin;
    resultaat["xmax"] = xmax;
    resultaat["ymin"] = ymin;
    resultaat["ymax"] = ymax;

    return resultaat;
}

// Xmin, Xmax, Ymin en Ymax berekenen
map<string, double> calculateMinMax(const Lines2D &lines) {
    map<string, double> resultaat;

    // initialiseren
    double xmin = numeric_limits<double>::max();
    double xmax = numeric_limits<double>::lowest();
    double ymin = numeric_limits<double>::max();
    double ymax = numeric_limits<double>::lowest();

    // Doorloop alle lijnen om de minimale en maximale waarden te vinden
    for (const Line2D &line : lines) {
        if (line.p1.x < xmin) xmin = line.p1.x;
        if (line.p1.x > xmax) xmax = line.p1.x;
        if (line.p2.x < xmin) xmin = line.p2.x;
        if (line.p2.x > xmax) xmax = line.p2.x;
        if (line.p1.y < ymin) ymin = line.p1.y;
        if (line.p1.y > ymax) ymax = line.p1.y;
        if (line.p2.y < ymin) ymin = line.p2.y;
        if (line.p2.y > ymax) ymax = line.p2.y;
    }

    // resultaten opslaan in de map
    resultaat["xmin"] = xmin;
    resultaat["xmax"] = xmax;
    resultaat["ymin"] = ymin;
    resultaat["ymax"] = ymax;

    return resultaat;
}

map<string, double> calculate(const Triangles &triangles, const int size) {
    map<string, double> resultaten;

    auto minMax = calculateMinMax(triangles);

    double xmin, xmax, ymin, ymax;

    double Xrange = minMax["xmax"] - minMax["xmin"];
    double Yrange = minMax["ymax"] - minMax["ymin"];

    //std::cout << "Xrange: " << Xrange << ", Yrange: " << Yrange << std::endl;

    // Set a minimum Xrange threshold to avoid near-zero values
    if (Xrange < 1e-6) {
        Xrange = 1.0; // Set a default value if Xrange is too small
    }
    if (Yrange < 1e-6) {
        Yrange = 1.0; // Set a default value if Yrange is too small
    }

    //std::cout << "Adjusted Xrange: " << Xrange << ", Adjusted Yrange: " << Yrange << std::endl;

    resultaten["Xmax"] = minMax["xmax"];
    resultaten["Xmin"] = minMax["xmin"];
    resultaten["Ymax"] = minMax["ymax"];
    resultaten["Ymin"] = minMax["ymin"];
    resultaten["Xrange"] = Xrange;
    resultaten["Yrange"] = Yrange;

    // Bereken de grootte van de image
    double imageX = size * (Xrange / max(Xrange, Yrange));
    double imageY = size * (Yrange / max(Xrange, Yrange));

    //std::cout << "Calculated image size: " << imageX << " x " << imageY << std::endl;

    // Ensure the calculated image size is not too small
    if (imageX < 1.0) imageX = 1.0;
    if (imageY < 1.0) imageY = 1.0;

    //std::cout << "Final image size: " << imageX << " x " << imageY << std::endl;

    resultaten["imageX"] = imageX;
    resultaten["imageY"] = imageY;

    // Bereken de schaalfactor
    double d = 0.95 * (imageX / Xrange);
    resultaten["d"] = d;

    // Bereken de verschuiving
    double DCx = d * ((minMax["xmin"] + minMax["xmax"]) / 2);
    double DCy = d * ((minMax["ymin"] + minMax["ymax"]) / 2);
    resultaten["DCx"] = DCx;
    resultaten["DCy"] = DCy;

    double dx = (imageX / 2) - DCx;
    double dy = (imageY / 2) - DCy;
    resultaten["dx"] = dx;
    resultaten["dy"] = dy;

    return resultaten;
}

// Bereken de grootte van de afbeelding, schaal de lijntekening en verschuif deze
map<string, double> calculate(const Lines2D &lines, const int size) {
    map<string, double> resultaten;

    auto minMax = calculateMinMax(lines);

    double xmin, xmax, ymin, ymax;

    double Xrange = minMax["xmax"] - minMax["xmin"];
    double Yrange = minMax["ymax"] - minMax["ymin"];

    //std::cout << "Xrange: " << Xrange << ", Yrange: " << Yrange << std::endl;

    // Set a minimum Xrange threshold to avoid near-zero values
    if (Xrange < 1e-6) {
        Xrange = 1.0; // Set a default value if Xrange is too small
    }
    if (Yrange < 1e-6) {
        Yrange = 1.0; // Set a default value if Yrange is too small
    }

    //std::cout << "Adjusted Xrange: " << Xrange << ", Adjusted Yrange: " << Yrange << std::endl;

    resultaten["Xmax"] = minMax["xmax"];
    resultaten["Xmin"] = minMax["xmin"];
    resultaten["Ymax"] = minMax["ymax"];
    resultaten["Ymin"] = minMax["ymin"];
    resultaten["Xrange"] = Xrange;
    resultaten["Yrange"] = Yrange;

    // Bereken de grootte van de image
    double imageX = size * (Xrange / max(Xrange, Yrange));
    double imageY = size * (Yrange / max(Xrange, Yrange));

    //std::cout << "Calculated image size: " << imageX << " x " << imageY << std::endl;

    // Ensure the calculated image size is not too small
    if (imageX < 1.0) imageX = 1.0;
    if (imageY < 1.0) imageY = 1.0;

    //std::cout << "Final image size: " << imageX << " x " << imageY << std::endl;

    resultaten["imageX"] = imageX;
    resultaten["imageY"] = imageY;

    // Bereken de schaalfactor
    double d = 0.95 * (imageX / Xrange);
    resultaten["d"] = d;

    // Bereken de verschuiving
    double DCx = d * ((minMax["xmin"] + minMax["xmax"]) / 2);
    double DCy = d * ((minMax["ymin"] + minMax["ymax"]) / 2);
    resultaten["DCx"] = DCx;
    resultaten["DCy"] = DCy;

    double dx = (imageX / 2) - DCx;
    double dy = (imageY / 2) - DCy;
    resultaten["dx"] = dx;
    resultaten["dy"] = dy;

    return resultaten;
}

img::EasyImage draw2DLines(Lines2D &lines, int size, const Color &backgroundColor) {

    //std::cout << "Starting to draw 2D lines" << std::endl;

    // Bereken de benodigde parameters voor de tekening
    auto resultaten = calculate(lines, size);

    double imageX = resultaten["imageX"];
    double imageY = resultaten["imageY"];

    // Maak een afbeelding met de berekende afmetingen
    img::EasyImage image(imageX, imageY, {static_cast<uint8_t>(backgroundColor.red * 255), static_cast<uint8_t>(backgroundColor.green * 255), static_cast<uint8_t>(backgroundColor.blue * 255)});

    /*
    // Fout aangezien dit de image pixel per pixel inkleurt, onnodig
    for (unsigned int i = 0; i < image.get_width(); i++) {
        for (unsigned int j = 0; j < image.get_height(); j++) {
            image(i,j).red = backgroundColor.red * 255;
            image(i,j).green = backgroundColor.green * 255;
            image(i,j).blue = backgroundColor.blue * 255;
        }
    }
    */

    double d = resultaten["d"];

    // Maak een kopie van de lijnen om te schalen en verschuiven
    //Lines2D scaledLines = lines;

    for (auto &line : lines) {
        line.p1.x *= d;
        line.p1.y *= d;
        line.p2.x *= d;
        line.p2.y *= d;
    }

    double dx = resultaten["dx"];
    double dy = resultaten["dy"];


    //std::cout << "Image size: " << imageX << " x " << imageY << std::endl;
    //std::cout << "Scaling factor: " << d << ", Translation: (" << dx << ", " << dy << ")" << std::endl;

    for (auto &line : lines) {
        line.p1.x += dx;
        line.p1.y += dy;
        line.p2.x += dx;
        line.p2.y += dy;

        //std::cout << "Drawing line from (" << line.p1.x << ", " << line.p1.y << ") to (" << line.p2.x << ", " << line.p2.y << ")" << std::endl;

        // Ensure the line coordinates are within the image bounds
        if (line.p1.x < 0 || line.p1.x >= imageX || line.p1.y < 0 || line.p1.y >= imageY ||
            line.p2.x < 0 || line.p2.x >= imageX || line.p2.y < 0 || line.p2.y >= imageY) {
            std::cerr << "Line coordinates out of bounds: (" << line.p1.x << ", " << line.p1.y << ") to (" << line.p2.x << ", " << line.p2.y << ") in image of size (" << imageX << ", " << imageY << ")" << std::endl;
        }
    }

    // De lijnen tekenen op de afbeelding
    for (const auto &line : lines) {
        img::Color color(line.color.red * 255, line.color.green * 255, line.color.blue * 255);
        image.draw_line(std::lround(line.p1.x), std::lround(line.p1.y), std::lround(line.p2.x), std::lround(line.p2.y), color);

        //std::cout << "Drew line from (" << line.p1.x << ", " << line.p1.y << ") to (" << line.p2.x << ", " << line.p2.y << ")" << std::endl;

    }

    //std::cout << "Finished drawing 2D lines" << std::endl;
    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system, const Color &lineColor) {
    Lines2D lines;
    std::stack<std::pair<Point2D, double>> positionStack;

    Point2D currentPosition(0.0, 0.0);
    double currentAngle = l_system.get_starting_angle() * Pi / 180.0;

    std::string instructions = l_system.get_initiator();
    for (int i = 0; i < l_system.get_nr_iterations(); ++i) {
        std::string newInstructions;
        for (char c : instructions) {
            if (l_system.get_alphabet().count(c)) {
                newInstructions += l_system.get_replacement(c);
            } else {
                newInstructions += c;
            }
        }
        instructions = newInstructions;
    }

    for (char c : instructions) {
        if (l_system.get_alphabet().count(c)) {
            if (l_system.draw(c)) {
                Point2D newPosition = currentPosition;
                newPosition.x += std::cos(currentAngle);
                newPosition.y += std::sin(currentAngle);
                lines.push_back(Line2D(currentPosition, newPosition, 0, 0, lineColor));
                currentPosition = newPosition;
            } else {
                // Move forward without drawing
                currentPosition.x += std::cos(currentAngle);
                currentPosition.y += std::sin(currentAngle);
            }
        } else if (c == '+') {
            currentAngle += l_system.get_angle() * Pi / 180.0;
        } else if (c == '-') {
            currentAngle -= l_system.get_angle() * Pi / 180.0;
        } else if (c == '(') {
            positionStack.push({currentPosition, currentAngle});
        } else if (c == ')') {
            if (!positionStack.empty()) {
                currentPosition = positionStack.top().first;
                currentAngle = positionStack.top().second;
                positionStack.pop();
            }
        }
        // Else: completely unknown symbol, skip
    }

    return lines;
}


    /*
    // Verwerk de instructies om de lijnen te tekenen
    for (char c : instructions) {
        switch (c) {
            case '+':
                currentAngle += l_system.get_angle() * Pi / 180.0;
                break;
            case '-':
                currentAngle -= l_system.get_angle() * Pi / 180.0;
                break;
            case '(':
                positionStack.push({currentPosition, currentAngle});
                break;
            case ')':
                std::tie(currentPosition, currentAngle) = positionStack.top();
                positionStack.pop();
                break;
            default:
                if (l_system.draw(c)) {
                    Point2D newPosition(
                            currentPosition.x + cos(currentAngle),
                            currentPosition.y + sin(currentAngle)
                    );
                    lines.push_back(Line2D(currentPosition, newPosition, lineColor));
                    currentPosition = newPosition;
                } else {
                    currentPosition.x += cos(currentAngle);
                    currentPosition.y += sin(currentAngle);
                }
                break;
        }
    }
    */


// Generate 3D figures based on the configuration
Figures3D generateFigures(const ini::Configuration &configuration) {
    Figures3D figures;
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    for (int i = 0; i < nrFigures; ++i) {
        std::string type = configuration["Figure" + std::to_string(i)]["type"].as_string_or_die();
        //std::cout << "Processing figure " << i << " of type: " << type << std::endl;
        Figure figure;

        // Read the color
        ini::DoubleTuple colorData = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
        Color color(colorData[0], colorData[1], colorData[2]);

        if (type == "Cube") {
            figure = Figure::createCube();
        } else if (type == "Tetrahedron") {
            figure = Figure::createTetrahedron();
        } else if (type == "Octahedron") {
            figure = Figure::createOctahedron();
        } else if (type == "Icosahedron") {
            figure = Figure::createIcosahedron();
        } else if (type == "Dodecahedron") {
            figure = Figure::createDodecahedron();
        } else if (type == "Cylinder") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            figure = Figure::createCylinder(n, height);
        } else if (type == "Cone") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            figure = Figure::createCone(n, height);
        } else if (type == "Sphere") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            figure = Figure::createSphere(n);
        } else if (type == "Torus") {
            double R = configuration["Figure" + std::to_string(i)]["R"].as_double_or_die();
            double r = configuration["Figure" + std::to_string(i)]["r"].as_double_or_die();
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            int m = configuration["Figure" + std::to_string(i)]["m"].as_int_or_die();
            figure = Figure::createTorus(R, r, n, m);
        } else if (type == "3DLSystem") {
            std::string inputFile = configuration["Figure" + std::to_string(i)]["inputfile"].as_string_or_die();
            figure = Figure::generate3DLSystem(inputFile, color);
        } else if (type == "LineDrawing") {

            //std::cout << "Zzzzzzzzzzzzzzzzzzzz " << std::endl;

            int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"].as_int_or_die();

            // Read points
            for (int j = 0; j < nrPoints; ++j) {
                std::vector<double> pointData = configuration["Figure" + std::to_string(i)]["point" + std::to_string(
                        j)].as_double_tuple_or_die();
                Vector3D point = Vector3D::point(pointData[0], pointData[1], pointData[2]);
                figure.points.push_back(point);
            }

            // Read lines
            for (int j = 0; j < nrLines; ++j) {
                std::vector<int> lineData = configuration["Figure" + std::to_string(i)]["line" + std::to_string(
                        j)].as_int_tuple_or_die();
                Face face;
                face.point_indexes = lineData;
                figure.faces.push_back(face);
            }
        }

        // Apply the color to the figure
        figure.color = color;

        // Apply transformations
        double scale = configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die();
        double rotateXAngle = configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die();
        double rotateYAngle = configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die();
        double rotateZAngle = configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die();
        std::vector<double> center = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();

        Matrix scaleMatrix = scaleFigure(scale);
        Matrix rotateXMatrix = rotateX(rotateXAngle);
        Matrix rotateYMatrix = rotateY(rotateYAngle);
        Matrix rotateZMatrix = rotateZ(rotateZAngle);
        Matrix translateMatrix = translate(Vector3D::point(center[0], center[1], center[2]));

        Matrix transformationMatrix = scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * translateMatrix;
        applyTransformation(figure, transformationMatrix);

        figures.push_back(figure);
        //std::cout << "Figure " << i << " (" << type << ") created successfully." << std::endl;
    }
    return figures;
}

//figure.color = Color(colorData[0], colorData[1], colorData[2]);
