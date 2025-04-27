#include "draw_triangulated_wireframes.h"
#include "easy_image.h"
#include "ZBuffer.h"
#include "vector3d.h"
#include "Line2D.h"
#include "Projection.h"
#include "Figure.h"
#include "easy_image.h"
#include "ini_configuration.h"
#include "lineDrawer.h"
#include "Transformations.h"
#include "Projection.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>

img::EasyImage draw_triangulated_wireframes(const ini::Configuration &configuration) {
    int size = configuration["General"]["size"].as_int_or_die();
    if (size <= 0) {
        throw std::runtime_error("Image size must be greater than 0.");
    }

    ini::DoubleTuple bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    ini::DoubleTuple eyeCoords = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye = Vector3D::point(eyeCoords[0], eyeCoords[1], eyeCoords[2]);

    Figures3D figures = generateFigures(configuration);

    // Triangulate all faces of all figures
    for (Figure &figure : figures) {
        std::vector<Face> new_faces;
        for (const Face &face : figure.faces) {
            for (size_t i = 1; i < face.point_indexes.size() - 1; ++i) {
                Face triangle;
                triangle.point_indexes.push_back(face.point_indexes[0]);
                triangle.point_indexes.push_back(face.point_indexes[i]);
                triangle.point_indexes.push_back(face.point_indexes[i + 1]);

                /*std::cout  << "Figure color: "
                        << figure.color.red << ", "
                        << figure.color.green << ", "
                        << figure.color.blue << std::endl;*/


                // triangle.color = face.color; // (Face has no color)
                new_faces.push_back(triangle);
            }
        }
        figure.faces = std::move(new_faces);

        // Make sure color isn't getting lost
        if (figure.color.red == 0 && figure.color.green == 0 && figure.color.blue == 0) {
            std::cerr << "Warning: figure has no color!" << std::endl;
            figure.color = {1.0, 0.0, 0.0};
        }
    }



    // Apply eye transformation
    Matrix eyeTransform = eyePointTrans(eye);
    applyTransformation(figures, eyeTransform);

    Triangles triangles = ZBuffer::doProjectTriangle(figures, 1);

    std::map<std::string, double> resultaten = calculate(triangles, size);
    int imageX = lround(resultaten["imageX"]);
    int imageY = lround(resultaten["imageY"]);

    img::EasyImage image(imageX, imageY, img::Color(
        (int)(bgColor[0] * 255),
        (int)(bgColor[1] * 255),
        (int)(bgColor[2] * 255)
    ));

    ZBuffer zbuffer(imageX, imageY);

    double d = resultaten["d"];
    double dx = resultaten["dx"];
    double dy = resultaten["dy"];

    for (auto &triangle : triangles) {
        triangle.p1.x = triangle.p1.x * d + dx;
        triangle.p1.y = triangle.p1.y * d + dy;

        triangle.p2.x = triangle.p2.x * d + dx;
        triangle.p2.y = triangle.p2.y * d + dy;

        triangle.p3.x = triangle.p3.x * d + dx;
        triangle.p3.y = triangle.p3.y * d + dy;

        /*
        std::cout << "Triangle: ("
          << triangle.p1.x << "," << triangle.p1.y << ") -> ("
          << triangle.p2.x << "," << triangle.p2.y << ") -> ("
          << triangle.p3.x << "," << triangle.p3.y << ") | Color: ("
          << triangle.color.red << ", "
          << triangle.color.green << ", "
          << triangle.color.blue << ")\n";*/

        zbuffer.draw_zbuf_triangle(image, triangle, 1);
    }

    //std::cout << "Finished drawing, image size: " << image.get_width() << "x" << image.get_height() << std::endl;

    return image;
}
