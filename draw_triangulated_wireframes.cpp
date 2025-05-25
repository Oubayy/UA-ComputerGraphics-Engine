#include "draw_triangulated_wireframes.h"
#include "easy_image.h"
#include "ZBuffer.h"
#include "Projection.h"
#include "ini_configuration.h"
#include "lineDrawer.h"
#include "Transformations.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

img::EasyImage draw_triangulated_wireframes(const ini::Configuration &configuration) {
    int size = configuration["General"]["size"].as_int_or_die();
    if (size <= 0) throw std::runtime_error("Image size must be > 0.");

    ini::DoubleTuple bgColorTuple = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color backgroundColor(bgColorTuple[0], bgColorTuple[1], bgColorTuple[2]);

    ini::DoubleTuple eyeCoords = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyePoint = Vector3D::point(eyeCoords[0], eyeCoords[1], eyeCoords[2]);
    Vector3D eyePositionInEyeSpace = Vector3D::point(0,0,0);

    Figures3D figures = generateFigures(configuration);

    for (Figure &figure : figures) {
        if (figure.faces.empty()) continue;
        std::vector<Face> new_triangulated_faces;
        for (const Face &face : figure.faces) {
            if (face.point_indexes.size() >= 3) {
                for (size_t i = 1; i < face.point_indexes.size() - 1; ++i) {
                    new_triangulated_faces.push_back(Face{{face.point_indexes[0], face.point_indexes[i], face.point_indexes[i + 1]}});
                }
            }
        }
        figure.faces = std::move(new_triangulated_faces);
    }

    Matrix eyeTransform = eyePointTrans(eyePoint);
    applyTransformation(figures, eyeTransform); // Figures are now in eye space

    Lights3D lights_list;
    std::string type = configuration["General"]["type"].as_string_or_die();

    if (type == "LightedZBuffering") {
        lights_list = generateLights(configuration, eyeTransform); // Pass eyeTransform
    } else {
        DirectionalLight* defaultLight = new DirectionalLight();
        defaultLight->ambientLight = Color(1.0, 1.0, 1.0);
        // Default light's direction is (0,0,-1) in world.
        defaultLight->direction = Vector3D::vector(0,0,-1.0) * eyeTransform;
        defaultLight->direction.normalise();
        lights_list.push_back(defaultLight);
    }

    const double d_projection_parameter = 1.0;
    Triangles projected_triangles = ZBuffer::doProjectTriangle(figures, d_projection_parameter);

    if (projected_triangles.empty()) {
        for (Light* l : lights_list) delete l;
        return img::EasyImage(size, size, img::Color(static_cast<uint8_t>(backgroundColor.red*255), static_cast<uint8_t>(backgroundColor.green*255), static_cast<uint8_t>(backgroundColor.blue*255)));
    }

    std::map<std::string, double> screen_params = calculate(projected_triangles, size);
    int imageX = lround(screen_params["imageX"]);
    int imageY = lround(screen_params["imageY"]);
    double d_scaling = screen_params["d"];
    double dx_translation = screen_params["dx"];
    double dy_translation = screen_params["dy"];

    img::EasyImage image(imageX, imageY, img::Color(static_cast<uint8_t>(backgroundColor.red * 255), static_cast<uint8_t>(backgroundColor.green * 255), static_cast<uint8_t>(backgroundColor.blue * 255)));
    ZBuffer zbuffer(imageX, imageY);

    for (auto &triangle_to_draw : projected_triangles) {
        triangle_to_draw.p1.x = triangle_to_draw.p1.x * d_scaling + dx_translation;
        triangle_to_draw.p1.y = triangle_to_draw.p1.y * d_scaling + dy_translation;
        triangle_to_draw.p2.x = triangle_to_draw.p2.x * d_scaling + dx_translation;
        triangle_to_draw.p2.y = triangle_to_draw.p2.y * d_scaling + dy_translation;
        triangle_to_draw.p3.x = triangle_to_draw.p3.x * d_scaling + dx_translation;
        triangle_to_draw.p3.y = triangle_to_draw.p3.y * d_scaling + dy_translation;

        zbuffer.draw_zbuf_triangle(image, triangle_to_draw, lights_list, eyePositionInEyeSpace, d_projection_parameter, d_scaling, dx_translation, dy_translation);
    }

    for (Light* l : lights_list) delete l;
    return image;
}