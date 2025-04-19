#include "lijntekeningen3D.h"
#include "lineDrawer.h"
#include "Transformations.h"
#include "Projection.h"
#include "vector3d.h"

img::EasyImage lijntekeningen3D::generate3DImage(const ini::Configuration &configuration) {
    int size = configuration["General"]["size"].as_int_or_die();
    ini::DoubleTuple bgTuple = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color background(bgTuple[0], bgTuple[1], bgTuple[2]);

    Figures3D figures = generateFigures(configuration);

    // Eye point
    ini::DoubleTuple eye = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyePoint = Vector3D::point(eye[0], eye[1], eye[2]);

    // Transform to eye space
    Matrix eyeTrans = eyePointTrans(eyePoint);
    applyTransformation(figures, eyeTrans);

    // Project to 2D lines
    Lines2D lines = doProjection(figures);

    // Draw lines
    return draw2DLines(lines, size, background);
}
