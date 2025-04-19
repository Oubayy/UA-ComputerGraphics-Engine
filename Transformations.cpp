#include "Transformations.h"
#include <cmath>

#ifndef Pi
#define Pi 3.14159265358979323846
#endif

Matrix scaleFigure(const double scale) {
    Matrix m;
    m(1, 1) = scale;
    m(2, 2) = scale;
    m(3, 3) = scale;
    return m;
}

Matrix rotateX(const double angle) {
    Matrix m;
    double rad = angle * Pi / 180.0;
    m(2, 2) = cos(rad);
    m(2, 3) = sin(rad);
    m(3, 2) = -sin(rad);
    m(3, 3) = cos(rad);
    return m;
}

Matrix rotateY(const double angle) {
    Matrix m;
    double rad = angle * Pi / 180.0;
    m(1, 1) = cos(rad);
    m(1, 3) = -sin(rad);
    m(3, 1) = sin(rad);
    m(3, 3) = cos(rad);
    return m;
}

Matrix rotateZ(const double angle) {
    Matrix m;
    double rad = angle * Pi / 180.0;
    m(1, 1) = cos(rad);
    m(1, 2) = sin(rad);
    m(2, 1) = -sin(rad);
    m(2, 2) = cos(rad);
    return m;
}

Matrix translate(const Vector3D &vector) {
    Matrix m;
    m(4, 1) = vector.x;
    m(4, 2) = vector.y;
    m(4, 3) = vector.z;
    return m;
}

void applyTransformation(Figure &fig, const Matrix &m) {
    for (auto &point : fig.points) {
        point *= m;
    }
    //std::cout << "Applied transformation to figure" << std::endl;
}

void applyTransformation(Figures3D &figs, const Matrix &m) {
    for (auto &fig : figs) {
        applyTransformation(fig, m);
    }
    //std::cout << "Applied transformation to all figures" << std::endl;
}

// De carthesische coördinaten omzetten naar poolcoördinaten
void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = point.length();
    theta = atan2(point.y, point.x);
    phi = acos(point.z / r);
}

Matrix eyePointTrans(const Vector3D &eyepoint) {
    double theta, phi, r;

    // De carthesische coördinaten omzetten naar poolcoördinaten
    toPolar(eyepoint, theta, phi, r);

    // Debug statements
    //std::cout << "Eye coordinates: " << eyepoint.x << ", " << eyepoint.y << ", " << eyepoint.z << std::endl;
    //std::cout << "Theta: " << theta << ", Phi: " << phi << ", R: " << r << std::endl;

    // eyePointTransformationMatrix
    Matrix tM;
    tM(1, 1) = -sin(theta);
    tM(1, 2) = -cos(theta) * cos(phi);
    tM(1, 3) = cos(theta) * sin(phi);
    tM(2, 1) = cos(theta);
    tM(2, 2) = -sin(theta) * cos(phi);
    tM(2, 3) = sin(theta) * sin(phi);
    tM(3, 2) = sin(phi);
    tM(3, 3) = cos(phi);
    tM(4, 3) = -r;

    // Debug statement
    //std::cout << "Transformation Matrix: " << std::endl << tM << std::endl;

    return tM;
}