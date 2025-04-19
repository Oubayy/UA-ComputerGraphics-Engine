//
// Created by Oubay on 27/07/2024.
//

#ifndef ENGINE_TRANSFORMATIONS_H
#define ENGINE_TRANSFORMATIONS_H


#include "vector3d.h"
#include "Figure.h"

// Function to create a scaling matrix
Matrix scaleFigure(const double scale);

// Function to create a rotation matrix around the X-axis
Matrix rotateX(const double angle);

// Function to create a rotation matrix around the Y-axis
Matrix rotateY(const double angle);

// Function to create a rotation matrix around the Z-axis
Matrix rotateZ(const double angle);

// Function to create a translation matrix
Matrix translate(const Vector3D &vector);

// Function to apply transformations to a figure
void applyTransformation(Figure &fig, const Matrix &m);

void applyTransformation(Figures3D &figs, const Matrix &m);


void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
Matrix eyePointTrans(const Vector3D &eyepoint);

#endif //ENGINE_TRANSFORMATIONS_H
