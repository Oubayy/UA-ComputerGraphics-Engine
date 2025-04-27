#ifndef ENGINE_FRACTALS3D_H
#define ENGINE_FRACTALS3D_H

#include "Figure.h" // Includes vector3d.h and Line2D.h (for Figures3D)

// Function declarations for generating fractals
// Note: actualScaleFactor is the factor applied (e.g., 0.5), not the input 'fractalScale' (e.g., 2)
Figure generateMengerSponge(int nrIterations);
Figure generateFractalTetrahedron(int nrIterations, double actualScaleFactor);
Figure generateFractalIcosahedron(int nrIterations, double actualScaleFactor);
Figure generateFractalCube(int nrIterations, double actualScaleFactor);
Figure generateFractalOctahedron(int nrIterations, double actualScaleFactor);
Figure generateFractalDodecahedron(int nrIterations, double actualScaleFactor);
Figure generateFractalBuckyBall(int nrIterations, double actualScaleFactor);

#endif //ENGINE_FRACTALS3D_H