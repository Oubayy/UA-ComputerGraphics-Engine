#ifndef ENGINE_FRACTALS3D_H
#define ENGINE_FRACTALS3D_H

#include "Figure.h"
#include "vector3d.h"
#include "easy_image.h"
#include "ini_configuration.h"
#include <list>

Figure generateMengerSponge(int nrIterations);
Figure generateFractalTetrahedron(int nrIterations);
Figure generateFractalIcosahedron(int nrIterations);
Figure generateFractalCube(int nrIterations);
Figure generateFractalOctahedron(int nrIterations);

#endif //ENGINE_FRACTALS3D_H
