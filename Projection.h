#ifndef ENGINE_PROJECTION_H
#define ENGINE_PROJECTION_H


#include "vector3d.h"
#include "Line2D.h"
#include "Figure.h"

Point2D doProjection(const Vector3D &point, const double d);
//Lines2D doProjection(const Figures3D &figures, const double d);
Lines2D doProjection(const Figures3D &figures);

#endif //ENGINE_PROJECTION_H
