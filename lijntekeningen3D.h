#ifndef ENGINE_LIJNTEKENINGEN3D_H
#define ENGINE_LIJNTEKENINGEN3D_H

#include "ini_configuration.h"
#include "easy_image.h"

namespace lijntekeningen3D {
    img::EasyImage generate3DImage(const ini::Configuration &configuration);
}

#endif // ENGINE_LIJNTEKENINGEN3D_H
