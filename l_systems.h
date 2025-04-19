#ifndef L_SYSTEMS_H
#define L_SYSTEMS_H

#include "ini_configuration.h"
#include "easy_image.h"

namespace l_systems {
    img::EasyImage generateLSystem(const ini::Configuration &configuration);
}

#endif // L_SYSTEMS_H
