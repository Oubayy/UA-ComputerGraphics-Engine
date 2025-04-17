#ifndef INTRO_EXERCISES_H
#define INTRO_EXERCISES_H

#include "easy_image.h"
#include "ini_configuration.h"

namespace intro {
    // Opdracht 1: ColorRectangle
    img::EasyImage generateColorRectangle(const ini::Configuration &configuration);

    // Opdracht 2: Blocks
    img::EasyImage generateBlocks(const ini::Configuration &configuration);

    // 3) Lines (opdracht 3: QuarterCircle, Eye, Diamond)
    img::EasyImage generateLines(const ini::Configuration &configuration);

}

#endif // INTRO_EXERCISES_H


// y