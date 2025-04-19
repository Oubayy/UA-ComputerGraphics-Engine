#include "l_systems.h"
#include <fstream>
#include "l_parser.h"
#include "lineDrawer.h"
#include "Line2D.h"
#include "Figure.h"

img::EasyImage l_systems::generateLSystem(const ini::Configuration &configuration) {
    int size = configuration["General"]["size"].as_int_or_die();
    ini::DoubleTuple bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    ini::DoubleTuple lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
    std::string inputFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();

    LParser::LSystem2D lSystem;
    std::ifstream inputStream(inputFile);
    inputStream >> lSystem;
    inputStream.close();

    Color lineColorObj(lineColor[0], lineColor[1], lineColor[2]);
    Lines2D lines = drawLSystem(lSystem, lineColorObj);

    Color bgColorObj(bgColor[0], bgColor[1], bgColor[2]);
    return draw2DLines(lines, size, bgColorObj);
}
