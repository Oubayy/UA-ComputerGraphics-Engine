#ifndef LINEDRAWER_H
#define LINEDRAWER_H

#include "easy_image.h"
#include "Line2D.h" // Includes Point2D, Color, Lines2D, Triangle, Triangles
#include "Figure.h" // Includes Figures3D
#include "l_parser.h"
#include "ini_configuration.h"
#include <map>
#include <vector>
#include <string>

/**
 * Bepaalt de minimale en maximale waarden voor X en Y in een verzameling lijnen.
 * @param lines De verzameling lijnen.
 * @return Een map met de minimale en maximale waarden voor X en Y.
 */
std::map<std::string, double> calculateMinMax(const Lines2D &lines);
std::map<std::string, double> calculateMinMax(const Triangles &triangles);

/**
 * Berekent de grootte van de afbeelding, schaalt de lijntekening en verschuift deze.
 * @param lines De verzameling lijnen.
 * @param size De gewenste grootte van de afbeelding.
 * @return Een map met de berekende waarden (imageX, imageY, d, dx, dy, etc.).
 */
std::map<std::string, double> calculate(const Lines2D &lines, const int size);
std::map<std::string, double> calculate(const Triangles &triangles, const int size);

/**
 * Tekent een verzameling lijnen op een afbeelding. Past scaling/translation toe.
 * @param lines De verzameling lijnen (modified in place).
 * @param size De gewenste grootte van de afbeelding (input for calculation).
 * @param backgroundColor De achtergrondkleur van de afbeelding.
 * @return Een afbeelding met de getekende lijnen.
 */
img::EasyImage draw2DLines(Lines2D &lines, int size, const Color &backgroundColor);

/**
 * Genereert een verzameling lijnen volgens een L-System.
 * @param l_system Het L-System.
 * @param lineColor De kleur van de lijnen.
 * @return Een verzameling lijnen.
 */
Lines2D drawLSystem(const LParser::LSystem2D &l_system, const Color &lineColor);

/**
 * Genereert een verzameling figuren op basis van de configuratie, including fractals.
 * @param configuration De configuratie.
 * @return Een vector van figuren.
 */
Figures3D generateFigures(const ini::Configuration &configuration);


#endif // LINEDRAWER_H