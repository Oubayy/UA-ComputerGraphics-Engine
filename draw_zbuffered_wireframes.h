#ifndef ENGINE_DRAW_ZBUFFERED_WIREFRAMES_H
#define ENGINE_DRAW_ZBUFFERED_WIREFRAMES_H

/**
 * @file draw_zbuffered_wireframes.h
 * @brief Contains the declaration for drawing 3D wireframes using Z-buffering.
 */

#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"

/**
 * @brief Generates an image of 3D wireframes using Z-buffering for hidden line removal.
 *
 * Reads configuration, generates 3D figures, applies transformations (including eye point),
 * projects figures to 2D lines (storing depth), calculates scaling/translation,
 * initializes a Z-buffer, performs line clipping (Cohen-Sutherland),
 * and draws visible line segments using the Z-buffer algorithm.
 *
 * @param configuration The INI configuration containing settings for the scene.
 * @return An img::EasyImage object representing the rendered wireframe image.
 */
img::EasyImage draw_zbuffered_wireframes(const ini::Configuration &configuration);

#endif //ENGINE_DRAW_ZBUFFERED_WIREFRAMES_H