#include "intro_exercises.h"
#include <cmath>    // std::lround, std::floor
#include <string>

// Helper: zet [0..1]-tuple om naar rgb [0..255]
static img::Color tupleToColor(const ini::DoubleTuple &t) {
    return {
        uint8_t(std::lround(t[0] * 255)),
        uint8_t(std::lround(t[1] * 255)),
        uint8_t(std::lround(t[2] * 255))
    };
}

namespace intro {

// ——————————————————————————————
// Opdracht 1: ColorRectangle
// ——————————————————————————————
img::EasyImage generateColorRectangle(const ini::Configuration &configuration) {
    int W = configuration["ImageProperties"]["width" ].as_int_or_die();
    int H = configuration["ImageProperties"]["height"].as_int_or_die();

    img::EasyImage image(W, H);
    
    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            uint8_t r = x % 256;
            uint8_t g = y % 256;
            uint8_t b = (r + g) % 256;
            image(x,y) = { r, g, b };
        }
    }
    return image;
}

// ——————————————————————————————
// Opdracht 2: Blocks
// ——————————————————————————————
img::EasyImage generateBlocks(const ini::Configuration &configuration) {
    int W = configuration["ImageProperties"]["width" ].as_int_or_die();
    int H = configuration["ImageProperties"]["height"].as_int_or_die();

    double nX = configuration["BlockProperties"]["nrXBlocks"].as_double_or_die();
    double nY = configuration["BlockProperties"]["nrYBlocks"].as_double_or_die();

    auto whiteT = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
    auto blackT = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
    bool invert = false;
    configuration["BlockProperties"]["invertColors"].as_bool_if_exists(invert);
    if (invert) std::swap(whiteT, blackT);

    img::Color white = tupleToColor(whiteT);
    img::Color black = tupleToColor(blackT);
    img::EasyImage image(W, H);

    double Bw = W / nX;
    double Bh = H / nY;

    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            int bx = std::floor(x / Bw);
            int by = std::floor(y / Bh);
            image(x,y) = ((bx + by) % 2 == 0 ? white : black);
        }
    }
    return image;
}

// —————————————————————————————————————————————
// Opdracht 3: Lines (QuarterCircle, Eye, Diamond)
// —————————————————————————————————————————————
img::EasyImage generateLines(const ini::Configuration &configuration) {
    int W      = configuration["ImageProperties"]["width" ].as_int_or_die();
    int H      = configuration["ImageProperties"]["height"].as_int_or_die();
    std::string fig = configuration["LineProperties"]["figure"].as_string_or_die();
    int M      = configuration["LineProperties"]["nrLines"].as_int_or_die();
    auto bgT   = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
    auto lnT   = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();

    img::Color bg = tupleToColor(bgT);
    img::Color ln = tupleToColor(lnT);
    img::EasyImage image(W, H, bg);

    // we gebruiken de werkelijke pixel-range [0..W-1], [0..H-1]
    double ws = double(W - 1) / double(M - 1);
    double hs = double(H - 1) / double(M - 1);

    if (fig == "QuarterCircle") {
        // basislijn onderaan
        image.draw_line(0, H - 1, W - 1, H - 1, ln);
        for (int i = 0; i < M; ++i) {
            int x = std::lround(i * ws);
            int y = std::lround(i * hs);
            image.draw_line(0, y, x, H - 1, ln);
        }
    }
    else if (fig == "Eye") {
        // linksonder → rechtsboven
        for (int i = 0; i < M; ++i) {
            int x = std::lround(i * ws);
            int y = std::lround(i * hs);
            image.draw_line(0, y, x, H - 1, ln);
        }
        // rechtsboven → linksonder
        for (int i = 0; i < M; ++i) {
            int x = std::lround(i * ws);
            int y = std::lround(i * hs);
            image.draw_line(W - 1, H - 1 - y, W - 1 - x, 0, ln);
        }
    }
    else if (fig == "Diamond") {
        auto drawFan = [&](int x1, int y1, int x2, int y2){
            for (int i = 0; i < M; ++i) {
                double t = double(i) / double(M - 1);
                int sx = std::lround(x1 + (x2 - x1) * t);
                int sy = std::lround(y1 + (y2 - y1) * (1 - t));
                image.draw_line(sx, y1, x1, sy, ln);
            }
        };
        int cx = (W - 1) / 2;
        int cy = (H - 1) / 2;
        drawFan(cx, cy, W - 1, H - 1);  // naar rechtsonder
        drawFan(cx, cy, W - 1, 0);      // naar rechtsboven
        drawFan(cx, cy, 0, H - 1);      // naar linksonder
        drawFan(cx, cy, 0, 0);          // naar linksboven
    }

    return image;
}

} // namespace intro
