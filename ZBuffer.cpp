#include "ZBuffer.h"
#include "Figure.h"
#include "Projection.h"   // Should be fine
#include "vector3d.h"
#include "easy_image.h"
#include "Line2D.h"
#include "Light.h"        // Included in ZBuffer.h
#include <limits>
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

ZBuffer::ZBuffer(const int width, const int height) {
    double negInf = std::numeric_limits<double>::lowest();
    if (width <= 0 || height <= 0) {
        int valid_width = std::max(0, width);
        int valid_height = std::max(0, height);
        this->resize(valid_width);
        if (valid_width > 0) {
            for (int i = 0; i < valid_width; ++i) {
                this->at(i).resize(valid_height, negInf);
            }
        }
        std::cerr << "Warning: ZBuffer created with non-positive dimensions (" << width << "x" << height << "). Size set to " << valid_width << "x" << valid_height << "." << std::endl;
        return;
    }
    this->resize(width);
    for (int i = 0; i < width; ++i) {
        this->at(i).resize(height, negInf);
    }
}

Triangles ZBuffer::doProjectTriangle(const Figures3D &figures, const double &d_projection_plane) {
    Triangles triangles;
    for (const auto &figure : figures) {
        if (figure.faces.empty() || figure.points.empty()) continue;

        FigurePropertiesForTriangle currentFigMaterial;
        currentFigMaterial.ambientReflection = figure.ambientReflection;
        currentFigMaterial.diffuseReflection = figure.diffuseReflection;
        currentFigMaterial.specularReflection = figure.specularReflection;
        currentFigMaterial.reflectionCoefficient = figure.reflectionCoefficient;

        for (const auto &face : figure.faces) {
            if (face.point_indexes.size() < 3) continue;

            const Vector3D &P1_3D_eye = figure.points[face.point_indexes[0]];
            const Vector3D &P2_3D_eye = figure.points[face.point_indexes[1]];
            const Vector3D &P3_3D_eye = figure.points[face.point_indexes[2]];

            Point2D p1_2D, p2_2D, p3_2D;
            double z1_val, z2_val, z3_val;

            auto project_point = [&](const Vector3D& p3d_eye, Point2D& p2d, double& z_val) {
                if (p3d_eye.z >= 0) {
                    p2d.x = (p3d_eye.x > 0) ? 1e6 : -1e6; p2d.y = (p3d_eye.y > 0) ? 1e6 : -1e6;
                    z_val = std::numeric_limits<double>::infinity();
                } else {
                    p2d.x = d_projection_plane * p3d_eye.x / -p3d_eye.z;
                    p2d.y = d_projection_plane * p3d_eye.y / -p3d_eye.z;
                    z_val = -p3d_eye.z;
                }
            };
            project_point(P1_3D_eye, p1_2D, z1_val);
            project_point(P2_3D_eye, p2_2D, z2_val);
            project_point(P3_3D_eye, p3_2D, z3_val);

            if (z1_val != std::numeric_limits<double>::infinity() &&
                z2_val != std::numeric_limits<double>::infinity() &&
                z3_val != std::numeric_limits<double>::infinity()) {
                triangles.emplace_back(
                    p1_2D, p2_2D, p3_2D,
                    z1_val, z2_val, z3_val,
                    P1_3D_eye, P2_3D_eye, P3_3D_eye,
                    currentFigMaterial
                );
            }
        }
    }
    return triangles;
}

void ZBuffer::draw_zbuf_triangle(
    img::EasyImage &image,
    const Triangle& triangle,
    const Lights3D& lights,
    const Vector3D& eyePosition_eyeSpace, // Should be (0,0,0)
    double d_projection_param,     // 'd' from x' = d*x/-z
    double image_scale_factor,     // 'd' from calculate() screen scaling
    double image_dx_translation,   // 'dx' from calculate()
    double image_dy_translation)   // 'dy' from calculate()
{
    const auto& mat = triangle.material;
    Vector3D N = Vector3D::cross(triangle.p2_eye_space - triangle.p1_eye_space, triangle.p3_eye_space - triangle.p1_eye_space);
    N.normalise();

    Vector3D p_on_triangle_eye = (triangle.p1_eye_space + triangle.p2_eye_space + triangle.p3_eye_space) / 3.0;
    if (N.dot(eyePosition_eyeSpace - p_on_triangle_eye) < 0) {
         N = -N;
    }

    Point2D p1_screen = triangle.p1; Point2D p2_screen = triangle.p2; Point2D p3_screen = triangle.p3;
    double inv_z1_eye = (triangle.z1_eye > 1e-9) ? 1.0 / triangle.z1_eye : std::numeric_limits<double>::lowest();
    double inv_z2_eye = (triangle.z2_eye > 1e-9) ? 1.0 / triangle.z2_eye : std::numeric_limits<double>::lowest();
    double inv_z3_eye = (triangle.z3_eye > 1e-9) ? 1.0 / triangle.z3_eye : std::numeric_limits<double>::lowest();

    double G_x_screen = (p1_screen.x + p2_screen.x + p3_screen.x) / 3.0;
    double G_y_screen = (p1_screen.y + p2_screen.y + p3_screen.y) / 3.0;
    double G_inv_z_eye = (inv_z1_eye + inv_z2_eye + inv_z3_eye) / 3.0;

    Vector3D vA_s = Vector3D::point(p1_screen.x, p1_screen.y, inv_z1_eye);
    Vector3D vB_s = Vector3D::point(p2_screen.x, p2_screen.y, inv_z2_eye);
    Vector3D vC_s = Vector3D::point(p3_screen.x, p3_screen.y, inv_z3_eye);
    Vector3D U_s = vB_s - vA_s; Vector3D V_s = vC_s - vA_s;
    Vector3D W_s = Vector3D::cross(U_s, V_s);
    double dz_inv_dx_s, dz_inv_dy_s;
    if (std::abs(W_s.z) < 1e-9) return;
    dz_inv_dx_s = -W_s.x / W_s.z; dz_inv_dy_s = -W_s.y / W_s.z;

    long y_min = lround(std::min({p1_screen.y, p2_screen.y, p3_screen.y}) + 0.5);
    long y_max = lround(std::max({p1_screen.y, p2_screen.y, p3_screen.y}) - 0.5);
    long x_min_overall = lround(std::min({p1_screen.x, p2_screen.x, p3_screen.x}) + 0.5);
    long x_max_overall = lround(std::max({p1_screen.x, p2_screen.x, p3_screen.x}) - 0.5);

    int img_w = image.get_width(); int img_h = image.get_height();
    y_min = std::max(0L, y_min); y_max = std::min((long)img_h - 1, y_max);
    x_min_overall = std::max(0L, x_min_overall); x_max_overall = std::min((long)img_w - 1, x_max_overall);

    for (long y_px = y_min; y_px <= y_max; ++y_px) {
        double x_int[3]; int n_int = 0;
        auto check_edge = [&](Point2D e1, Point2D e2){
            if((e1.y <= y_px && e2.y > y_px) || (e2.y <= y_px && e1.y > y_px)){
                if(std::abs(e2.y - e1.y)>1e-9) x_int[n_int++] = e1.x + (e2.x-e1.x)*(y_px-e1.y)/(e2.y-e1.y);
            }
        };
        check_edge(p1_screen, p2_screen); check_edge(p2_screen, p3_screen); check_edge(p3_screen, p1_screen);
        if (n_int < 2) continue;
        long x_L = lround(std::min(x_int[0], x_int[1]) + 0.5);
        long x_R = lround(std::max(x_int[0], x_int[1]) - 0.5);
        if (n_int == 3) { x_L = lround(std::min({x_int[0],x_int[1],x_int[2]}) + 0.5); x_R = lround(std::max({x_int[0],x_int[1],x_int[2]}) - 0.5); }
        x_L = std::max(x_min_overall, x_L); x_R = std::min(x_max_overall, x_R);

        for (long x_px = x_L; x_px <= x_R; ++x_px) {
            double current_inv_Z = G_inv_z_eye + (x_px - G_x_screen) * dz_inv_dx_s + (y_px - G_y_screen) * dz_inv_dy_s;
            if (x_px >=0 && x_px < img_w && y_px >=0 && y_px < img_h && current_inv_Z > this->at(x_px)[y_px]) {
                this->at(x_px)[y_px] = current_inv_Z;
                double Z_eye = 1.0 / current_inv_Z;
                double xp = (x_px - image_dx_translation) / image_scale_factor;
                double yp = (y_px - image_dy_translation) / image_scale_factor;
                Vector3D P_surf_eye = Vector3D::point(xp*Z_eye/d_projection_param, yp*Z_eye/d_projection_param, -Z_eye);
                Vector3D V_vec = Vector3D::normalise(eyePosition_eyeSpace - P_surf_eye);
                Color totalColor(0,0,0);

                for (const Light* light_ptr : lights) {
                    const Light& light = *light_ptr; Color lightCont(0,0,0);
                    lightCont.red += light.ambientLight.red * mat.ambientReflection.red;
                    lightCont.green += light.ambientLight.green * mat.ambientReflection.green;
                    lightCont.blue += light.ambientLight.blue * mat.ambientReflection.blue;
                    Vector3D L_vec = light.getLightVector(P_surf_eye); // Normalized
                    double spot = light.getAttenuation(P_surf_eye, L_vec);
                    if (spot <= 0.0) { totalColor.red+=lightCont.red; totalColor.green+=lightCont.green; totalColor.blue+=lightCont.blue; continue; }
                    double NdotL = std::max(0.0, N.dot(L_vec));
                    if (NdotL > 0.0) {
                        lightCont.red   += light.diffuseLight.red   * mat.diffuseReflection.red   * NdotL * spot;
                        lightCont.green += light.diffuseLight.green * mat.diffuseReflection.green * NdotL * spot;
                        lightCont.blue  += light.diffuseLight.blue  * mat.diffuseReflection.blue  * NdotL * spot;
                        Vector3D R_vec = Vector3D::normalise((2.0 * NdotL * N) - L_vec);
                        double RdotV = std::max(0.0, R_vec.dot(V_vec));
                        if (RdotV > 0.0 && mat.reflectionCoefficient > 0) {
                            double specF = std::pow(RdotV, mat.reflectionCoefficient);
                            lightCont.red   += light.specularLight.red   * mat.specularReflection.red   * specF * spot;
                            lightCont.green += light.specularLight.green * mat.specularReflection.green * specF * spot;
                            lightCont.blue  += light.specularLight.blue  * mat.specularReflection.blue  * specF * spot;
                        }
                    }
                    totalColor.red+=lightCont.red; totalColor.green+=lightCont.green; totalColor.blue+=lightCont.blue;
                }
                totalColor.red=min(1.0,max(0.0,totalColor.red)); totalColor.green=min(1.0,max(0.0,totalColor.green)); totalColor.blue=min(1.0,max(0.0,totalColor.blue));
                image(x_px,y_px) = img::Color(round(totalColor.red*255), round(totalColor.green*255), round(totalColor.blue*255));
            }
        }
    }
}

void ZBuffer::draw_zbuf_line(img::EasyImage &image, Line2D line) {
    double x0_d = line.p1.x; double y0_d = line.p1.y;
    double x1_d = line.p2.x; double y1_d = line.p2.y;
    double inv_z0 = (line.z1 > 1e-9 && !std::isinf(line.z1)) ? 1.0 / line.z1 : std::numeric_limits<double>::lowest();
    double inv_z1 = (line.z2 > 1e-9 && !std::isinf(line.z2)) ? 1.0 / line.z2 : std::numeric_limits<double>::lowest();
    img::Color color_img(round(line.color.red*255), round(line.color.green*255), round(line.color.blue*255));
    int x0_i=round(x0_d), y0_i=round(y0_d), x1_i=round(x1_d), y1_i=round(y1_d);

    if ((x0_i > x1_i) || ((x0_i == x1_i) && (y0_i > y1_i))) {
        std::swap(x0_d, x1_d);
        std::swap(y0_d, y1_d);
        std::swap(inv_z0, inv_z1);
        x0_i = round(x0_d); y0_i = round(y0_d);
        x1_i = round(x1_d); y1_i = round(y1_d);
    }

    int imgW = image.get_width(), imgH = image.get_height();
    int dx_i = x1_i - x0_i, dy_i = y1_i - y0_i;

    if(dx_i == 0 && dy_i == 0){
        if(x0_i >= 0 && x0_i < imgW && y0_i >= 0 && y0_i < imgH && inv_z0 > this->at(x0_i)[y0_i]){
            image(x0_i,y0_i) = color_img;
            this->at(x0_i)[y0_i] = inv_z0;
        }
    } else if(dx_i == 0){
        for(int y = y0_i; y <= y1_i; ++y){
            if(x0_i >= 0 && x0_i < imgW && y >= 0 && y < imgH){
                double t = (dy_i == 0) ? 0.0 : static_cast<double>(y - y0_i) / dy_i;
                double cur_invZ = inv_z0 + t * (inv_z1 - inv_z0);
                if(cur_invZ > this->at(x0_i)[y]){
                    image(x0_i,y) = color_img;
                    this->at(x0_i)[y] = cur_invZ;
                }
            }
        }
    } else if(dy_i == 0){
        for(int x = x0_i; x <= x1_i; ++x){
            if(x >= 0 && x < imgW && y0_i >= 0 && y0_i < imgH){
                double t = (dx_i == 0) ? 0.0 : static_cast<double>(x - x0_i) / dx_i;
                double cur_invZ = inv_z0 + t * (inv_z1 - inv_z0);
                if(cur_invZ > this->at(x)[y0_i]){
                    image(x,y0_i) = color_img;
                    this->at(x)[y0_i] = cur_invZ;
                }
            }
        }
    } else {
        double m = static_cast<double>(dy_i) / dx_i;
        if(-1.0 <= m && m <= 1.0){
            for(int x = x0_i; x <= x1_i; ++x){
                int y = round(y0_d + m * (x - x0_i));
                if(x >= 0 && x < imgW && y >= 0 && y < imgH){
                    double t = (dx_i == 0) ? 0.5 : static_cast<double>(x - x0_i) / dx_i;
                    double cur_invZ = inv_z0 + t * (inv_z1 - inv_z0);
                    if(cur_invZ > this->at(x)[y]){
                        image(x,y) = color_img;
                        this->at(x)[y] = cur_invZ;
                    }
                }
            }
        } else {
            for(int y = y0_i; ; (y0_i <= y1_i ? y++ : y--)){
                int x = round(x0_d + (y - y0_d) / m);
                if(x >= 0 && x < imgW && y >= 0 && y < imgH){
                    double t = (dy_i == 0) ? 0.5 : static_cast<double>(y - y0_i) / dy_i;
                    double cur_invZ = inv_z0 + t * (inv_z1 - inv_z0);
                    if(cur_invZ > this->at(x)[y]){
                        image(x,y) = color_img;
                        this->at(x)[y] = cur_invZ;
                    }
                }
                if(y == y1_i) break;
            }
        }
    }
}