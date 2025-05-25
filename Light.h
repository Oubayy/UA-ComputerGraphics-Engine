// Light.h
#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H

#include "vector3d.h"
#include "Line2D.h"
#include <list>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Light {
public:
    Color ambientLight;
    Color diffuseLight;
    Color specularLight;

    Light() : ambientLight(0,0,0), diffuseLight(0,0,0), specularLight(0,0,0) {}
    virtual ~Light() = default;

    virtual Vector3D getLightVector(const Vector3D& pointOnSurface_eye) const = 0;

    virtual double getAttenuation(const Vector3D& pointOnSurface_eye, const Vector3D& lightVectorToPointFromSurface_eye) const {
        return 1.0;
    }
};

class DirectionalLight : public Light {
public:
    // This 'direction' will store Ld (vector from light to origin) ALREADY TRANSFORMED TO EYE SPACE.
    Vector3D direction;

    DirectionalLight() : direction(Vector3D::vector(0,0,-1.0)) {} // Default Ld (world), will be transformed

    Vector3D getLightVector(const Vector3D& pointOnSurface_eye) const override {
        // For directional light, L is constant and points TOWARDS the light source.
        // this->direction is Ld in eye space. L_eye = -Ld_eye (normalized).
        return Vector3D::normalise(-this->direction);
    }
};

class PointLight : public Light {
public:
    // This 'location' will be ALREADY TRANSFORMED TO EYE SPACE.
    Vector3D location;
    double spotAngleDegrees;
    // This 'spotDirection' will be ALREADY TRANSFORMED/CALCULATED IN EYE SPACE.
    Vector3D spotDirection;

    PointLight() :
        location(Vector3D::point(0,0,0)),
        spotAngleDegrees(181.0),
        spotDirection(Vector3D::vector(0,0,-1.0))
    {}

    Vector3D getLightVector(const Vector3D& pointOnSurface_eye) const override {
        // Both location and pointOnSurface_eye are in eye space.
        return Vector3D::normalise(this->location - pointOnSurface_eye);
    }

    double getAttenuation(const Vector3D& pointOnSurface_eye, const Vector3D& lightVectorToPointFromSurface_eye) const override {
        if (spotAngleDegrees <= 0.0 || spotAngleDegrees >= 90.0) {
            return 1.0;
        }
        Vector3D vectorFromLightToSurface_eye = -lightVectorToPointFromSurface_eye;
        double cosAngle = this->spotDirection.dot(vectorFromLightToSurface_eye); // Both in eye space, spotDirection normalized
        double cosCutoff = std::cos(spotAngleDegrees * M_PI / 180.0);
        if (cosAngle >= cosCutoff) return 1.0;
        else return 0.0;
    }
};

typedef std::list<Light*> Lights3D;

#endif // ENGINE_LIGHT_H