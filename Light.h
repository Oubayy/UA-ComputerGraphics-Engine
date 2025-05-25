// Light.h
#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H

#include "vector3d.h"
#include "Line2D.h" // For Color
#include <list>
#include <cmath> // For M_PI, cos, std::pow, std::cos

#ifndef M_PI // Define M_PI if not already defined by cmath
#define M_PI 3.14159265358979323846
#endif

class Light {
public:
    Color ambientLight;    // Ia
    Color diffuseLight;    // Id
    Color specularLight;   // Is

    Light() : ambientLight(0,0,0), diffuseLight(0,0,0), specularLight(0,0,0) {}
    virtual ~Light() = default;

    virtual Vector3D getLightVector(const Vector3D& pointOnSurface) const = 0;

    virtual double getAttenuation(const Vector3D& pointOnSurface, const Vector3D& lightVectorToPointFromSurface) const {
        return 1.0;
    }
};

class DirectionalLight : public Light {
public:
    Vector3D direction; // Ld (points FROM light TO origin)

    DirectionalLight() : direction(Vector3D::vector(0,0,-1.0)) {}

    Vector3D getLightVector(const Vector3D& pointOnSurface) const override {
        // L points TOWARDS the light source
        return Vector3D::normalise(-direction);
    }
};

class PointLight : public Light {
public:
    Vector3D location;
    double spotAngleDegrees;
    Vector3D spotDirection;  // Normalized direction the spotlight is pointing.

    PointLight() :
        location(Vector3D::point(0,0,0)),
        spotAngleDegrees(181.0), // Default: omni
        spotDirection(Vector3D::vector(0,0,-1.0))
    {}

    Vector3D getLightVector(const Vector3D& pointOnSurface) const override {
        return Vector3D::normalise(location - pointOnSurface);
    }

    double getAttenuation(const Vector3D& pointOnSurface, const Vector3D& lightVectorToPointFromSurface) const override {
        if (spotAngleDegrees <= 0.0 || spotAngleDegrees >= 90.0) { // Per spec PDF: angle between 0 and 90 inclusive. ">90" means omni.
            return 1.0; // Omni-directional
        }

        // lightVectorToPointFromSurface is L (from surface to light).
        // We need vector from light TO surface for spotlight calculation.
        Vector3D vectorFromLightToSurface = -lightVectorToPointFromSurface; // This is already normalized if L was.

        // spotDirection should be normalized when set.
        double cosAngle = spotDirection.dot(vectorFromLightToSurface);
        double cosCutoff = std::cos(spotAngleDegrees * M_PI / 180.0);

        if (cosAngle >= cosCutoff) {
            // Optional: smooth falloff (e.g., cosine weighted or phong-like)
            // For now, hard cutoff as per basic requirement
            return 1.0; // Point is inside the cone
        } else {
            return 0.0; // Point is outside the cone
        }
    }
};

typedef std::list<Light*> Lights3D;

#endif // ENGINE_LIGHT_H