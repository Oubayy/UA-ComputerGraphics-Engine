#include "fractals3D.h"
#include "Transformations.h"

static void fractalize(Figure &base, Figures3D &result, int iterations, double scaleFactor) {
    if (iterations == 0) {
        result.push_back(base);
        return;
    }

    for (const auto &point : base.points) {
        Figure copy = base;

        Matrix scale = scaleFigure(scaleFactor);
        applyTransformation(copy, scale);

        Vector3D translationVector = point - Vector3D::point(0, 0, 0); // 🔥 Fix here!
        Matrix translation = translate(translationVector);
        applyTransformation(copy, translation);

        fractalize(copy, result, iterations - 1, scaleFactor);
    }
}


Figure generateMengerSponge(int nrIterations) {
    Figure cube = Figure::createCube();
    Figures3D fractalResult;
    fractalize(cube, fractalResult, nrIterations, 1.0 / 3.0);

    // Merge all Figures into one big Figure
    Figure merged;
    for (const auto &fig : fractalResult) {
        int pointOffset = merged.points.size();
        for (const auto &p : fig.points)
            merged.points.push_back(p);

        for (const auto &f : fig.faces) {
            Face newFace;
            for (int idx : f.point_indexes)
                newFace.point_indexes.push_back(idx + pointOffset);
            merged.faces.push_back(newFace);
        }
    }

    return merged;
}

Figure generateFractalTetrahedron(int nrIterations) {
    Figure tetra = Figure::createTetrahedron();
    Figures3D fractalResult;
    fractalize(tetra, fractalResult, nrIterations, 0.5);

    Figure merged;
    for (const auto &fig : fractalResult) {
        int pointOffset = merged.points.size();
        for (const auto &p : fig.points)
            merged.points.push_back(p);

        for (const auto &f : fig.faces) {
            Face newFace;
            for (int idx : f.point_indexes)
                newFace.point_indexes.push_back(idx + pointOffset);
            merged.faces.push_back(newFace);
        }
    }

    return merged;
}

Figure generateFractalIcosahedron(int nrIterations) {
    Figure icosahedron = Figure::createIcosahedron();
    Figures3D fractalResult;
    fractalize(icosahedron, fractalResult, nrIterations, 0.5);

    Figure merged;
    for (const auto &fig : fractalResult) {
        int pointOffset = merged.points.size();
        for (const auto &p : fig.points)
            merged.points.push_back(p);

        for (const auto &f : fig.faces) {
            Face newFace;
            for (int idx : f.point_indexes)
                newFace.point_indexes.push_back(idx + pointOffset);
            merged.faces.push_back(newFace);
        }
    }

    return merged;
}

Figure generateFractalCube(int nrIterations) {
    Figure cube = Figure::createCube();
    Figures3D fractalResult;
    fractalize(cube, fractalResult, nrIterations, 0.5);

    Figure merged;
    for (const auto &fig : fractalResult) {
        int pointOffset = merged.points.size();
        for (const auto &p : fig.points)
            merged.points.push_back(p);

        for (const auto &f : fig.faces) {
            Face newFace;
            for (int idx : f.point_indexes)
                newFace.point_indexes.push_back(idx + pointOffset);
            merged.faces.push_back(newFace);
        }
    }

    return merged;
}

Figure generateFractalOctahedron(int nrIterations) {
    Figure octa = Figure::createOctahedron();
    Figures3D fractalResult;
    fractalize(octa, fractalResult, nrIterations, 0.5);

    Figure merged;
    for (const auto &fig : fractalResult) {
        int pointOffset = merged.points.size();
        for (const auto &p : fig.points)
            merged.points.push_back(p);

        for (const auto &f : fig.faces) {
            Face newFace;
            for (int idx : f.point_indexes)
                newFace.point_indexes.push_back(idx + pointOffset);
            merged.faces.push_back(newFace);
        }
    }

    return merged;
}
