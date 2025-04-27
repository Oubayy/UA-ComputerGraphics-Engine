#include "fractals3D.h"
#include "Transformations.h" // For scaleFigure, translate, applyTransformation
#include "Figure.h"          // For Figure class definition and center()

// Helper function to merge a list of figures into a single figure
static Figure mergeFigures(const Figures3D& figures) {
    Figure merged;
    if (figures.empty()) return merged;

    // Preserve the color of the first figure in the list (assuming they might be the same)
    merged.color = figures.front().color;

    int currentPointOffset = 0;
    for (const auto &fig : figures) {
        // Append points from the current figure
        for (const auto &p : fig.points) {
            merged.points.push_back(p);
        }
        // Append faces from the current figure, adjusting indices
        for (const auto &f : fig.faces) {
            Face newFace;
            newFace.point_indexes.reserve(f.point_indexes.size()); // Optimize vector allocation
            for (int oldIndex : f.point_indexes) {
                newFace.point_indexes.push_back(oldIndex + currentPointOffset);
            }
            merged.faces.push_back(newFace);
        }
        // Update the offset for the next figure's points
        currentPointOffset += fig.points.size();
    }
    return merged;
}


// The recursive fractal generation function (Sierpinski-style replacement)
static void fractalize(const Figure &baseFigure, Figures3D &resultingFigures, int iterationsRemaining, double scaleFactor) {
    // Base case: If no more iterations, add the current figure to the results list
    if (iterationsRemaining <= 0) {
        resultingFigures.push_back(baseFigure);
        return;
    }

    // Recursive step: For each vertex in the baseFigure...
    for (const Vector3D &vertex : baseFigure.points) {
        // 1. Create a copy of the base figure
        Figure copy = baseFigure;

        // 2. Scale the copy
        Matrix scaleMatrix = scaleFigure(scaleFactor);
        applyTransformation(copy, scaleMatrix); // Apply scaling to the copy

        // 3. Calculate the translation vector
        // We want to move the scaled copy so its center aligns with the original vertex
        Vector3D scaledCenter = copy.center(); // Get the center *after* scaling
        Vector3D translationVector = vertex - scaledCenter; // Vector from scaled center to original vertex

        // 4. Translate the scaled copy
        Matrix translationMatrix = translate(translationVector);
        applyTransformation(copy, translationMatrix); // Apply translation

        // 5. Recursively call fractalize on the transformed copy
        // Decrease the iteration count for the recursive call
        fractalize(copy, resultingFigures, iterationsRemaining - 1, scaleFactor);
    }
    // Note: The original baseFigure is implicitly replaced by the collection of smaller figures
    // generated in the recursive calls and added to resultingFigures.
}


// --- Implementations of the public generator functions ---

// Menger Sponge - Special case, potentially needs its own algorithm
// For now, uses the 'fractalize' method which isn't strictly correct for Menger
Figure generateMengerSponge(int nrIterations) {
    Figure baseCube = Figure::createCube();
    Figures3D fractalResult;
    double actualScaleFactor = 1.0 / 3.0; // Menger uses 1/3 scale

    // Warning: This `fractalize` creates replacements at vertices, not the hole-punching Menger algorithm.
    // If a true Menger Sponge is required, this needs a different implementation.
    fractalize(baseCube, fractalResult, nrIterations, actualScaleFactor);

    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = baseCube.color; // Assign color explicitly
    return finalFigure;
    // TODO: Implement true Menger Sponge algorithm if needed.
}

Figure generateFractalTetrahedron(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createTetrahedron();
    Figures3D fractalResult;
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}

Figure generateFractalIcosahedron(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createIcosahedron();
    Figures3D fractalResult;
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}

Figure generateFractalCube(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createCube();
    Figures3D fractalResult;
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}

Figure generateFractalOctahedron(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createOctahedron();
    Figures3D fractalResult;
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}

Figure generateFractalDodecahedron(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createDodecahedron();
    Figures3D fractalResult;
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}

Figure generateFractalBuckyBall(int nrIterations, double actualScaleFactor) {
    Figure base = Figure::createBuckyBall(); // Use the (placeholder) Buckyball
    Figures3D fractalResult;
    // Use the same fractalize logic
    fractalize(base, fractalResult, nrIterations, actualScaleFactor);
    Figure finalFigure = mergeFigures(fractalResult);
    finalFigure.color = base.color; // Assign color
    return finalFigure;
}