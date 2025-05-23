#include "fractals3D.h"
#include "Transformations.h" // For scaleFigure, translate, applyTransformation
#include "Figure.h"          // For Figure class definition and center()
#include <vector>
#include <list>              // Figures3D (std::list<Figure>)
#include <iostream>

// Helper function to merge a list of figures into a single figure
static Figure mergeFigures(const Figures3D& figures) {
    Figure merged;
    if (figures.empty()) {
        return merged;
    }

    // Preserve the color of the first figure in the list.
    // This color will likely be overridden by the color specified in the .ini file
    // for the entire FigureX by the generateFigures function in lineDrawer.cpp.
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
            if (!f.point_indexes.empty()) { // Ensure face has points
                newFace.point_indexes.reserve(f.point_indexes.size());
                for (int oldIndex : f.point_indexes) {
                    newFace.point_indexes.push_back(oldIndex + currentPointOffset);
                }
            }
            merged.faces.push_back(newFace);
        }
        // Update the offset for the next figure's points
        currentPointOffset += fig.points.size();
    }
    return merged;
}

// Fractalize function with vertex-to-vertex alignment
static void fractalize(const Figure &parentFig, // Figure for current recursion level
                       Figures3D &resultingFigures,
                       int iterationsRemaining,
                       double perIterationScaleFactor) // This is 1.0 / ini_fractalScale
{
    if (iterationsRemaining <= 0) {
        resultingFigures.push_back(parentFig);
        return;
    }

    if (parentFig.points.empty()) { // Safety check
        if (iterationsRemaining > 0) { // Only warn if we expected to generate children
            std::cerr << "Warning: Attempting to fractalize a figure with no points." << std::endl;
        }
        // If parent has no points, just add it if it's the final step, otherwise stop.
        // This case (parentFig with no points but iterationsRemaining > 0) should ideally not happen
        // if starting from valid base shapes.
        if (iterationsRemaining == 0) resultingFigures.push_back(parentFig);
        return;
    }

    // For each vertex of the parentFig, create a shrunken child and align it.
    for (unsigned int i = 0; i < parentFig.points.size(); ++i) {
        const Vector3D& parentPlacementVertex = parentFig.points[i];

        Figure childFig = parentFig; // Template for the child is the parent itself

        Matrix SMatrix = scaleFigure(perIterationScaleFactor);
        applyTransformation(childFig, SMatrix); // childFig is now shrunken

        // Align the i-th vertex of the shrunken childFig
        // with the i-th vertex of the original parentFig (which is parentPlacementVertex).
        if (i < childFig.points.size()) { // Ensure index 'i' is valid for the (shrunken) child
            const Vector3D& shrunkenChildVertexToAlign = childFig.points[i];
            Vector3D translationVector = parentPlacementVertex - shrunkenChildVertexToAlign;
            Matrix TMatrix = translate(translationVector);
            applyTransformation(childFig, TMatrix);
        } else {
            // Fallback: if index 'i' is somehow out of bounds for the child's points
            // (should not happen if parentFig and childFig are congruent before scaling).
            // Default to centering the shrunken child on the parent's vertex.
            std::cerr << "Warning: Point index " << i << " out of bounds for child figure in fractalize. Using center placement." << std::endl;
            Vector3D shrunkenChildCenter = childFig.center();
            Vector3D translationVector = parentPlacementVertex - shrunkenChildCenter;
            Matrix TMatrix = translate(translationVector);
            applyTransformation(childFig, TMatrix);
        }

        fractalize(childFig, resultingFigures, iterationsRemaining - 1, perIterationScaleFactor);
    }
}



Figure generateFractalTetrahedron(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createTetrahedron();
    // The color will be applied by generateFigures from the .ini file settings
    // baseUnitShape.color = Color(1,0,0); // Default for tetrahedron if needed here

    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }

    return mergeFigures(fractalComponents);
}

Figure generateFractalCube(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createCube();
    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }
    return mergeFigures(fractalComponents);
}

Figure generateFractalIcosahedron(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createIcosahedron();
    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }
    return mergeFigures(fractalComponents);
}

Figure generateFractalOctahedron(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createOctahedron();
    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }
    return mergeFigures(fractalComponents);
}

Figure generateFractalDodecahedron(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createDodecahedron();
    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }
    return mergeFigures(fractalComponents);
}

Figure generateFractalBuckyBall(int nrIterations, double actualScaleFactor) {
    Figure baseUnitShape = Figure::createBuckyBall(); // Assuming createBuckyBall() is implemented
    Figures3D fractalComponents;
    if (nrIterations == 0) {
        fractalComponents.push_back(baseUnitShape);
    } else {
        fractalize(baseUnitShape, fractalComponents, nrIterations, actualScaleFactor);
    }
    return mergeFigures(fractalComponents);
}

Figure generateMengerSponge(int nrIterations) {
    if (nrIterations < 0) nrIterations = 0;

    // For nrIterations = 0, return a simple cube.
    if (nrIterations == 0) {
        Figures3D components;
        components.push_back(Figure::createCube());
        return mergeFigures(components); // or just return Figure::createCube();
    }

    std::list<Figure> currentGenerationCubes;
    currentGenerationCubes.push_back(Figure::createCube()); // Start with one unit cube

    double currentGlobalScale = 1.0; // Represents the scale of cubes in currentGenerationCubes relative to initial unit cube

    for (int iter = 0; iter < nrIterations; ++iter) {
        std::list<Figure> nextGenerationCubes;
        double childGlobalScale = currentGlobalScale / 3.0;
        Matrix childShrinkRelativeToUnit = scaleFigure(childGlobalScale);

        for (const Figure& parentCube : currentGenerationCubes) {
            Vector3D parentCenter = parentCube.center();

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        // Menger sponge condition: remove central cube and face-center cubes
                        if (std::abs(dx) + std::abs(dy) + std::abs(dz) >= 2) {
                            Figure newChildTemplate = Figure::createCube(); // Fresh unit cube

                            // Scale it to the correct absolute size for this iteration
                            applyTransformation(newChildTemplate, childShrinkRelativeToUnit);

                            // Calculate child's center position
                            // Displacement from parent's center, scaled by parent's effective radius contribution
                            Vector3D displacementFromParentCenter = Vector3D::vector(
                                dx * (currentGlobalScale * (2.0 / 3.0)),
                                dy * (currentGlobalScale * (2.0 / 3.0)),
                                dz * (currentGlobalScale * (2.0 / 3.0))
                            );
                            Vector3D childTargetCenter = parentCenter + displacementFromParentCenter;

                            // Translate the correctly-sized child to its target center
                            Vector3D translationForChild = childTargetCenter - newChildTemplate.center(); // newChildTemplate.center() is origin
                            applyTransformation(newChildTemplate, translate(translationForChild));

                            nextGenerationCubes.push_back(newChildTemplate);
                        }
                    }
                }
            }
        }
        currentGenerationCubes = nextGenerationCubes;
        currentGlobalScale = childGlobalScale; // Update for next iteration
    }

    Figures3D resultFigures;
    for(const auto& cube : currentGenerationCubes) {
        resultFigures.push_back(cube);
    }
    return mergeFigures(resultFigures);
}