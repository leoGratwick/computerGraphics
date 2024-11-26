#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Utils.h>
#include <iostream>
#include <TextureMap.h>
#include <vector>
#include <limits>
#include <bits/stdc++.h>
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>
#include "ModelTriangle.h"
#include "Helpers.h"
#include "RayTriangleIntersection.h"
#include "ModelLoader.h"
#include "Camera.h"
#include "LightSource.h"
#include "Model.h"
#include "World.h"

#define WIDTH 900
#define HEIGHT 900
#define INFINITY std::numeric_limits<float>::infinity()
#define FRAME_RATE 60

// GLOBAL VARS

// window
bool windowOpen = true;
int frameNumber = 0;

// world
glm::vec3 worldOrign(0, 0, 0);
World world = World();


// camera
bool lookAt = false;
bool record = false;

// for rotation
float theta = 0;
bool moveLight = false;

// depth buffer
std::vector<std::vector<float> > zDepth(WIDTH, std::vector<float>(HEIGHT, -INFINITY));


void drawCanvasPoint(CanvasPoint point, Colour col, DrawingWindow &window) {
    if (pointInCanvas(point, window)) {
        // round coordinates to pixel values
        int x = std::round(point.x);
        int y = std::round(point.y);

        uint32_t colour = colourToInt(col);

        // check depth buffer
        if (point.depth < 0.00001f && point.depth > -0.00001f) {
            window.setPixelColour(x, y, colour);
            zDepth[x][y] = 0;
        }
        // check if point is closer than current point on that pixel
        else if ((zDepth[x][y] < 1 / point.depth) && !(zDepth[x][y] == 0)) {
            window.setPixelColour(x, y, colour);
            zDepth[x][y] = 1 / point.depth;
        }
    }
}

void drawCanvasPointUint(CanvasPoint point, uint32_t col, DrawingWindow &window) {
    if (pointInCanvas(point, window)) {
        // round coordinates to pixel values
        int x = std::round(point.x);
        int y = std::round(point.y);

        // check depth buffer
        if (fabs(point.depth) < 1e-6f) {
            return; // Skip this point
        }

        float inverseDepth = 1.0f / point.depth;

        if (zDepth[x][y] == 0 || inverseDepth < zDepth[x][y]) {
            window.setPixelColour(x, y, col);
            zDepth[x][y] = inverseDepth; // Update the depth buffer
        }
    }
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    // if the line is a point
    if (from.x == to.x and to.y == from.y) {
        drawCanvasPoint(from, colour, window);
        return;
    }

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;

    // number of points to be drawn, multiply by 2 to avoid gaps
    int steps = std::max(std::abs(xDiff), std::abs(yDiff)) * 2;
    float xStep = xDiff / steps;
    float yStep = yDiff / steps;

    // interpolate point depths
    std::vector<float> lerpDepths = lerpSingleFloats(from.depth, to.depth, steps);

    // draw each point on the line
    for (int i = 0; i < steps; i++) {
        CanvasPoint point = CanvasPoint(from.x + i * xStep, from.y + i * yStep);
        point.depth = lerpDepths.at(i);
        drawCanvasPoint(point, colour, window);
    }
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v0(), triangle.v2(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
}

void fillFlatBottomTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // outline triangle to avoid gaps
    drawTriangle(window, triangle, colour);
    // the v0 is at the top of the triangle
    CanvasPoint v0 = triangle.v0();
    CanvasPoint v1 = triangle.v1();
    CanvasPoint v2 = triangle.v2();

    // calculate x slope dx/dy
    float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
    float slope2 = (v2.x - v0.x) / (v2.y - v0.y);

    float currentX1 = v0.x;
    float currentX2 = v0.x;

    // interpolate depths
    std::vector<float> depths1 = lerpSingleFloats(v0.depth, v1.depth, std::round(v1.y) - std::round(v0.y));
    std::vector<float> depths2 = lerpSingleFloats(v0.depth, v2.depth, std::round(v1.y) - std::round(v0.y));

    int depthCount = 0;

    // draw a line from currentX1 to currentX2 for each y value
    for (int y = std::round(v0.y); y < std::round(v1.y); y++) {
        CanvasPoint from = CanvasPoint(currentX1, y);
        from.depth = depths1.at(depthCount);

        CanvasPoint to = CanvasPoint(currentX2, y);
        to.depth = depths2.at(depthCount);

        drawLine(window, from, to, colour);

        // update vars
        currentX1 += slope1;
        currentX2 += slope2;
        depthCount += 1;
    }
}

void fillFlatTopTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // draw triangle outline to avoid gaps
    drawTriangle(window, triangle, colour);

    // v2 is the bottom vertex
    CanvasPoint v0 = triangle.v0();
    CanvasPoint v1 = triangle.v1();
    CanvasPoint v2 = triangle.v2();

    // calculate x slope dx/dy
    float slope1 = (v0.x - v2.x) / (v0.y - v2.y);
    float slope2 = (v1.x - v2.x) / (v1.y - v2.y);

    float currentX1 = v2.x;
    float currentX2 = v2.x;

    // interpolate depths
    std::vector<float> depths1 = lerpSingleFloats(v2.depth, v0.depth, std::round(v2.y) - std::round(v1.y));
    std::vector<float> depths2 = lerpSingleFloats(v2.depth, v1.depth, std::round(v2.y) - std::round(v1.y));

    int depthCount = 0;

    // draw a line from currentX1 to currentX2 for each y value
    for (int y = std::round(v2.y); y > std::round(v1.y); y--) {
        CanvasPoint from = CanvasPoint(currentX1, y);
        from.depth = depths1.at(depthCount);

        CanvasPoint to = CanvasPoint(currentX2, y);
        to.depth = depths2.at(depthCount);

        drawLine(window, from, to, colour);

        // update vars
        currentX1 -= slope1;
        currentX2 -= slope2;
        depthCount += 1;
    }
}

void fillTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    std::vector<CanvasPoint> verts = {triangle.v0(), triangle.v1(), triangle.v2()};

    // draw triangle outline to avoid gaps
    drawTriangle(window, triangle, colour);

    // sort verticies by y coordinate
    std::sort(verts.begin(), verts.end(), [](const CanvasPoint &a, const CanvasPoint &b) {
        return a.y < b.y; // Sort by y
    });

    CanvasPoint v0 = verts[0];;
    CanvasPoint v1 = verts[1];
    CanvasPoint v2 = verts[2];

    // flat bottom triangle
    if (v2.y == v1.y) {
        CanvasTriangle tri = CanvasTriangle(v0, v1, v2);
        fillFlatBottomTriangle(window, tri, colour);
    }
    // flat top triangle
    else if (v0.y == v1.y) {
        CanvasTriangle tri = CanvasTriangle(v0, v1, v2);
        fillFlatTopTriangle(window, tri, colour);
    }
    // all other triangles
    else {
        // calculate the extra point to split the triangle into a flat bottom and flat top triangle
        float m = (v2.y - v0.y) / (v2.x - v0.x);
        float x = (v1.y - v0.y) / m + v0.x;
        CanvasPoint v3 = CanvasPoint(x, v1.y);

        // calculate depth of extra point
        float sp = splitPercent(v0, v2, v3);
        v3.depth = v0.depth + (v2.depth - v0.depth) * sp;

        // draw flat-bottomed triangle
        CanvasTriangle flatBottomTri = CanvasTriangle(v0, v1, v3);
        fillFlatBottomTriangle(window, flatBottomTri, colour);

        // draw flat-topped triangle
        CanvasTriangle flatTopTri = CanvasTriangle(v1, v3, v2);
        fillFlatTopTriangle(window, flatTopTri, colour);
    }
}

void randomTriangle(DrawingWindow &window) {
    // draw a random triangle with white outline

    CanvasPoint v0 = CanvasPoint(rand() % (window.width - 1), rand() % (window.height - 1));
    CanvasPoint v1 = CanvasPoint(rand() % (window.width - 1), rand() % (window.height - 1));
    CanvasPoint v2 = CanvasPoint(rand() % (window.width - 1), rand() % (window.height - 1));
    Colour randomCol = Colour(rand() % 255, rand() % 255, rand() % 255);

    CanvasTriangle randomTriangle = CanvasTriangle(v0, v1, v2);

    fillTriangle(window, randomTriangle, randomCol);
    drawTriangle(window, randomTriangle, Colour(255, 255, 255));
    window.renderFrame();
}

void textureFlatBottomTriangle(DrawingWindow &window, CanvasTriangle drawingTriangle, TextureMap textureMap,
                               CanvasTriangle textureTriangle) {
    // v0 is the top of the triangle
    CanvasPoint v0 = drawingTriangle.v0();
    CanvasPoint v1 = drawingTriangle.v1();
    CanvasPoint v2 = drawingTriangle.v2();

    // calculate x slope dx/dy
    float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
    float slope2 = (v2.x - v0.x) / (v2.y - v0.y);

    float currentX1 = v0.x;
    float currentX2 = v0.x;

    int numPixels;
    float ratio1;
    float ratio2;
    CanvasPoint textureFrom;
    CanvasPoint textureTo;

    // iterate through each horizontal row of the canvas triangle
    for (int y = v0.y; y <= v1.y; y++) {
        // for first iteration draw 1 pixel
        if (currentX1 == currentX2) {
            uint32_t col = getColourFromTexture(textureTriangle.v0(), textureMap);
            drawCanvasPointUint(CanvasPoint(currentX1, y), col, window);

            //update X1 and X2
            currentX1 += slope1;
            currentX2 += slope2;
            continue;
        }

        CanvasPoint from = CanvasPoint(currentX1, y);
        CanvasPoint to = CanvasPoint(currentX2, y);

        // number of pixels between currentX1 and currentX2
        numPixels = std::abs(std::round(currentX1) - std::round(currentX2)) + 1;

        // ratio of currentX1 and currentx2 along v0v1 and v0v2 respectively
        ratio1 = splitPercent(v0, v1, from);
        ratio2 = splitPercent(v0, v2, to);

        //calculate the points equivalent to from and to in the texture
        textureFrom = getPointAlongLine(textureTriangle.v0(), textureTriangle.v1(), ratio1);
        textureTo = getPointAlongLine(textureTriangle.v0(), textureTriangle.v2(), ratio2);

        // vector of numPixel points along a line on the texture
        std::vector<CanvasPoint> texturePixels = lerpCanvasPoints(textureFrom, textureTo, numPixels);
        std::vector<uint32_t> texturePixelsColours = {};

        // find colour of each point on the texture
        for (int i = 0; i < numPixels; i++) {
            uint32_t colour = getColourFromTexture(texturePixels.at(i), textureMap);
            texturePixelsColours.push_back(colour);
        }

        // draw each pixel in the row with the colour from texture pixels colours
        int j = 0;
        for (int x = std::round(currentX1); x <= std::round(currentX2); x++) {
            drawCanvasPointUint(CanvasPoint(x, y), texturePixelsColours[j], window);
            j += 1;
        }


        // update X1 and X2
        currentX1 += slope1;
        currentX2 += slope2;
    }
}

void textureFlatTopTriangle(DrawingWindow &window, CanvasTriangle drawingTriangle, TextureMap textureMap,
                            CanvasTriangle textureTriangle) {
    // v2 is the bottom of the triangle
    CanvasPoint v0 = drawingTriangle.v0();
    CanvasPoint v1 = drawingTriangle.v1();
    CanvasPoint v2 = drawingTriangle.v2();

    // calculate x slope dx/dy
    float slope1 = (v0.x - v2.x) / (v0.y - v2.y);
    float slope2 = (v1.x - v2.x) / (v1.y - v2.y);

    float currentX1 = v2.x;
    float currentX2 = v2.x;


    int numPixels;
    float ratio1;
    float ratio2;
    CanvasPoint textureFrom;
    CanvasPoint textureTo;

    // iterate through each horizontal row of the canvas triangle
    for (int y = v2.y; y >= v1.y; y--) {
        // for first itteration draw 1 pixel
        if (currentX1 == currentX2) {
            uint32_t col = getColourFromTexture(textureTriangle.v0(), textureMap);
            drawCanvasPointUint(CanvasPoint(currentX1, y), col, window);

            currentX1 -= slope1;
            currentX2 -= slope2;
            continue;
        }

        CanvasPoint from = CanvasPoint(currentX1, y);
        CanvasPoint to = CanvasPoint(currentX2, y);

        // number of pixels between currentX1 and currentX2
        numPixels = std::abs(std::round(currentX1) - std::round(currentX2)) + 1;

        // ratio of currentX1 and currentx2 along v0v1 and v0v2 respectively
        ratio1 = splitPercent(v2, v0, from);
        ratio2 = splitPercent(v2, v1, to);

        //calculate the points equivalent to from and to in the texture

        textureFrom = getPointAlongLine(textureTriangle.v2(), textureTriangle.v0(), ratio1);
        textureTo = getPointAlongLine(textureTriangle.v2(), textureTriangle.v1(), ratio2);

        // vector of numPixel points along a line on the texture
        std::vector<CanvasPoint> texturePixels = lerpCanvasPoints(textureFrom, textureTo, numPixels);
        std::vector<uint32_t> texturePixelsColours = {};

        // find colour of each point on the texture
        for (int i = 0; i < numPixels; i++) {
            uint32_t colour = getColourFromTexture(texturePixels.at(i), textureMap);

            texturePixelsColours.push_back(colour);
        }

        // draw each pixel in the row with the colour from texture pixels colours
        int j = 0;
        for (int x = std::round(currentX1); x <= std::round(currentX2); x++) {
            drawCanvasPointUint(CanvasPoint(x, y), texturePixelsColours[j], window);
            j += 1;
        }


        // update X1 and X2
        currentX1 -= slope1;
        currentX2 -= slope2;
    }
}

void textureTriangle(DrawingWindow &window, CanvasTriangle drawingTriangle, TextureMap textureMap,
                     CanvasTriangle textureTriangle) {
    // makes sure vertices are sorted the same way for both triangles
    std::vector<std::vector<CanvasPoint> > verts = {
        {drawingTriangle.v0(), textureTriangle.v0()},
        {drawingTriangle.v1(), textureTriangle.v1()},
        {drawingTriangle.v2(), textureTriangle.v2()}
    };

    // sort canvas vertices by drawingTriangle y coordinate
    std::sort(verts.begin(), verts.end(), [](const std::vector<CanvasPoint> &a, const std::vector<CanvasPoint> &b) {
        return a.at(0).y < b.at(0).y; // Sort by y
    });


    CanvasPoint vD0 = verts.at(0).at(0);
    CanvasPoint vD1 = verts.at(1).at(0);
    CanvasPoint vD2 = verts.at(2).at(0);

    CanvasPoint vT0 = verts.at(0).at(1);
    CanvasPoint vT1 = verts.at(1).at(1);
    CanvasPoint vT2 = verts.at(2).at(1);


    // flat bottom triangle
    if (vD2.y == vD1.y) {
        CanvasTriangle tri = CanvasTriangle(vD0, vD1, vD2);
        textureFlatBottomTriangle(window, tri, textureMap, textureTriangle);
    }
    // flat top triangle
    else if (vD0.y == vD1.y) {
        CanvasTriangle tri = CanvasTriangle(vD0, vD1, vD2);
        textureFlatBottomTriangle(window, tri, textureMap, textureTriangle);
    }
    // all other triangles
    else {
        // calculate the point which splits drawingTriangle into a flat top and bottomed triangle
        float m = (vD2.y - vD0.y) / (vD2.x - vD0.x);
        float x = (vD1.y - vD0.y) / m + vD0.x;
        CanvasPoint vD3 = CanvasPoint(x, vD1.y);

        // calculate equivalent point to split texture at
        float sp = splitPercent(vD0, vD2, vD3);
        CanvasPoint vT3 = getPointAlongLine(vT0, vT2, sp);

        //draw flat bottom triangle
        CanvasTriangle flatBottomTri = CanvasTriangle(vD0, vD1, vD3);
        CanvasTriangle topTextTri = CanvasTriangle(vT0, vT1, vT3);
        textureFlatBottomTriangle(window, flatBottomTri, textureMap, topTextTri);

        // draw flat top triangle
        CanvasTriangle flatTopTri = CanvasTriangle(vD1, vD3, vD2);
        CanvasTriangle bottomTextTri = CanvasTriangle(vT1, vT3, vT2);
        textureFlatTopTriangle(window, flatTopTri, textureMap, bottomTextTri);
    }
}

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 point, float focalLength, DrawingWindow &window) {
    // in camera coord system
    CanvasPoint out;

    // apply camera orientation
    glm::vec3 adjCoord = point * world.camera.orientation;
    point = adjCoord;

    // adjust for canvas size otherwise (0,0,0) will be in the top left corner
    float heightShift = window.height / 2;
    float widthShift = window.width / 2;

    // if behind camera
    if (point.z > -0.00001f) {
        return {INFINITY, INFINITY};
    }

    // calculate projected coordinates
    out.x = -(point.x * (focalLength / point.z)) * world.screenScale + widthShift;
    out.y = (point.y * (focalLength / point.z)) * world.screenScale + heightShift;

    // calculate depth
    if (point.z < 0.0001f && point.z > -0.0001f) {
        out.depth = 0;
    } else {
        out.depth = 1 / point.z;
    }


    return out;
}


RayTriangleIntersection getClosestValidIntersection(glm::vec3 rayDirection, std::vector<ModelTriangle> modelTriangles,
                                                    glm::vec3 rayOrigin) {
    RayTriangleIntersection closestIntersection = RayTriangleIntersection();
    closestIntersection.distanceFromCamera = INFINITY;
    float epsilion = 1e-6f;

    rayDirection = normalize(rayDirection);

    // iterate through model triangles and look for an intersection

    for (int i = 0; i < modelTriangles.size(); i++) {
        ModelTriangle triangle = modelTriangles[i];

        // Calculate edge vectors of the triangle
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];

        // Compute the vector from the ray's origin to the first vertex of the triangle
        glm::vec3 SPVector = rayOrigin - triangle.vertices[0];

        // calculate matrix used to solve for intersection points
        glm::mat3 DEMatrix(-rayDirection, e0, e1);

        // Check if matrix is invertible
        float det = determinant(DEMatrix);
        if (std::abs(det) < epsilion) continue;


        // calculate possible solution
        glm::vec3 possibleSolution = inverse(DEMatrix) * SPVector;


        float distance = possibleSolution.x; // distance
        float u = possibleSolution.y; // barycentric coordinate
        float v = possibleSolution.z; // barycentric coordinate

        // check for valid intersection
        if (distance >= 0.0f && // intersection is in front of ray origin
            u >= -epsilion && u <= 1.0f + epsilion && // epsilon for floating point error
            v >= -epsilion && v <= 1.0f + epsilion && // epsilon for floating point error
            (u + v) <= 1.0f + epsilion) {
            if (distance < closestIntersection.distanceFromCamera) {
                // new closest solution found
                closestIntersection = RayTriangleIntersection(
                    rayOrigin + distance * rayDirection,
                    distance,
                    triangle,
                    i,
                    glm::vec2(u, v)
                );
            }
        }
    }


    return closestIntersection;
}

float getProximityBrightness(float distanceFromLight, float intensity) {
    float epsilon = 0.0000001;
    if (distanceFromLight < epsilon) {
        return 1.0f;
    }

    float proximityBrightness = intensity / (4 * M_PI * distanceFromLight * distanceFromLight);
    // cap at 0-1
    proximityBrightness = std::min(proximityBrightness, 1.0f);
    proximityBrightness = std::max(0.0f, proximityBrightness);
    return proximityBrightness;
}

float getSpecularBrightness(glm::vec3 pointToCameraDir, glm::vec3 pointToLightDir, glm::vec3 triangleNormal,
                            float specularExponent, float specularIntensity) {
    glm::vec3 reflectionDir = reflect(-pointToLightDir, triangleNormal);

    float specularBrightness = std::min(1.0f, dot(pointToCameraDir, reflectionDir));

    // add "shine"
    specularBrightness = std::pow(specularBrightness, specularExponent);

    return specularIntensity * std::max(specularBrightness, 0.0f);
}

void draw3DPoint(glm::vec3 point, Colour colour, DrawingWindow &window) {
    // get position of point in camera space and project it onto the canvas
    glm::vec3 pointCameraSpace = changeCoordSystem(worldOrign, world.camera.position, point);
    CanvasPoint canp1 = projectVertexOntoCanvasPoint(pointCameraSpace, world.camera.focalLength, window);

    // draw 5 by 5 square around projected point
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            CanvasPoint canp = CanvasPoint(canp1.x + i, canp1.y + j);
            drawCanvasPoint(canp, colour, window);
        }
    }
}

void drawRasterizedModel(Model model, DrawingWindow &window) {
    // initialise depth buffer
    for (auto &row: zDepth) {
        for (auto &elem: row) {
            elem = -INFINITY;
        }
    }

    // move triangle coordinates from the model coordinate system to the camera coordinate system
    model.changeCoordSystems(world.camera.position);

    // loop through each triangle in the model
    for (auto &tri: model.triangles) {
        Colour colour = tri.colour;

        // project the vertices onto the canvas
        CanvasPoint canp1 = projectVertexOntoCanvasPoint(tri.vertices[0], world.camera.focalLength, window);
        CanvasPoint canp2 = projectVertexOntoCanvasPoint(tri.vertices[1], world.camera.focalLength, window);
        CanvasPoint canp3 = projectVertexOntoCanvasPoint(tri.vertices[2], world.camera.focalLength, window);

        // dont render behind camera
        if (canp1.x == INFINITY || canp2.x == INFINITY || canp3.x == INFINITY || std::isnan(canp1.x) ||
            std::isnan(canp2.x) || std::isnan(canp3.x)) { continue; }

        // draw the triangles created by projected vertices
        if (world.renderMode == WireFrame) {
            drawTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
        } else {
            fillTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
            // window.renderFrame();
        }
    }
}

void drawRayTracedModel(Model model,
                        DrawingWindow &window) {
    // use rasterizer to get depth buffer for optimisation
    drawRasterizedModel(model, window);

    // change model coordinate system and record triangle normals for aoi brightness
    model.changeCoordSystems(world.camera.position);

    float z = -world.camera.focalLength;
    CanvasPoint canvasCentre((window.width / 2), (window.height / 2));

    // loop through each pixel on the screen
    for (int i = 0; i < window.width; i++) {
        for (int j = 0; j < window.height; j++) {
            // if nothing is drawn on the pixel move to next pixel
            if (zDepth[i][j] == -INFINITY) {
                continue;
            }

            CanvasPoint pixelLocation = CanvasPoint(i, j);

            // calculate pixel coordinates in camera space
            float x = (pixelLocation.x - canvasCentre.x) / world.screenScale;
            float y = -(pixelLocation.y - canvasCentre.y) / world.screenScale;

            glm::vec3 canvasPointLocationCameraSpace = glm::vec3(x, y, z);

            // transform to worldSpace by undoing camera orientation and changing coordinate system
            glm::vec3 canvasPointWorldSpace = changeCoordSystem(world.camera.position, worldOrign,
                                                                canvasPointLocationCameraSpace * transpose(
                                                                    world.camera.orientation));


            // direction of ray from camera to pixel
            glm::vec3 rayDir = normalize(canvasPointWorldSpace - world.camera.position);

            // find the closest intersection with a model triangle
            RayTriangleIntersection intersection = getClosestValidIntersection(rayDir, model.triangles, worldOrign);

            // mirror for blue box in cornell box
            // if (intersection.triangleIndex == 31 or intersection.triangleIndex == 26) {
            //     glm::vec3 reflectionDir;
            //     reflectionDir = rayDir - 2 * dot(rayDir, intersection.intersectedTriangle.normal) * intersection.
            //                     intersectedTriangle.normal;
            //     reflectionDir = normalize(reflectionDir);
            //
            //     glm::vec3 mirrorOrigin = intersection.intersectionPoint - rayDir * 0.001f;
            //
            //     RayTriangleIntersection mirrorReflection = getClosestValidIntersection(
            //         reflectionDir, model.triangles, mirrorOrigin);
            //     intersection = mirrorReflection;
            // }


            // if there is a valid intersection
            if (intersection.distanceFromCamera < INFINITY) {
                // change light source position from world coordinates to camera coordinates
                glm::vec3 light =
                        changeCoordSystem(glm::vec3(0, 0, 0), world.camera.position, world.lights[0].position);

                float epsilon = 0.001;

                // direction vector going from the point to the light
                glm::vec3 pointToLightRayDir = normalize(light - intersection.intersectionPoint);

                // move the origin of the shadow ray slightly towards the light to avoid precision errors
                glm::vec3 shadowRayOrigin = intersection.intersectionPoint + pointToLightRayDir * epsilon;

                float distanceFromLight = length(shadowRayOrigin - light);


                // send a ray from the current point to the light source to see if the point can see the light
                RayTriangleIntersection shadowRayIntersection =
                        getClosestValidIntersection(pointToLightRayDir, model.triangles, shadowRayOrigin);


                // soft shadows
                int lightsNotSeen = 0;
                std::vector<glm::vec3> lightPoints = world.lights[0].getLightPoints(-pointToLightRayDir);
                bool inCastShadow = false;


                for (auto &lightPoint: lightPoints) {
                    // light in camera space
                    glm::vec3 lightPointCameraSpace = changeCoordSystem(glm::vec3(0, 0, 0), world.camera.position,
                                                                        lightPoint);
                    // direction vector going from the point to the light
                    glm::vec3 pointToLightRaySoftDir =
                            normalize(lightPointCameraSpace - intersection.intersectionPoint);

                    // send a ray from the current point to the light source to see if the point can see the light
                    RayTriangleIntersection shadowRayIntersection =
                            getClosestValidIntersection(pointToLightRaySoftDir, model.triangles, shadowRayOrigin);

                    bool lightPointNotSeen = !(shadowRayIntersection.triangleIndex == intersection.triangleIndex or
                                               shadowRayIntersection.
                                               distanceFromCamera == INFINITY) &&
                                             shadowRayIntersection.
                                             distanceFromCamera <= distanceFromLight;

                    if (lightPointNotSeen) {
                        inCastShadow = true;
                        lightsNotSeen++;
                    }
                }


                // BRIGHTNESS CALCULATIONS
                // get normal to use in brightness calculations
                glm::vec3 normal;


                if (world.phongShading) {
                    ModelTriangle triangle = intersection.intersectedTriangle;

                    // get vertex normals of triangle
                    glm::vec3 vert1Normal = model.vertexNormalMap[triangle.verticeIndecies[0]];
                    glm::vec3 vert2Normal = model.vertexNormalMap[triangle.verticeIndecies[1]];
                    glm::vec3 vert3Normal = model.vertexNormalMap[triangle.verticeIndecies[2]];

                    // barycentric coordinates of point on triangle
                    float u = intersection.baryCoords[0];
                    float v = intersection.baryCoords[1];

                    // calculate point normal
                    glm::vec3 a = u * (vert2Normal - vert1Normal);
                    glm::vec3 b = v * (vert3Normal - vert1Normal);
                    normal = vert1Normal + a + b;
                    normal = normalize(normal);
                } else {
                    normal = intersection.intersectedTriangle.normal;
                }

                // calculate proximity brightness
                float proximityBrightness = getProximityBrightness(distanceFromLight, world.lights[0].intensity);

                // angle of incidence brightness
                float aoiBrightness;
                aoiBrightness = dot(pointToLightRayDir, normal);
                // cap between 0-1
                aoiBrightness = std::min(aoiBrightness, 1.0f);
                aoiBrightness = std::max(0.0f, aoiBrightness);

                //specular brightness
                float specularBrightness = getSpecularBrightness(-rayDir, pointToLightRayDir,
                                                                 normal, model.specularExponent,
                                                                 model.specularIntensity);


                // combine proximity brightness and aoi brightness
                float brightness = proximityBrightness * aoiBrightness;

                specularBrightness = specularBrightness * brightness;

                // apply ambient light threshold
                brightness = std::min(1.0f, brightness);
                brightness = std::max(world.ambientLightThresh, brightness);

                // if darkest brightness (ambient light threshold) then no specular highlight
                if (brightness == world.ambientLightThresh) {
                    specularBrightness = 0.0f;
                }

                // calculater the brightness of the cast shadows
                float castShadowBrightness = world.shadowDarkness * aoiBrightness;
                castShadowBrightness = std::min(1.0f, castShadowBrightness);
                castShadowBrightness = std::max(world.ambientLightThresh, castShadowBrightness);

                // get the colour of the model triangle
                Colour pixelColour;

                if (model.hasTexture) {
                    // get testure point using barycentric coordinates of intersection
                    float u = intersection.baryCoords[0];
                    float v = intersection.baryCoords[1];

                    std::array<TexturePoint, 3> texturePoints = intersection.intersectedTriangle.texturePoints;

                    glm::vec2 texturep1 = {texturePoints.at(0).x, texturePoints.at(0).y};
                    glm::vec2 texturep2 = {texturePoints.at(1).x, texturePoints.at(1).y};
                    glm::vec2 texturep3 = {texturePoints.at(2).x, texturePoints.at(2).y};

                    CanvasPoint pointOnTextureMap;

                    // calculate point on texture
                    glm::vec2 a = u * (texturep2 - texturep1);
                    glm::vec2 b = v * (texturep3 - texturep1);
                    glm::vec2 texturePoint = texturep1 + a + b;

                    // convert from vec2 to canvas point
                    pointOnTextureMap.x = texturePoint.x;
                    pointOnTextureMap.y = texturePoint.y;
                    uint32_t colour = getColourFromTexture(pointOnTextureMap, model.textureMap);
                    pixelColour = intToColour(colour);
                } else {
                    pixelColour = intersection.intersectedTriangle.colour;
                }

                // apply brightness
                if (inCastShadow) {
                    // calculate soft shadow, the darker, the more lights the point has seen the smaller shadow softness
                    // if the point has not seen any lights softness will be 0
                    float shadowSoftness =
                            1.0f - static_cast<float>(lightsNotSeen) / static_cast<float>(lightPoints.size());

                    float brighDiff = brightness - castShadowBrightness;
                    brightness = castShadowBrightness + shadowSoftness * brighDiff;
                }

                pixelColour.red = std::min(pixelColour.red * brightness + 255 * specularBrightness,
                                           255.0f);
                pixelColour.green = std::min(
                    pixelColour.green * brightness + 255 * specularBrightness,
                    255.0f);
                pixelColour.blue = std::min(pixelColour.blue * brightness + 255 * specularBrightness,
                                            255.0f);


                drawCanvasPoint(pixelLocation, pixelColour, window);
            } else {
                // corrects outline differnce between rasterizer and ray tracer
                drawCanvasPoint(pixelLocation, Colour(0, 0, 0), window);
            }
        }
    }
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
    // handle key presses
    if (event.type == SDL_KEYDOWN) {
        // camera and light source movement
        if (event.key.keysym.sym == SDLK_a) {
            // left
            glm::vec3 right = glm::vec3(world.camera.orientation[0][0], world.camera.orientation[1][0],
                                        world.camera.orientation[2][0]);
            if (moveLight) {
                world.lights[0].moveLeft();
            } else {
                world.camera.moveLeft();
            }
        } else if (event.key.keysym.sym == SDLK_d) {
            // right
            if (moveLight) {
                world.lights[0].moveRight();
            } else {
                world.camera.moveRight();
            }
        } else if (event.key.keysym.sym == SDLK_SPACE) {
            // up
            if (moveLight) {
                world.lights[0].moveUp();
            } else {
                world.camera.moveUp();
            }
        } else if (event.key.keysym.sym == SDLK_LSHIFT) {
            // down
            if (moveLight) {
                world.lights[0].moveDown();
            } else {
                world.camera.moveDown();
            }
        } else if (event.key.keysym.sym == SDLK_w) {
            // forward
            if (moveLight) {
                world.lights[0].moveForward();
            } else {
                world.camera.moveForward();
            }
        } else if (event.key.keysym.sym == SDLK_s) {
            // backward
            if (moveLight) {
                world.lights[0].moveBackward();
            } else {
                world.camera.moveBackward();
            }
        } else if (event.key.keysym.sym == SDLK_i) {
            // rotate on x axis
            world.camera.rotateX(0.05f);
        } else if (event.key.keysym.sym == SDLK_k) {
            // rotate on x axis
            world.camera.rotateX(-0.05f);
        } else if (event.key.keysym.sym == SDLK_l) {
            // rotate on y axis
            world.camera.rotateY(-0.05f);
        } else if (event.key.keysym.sym == SDLK_j) {
            // rotate on y axis
            world.camera.rotateY(0.05f);
        } else if (event.key.keysym.sym == SDLK_h) {
            // spin around model
            lookAt = true;
            if (theta + 0.05 > M_PI * 2) {
                theta = 0.0;
            } else theta += 0.05;

            float r = 4.0;
            glm::vec3 centre(0, 0, 0);
            centre = glm::vec3(0.0f, 1.0f, 0.0f);
            world.camera.position.x = centre.x + std::cos(theta) * r;
            world.camera.position.z = centre.z + std::sin(theta) * r;
        } else if (event.key.keysym.sym == SDLK_v) {
            // toggle lookAt
            lookAt = !lookAt;
            // camera.lookAt(glm::vec3(0, 1, 1));
        } else if (event.key.keysym.sym == SDLK_m) {
            // toggle move light
            moveLight = !moveLight;
            if (moveLight) {
                std::cout << "You are now moving the light" << std::endl;
            } else {
                std::cout << "You are now moving the camera" << std::endl;
            }
        } else if (event.key.keysym.sym == SDLK_p) {
            // print world info
            world.printWorldVars();
        } else if (event.key.keysym.sym == SDLK_UP) {
            // increase shinyness
            world.models[world.currentModelIndex].specularExponent *= 2;
            std::cout << "specular exponent: " << world.models[0].specularExponent << std::endl;
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            // decrease shininess
            if (world.models[world.currentModelIndex].specularExponent > 1) {
                world.models[world.currentModelIndex].specularExponent /= 2;
            }
            std::cout << "specular exponent: " << world.models[0].specularExponent << std::endl;
        } else if (event.key.keysym.sym == SDLK_LEFT) {
            // toggle phong shading
            world.phongShading = !world.phongShading;
            if (world.phongShading) {
                std::cout << "Phong shading on" << std::endl;
            } else {
                std::cout << "Phong shading off" << std::endl;
            }
        } else if (event.key.keysym.sym == SDLK_RIGHT) {
            // toggle soft shading
            if (world.lights[0].size == 0.0f) {
                world.lights[0].size = 0.1;
                std::cout << "soft shading on" << std::endl;
            } else {
                world.lights[0].size = 0.0f;
                std::cout << "soft shading off" << std::endl;
            }
        } else if (event.key.keysym.sym == SDLK_r) {
            record = !record;
            if (record) {
                std::cout << "recording" << std::endl;
            } else {
                std::cout << "recording off" << std::endl;
            }
        }


        // Change Render Type
        else if (event.key.keysym.sym == SDLK_6) {
            world.renderMode = WireFrame;
            std::cout << "switched to wireframe mode" << std::endl;
        } else if (event.key.keysym.sym == SDLK_7) {
            world.renderMode = Rasterized;
            std::cout << "switched to rasterized mode" << std::endl;
        } else if (event.key.keysym.sym == SDLK_8) {
            world.renderMode = RayTraced;
            // record = true;
            std::cout << "switched to raytraced mode" << std::endl;
        }

        // change model drawn
        else if (event.key.keysym.sym == SDLK_1) {
            world.currentModelIndex = 0;
        } else if (event.key.keysym.sym == SDLK_2) {
            world.currentModelIndex = 1;
        } else if (event.key.keysym.sym == SDLK_3) {
            world.currentModelIndex = 2;
        }
    }
    // handle mouse clicks
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    // close window
    else if (event.type == SDL_WINDOWEVENT) {
        if (event.window.event == SDL_WINDOWEVENT_CLOSE) {
            windowOpen = false;
        }
    }
}

void draw(DrawingWindow &window, Model model) {
    // render model

    if (world.renderMode == WireFrame or world.renderMode == Rasterized) {
        drawRasterizedModel(model, window);
    } else if (world.renderMode == RayTraced) {
        drawRayTracedModel(model, window);
    }
}


int main(int argc, char *argv[]) {
    srand(time(nullptr));
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;


    // Create models by parsing files

    Model boxModel = parseObj("cornell-box.obj", 0.35f, "cornell-box.mtl");
    world.models.push_back(boxModel);
    //
    //
    // Model sphereModel = parseObj("sphere.obj", 1.0f);
    // sphereModel.origin = glm::vec3(0.1f, -1.3f, -1.0f);
    // world.models.push_back(sphereModel);
    //
    //
    // Model logoModel = parseObj("logo.obj", 1.0f / 450.0f);
    // logoModel.origin = glm::vec3(-0.7f, -0.5f, 1.7f);
    // logoModel.center = glm::vec3(0.0f, 0.1f, 1.7f);
    // logoModel.specularIntensity = 0.1f;
    // world.models.push_back(logoModel);


    // world.lights[0].size = 0.3f;
    world.currentModelIndex = 0;
    // world.shadowDarkness = 0.1f;


    // world.lights[0].speed = 0.02f;
    // world.lights[0].intensity = 250.0f;
    // world.camera.position = glm::vec3(0.825764, -0.312230, 0.350722);

    // world.camera.speed = 0.05f;
    world.ambientLightThresh = 0.1f;
    world.renderMode = Rasterized;

    while (windowOpen) {
        Uint32 frameStart = SDL_GetTicks();
        //poll for events - otherwise the window will freeze
        if (window.pollForInputEvents(event)) handleEvent(event, window);

        window.clearPixels();

        if (lookAt) {
            world.camera.lookAt(world.getCurrentModel().center);
        }

        //draw all models in world
        // for (auto model: world.models) {
        //     draw(window, model);
        // }

        draw(window, world.getCurrentModel());

        // draw light source
        // draw3DPoint(world.lights[0].position, Colour(255, 255, 255), window);
        // draw3DPoint(world.getCurrentModel().center, Colour(255, 0, 0), window);

        // render the frame so it gets shown on screen
        window.renderFrame();

        // lookAt = true;
        // if (theta + 0.05 > M_PI * 2) {
        //     theta = 0.0;
        // } else theta += 0.05;
        //
        // float r = 4.0;
        // glm::vec3 centre(0, 0, 0);
        // centre = glm::vec3(0.0f, 1.0f, 0.0f);
        // world.camera.position.x = centre.x + std::cos(theta) * r;
        // world.camera.position.z = centre.z + std::sin(theta) * r;

        if (record) {
            // if (frameNumber < 33) {
            //     world.camera.moveForward();
            // } else if (frameNumber < 60) {
            //     world.camera.moveRight();
            // } else if (frameNumber < 100) {
            //     world.camera.moveLeft();
            // } else if (frameNumber < 120) {
            //     world.camera.moveRight();
            //     world.camera.moveUp();
            // } else if (frameNumber < 140) {
            //     // world.camera.moveRight();
            //     world.camera.moveDown();
            // } else {
            //     windowOpen = false;
            // }

            lookAt = true;
            if (frameNumber < 100) {
                world.camera.moveRight();
                world.lights[0].moveBackward();
            } else if (frameNumber < 200) {
                world.lights[0].moveBackward();

                if (frameNumber % 3 == 0) {
                    world.camera.moveRight();
                }
            } else {
                windowOpen = false;
            }

            window.savePPM("HackspaceLogoAnimation/output_" + std::to_string(frameNumber) + ".ppm");
            frameNumber++;
        }


        // Cap the frame rate
        Uint32 frameTime = SDL_GetTicks() - frameStart; // Calculate frame duration
        if (frameTime < FRAME_RATE) {
            SDL_Delay(FRAME_RATE - frameTime); // Delay to maintain 60 FPS
        }
    }

    return 0;
}
