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
#include <tuple>

#define WIDTH 900
#define HEIGHT 900

// GLOBAL VARS

// window
bool windowOpen = true;

// world
glm::vec3 worldOrign(0, 0, 0);


// camera
Camera camera = Camera();
bool lookAt = false;

// for rotation
float theta = 0;

// lighting

LightSource lightSource = LightSource();

float ambientLightThresh = 0.05f;
float shadowDarkness = 0.4f;
float specularExponent = 16.0f;
bool moveLight = false;
bool phongShading = true;

// depth buffer
std::vector<std::vector<float> > zDepth(WIDTH, std::vector<float>(HEIGHT, -10000));

// Render Mode
enum Mode {
    WireFrame,
    Rasterized,
    RayTraced
};

Mode renderMode = WireFrame;


void drawCanvasPoint(CanvasPoint point, Colour col, DrawingWindow &window) {
    if (pointInCanvas(point, window)) {
        // round coordinates to pixel values
        int x = std::round(point.x);
        int y = std::round(point.y);

        uint32_t colour = colourToInt(col);

        // check depth buffer
        if (point.depth == 0) {
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
        if (point.depth == 0) {
            window.setPixelColour(x, y, col);
            zDepth[x][y] = 0;
        }
        // check if point is closer than current point on that pixel
        else if ((zDepth[x][y] < 1 / point.depth) && !(zDepth[x][y] == 0)) {
            window.setPixelColour(x, y, col);
            zDepth[x][y] = 1 / point.depth;
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

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 point, float focalLength, float scale, DrawingWindow &window) {
    // in camera coord system
    CanvasPoint out;

    // apply camera orientation
    glm::vec3 adjCoord = point * camera.orientation;
    point = adjCoord;

    // adjust for canvas size otherwise (0,0,0) will be in the top left corner
    float heightShift = window.height / 2;
    float widthShift = window.width / 2;

    // calculate projected coordinates
    out.x = -(point.x * (focalLength / point.z)) * scale + widthShift;
    out.y = (point.y * (focalLength / point.z)) * scale + heightShift;

    // calcualte depth
    if (point.z == 0) {
        out.depth = 0;
    } else {
        out.depth = 1 / point.z;
    }


    return out;
}


RayTriangleIntersection getClosestValidIntersection(glm::vec3 rayDirection, std::vector<ModelTriangle> modelTriangles,
                                                    glm::vec3 rayOrigin) {
    glm::vec3 closestSolution;
    ModelTriangle closestTriangle;
    size_t closestTriangleIndex = 0;
    bool validIntersection = false;
    int validIntersections = 0;
    float epsilion = 1e-6f;

    rayDirection = normalize(rayDirection);

    // iterate through model triangles and look for an intersection
    for (int i = 0; i < modelTriangles.size(); i++) {
        validIntersection = false;

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


        float t = possibleSolution.x; // distance
        float u = possibleSolution.y; // barycentric coordinate
        float v = possibleSolution.z; // barycentric coordinate

        // check for valid intersection
        if (t >= 0.0f && // intersection is in front of ray origin
            u >= -epsilion && u <= 1.0f + epsilion && // epsilon for floating point comparison
            v >= -epsilion && v <= 1.0f + epsilion && // epsilon for floating point comparison
            (u + v) <= 1.0f + epsilion) {
            validIntersection = true;
            validIntersections++;
        }


        if (validIntersection) {
            if (validIntersections == 1) {
                // for first valid intersection
                closestSolution = possibleSolution;
                closestTriangle = triangle;
                closestTriangleIndex = i;
            } else {
                // if distance of possible solution is shorter than current closest solution
                if (t < closestSolution.x) {
                    // new closest solution found
                    closestSolution = possibleSolution;
                    closestTriangle = triangle;
                    closestTriangleIndex = i;
                }
            }
        }
    }

    // no valid intersections found - return will have distanceFromCamera = negative infinity
    if (validIntersections == 0) {
        RayTriangleIntersection invalid = RayTriangleIntersection();
        invalid.distanceFromCamera = -std::numeric_limits<float>::infinity();

        return invalid;
    }

    // get the distance and barycentric coords of the closest solution
    float t = closestSolution.x; // distance
    float u = closestSolution.y; // barycentric coordinate
    float v = closestSolution.z; // barycentric coordinate

    // calculate x,y,z coordinates of closest point using barycentric coordinates or triangle
    glm::vec3 a = u * (closestTriangle.vertices[1] - closestTriangle.vertices[0]);
    glm::vec3 b = v * (closestTriangle.vertices[2] - closestTriangle.vertices[0]);
    glm::vec3 closestPoint = closestTriangle.vertices[0] + a + b;
    glm::vec2 baryCoords = {u, v};


    return RayTriangleIntersection(
        closestPoint, t, closestTriangle, closestTriangleIndex, baryCoords);
}

float getProximityBrightness(float distanceFromLight, float intensity) {
    float epsilon = 0.0000001;
    if (distanceFromLight < epsilon) {
        return 1.0f;
    } else {
        float proximityBrightness = intensity / (4 * M_PI * distanceFromLight * distanceFromLight);
        // cap at 0-1
        proximityBrightness = std::min(proximityBrightness, 1.0f);
        proximityBrightness = std::max(0.0f, proximityBrightness);
        return proximityBrightness;
    }
}

float getSpecularBrightness(glm::vec3 pointToCameraDir, glm::vec3 pointToLightDir, glm::vec3 triangleNormal) {
    glm::vec3 reflectionDir = reflect(pointToLightDir, triangleNormal);

    float specularBrightness = std::min(1.0f, dot(pointToCameraDir, reflectionDir));
    specularBrightness = std::pow(specularBrightness, specularExponent);

    return std::max(specularBrightness, 0.0f);
}

void draw3DPoint(glm::vec3 point, Colour colour, float screenScale, DrawingWindow &window) {
    glm::vec3 pointCameraSpace = changeCoordSystem(worldOrign, camera.position, point);
    CanvasPoint canp1 = projectVertexOntoCanvasPoint(pointCameraSpace, camera.focalLength, screenScale, window);

    // draw 3 by 3 square around projected point
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            CanvasPoint canp = CanvasPoint(canp1.x + i, canp1.y + j);
            drawCanvasPoint(canp, colour, window);
        }
    }
}

void drawWireframeModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale, DrawingWindow &window) {
    // initialise depth buffer
    for (auto &row: zDepth) {
        for (auto &elem: row) {
            elem = -std::numeric_limits<float>::infinity();
        }
    }

    // for each triangle in the model
    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        Colour colour = tri.colour;

        // move triangle coordinates from the model coordinate system to the camera coordinate system
        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[1]);
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[2]);

        // project the verticies onto the canvas
        CanvasPoint canp1 = projectVertexOntoCanvasPoint(vert1, camera.focalLength, scale, window);
        CanvasPoint canp2 = projectVertexOntoCanvasPoint(vert2, camera.focalLength, scale, window);
        CanvasPoint canp3 = projectVertexOntoCanvasPoint(vert3, camera.focalLength, scale, window);

        // draw the triangle created by projected verticies
        drawTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
    }
}

void drawRasterizedModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale, DrawingWindow &window) {
    // initialise depth buffer
    for (auto &row: zDepth) {
        for (auto &elem: row) {
            elem = -std::numeric_limits<float>::infinity();
        }
    }
    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        Colour colour = tri.colour;
        // std::cout << glm::to_string(tri.vertices[0]) << glm::to_string(tri.vertices[1]) << glm::to_string(tri.vertices[2]) << std::endl;
        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[1]);
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera.position, tri.vertices[2]);


        CanvasPoint canp1 = projectVertexOntoCanvasPoint(vert1, camera.focalLength, scale, window);
        CanvasPoint canp2 = projectVertexOntoCanvasPoint(vert2, camera.focalLength, scale, window);
        CanvasPoint canp3 = projectVertexOntoCanvasPoint(vert3, camera.focalLength, scale, window);

        // if (pointInCanvas(canp1, window) and pointInCanvas(canp2, window) and pointInCanvas(canp3, window)) {
        fillTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
    }
}

void drawRayTracedModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale,
                        std::map<int, glm::vec3> vertexNormals, DrawingWindow &window) {
    // change model coordinate system and record triangle normals for aoi brightness

    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        // std::cout << glm::to_string(tri.vertices[0]) << glm::to_string(tri.vertices[1]) << glm::to_string(tri.vertices[2]) << std::endl;
        glm::vec3 modelVert1 = tri.vertices[0];
        glm::vec3 modelVert2 = tri.vertices[1];
        glm::vec3 modelVert3 = tri.vertices[2];

        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera.position, modelVert1);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera.position, modelVert2);;
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera.position, modelVert3);

        model[i].vertices[0] = vert1;
        model[i].vertices[1] = vert2;
        model[i].vertices[2] = vert3;

        model[i].normal = triangleNormal(vert1, vert2, vert3);
    }

    float z = -camera.focalLength;
    CanvasPoint canvasCentre((window.width / 2), (window.height / 2));

    for (int i = 0; i < window.width; i++) {
        for (int j = 0; j < window.height; j++) {
            CanvasPoint pixelLocation = CanvasPoint(i, j);

            // calculate pixel coordinates in camera space
            float x = (pixelLocation.x - canvasCentre.x) / scale;
            float y = -(pixelLocation.y - canvasCentre.y) / scale;

            glm::vec3 canvasPointLocationCameraSpace = glm::vec3(x, y, z);

            // transform to worldSpace by undoing camera orientation and changing coordinate system
            glm::vec3 canvasPointWorldSpace = changeCoordSystem(camera.position, worldOrign,
                                                                canvasPointLocationCameraSpace * transpose(
                                                                    camera.orientation));


            // direction of ray from camera to pixel
            glm::vec3 rayDir = normalize(canvasPointWorldSpace - camera.position);

            //find the closest intersection with a model triangle
            RayTriangleIntersection intersection = getClosestValidIntersection(rayDir, model, worldOrign);


            // if not invalid intersection
            if (intersection.distanceFromCamera > 0) {
                // change light source position from world coordinates to camera coordinates
                glm::vec3 light = changeCoordSystem(glm::vec3(0, 0, 0), camera.position, lightSource.position);
                // light = lightSource;

                float epsilon = 0.001;
                glm::vec3 pointToLightRayDir = normalize(light - intersection.intersectionPoint);
                glm::vec3 shadowRayOrigin = intersection.intersectionPoint + pointToLightRayDir * epsilon;
                RayTriangleIntersection shadowRayIntersection =
                        getClosestValidIntersection(pointToLightRayDir, model, shadowRayOrigin);

                // calculate proximity brightness
                float distanceFromLight = length(shadowRayOrigin - light);
                float proximityBrightness = getProximityBrightness(distanceFromLight, lightSource.intensity);


                // get normal of the triangle that the point is on
                glm::vec3 normal;

                // angle of incidence brightness
                float aoiBrightness;
                if (phongShading) {
                    ModelTriangle triangle = intersection.intersectedTriangle;

                    // get vertex normals of triangle
                    glm::vec3 vert1Normal = vertexNormals[triangle.verticeIndecies[0]];
                    glm::vec3 vert2Normal = vertexNormals[triangle.verticeIndecies[1]];
                    glm::vec3 vert3Normal = vertexNormals[triangle.verticeIndecies[2]];

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


                aoiBrightness = dot(pointToLightRayDir, normal);


                // cap between 0-1
                aoiBrightness = std::min(aoiBrightness, 1.0f);
                aoiBrightness = std::max(0.0f, aoiBrightness);

                //specular brightness
                float specularBrightness = getSpecularBrightness(-rayDir, pointToLightRayDir,
                                                                 normal);


                //combine proximity brightness and aoi brightness
                float brightness = proximityBrightness * aoiBrightness;

                // apply ambient light threshold
                brightness = std::min(1.0f, brightness);
                brightness = std::max(ambientLightThresh, brightness);

                // if darkest brightness then no specular highlight
                if (brightness == ambientLightThresh) {
                    specularBrightness = 0.0f;
                }

                float castShadowBrightness = shadowDarkness * aoiBrightness * proximityBrightness;
                castShadowBrightness = std::min(1.0f, castShadowBrightness);
                castShadowBrightness = std::max(ambientLightThresh, castShadowBrightness);

                Colour pixelColour = intersection.intersectedTriangle.colour;

                float directBrightness;
                // if in shadow
                if (!(shadowRayIntersection.triangleIndex == intersection.triangleIndex or shadowRayIntersection.
                      distanceFromCamera == -std::numeric_limits<float>::infinity()) && shadowRayIntersection.
                    distanceFromCamera <= distanceFromLight) {
                    // in shadow
                    pixelColour.red *= castShadowBrightness;
                    pixelColour.green *= castShadowBrightness;
                    pixelColour.blue *= castShadowBrightness;


                    // pixelColour = shadowRayIntersection.intersectedTriangle.colour;
                } else {
                    pixelColour.red = std::min(pixelColour.red * brightness * (1 + specularBrightness),
                                               255.0f);
                    pixelColour.green = std::min(
                        pixelColour.green * brightness * (1 + specularBrightness),
                        255.0f);
                    pixelColour.blue = std::min(pixelColour.blue * brightness * (1 + specularBrightness),
                                                255.0f);
                }


                drawCanvasPoint(pixelLocation, pixelColour, window);
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
            glm::vec3 right = glm::vec3(camera.orientation[0][0], camera.orientation[1][0], camera.orientation[2][0]);
            if (moveLight) {
                lightSource.moveLeft();
            } else {
                camera.moveLeft();
            }
        } else if (event.key.keysym.sym == SDLK_d) {
            // right
            if (moveLight) {
                lightSource.moveRight();
            } else {
                camera.moveRight();
            }
        } else if (event.key.keysym.sym == SDLK_SPACE) {
            // up
            if (moveLight) {
                lightSource.moveUp();
            } else {
                camera.moveUp();
            }
        } else if (event.key.keysym.sym == SDLK_LSHIFT) {
            // down
            if (moveLight) {
                lightSource.moveDown();
            } else {
                camera.moveDown();
            }
        } else if (event.key.keysym.sym == SDLK_w) {
            // forward
            if (moveLight) {
                lightSource.moveForward();
            } else {
                camera.moveForward();
            }
        } else if (event.key.keysym.sym == SDLK_s) {
            // backward
            if (moveLight) {
                lightSource.moveBackward();
            } else {
                camera.moveBackward();
            }
        } else if (event.key.keysym.sym == SDLK_l) {
            // rotate on x axis
            camera.rotateX(0.1f);
        } else if (event.key.keysym.sym == SDLK_k) {
            // rotate on y axis
            camera.rotateY(0.1f);
        } else if (event.key.keysym.sym == SDLK_o) {
            // rotate on z axis
            camera.rotateZ(0.1f);
        } else if (event.key.keysym.sym == SDLK_h) {
            // spin around model
            lookAt = true;
            if (theta + 0.05 > M_PI * 2) {
                theta = 0.0;
            } else theta += 0.05;

            float r = 4.0;
            glm::vec3 centre(0, 0, 0);
            centre = glm::vec3(0.0f, 1.0f, 0.0f);
            camera.position.x = centre.x + std::cos(theta) * r;
            camera.position.z = centre.z + std::sin(theta) * r;
        } else if (event.key.keysym.sym == SDLK_v) {
            // toggle lookAt
            lookAt = !lookAt;
            // camera.lookAt(glm::vec3(0, 1, 1));
        } else if (event.key.keysym.sym == SDLK_m) {
            // toggle move lightp
            moveLight = !moveLight;
        } else if (event.key.keysym.sym == SDLK_p) {
            // print world info
            std::cout << "light: " << to_string(lightSource.position) << std::endl;
            std::cout << "camera: " << to_string(camera.position) << std::endl;
            std::cout << "camera orientation: " << to_string(camera.orientation) << std::endl;
            std::cout << "specular exponent: " << specularExponent << std::endl;
        } else if (event.key.keysym.sym == SDLK_UP) {
            // increase shinyness
            specularExponent *= 2;
            std::cout << "specular exponent: " << specularExponent << std::endl;
        } else if (event.key.keysym.sym == SDLK_DOWN) {
            // decrease shinyness
            specularExponent /= 2;
            std::cout << "specular exponent: " << specularExponent << std::endl;
        } else if (event.key.keysym.sym == SDLK_LEFT) {
            // toggle phong shading
            phongShading = !phongShading;
            std::cout << "phong shading: " << phongShading << std::endl;
        }


        // Change Render Type
        else if (event.key.keysym.sym == SDLK_6) {
            renderMode = WireFrame;
        } else if (event.key.keysym.sym == SDLK_7) {
            renderMode = Rasterized;
        } else if (event.key.keysym.sym == SDLK_8) {
            renderMode = RayTraced;
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

void draw(DrawingWindow &window, std::vector<ModelTriangle> model, glm::vec3 modelOrigin,
          std::map<int, glm::vec3> vertexNormals, float scale) {
    // render model

    if (renderMode == WireFrame) {
        drawWireframeModel(model, modelOrigin, scale, window);
    } else if (renderMode == Rasterized) {
        drawRasterizedModel(model, modelOrigin, scale, window);
    } else if (renderMode == RayTraced) {
        drawRayTracedModel(model, modelOrigin, scale, vertexNormals, window);
    }
}


int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;

    // parse obj file to create model
    // auto [boxModelTriangles, boxModelVertNormals] = parseObj("cornell-box.obj", 0.35f, "cornell-box.mtl");
    auto [bunnyModelTriangles, bunnyModelVertNormals] = parseObj("stanford-bunny.obj", 10.0f);
    // auto [sphereModelTriangles, sphereModelVertNormals] = parseObj("sphere.obj", 1.0f);

    // put light in ceiling
    lightSource.position = glm::vec3(0.0f, 0.7f, 0.0f);
    // lightSource.position = glm::vec3(0.0f, 0.2f, 0.6f);

    // for bunny
    lightSource.position = glm::vec3(-0.500000, 2.000000, 1.200000);


    while (windowOpen) {
        //poll for events - otherwise the window will freeze
        if (window.pollForInputEvents(event)) handleEvent(event, window);

        window.clearPixels();

        // draw box model
        glm::vec3 modelOrigin = glm::vec3(0, 0, 0);
        float screenScale = 500;
        glm::vec3 modelCentre = glm::vec3(0.0f, 1.0f, 0.0f);

        if (lookAt) {
            camera.lookAt(modelCentre);
        }

        draw3DPoint(lightSource.position, Colour(155, 155, 0), screenScale, window);
        draw3DPoint(modelOrigin, Colour(255, 0, 0), screenScale, window);
        draw3DPoint(modelCentre, Colour(0, 255, 0), screenScale, window);

        draw(window, bunnyModelTriangles, modelOrigin, bunnyModelVertNormals, screenScale);

        // render the frame so it gets shown on screen
        window.renderFrame();
    }

    return 0;
}
