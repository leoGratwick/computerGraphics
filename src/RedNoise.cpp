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
#include <filesystem>
#include <cerrno>
#include "Helpers.h"
#include "RayTriangleIntersection.h"

#define WIDTH 900
#define HEIGHT 900

// GLOBAL VARS

bool windowOpen = true;

// camera
glm::vec3 camera(0, 0, 4);
glm::vec3 camDir(1, 1, 1);
glm::mat3 camOr(1, 0, 0, 0, 1, 0, 0, 0, 1);
float cameraSpeed = 0.1;
float focalLength = 2.0;

// lighting
// glm::vec3 lightSource(0, 0.7, 1);
glm::vec3 lightSource(0, 0.9, 0);
float ambientLightThresh = 0.2;
float shadowDarkness = 0.3;
bool moveLight = true;
float specularExponent = 8;


float theta = 0;
float mouseSensitivity = 0.1;
bool mouseMovement = false;

std::vector<std::vector<float> > zDepth(WIDTH, std::vector<float>(HEIGHT, -10000));


bool lookAt = false;

enum Mode {
    WireFrame,
    Rasterized,
    RayTraced
};


Mode renderMode = RayTraced;

void processMouseMovement(int xRot, int yRot) {
    // Convert offsets to radians and apply sensitivity
    float xAngle = glm::radians(xRot * mouseSensitivity);
    float yAngle = glm::radians(yRot * mouseSensitivity);

    // Rotate around the up axis (Y) for horizontal mouse movement
    camOr = rotateOrientation("y", xAngle, camOr);

    // Rotate around the right axis (X) for vertical mouse movement
    glm::vec3 right = camOr[0]; // Assuming first column is the right direction
    camOr = rotateOrientation("x", yAngle, camOr);

    // prevent drift by making forward, up and right are orthogonal to eachoter
    camOr = orthonormalize(camOr);
}

bool pointInCanvas(CanvasPoint point, DrawingWindow &window) {
    return ((point.x < window.width - 1) && (point.y < window.height - 1) && (point.y >= 0) && (point.x >= 0));
}

std::map<std::string, Colour> createPalette(std::string filename) {
    std::map<std::string, Colour> palette;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            if (line.substr(0, 6) == "newmtl") {
                std::vector<std::string> tokens = splitByDelimiter(line, ' ');
                std::string name = tokens.at(1);

                getline(file, line);
                std::vector<std::string> colourVals = splitByDelimiter(line, ' ');

                // convert from 0-1 to 0-225 for rgb values
                int r = std::round(std::stof(colourVals.at(1)) * 225);
                int g = std::round(std::stof(colourVals.at(2)) * 225);
                int b = std::round(std::stof(colourVals.at(3)) * 225);
                palette[name] = Colour(r, g, b);
                // std::cout << name << " added to palette as " << Colour(r,g,b)<< std::endl;
            }
        }
        file.close();
    } else {
        std::cerr << "Could not open file " << filename << std::endl;
    }
    return palette;
}

std::vector<ModelTriangle> parseObj(std::string objFilename, std::string mtlFilename, float scale) {
    std::map<std::string, Colour> palette = createPalette(mtlFilename);
    std::ifstream file(objFilename);
    std::vector<ModelTriangle> tris = {};
    std::string line;
    ModelTriangle tri;
    std::vector<glm::vec3> vertices = {};
    std::vector<std::vector<std::string> > faces = {};
    std::string col;


    if (file.is_open()) {
        while (std::getline(file, line)) {
            // std::cout << line << std::endl;

            // add verticies
            if (line[0] == 'v') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    float x = std::stof(nums.at(1)) * scale;
                    float y = std::stof(nums.at(2)) * scale;
                    float z = std::stof(nums.at(3)) * scale;
                    // std::cout << x << y << z << std::endl;
                    glm::vec3 vert(x, y, z);
                    // std::cout << line << std::endl;
                    // std::cout << vert.x << vert.y << vert.z << std::endl;
                    vertices.push_back(vert);
                } catch (...) {
                    std::cout << "couldn't convert from string to float - vert" << std::endl;
                    return tris;
                }
            }
            // add faces
            else if (line[0] == 'f') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    nums.at(1).pop_back();
                    nums.at(2).pop_back();
                    nums.at(3).pop_back();
                    std::vector<std::string> face = {};
                    face.push_back(nums.at(1));
                    face.push_back(nums.at(2));
                    face.push_back(nums.at(3));
                    face.push_back(col);
                    faces.push_back(face);
                } catch (...) {
                    std::cout << "couldn't convert from string to float - face" << std::endl;
                    return tris;
                }
            } else if (line.substr(0, 6) == "usemtl") {
                // changing current coulur
                std::vector<std::string> tokens = splitByDelimiter(line, ' ');
                std::string name = tokens.at(1);
                col = name;
            }
        }

        file.close();

        // add every face along with its colour
        for (int i = 0; i < faces.size(); i++) {
            int v1loc = std::stoi(faces.at(i)[0]) - 1;
            glm::vec3 vert1(vertices.at(v1loc));
            int v2loc = std::stoi(faces.at(i)[1]) - 1;
            glm::vec3 vert2(vertices.at(v2loc));
            int v3loc = std::stoi(faces.at(i)[2]) - 1;
            glm::vec3 vert3(vertices.at(v3loc));
            std::string colour = faces.at(i)[3];
            ModelTriangle tri = ModelTriangle(vert1, vert2, vert3, palette[colour]);
            // std::cout << tri << palette[colour] << std::endl;
            tris.push_back(tri);
        }
    } else {
        std::cerr << "Could not open file " << objFilename << std::endl;
    }
    return tris;
}

void drawCanvasPoint(CanvasPoint point, Colour col, DrawingWindow &window) {
    if (pointInCanvas(point, window)) {
        int x = std::round(point.x);
        int y = std::round(point.y);

        if (point.depth == 0) {
            uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue);
            window.setPixelColour(x, y, colour);
            zDepth[x][y] = 0;
        } else if (zDepth[x][y] < 1 / point.depth && !zDepth[x][y] == 0) {
            uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue);
            window.setPixelColour(x, y, colour);
            zDepth[x][y] = 1 / point.depth;
        }
    }
}

void drawCanvasPointUint(CanvasPoint point, uint32_t col, DrawingWindow &window) {
    if (pointInCanvas(point, window)) {
        int x = std::round(point.x);
        int y = std::round(point.y);

        if (point.depth == 0) {
            window.setPixelColour(x, y, col);
            zDepth[x][y] = 0;
        } else if (zDepth[x][y] < 1 / point.depth && !zDepth[x][y] == 0) {
            window.setPixelColour(x, y, col);
            zDepth[x][y] = 1 / point.depth;
        }
    }
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    // std::cout << "drawing line from: " << from.x << ", " << from.y << " to: " << to.x << ", " << to.y << std::endl;
    if (from.x == to.x and to.y == from.y) {
        drawCanvasPoint(from, colour, window);
        return;
    }
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    int steps = std::max(std::abs(xDiff), std::abs(yDiff)) * 2;
    float xStep = xDiff / steps;
    float yStep = yDiff / steps;
    // std::vector<float> lerpX = interpolateSingleFloats(from.x, to.x, range);
    // std::vector<float> lerpY = interpolateSingleFloats(from.y, to.y, range);
    std::vector<float> lerpDepths = interpolateSingleFloats(from.depth, to.depth, steps);

    for (int i = 0; i < steps; i++) {
        CanvasPoint point = CanvasPoint(from.x + i * xStep, from.y + i * yStep);
        point.depth = lerpDepths.at(i);
        // std::cout << "drawing point: " << point << std::endl;
        drawCanvasPoint(point, colour, window);
    }
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v0(), triangle.v2(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
}

void fillFlatBottomTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // std::cout << "drawing flat bottom triangle" << std::endl;
    // the first vertex is at the top of the triangle
    drawTriangle(window, triangle, colour);
    drawTriangle(window, triangle, colour);

    CanvasPoint v0 = triangle.v0();
    CanvasPoint v1 = triangle.v1();
    CanvasPoint v2 = triangle.v2();


    float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
    float slope2 = (v2.x - v0.x) / (v2.y - v0.y);

    std::string dir1;
    std::string dir2;


    float currentX1 = v0.x;
    float currentX2 = v0.x;

    std::vector<float> depths1 = interpolateSingleFloats(v0.depth, v1.depth, std::round(v1.y) - std::round(v0.y));
    std::vector<float> depths2 = interpolateSingleFloats(v0.depth, v2.depth, std::round(v1.y) - std::round(v0.y));

    int depthCount = 0;

    for (int y = std::round(v0.y); y < std::round(v1.y); y++) {
        CanvasPoint from = CanvasPoint(currentX1, y);
        from.depth = depths1.at(depthCount);

        CanvasPoint to = CanvasPoint(currentX2, y);
        to.depth = depths2.at(depthCount);


        drawLine(window, from, to, colour);
        currentX1 += slope1;
        currentX2 += slope2;
        depthCount += 1;
    }
}

void fillFlatTopTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // std::cout << "drawing flat top filled triangle" << std::endl;
    drawTriangle(window, triangle, colour);
    // v2 is the bottom vertex
    CanvasPoint v0 = triangle.v0();
    CanvasPoint v1 = triangle.v1();
    CanvasPoint v2 = triangle.v2();


    float slope1 = (v0.x - v2.x) / (v0.y - v2.y);
    float slope2 = (v1.x - v2.x) / (v1.y - v2.y);

    float currentX1 = v2.x;
    float currentX2 = v2.x;

    std::vector<float> depths1 = interpolateSingleFloats(v2.depth, v0.depth, std::round(v2.y) - std::round(v1.y));
    std::vector<float> depths2 = interpolateSingleFloats(v2.depth, v1.depth, std::round(v2.y) - std::round(v1.y));

    int depthCount = 0;

    for (int y = std::round(v2.y); y > std::round(v1.y); y--) {
        CanvasPoint from = CanvasPoint(currentX1, y);
        from.depth = depths1.at(depthCount);

        CanvasPoint to = CanvasPoint(currentX2, y);
        to.depth = depths2.at(depthCount);
        // std::cout << "from: "<< from << " to: " << to << std::endl;
        drawLine(window, from, to, colour);
        currentX1 -= slope1;
        currentX2 -= slope2;
        depthCount += 1;
    }
}

void fillTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    std::vector<CanvasPoint> verts = {triangle.v0(), triangle.v1(), triangle.v2()};

    drawTriangle(window, triangle, colour);

    // sort verticies by y coordinate
    std::sort(verts.begin(), verts.end(), [](const CanvasPoint &a, const CanvasPoint &b) {
        return a.y < b.y; // Sort by y
    });

    CanvasPoint v0 = verts[0];;
    // std::cout << v0 << std::endl;
    CanvasPoint v1 = verts[1];
    // std::cout << v1 << std::endl;
    CanvasPoint v2 = verts[2];
    // std::cout << v2 << std::endl;

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
        float m = (v2.y - v0.y) / (v2.x - v0.x);
        float x = (v1.y - v0.y) / m + v0.x;
        CanvasPoint v3 = CanvasPoint(x, v1.y);
        float sp = splitPercent(v0, v2, v3);
        v3.depth = v0.depth + (v2.depth - v0.depth) * sp;
        // std::cout << v3 << std::endl;
        CanvasTriangle flatBottomTri = CanvasTriangle(v0, v1, v3);
        fillFlatBottomTriangle(window, flatBottomTri, colour);

        CanvasTriangle flatTopTri = CanvasTriangle(v1, v3, v2);
        fillFlatTopTriangle(window, flatTopTri, colour);
    }
}

void randomTriangle(DrawingWindow &window) {
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
    std::cout << "drawing textured flat bottom triangle" << std::endl;

    CanvasPoint v0 = drawingTriangle.v0();
    CanvasPoint v1 = drawingTriangle.v1();
    CanvasPoint v2 = drawingTriangle.v2();

    float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
    float slope2 = (v2.x - v0.x) / (v2.y - v0.y);

    float currentX1 = v0.x;
    float currentX2 = v0.x;

    // number of pixels between currentX1 and currentX2
    int numPixels;
    float ratio1;
    float ratio2;
    CanvasPoint textureFrom;
    CanvasPoint textureTo;


    for (int y = v0.y; y <= v1.y; y++) {
        if (currentX1 == currentX2) {
            uint32_t col = getColourFromTexture(textureTriangle.v0(), textureMap);
            drawCanvasPointUint(CanvasPoint(currentX1, y), col, window);

            currentX1 += slope1;
            currentX2 += slope2;
            continue;
        }

        CanvasPoint from = CanvasPoint(currentX1, y);
        CanvasPoint to = CanvasPoint(currentX2, y);

        numPixels = std::abs(std::round(currentX1) - std::round(currentX2)) + 1;

        ratio1 = splitPercent(v0, v1, from);
        ratio2 = splitPercent(v0, v2, to);

        textureFrom = getPointAlongLine(textureTriangle.v0(), textureTriangle.v1(), ratio1);
        textureTo = getPointAlongLine(textureTriangle.v0(), textureTriangle.v2(), ratio2);

        // vector of numPixel points along a line on the texture
        std::vector<CanvasPoint> texturePixels = lerpCanvasPoints(textureFrom, textureTo, numPixels);
        std::vector<uint32_t> texturePixelsColours = {};
        // std::cout << texturePixels.at(0) << std::endl;
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


        // drawLine(window, from, to, colour);
        currentX1 += slope1;
        currentX2 += slope2;
    }
}

void textureFlatTopTriangle(DrawingWindow &window, CanvasTriangle drawingTriangle, TextureMap textureMap,
                            CanvasTriangle textureTriangle) {
    std::cout << "drawing textured flat top triangle" << std::endl;


    CanvasPoint v0 = drawingTriangle.v0();
    CanvasPoint v1 = drawingTriangle.v1();
    CanvasPoint v2 = drawingTriangle.v2();

    float slope1 = (v0.x - v2.x) / (v0.y - v2.y);
    float slope2 = (v1.x - v2.x) / (v1.y - v2.y);

    float currentX1 = v2.x;
    float currentX2 = v2.x;


    // number of pixels between currentX1 and currentX2
    int numPixels;
    float ratio1;
    float ratio2;
    CanvasPoint textureFrom;
    CanvasPoint textureTo;


    for (int y = v2.y; y >= v1.y; y--) {
        if (currentX1 == currentX2) {
            uint32_t col = getColourFromTexture(textureTriangle.v0(), textureMap);
            drawCanvasPointUint(CanvasPoint(currentX1, y), col, window);

            currentX1 -= slope1;
            currentX2 -= slope2;
            continue;
        }

        CanvasPoint from = CanvasPoint(currentX1, y);
        CanvasPoint to = CanvasPoint(currentX2, y);

        numPixels = std::abs(std::round(currentX1) - std::round(currentX2)) + 1;

        ratio1 = splitPercent(v2, v0, from);
        ratio2 = splitPercent(v2, v1, to);

        textureFrom = getPointAlongLine(textureTriangle.v2(), textureTriangle.v0(), ratio1);
        textureTo = getPointAlongLine(textureTriangle.v2(), textureTriangle.v1(), ratio2);

        // vector of numPixel points along a line on the texture
        std::vector<CanvasPoint> texturePixels = lerpCanvasPoints(textureFrom, textureTo, numPixels);
        std::vector<uint32_t> texturePixelsColours = {};
        // std::cout << texturePixels.at(0) << std::endl;
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


        // drawLine(window, from, to, colour);
        currentX1 -= slope1;
        currentX2 -= slope2;
    }
}

void textureTriangle(DrawingWindow &window, CanvasTriangle drawingTriangle, TextureMap textureMap,
                     CanvasTriangle textureTriangle) {
    // makes sure verts are sorted the same way for both triangles
    std::vector<std::vector<CanvasPoint> > verts = {
        {drawingTriangle.v0(), textureTriangle.v0()},
        {drawingTriangle.v1(), textureTriangle.v1()},
        {drawingTriangle.v2(), textureTriangle.v2()}
    };
    // std::vector<CanvasPoint> vertsTx = {textureTriangle.v0(), textureTriangle.v1(), textureTriangle.v2()};

    // sort canvas verticies by drawingTriangle y coordinate
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
        float m = (vD2.y - vD0.y) / (vD2.x - vD0.x);
        float x = (vD1.y - vD0.y) / m + vD0.x;
        CanvasPoint vD3 = CanvasPoint(x, vD1.y);
        float sp = splitPercent(vD0, vD2, vD3);
        // calculate point to split texture at
        CanvasPoint vT3 = getPointAlongLine(vT0, vT2, sp);

        CanvasTriangle flatBottomTri = CanvasTriangle(vD0, vD1, vD3);
        CanvasTriangle topTextTri = CanvasTriangle(vT0, vT1, vT3);
        textureFlatBottomTriangle(window, flatBottomTri, textureMap, topTextTri);

        CanvasTriangle flatTopTri = CanvasTriangle(vD1, vD3, vD2);
        CanvasTriangle bottomTextTri = CanvasTriangle(vT1, vT3, vT2);
        textureFlatTopTriangle(window, flatTopTri, textureMap, bottomTextTri);
    }
}

glm::vec3 changeCoordSystem(glm::vec3 fromOrigin, glm::vec3 toOrigin, glm::vec3 point) {
    glm::vec3 transformationVector = fromOrigin - toOrigin;
    glm::vec3 outPoint = point + transformationVector;
    return outPoint;
}


CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 point, float focalLength, float scale, DrawingWindow &window) {
    // in camera coord system
    CanvasPoint out;
    glm::vec3 adjCoord = point * camOr;

    point = adjCoord;
    float heightShift = window.height / 2;
    float widthShift = window.width / 2;

    if (point.z == 0) {
        out.x = 100000000000000;
        out.y = 100000000000000;
        out.depth = 0;
    } else {
        out.x = -(point.x * (focalLength / point.z)) * scale + widthShift;
        out.y = (point.y * (focalLength / point.z)) * scale + heightShift;
        out.depth = 1 / point.z;

        // std::cout<< out << point.x << point.y << point.z << std::endl;
    }


    return out;
}

void camLookAt(glm::vec3 point) {
    glm::vec3 camDir = glm::normalize(point - camera);
    glm::vec3 up(0, 1, 0);
    glm::vec3 left = -glm::normalize(glm::cross(up, camDir));
    glm::vec3 newUp = glm::normalize(glm::cross(camDir, -left));

    glm::mat3 newCamOr = glm::mat3(
        left, // First column: right vector
        newUp, // Second column: up vector
        -camDir // Third column: negative direction vector (looking direction)
    );


    camOr = newCamOr;
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

    for (int i = 0; i < modelTriangles.size(); i++) {
        validIntersection = false;

        ModelTriangle triangle = modelTriangles.at(i);
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 SPVector = rayOrigin - triangle.vertices[0];

        glm::mat3 DEMatrix(-rayDirection, e0, e1);

        // Check if matrix is invertible
        float det = glm::determinant(DEMatrix);
        if (std::abs(det) < epsilion) continue;

        glm::vec3 possibleSolution = inverse(DEMatrix) * SPVector;


        float t = possibleSolution.x; // distance
        float u = possibleSolution.y; // barycentric coordinate
        float v = possibleSolution.z; // barycentric coordinate

        // check for valid intersection with numerical stability
        if (t >= 0.0f && // Intersection is in front of ray origin
            u >= -epsilion && u <= 1.0f + epsilion && // Added epsilon for floating point comparison
            v >= -epsilion && v <= 1.0f + epsilion && // Added epsilon for floating point comparison
            (u + v) <= 1.0f + epsilion) {
            validIntersection = true;
            validIntersections++;
        }


        if (validIntersection) {
            if (validIntersections == 1) {
                closestSolution = possibleSolution;
                closestTriangle = triangle;
                closestTriangleIndex = i;
            } else {
                // if distance of possible solution is shorter than current closest solution
                if (std::abs(t) < std::abs(closestSolution.x)) {
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

    float t = closestSolution.x; // distance
    float u = closestSolution.y; // barycentric coordinate
    float v = closestSolution.z; // barycentric coordinate


    glm::vec3 a = u * (closestTriangle.vertices[1] - closestTriangle.vertices[0]);
    glm::vec3 b = v * (closestTriangle.vertices[2] - closestTriangle.vertices[0]);

    glm::vec3 closestPoint = closestTriangle.vertices[0] + a + b;
    RayTriangleIntersection closestIntersection = RayTriangleIntersection(
        closestPoint, t, closestTriangle, closestTriangleIndex);

    return closestIntersection;
}

void drawWireframeModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale, DrawingWindow &window) {
    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        Colour colour = tri.colour;
        // std::cout << glm::to_string(tri.vertices[0]) << glm::to_string(tri.vertices[1]) << glm::to_string(tri.vertices[2]) << std::endl;
        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera, tri.vertices[1]);
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera, tri.vertices[2]);


        CanvasPoint canp1 = projectVertexOntoCanvasPoint(vert1, focalLength, scale, window);
        CanvasPoint canp2 = projectVertexOntoCanvasPoint(vert2, focalLength, scale, window);
        CanvasPoint canp3 = projectVertexOntoCanvasPoint(vert3, focalLength, scale, window);

        // if (pointInCanvas(canp1, window) and pointInCanvas(canp2, window) and pointInCanvas(canp3, window)) {
        drawTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
    }
}

void drawRasterizedModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale, DrawingWindow &window) {
    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        Colour colour = tri.colour;
        // std::cout << glm::to_string(tri.vertices[0]) << glm::to_string(tri.vertices[1]) << glm::to_string(tri.vertices[2]) << std::endl;
        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera, tri.vertices[1]);
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera, tri.vertices[2]);


        CanvasPoint canp1 = projectVertexOntoCanvasPoint(vert1, focalLength, scale, window);
        CanvasPoint canp2 = projectVertexOntoCanvasPoint(vert2, focalLength, scale, window);
        CanvasPoint canp3 = projectVertexOntoCanvasPoint(vert3, focalLength, scale, window);

        // if (pointInCanvas(canp1, window) and pointInCanvas(canp2, window) and pointInCanvas(canp3, window)) {
        fillTriangle(window, CanvasTriangle(canp1, canp2, canp3), colour);
    }
}

void drawRayTracedModel(std::vector<ModelTriangle> model, glm::vec3 modelOrigin, float scale, DrawingWindow &window) {
    // change model coordinate system
    std::vector<glm::vec3> triangleNormals;

    for (int i = 0; i < model.size(); i++) {
        ModelTriangle tri = model.at(i);
        // std::cout << glm::to_string(tri.vertices[0]) << glm::to_string(tri.vertices[1]) << glm::to_string(tri.vertices[2]) << std::endl;

        glm::vec3 vert1 = changeCoordSystem(modelOrigin, camera, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(modelOrigin, camera, tri.vertices[1]);;
        glm::vec3 vert3 = changeCoordSystem(modelOrigin, camera, tri.vertices[2]);

        model[i].vertices[0] = vert1;
        model[i].vertices[1] = vert2;
        model[i].vertices[2] = vert3;

        triangleNormals.push_back(triangleNormal(vert1, vert2, vert3));
    }

    float z = -focalLength;
    CanvasPoint canvasCentre((window.width / 2), (window.height / 2));

    for (int i = 0; i < window.width; i++) {
        for (int j = 0; j < window.height; j++) {
            CanvasPoint pixelLocation = CanvasPoint(i, j);

            // // calculate unit vector from camera to current pixel
            float x = (i - canvasCentre.x) / scale;
            float y = -(j - canvasCentre.y) / scale;

            glm::vec3 canvasPointLocation3D = glm::vec3(x, y, z);

            // apply camera orientation
            // direction of ray from camera to pixel
            glm::vec3 rayDir = normalize(camOr * canvasPointLocation3D - camera);


            //find the closest intersection with a model triangle
            RayTriangleIntersection intersection = getClosestValidIntersection(rayDir, model, camera);


            // if not invalid intersection
            if (intersection.distanceFromCamera > 0) {
                // change light source position from world coordinates to camera coordinates
                glm::vec3 light = changeCoordSystem(glm::vec3(0, 0, 0), camera, lightSource);
                // light = lightSource;

                float epsilon = 0.001;
                glm::vec3 pointToLightRayDir = normalize(light - intersection.intersectionPoint);
                glm::vec3 shadowRayOrigin = intersection.intersectionPoint + pointToLightRayDir * epsilon;
                RayTriangleIntersection shadowRayIntersection =
                        getClosestValidIntersection(pointToLightRayDir, model, shadowRayOrigin);

                // calculate proximity brightness
                float distanceFromLight = length(shadowRayOrigin - light);
                float brightness;
                if (distanceFromLight < 0.0000001) {
                    brightness = 1;
                } else {
                    brightness = 50.0f / (4 * M_PI * distanceFromLight * distanceFromLight);
                    // cap between 0-1
                    // brightness = std::min(brightness, 1.0f);
                    brightness = std::max(0.0f, brightness);
                }

                // angle of incidence brightness
                float aoiBrightness = dot(pointToLightRayDir, triangleNormals[intersection.triangleIndex]);
                // cap between 0-1
                // aoiBrightness = std::min(aoiBrightness, 1.0f);
                aoiBrightness = std::max(0.0f, aoiBrightness);

                //combine proximity brightness and aoi brighness


                //specular brightness
                glm::vec3 reflectionDir = reflect(-pointToLightRayDir, triangleNormals[intersection.triangleIndex]);
                glm::vec3 pointToCameraDir = -rayDir;
                float specularBrightness = std::pow(dot(pointToCameraDir, reflectionDir), specularExponent);
                specularBrightness = std::max(specularBrightness, 0.0f);

                brightness = brightness * aoiBrightness;
                // brightness = aoiBrightness;

                // apply ambient light threshold
                brightness = std::min(1.0f, brightness);
                brightness = std::max(ambientLightThresh, brightness);

                Colour pixelColour = intersection.intersectedTriangle.colour;
                pixelColour.red = std::min(pixelColour.red * brightness + specularBrightness * 100, 255.0f);
                pixelColour.green = std::min(pixelColour.green * brightness + specularBrightness * 100, 255.0f);
                pixelColour.blue = std::min(pixelColour.blue * brightness + specularBrightness * 100, 255.0f);

                // if in shadow
                if (!(shadowRayIntersection.triangleIndex == intersection.triangleIndex or shadowRayIntersection.
                      distanceFromCamera == -std::numeric_limits<float>::infinity())) {
                    // in shadow
                    // pixelColour.red *= shadowDarkness;
                    // pixelColour.green *= shadowDarkness;
                    // pixelColour.blue *= shadowDarkness;

                    pixelColour = shadowRayIntersection.intersectedTriangle.colour;
                }

                drawCanvasPoint(pixelLocation, pixelColour, window);
            }
        }
    }
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_a) {
            // left
            glm::vec3 right = glm::vec3(camOr[0][0], camOr[1][0], camOr[2][0]);
            if (moveLight) {
                lightSource -= cameraSpeed * right;
            } else {
                camera -= cameraSpeed * right;
            }
        } else if (event.key.keysym.sym == SDLK_d) {
            // right
            glm::vec3 right = glm::vec3(camOr[0][0], camOr[1][0], camOr[2][0]);
            if (moveLight) {
                lightSource += cameraSpeed * right;
            } else {
                camera += cameraSpeed * right;
            }
        } else if (event.key.keysym.sym == SDLK_SPACE) {
            // up
            glm::vec3 up = glm::vec3(camOr[0][1], camOr[1][1], camOr[2][1]);
            if (moveLight) {
                lightSource += cameraSpeed * up;
            } else {
                camera += cameraSpeed * up;
            }
        } else if (event.key.keysym.sym == SDLK_LSHIFT) {
            // down
            glm::vec3 up = glm::vec3(camOr[0][1], camOr[1][1], camOr[2][1]);
            if (moveLight) {
                lightSource -= cameraSpeed * up;
            } else {
                camera -= cameraSpeed * up;
            }
        } else if (event.key.keysym.sym == SDLK_w) {
            // forward
            glm::vec3 forward = glm::vec3(camOr[0][2], camOr[1][2], camOr[2][2]);
            if (moveLight) {
                lightSource -= cameraSpeed * forward;
            } else {
                camera -= cameraSpeed * forward;
            }
        } else if (event.key.keysym.sym == SDLK_s) {
            // backward
            glm::vec3 forward = glm::vec3(camOr[0][2], camOr[1][2], camOr[2][2]);
            if (moveLight) {
                lightSource += cameraSpeed * forward;
            } else {
                camera += cameraSpeed * forward;
            }
        } else if (event.key.keysym.sym == SDLK_l) {
            // rotate on x axis
            camOr = rotateOrientation("x", 0.1, camOr);
        } else if (event.key.keysym.sym == SDLK_k) {
            // rotate on y axis
            camOr = rotateOrientation("y", 0.1, camOr);;
        } else if (event.key.keysym.sym == SDLK_o) {
            // rotate on z axis
            camOr = rotateOrientation("z", 0.1, camOr);
        } else if (event.key.keysym.sym == SDLK_v) {
            // toggle lookAt
            lookAt = !lookAt;
        } else if (event.key.keysym.sym == SDLK_h) {
            // spin around model
            if (theta + 0.05 > M_PI * 2) {
                theta = 0.0;
            } else theta += 0.05;

            float r = 4.0;
            glm::vec3 centre(0, 0, 0);

            camera.x = centre.x + std::cos(theta) * r;
            camera.z = centre.z + std::sin(theta) * r;
        } else if (event.key.keysym.sym == SDLK_6) {
            renderMode = WireFrame;
        } else if (event.key.keysym.sym == SDLK_7) {
            renderMode = Rasterized;
        } else if (event.key.keysym.sym == SDLK_8) {
            renderMode = RayTraced;
        } else if (event.key.keysym.sym == SDLK_t) {
            // toggle mouse
            mouseMovement = !mouseMovement;
        } else if (event.key.keysym.sym == SDLK_m) {
            // toggle move light
            moveLight = !moveLight;
        } else if (event.key.keysym.sym == SDLK_p) {
            // print light source position
            std::cout << glm::to_string(lightSource) << std::endl;
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    } else if (event.type == SDL_WINDOWEVENT) {
        if (event.window.event == SDL_WINDOWEVENT_CLOSE) {
            std::cout << "Quiting" << std::endl;
            windowOpen = false;
        }
    } else if (event.type == SDL_MOUSEMOTION) {
        // Get the relative movement of the mouse
        // std::cout << "mouse moved" << std::endl;
        if (mouseMovement) {
            int xOffset = event.motion.xrel;
            int yOffset = event.motion.yrel;

            // Update camera orientation based on mouse movement
            // processMouseMovement(xOffset, yOffset);
        }
    }
}

void draw(DrawingWindow &window, std::vector<ModelTriangle> model) {
    window.clearPixels();


    // initialise depth buffer
    for (auto &row: zDepth) {
        for (auto &elem: row) {
            elem = -100000;
        }
    }

    // render box model
    glm::vec3 modelOrigin = glm::vec3(0, 0, 0);
    float scale = 500;

    if (lookAt) {
        camLookAt(modelOrigin);
    }

    if (renderMode == WireFrame) {
        drawWireframeModel(model, modelOrigin, scale, window);
    } else if (renderMode == Rasterized) {
        drawRasterizedModel(model, modelOrigin, scale, window);
    } else if (renderMode == RayTraced) {
        drawRayTracedModel(model, modelOrigin, scale, window);
    }
}


int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;


    std::vector<ModelTriangle> boxModel = parseObj("cornell-box.obj", "cornell-box.mtl", 0.35);


    while (windowOpen) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window);

        draw(window, boxModel);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }

    return 0;
}
