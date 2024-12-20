//
// Created by User on 9/10/2024.
//
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Utils.h>
#include <TextureMap.h>
#include <vector>
#include <bits/stdc++.h>
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>
#include <cerrno>

std::vector<float> lerpSingleFloats(float from, float to, int numberOfValues) {
    float gap = (to - from) / (numberOfValues - 1);
    std::vector<float> vec = {};
    for (int i = 0; i < numberOfValues; i++) {
        vec.push_back(from + gap * i);
    }

    return vec;
}

std::vector<glm::vec3> lerpThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    std::vector<glm::vec3> lerpd = {};
    std::vector<float> col1 = lerpSingleFloats(from[0], to[0], numberOfValues);
    std::vector<float> col2 = lerpSingleFloats(from[1], to[1], numberOfValues);
    std::vector<float> col3 = lerpSingleFloats(from[2], to[2], numberOfValues);

    for (int i = 0; i < numberOfValues; i++) {
        glm::vec3 add(col1.at(i), col2.at(i), col3.at(i));
        lerpd.push_back(add);
        // std::cout << glm::to_string(add) << std::endl;
    }
    return lerpd;
}

std::vector<CanvasPoint> lerpCanvasPoints(CanvasPoint from, CanvasPoint to, int numberOfValues) {
    std::vector<CanvasPoint> lerpd = {};
    std::vector<float> col1 = lerpSingleFloats(from.x, to.x, numberOfValues);
    std::vector<float> col2 = lerpSingleFloats(from.y, to.y, numberOfValues);

    for (int i = 0; i < numberOfValues; i++) {
        CanvasPoint add(col1.at(i), col2.at(i));
        lerpd.push_back(add);
    }
    return lerpd;
}

uint32_t colourToInt(Colour col) {
    uint32_t colour = (255 << 24) + (col.red << 16) + (col.green << 8) + col.blue;

    return colour;
}

Colour intToColour(uint32_t colour) {
    Colour col;
    col.red = (colour >> 16) & 0xFF;
    col.green = (colour >> 8) & 0xFF;
    col.blue = colour & 0xFF;

    return col;
}

float splitPercent(CanvasPoint lineStart, CanvasPoint lineEnd, CanvasPoint point) {
    float splitP = std::abs(point.x - lineStart.x) / std::abs(lineEnd.x - lineStart.x);
    return splitP;
}

bool pointInCanvas(CanvasPoint point, DrawingWindow &window) {
    return ((point.x < window.width - 1) && (point.y < window.height - 1) && (point.y >= 0) && (point.x >= 0));
}

CanvasPoint getPointAlongLine(CanvasPoint lineStart, CanvasPoint lineEnd, float ratio) {
    CanvasPoint point = CanvasPoint(lineStart.x + (lineEnd.x - lineStart.x) * ratio,
                                    lineStart.y + (lineEnd.y - lineStart.y) * ratio);
    return point;
}

std::vector<std::string> splitByDelimiter(std::string str, char delimiter) {
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream ss(str);
    std::string stringDelim = std::string(1, delimiter);

    while (std::getline(ss >> std::ws, token, delimiter)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

uint32_t getColourFromTexture(CanvasPoint point, TextureMap &textureMap) {
    int h = textureMap.height;
    int w = textureMap.width;

    if (std::round(point.x) >= w or std::round(point.y) >= h) {
        // std::cout << point ;
        std::stringstream ss;
        ss << "this point doesnt exist on the texture: " << point << " texture height " << h << " texture width " << w;
        throw std::runtime_error(ss.str());
    }
    int ind = std::round(point.y) * w + std::round(point.x);

    return textureMap.pixels[ind];
}


glm::mat3 rotateOrientation(std::string axis, float angle, glm::mat3 currentOr) {
    glm::mat3 rotationMatrix = glm::mat3();
    if (axis == "x") {
        rotationMatrix = glm::mat3(1, 0, 0, 0, std::cos(angle), std::sin(angle), 0, -std::sin(angle), std::cos(angle));
    } else if (axis == "y") {
        rotationMatrix = glm::mat3(std::cos(angle), 0, -std::sin(angle), 0, 1, 0, std::sin(angle), 0, std::cos(angle));
    } else if (axis == "z") {
        rotationMatrix = glm::mat3(std::cos(angle), std::sin(angle), 0, -std::sin(angle), std::cos(angle), 0, 0, 0, 1);
    }
    glm::mat3 newOr = rotationMatrix * currentOr;
    return newOr;
}

glm::vec3 triangleNormal(glm::vec3 vert1, glm::vec3 vert2, glm::vec3 vert3) {
    return normalize(glm::cross(vert2 - vert1, vert3 - vert1));
}

glm::vec3 changeCoordSystem(glm::vec3 fromOrigin, glm::vec3 toOrigin, glm::vec3 point) {
    glm::vec3 transformationVector = fromOrigin - toOrigin;
    glm::vec3 outPoint = point + transformationVector;
    return outPoint;
}

glm::vec3 addNoisetoDir(glm::vec3 direction, float intensity) {
    direction = glm::normalize(direction);
    glm::vec3 noise = glm::vec3(0, 0, 0);
    // random number between 1 and -1
    noise.x = static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f;
    noise.y = static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f;
    noise.z = static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f;

    noise = noise * intensity;

    return normalize(noise + direction);
}

void redNoise(DrawingWindow &window) {
    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            float red = rand() % 256;
            float green = 0.0;
            float blue = 0.0;
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}

void grayScale1D(DrawingWindow &window) {
    // get horizonal line of values
    std::vector<float> lineCols = lerpSingleFloats(255, 0, window.width);

    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            int col = int(lineCols.at(x));
            uint32_t colour = (225 << 24) + (col << 16) + (col << 8) + col;
            window.setPixelColour(x, y, colour);
        }
    }
}

void rainbow(DrawingWindow &window) {
    glm::vec3 topLeft(255, 0, 0); // red
    glm::vec3 topRight(0, 0, 255); // blue
    glm::vec3 bottomRight(0, 255, 0); // green
    glm::vec3 bottomLeft(255, 255, 0); // yellow

    // left
    std::vector<glm::vec3> left = lerpThreeElementValues(topLeft, bottomLeft, window.height);

    // right
    std::vector<glm::vec3> right = lerpThreeElementValues(topRight, bottomRight, window.height);


    // fill
    for (size_t y = 0; y < window.height; y++) {
        std::vector<glm::vec3> lineCols = lerpThreeElementValues(left.at(y), right.at(y), window.width);
        for (size_t x = 0; x < window.width; x++) {
            glm::vec3 col = lineCols.at(x);
            int colour = colourToInt(Colour(col.r, col.g, col.b));
            window.setPixelColour(x, y, colour);
        }
    }
}



