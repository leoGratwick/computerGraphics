//
// Created by User on 9/10/2024.
//
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Utils.h>
#include <iostream>
#include <TextureMap.h>
#include <vector>
#include <bits/stdc++.h>
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>
#include "ModelTriangle.h"
#include <filesystem>
#include <cerrno>

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    float gap = (to - from) / (numberOfValues - 1);
    std::vector<float> vec = {};
    for (int i = 0; i < numberOfValues; i++) {
        vec.push_back(from + gap * i);
    }

    return vec;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    std::vector<glm::vec3> lerpd = {};
    std::vector<float> col1 = interpolateSingleFloats(from[0], to[0], numberOfValues);
    std::vector<float> col2 = interpolateSingleFloats(from[1], to[1], numberOfValues);
    std::vector<float> col3 = interpolateSingleFloats(from[2], to[2], numberOfValues);

    for (int i = 0; i < numberOfValues; i++) {
        glm::vec3 add(col1.at(i), col2.at(i), col3.at(i));
        lerpd.push_back(add);
        // std::cout << glm::to_string(add) << std::endl;
    }
    return lerpd;
}

std::vector<CanvasPoint> lerpCanvasPoints(CanvasPoint from, CanvasPoint to, int numberOfValues) {
    std::vector<CanvasPoint> lerpd = {};
    std::vector<float> col1 = interpolateSingleFloats(from.x, to.x, numberOfValues);
    std::vector<float> col2 = interpolateSingleFloats(from.y, to.y, numberOfValues);

    for (int i = 0; i < numberOfValues; i++) {
        CanvasPoint add(col1.at(i), col2.at(i));
        lerpd.push_back(add);
    }
    return lerpd;
}

uint32_t vec3toColour(glm::vec3 vec) {
    float red = vec[0];
    float green = vec[1];
    float blue = vec[2];
    uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

    return colour;
}

float splitPercent(CanvasPoint lineStart, CanvasPoint lineEnd, CanvasPoint point) {
    float splitP = std::abs(point.x - lineStart.x) / std::abs(lineEnd.x - lineStart.x);
    return splitP;
}

bool pointInWindow(CanvasPoint point, DrawingWindow &window) {
    if (point.x > window.width || point.y > window.height) {
        return false;
    }
    return true;
}

CanvasPoint getPointAlongLine(CanvasPoint lineStart, CanvasPoint lineEnd, float ratio) {
    CanvasPoint point = CanvasPoint(lineStart.x + (lineEnd.x - lineStart.x) * ratio,
                                    lineStart.y + (lineEnd.y - lineStart.y) * ratio);
    return point;
}

std::vector<std::string> splitByDelimiter(std::string str, char delimiter) {
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream stream(str);
    while (std::getline(stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

uint32_t getColourFromTexture(CanvasPoint point, TextureMap &textureMap) {
    int h = textureMap.height;
    int w = textureMap.width;

    if (std::round(point.x) >= w or std::round(point.y) >= h) {
        // std::cout << point ;
        throw std::runtime_error(" this point doesnt exist on the texture");
    }
    int ind = std::round(point.y) * w + std::round(point.x);

    return textureMap.pixels[ind];
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
    std::vector<float> lineCols = interpolateSingleFloats(255, 0, window.width);

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
    std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, window.height);

    // right
    std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, window.height);


    // fill
    for (size_t y = 0; y < window.height; y++) {
        std::vector<glm::vec3> lineCols = interpolateThreeElementValues(left.at(y), right.at(y), window.width);
        for (size_t x = 0; x < window.width; x++) {
            int col = vec3toColour(lineCols.at(x));
            window.setPixelColour(x, y, col);
        }
    }
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

glm::mat3 orthonormalize(const glm::mat3 &mat) {
    glm::vec3 right = glm::normalize(mat[0]);
    glm::vec3 up = glm::normalize(mat[1]);
    glm::vec3 forward = glm::normalize(mat[2]);

    // Recalculate up to ensure orthogonality
    up = glm::normalize(up - glm::dot(up, forward) * forward);
    forward = glm::normalize(forward - glm::dot(forward, right) * right);

    return glm::mat3(right, up, forward);
}

glm::vec3 triangleNormal(glm::vec3 vert1, glm::vec3 vert2, glm::vec3 vert3) {
    return normalize(glm::cross(vert2 - vert1, vert3 - vert1));
}


