
#ifndef HELPERS_H
#define HELPERS_H
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

std::vector<float> lerpSingleFloats(float from, float to, int numberOfValues);

std::vector<glm::vec3> lerpThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);

std::vector<CanvasPoint> lerpCanvasPoints(CanvasPoint from, CanvasPoint to, int numberOfValues);

uint32_t colourToInt(Colour col);

float splitPercent(CanvasPoint lineStart, CanvasPoint lineEnd, CanvasPoint point);

bool pointInCanvas(CanvasPoint point, DrawingWindow &window);

CanvasPoint getPointAlongLine(CanvasPoint lineStart, CanvasPoint lineEnd, float ratio);

std::vector<std::string> splitByDelimiter(std::string str, char delimiter);

uint32_t getColourFromTexture(CanvasPoint point, TextureMap &textureMap);

glm::mat3 rotateOrientation(std::string axis, float angle, glm::mat3 currentOr);

glm::vec3 triangleNormal(glm::vec3 vert1, glm::vec3 vert2, glm::vec3 vert3);

void redNoise(DrawingWindow &window);

void grayScale1D(DrawingWindow &window);

void rainbow(DrawingWindow &window);


#endif //HELPERS_H
