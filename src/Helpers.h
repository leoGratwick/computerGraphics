//
// Created by User on 9/10/2024.
//

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

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) ;

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) ;
std::vector<CanvasPoint> lerpCanvasPoints(CanvasPoint from, CanvasPoint to, int numberOfValues);

uint32_t vec3toColour(glm::vec3 vec);

float splitPercent(CanvasPoint lineStart, CanvasPoint lineEnd, CanvasPoint point);

bool pointInWindow(CanvasPoint point, DrawingWindow& window);

CanvasPoint getPointAlongLine(CanvasPoint lineStart, CanvasPoint lineEnd, float ratio);

std::vector<std::string> splitByDelimiter(std::string str, char delimiter);

uint32_t getColourFromTexture(CanvasPoint point, TextureMap& textureMap);

void redNoise(DrawingWindow &window);

void grayScale1D(DrawingWindow &window);

void rainbow(DrawingWindow &window);


#endif //HELPERS_H
