//
// Created by User on 4/11/2024.
//

#include "LightSource.h"

#include <stdexcept>
#include <vector>

LightSource::LightSource()
    : position(0.0f, 0.0f, 0.0f), // Initialize position to (0, 0, 0)
      intensity(100.0f),
      speed(0.1f),
      size(0.0f) {
    // Set default speed
}

void LightSource::moveForward() {
    position.z -= speed;
}

void LightSource::moveBackward() {
    position.z += speed;
}

void LightSource::moveLeft() {
    position.x -= speed;
}

void LightSource::moveRight() {
    position.x += speed;
}

void LightSource::moveUp() {
    position.y += speed;
}

void LightSource::moveDown() {
    position.y -= speed;
}

std::vector<glm::vec3> LightSource::getLightPoints(glm::vec3 lightDir) {
    std::vector<glm::vec3> lightPoints = {position};

    // size = 0 means a singular light point

    if (size == 0) {
        return lightPoints;
    }

    // invalid direction vector
    if (length(lightDir) == 0.0f) {
        throw std::invalid_argument("light direction cannot be a zero vector");
    }

    lightDir = normalize(lightDir);

    // get orthogonal basis perpendicular to light direction

    // Choose an arbitrary vector not parallel to the light direction
    glm::vec3 arbitrary = (std::abs(lightDir[0]) < 0.9) ? glm::vec3(1, 0, 0) : glm::vec3(0, 1, 0);
    glm::vec3 w1 = normalize(cross(lightDir, arbitrary));
    glm::vec3 w2 = normalize(cross(lightDir, w1));

    // plane representing circle using orthogonal basis
    // circle(theta) = cos(theta)w1 + sin(theta)w2

    //generate N perpendicular vectors
    // number of angles
    int numberOfVectors = 20;
    //number of points along each perpendicular vector
    int pointsPerVector = 15;
    float spaceBetweenPoints = size / pointsPerVector;


    // for each angle calculate the  perpendicular vector using the circle plane and normalise
    for (int i = 0; i < numberOfVectors; i++) {
        float theta = (2 * M_PI * i) / numberOfVectors;
        glm::vec3 perpendicularVector = std::cos(theta) * w1 + std::sin(theta) * w2;
        perpendicularVector = glm::normalize(perpendicularVector);

        // for each perpendicular vector get M points along it making the max distance the size of the light
        for (int j = 0; j < pointsPerVector; j++) {
            glm::vec3 point = position + perpendicularVector * (spaceBetweenPoints * j);
            lightPoints.push_back(point);
        }
    }


    return lightPoints;
}
