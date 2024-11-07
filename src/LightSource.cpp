//
// Created by User on 4/11/2024.
//

#include "LightSource.h"

LightSource::LightSource()
    : position(0.0f, 0.0f, 0.0f), // Initialize position to (0, 0, 0)
      intensity(100.0f),
      speed(0.1f) {
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
