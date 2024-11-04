//
// Created by User on 4/11/2024.
//

#include "Camera.h"
#include <iostream>
#include "Helpers.h"

// Default constructor
Camera::Camera()
    : position(0.0f, 0.0f, 4.0f), // Initialize position to (0, 0, 4)
      orientation(glm::mat3(1.0f)), // Initialize orientation to the identity matrix
      speed(0.1f), // Set default speed
      focalLength(2.0f) {
    // Set default focal length
}

void Camera::printStatus() const {
    std::cout << "Camera Position: (" << position.x << ", " << position.y << ", " << position.z << ")\n";
    std::cout << "Camera Speed: " << speed << "\n";
    std::cout << "Focal Length: " << focalLength << "\n";
}

void Camera::moveForward() {
    glm::vec3 forward = glm::vec3(orientation[0][2], orientation[1][2], orientation[2][2]);
    position -= speed * forward;
}

void Camera::moveBackward() {
    glm::vec3 forward = glm::vec3(orientation[0][2], orientation[1][2], orientation[2][2]);
    position += speed * forward;
}

void Camera::moveLeft() {
    glm::vec3 right = glm::vec3(orientation[0][0], orientation[1][0], orientation[2][0]);
    position -= speed * right;
}

void Camera::moveRight() {
    glm::vec3 right = glm::vec3(orientation[0][0], orientation[1][0], orientation[2][0]);
    position += speed * right;
}

void Camera::moveUp() {
    glm::vec3 up = glm::vec3(orientation[0][1], orientation[1][1], orientation[2][1]);
    position += speed * up;
}

void Camera::moveDown() {
    glm::vec3 up = glm::vec3(orientation[0][1], orientation[1][1], orientation[2][1]);
    position -= speed * up;
}

void Camera::rotateX(float angle) {
    orientation = ("x", angle, orientation);
}

void Camera::rotateY(float angle) {
    orientation = ("y", angle, orientation);
}

void Camera::rotateZ(float angle) {
    orientation = ("z", angle, orientation);
}
