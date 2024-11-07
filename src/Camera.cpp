//
// Created by User on 4/11/2024.
//

#include "Camera.h"
#include "Helpers.h"
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

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
    glm::vec3 forward = orientation[2];
    position -= speed * forward;
}

void Camera::moveBackward() {
    glm::vec3 forward = orientation[2];
    position += speed * forward;
}

void Camera::moveLeft() {
    glm::vec3 left = orientation[0];
    position -= speed * left;
}

void Camera::moveRight() {
    glm::vec3 left = orientation[0];
    position += speed * left;
}

void Camera::moveUp() {
    glm::vec3 up = orientation[1];
    position += speed * up;
}

void Camera::moveDown() {
    glm::vec3 up = orientation[1];
    position -= speed * up;
}

void Camera::rotateX(float angle) {
    orientation = rotateOrientation("x", angle, orientation);
}

void Camera::rotateY(float angle) {
    orientation = rotateOrientation("y", angle, orientation);
}

void Camera::rotateZ(float angle) {
    orientation = rotateOrientation("z", angle, orientation);
}

void Camera::lookAt(glm::vec3 point) {
    // calculate new camera direction unit vectors
    glm::vec3 camDir = normalize(point - position);
    glm::vec3 up(0, 1, 0);
    glm::vec3 left = -normalize(glm::cross(up, camDir));
    glm::vec3 newUp = normalize(glm::cross(camDir, -left));

    // apply camera direction unit vectors to camera orientation
    glm::mat3 newCamOr = glm::mat3(
        left, // First column: right vector
        newUp, // Second column: up vector
        -camDir // Third column: negative direction vector (looking direction)
    );

    orientation = newCamOr;
}
