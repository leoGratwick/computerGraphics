//
// Created by User on 4/11/2024.
//

#ifndef CAMERA_H
#define CAMERA_H

#include "Helpers.h"
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>


struct Camera {
    glm::vec3 position;
    glm::mat3 orientation;
    float speed;
    float focalLength;

    // constructor
    Camera();

    void printStatus() const;

    void moveForward();

    void moveBackward();

    void moveLeft();

    void moveRight();

    void moveUp();

    void moveDown();

    void rotateX(float angle);

    void rotateY(float angle);

    void rotateZ(float angle);

    void lookAt(glm::vec3 point);
};


#endif //CAMERA_H
