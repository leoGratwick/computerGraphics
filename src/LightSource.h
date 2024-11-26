//
// Created by User on 4/11/2024.
//

#ifndef LIGHTSOURCE_H
#define LIGHTSOURCE_H

#include <vector>
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>


struct LightSource {
    glm::vec3 position;
    float intensity;
    float speed;
    float size;

    LightSource();

    void moveForward();

    void moveBackward();

    void moveLeft();

    void moveRight();

    void moveUp();

    void moveDown();

    std::vector<glm::vec3> getLightPoints(glm::vec3 lightDir);
};


#endif //LIGHTSOURCE_H
