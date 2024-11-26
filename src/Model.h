//
// Created by User on 12/11/2024.
//

#ifndef MODEL_H
#define MODEL_H
#include <vector>
#include <map>
#include "ModelTriangle.h"
#include "TextureMap.h"


struct Model {
    std::vector<ModelTriangle> triangles{};
    std::map<int, glm::vec3> vertexNormalMap{};
    glm::vec3 center{};
    glm::vec3 origin{};
    bool hasTexture;
    TextureMap textureMap{};
    float specularExponent;
    float specularIntensity;

    Model();

    Model(std::vector<ModelTriangle> triangles, std::map<int, glm::vec3> vertexNormalMap, glm::vec3 center,
          glm::vec3 origin);

    Model(std::vector<ModelTriangle> triangles, std::map<int, glm::vec3> vertexNormalMap);

    void changeCoordSystems(glm::vec3 toOrigin);
};


#endif //MODEL_H
