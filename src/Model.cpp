//
// Created by User on 12/11/2024.
//

#include "Model.h"
#include "Helpers.h"

// Default constructor
Model::Model() = default;

Model::Model(std::vector<ModelTriangle> triangles, std::map<int, glm::vec3> vertexNormalMap, glm::vec3 center,
             glm::vec3 origin) : triangles(triangles),
                                 vertexNormalMap(vertexNormalMap), center(center), origin(origin), hasTexture(false),
                                 specularExponent(16.0f),
                                 specularIntensity(0.5f) {
    // when model is initialised make sure the triangles have normals
    for (auto &tri: this->triangles) {
        // add a triangle normal
        tri.normal = triangleNormal(tri.vertices[0], tri.vertices[1], tri.vertices[2]);
    }
}

Model::Model(std::vector<ModelTriangle> triangles, std::map<int, glm::vec3> vertexNormalMap) : triangles(triangles),
    vertexNormalMap(vertexNormalMap), center(glm::vec3(0, 0, 0)), origin(glm::vec3(0, 0, 0)), hasTexture(false),
    specularExponent(16.0f),
    specularIntensity(0.5f) {
    // when model is initialised make sure the triangles have normals
    for (auto &tri: this->triangles) {
        // add a triangle normal
        tri.normal = triangleNormal(tri.vertices[0], tri.vertices[1], tri.vertices[2]);
    }
}


void Model::changeCoordSystems(glm::vec3 toOrigin) {
    for (auto &tri: triangles) {
        // calculate vertices in new coordinate system
        glm::vec3 vert1 = changeCoordSystem(origin, toOrigin, tri.vertices[0]);
        glm::vec3 vert2 = changeCoordSystem(origin, toOrigin, tri.vertices[1]);;
        glm::vec3 vert3 = changeCoordSystem(origin, toOrigin, tri.vertices[2]);

        // adjust vertices in model
        tri.vertices[0] = vert1;
        tri.vertices[1] = vert2;
        tri.vertices[2] = vert3;
    }

    center = changeCoordSystem(origin, toOrigin, center);
    origin = toOrigin;
}
