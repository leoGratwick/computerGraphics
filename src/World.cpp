//
// Created by User on 8/11/2024.
//

#include "World.h"

World::World() : lights({LightSource()}),
                 camera(Camera()),
                 ambientLightThresh(0.05f),
                 shadowDarkness(0.4f),
                 phongShading(false),
                 renderMode(WireFrame),
                 models(),
                 screenScale(500.0f),
                 currentModelIndex(0) {
}

void World::printWorldVars() {
    std::cout << "light: " << to_string(lights[0].position) << std::endl;
    std::cout << "camera: " << to_string(camera.position) << std::endl;
    std::cout << "camera orientation: " << to_string(camera.orientation) << std::endl;
    std::cout << "phong shading: " << phongShading << std::endl;

    std::cout << "model info:" << std::endl;

    for (auto &model: models) {
        std::cout << "specular exponent: " << model.specularExponent;
        std::cout << "specular intensity: " << model.specularIntensity << std::endl;
    }
}

Model World::getCurrentModel() {
    if (currentModelIndex < models.size()) {
        return models[currentModelIndex];
    }
    return models[0];
}
