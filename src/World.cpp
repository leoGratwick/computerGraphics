//
// Created by User on 8/11/2024.
//

#include "World.h"

World::World() : lights({LightSource()}),
                 camera(Camera()),
                 ambientLightThresh(0.05f),
                 shadowDarkness(0.4f),
                 specularExponent(16.0f),
                 specularIntensity(0.5f),
                 phongShading(false),
                 renderMode(WireFrame) {
}

void World::printWorldVars() {
    std::cout << "light: " << to_string(lights[0].position) << std::endl;
    std::cout << "camera: " << to_string(camera.position) << std::endl;
    std::cout << "camera orientation: " << to_string(camera.orientation) << std::endl;
    std::cout << "specular exponent: " << specularExponent << std::endl;
    std::cout << "specular intensity: " << specularIntensity << std::endl;
    std::cout << "phong shading: " << phongShading << std::endl;
}
