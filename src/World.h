//
// Created by User on 8/11/2024.
//

#ifndef WORLD_H
#define WORLD_H
#include <vector>

#include "Camera.h"
#include "LightSource.h"

enum Mode {
    WireFrame,
    Rasterized,
    RayTraced
};

struct World {
    std::vector<LightSource> lights;
    Camera camera;


    // lighting vars
    float ambientLightThresh;
    float shadowDarkness;
    float specularExponent;
    float specularIntensity;

    bool phongShading;

    Mode renderMode;

    World();

    void printWorldVars();
};


#endif //WORLD_H
