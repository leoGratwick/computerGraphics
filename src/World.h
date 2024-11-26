//
// Created by User on 8/11/2024.
//

#ifndef WORLD_H
#define WORLD_H
#include <vector>

#include "Camera.h"
#include "LightSource.h"
#include "Model.h"

enum Mode {
    WireFrame,
    Rasterized,
    RayTraced
};

struct World {
    std::vector<LightSource> lights;
    Camera camera;

    // models
    std::vector<Model> models;

    float screenScale;


    // lighting vars
    float ambientLightThresh;
    float shadowDarkness;


    bool phongShading;

    Mode renderMode;
    int currentModelIndex;

    World();

    void printWorldVars();

    Model getCurrentModel();
};


#endif //WORLD_H
