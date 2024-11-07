
#ifndef MODELLOADER_H
#define MODELLOADER_H

#include <Colour.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "ModelTriangle.h"
#include <tuple>

#include "Helpers.h"

std::map<std::string, Colour> createPalette(const std::string &filename);

std::unordered_map<int, glm::vec3> getVertexNormals(
    std::map<glm::vec3, std::vector<glm::vec3> > verticesMap);

std::tuple<std::vector<ModelTriangle>, std::map<int, glm::vec3> > parseObj(
    const std::string &objFilename, float scale, const std::string &mtlFilename = "");

#endif //MODELLOADER_H
