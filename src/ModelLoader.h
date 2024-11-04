
#ifndef MODELLOADER_H
#define MODELLOADER_H

#include <Colour.h>
#include <vector>
#include <map>
#include "ModelTriangle.h"
#include "Helpers.h"

std::map<std::string, Colour> createPalette(const std::string &filename);

std::vector<ModelTriangle> parseObj(const std::string &objFilename, const std::string &mtlFilename, float scale);

#endif //MODELLOADER_H
