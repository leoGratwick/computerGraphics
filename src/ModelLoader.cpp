#include <Colour.h>
#include <vector>
#include <map>
#include "ModelTriangle.h"
#include "Helpers.h"


std::map<std::string, Colour> createPalette(const std::string &filename) {
    std::map<std::string, Colour> palette;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            // when a new material is defined
            if (line.substr(0, 6) == "newmtl") {
                std::vector<std::string> tokens = splitByDelimiter(line, ' ');
                std::string colourName = tokens.at(1);

                getline(file, line);
                std::vector<std::string> colourVals = splitByDelimiter(line, ' ');

                // convert from 0-1 to 0-225 for rgb values
                int r = std::round(std::stof(colourVals.at(1)) * 225);
                int g = std::round(std::stof(colourVals.at(2)) * 225);
                int b = std::round(std::stof(colourVals.at(3)) * 225);

                // add colour to palette
                palette[colourName] = Colour(r, g, b);
            }
        }
        file.close();
    } else {
        std::cerr << "Could not open file " << filename << std::endl;
    }
    return palette;
}

std::vector<ModelTriangle> parseObj(const std::string &objFilename, const std::string &mtlFilename, float scale) {
    // create colour palette
    std::map<std::string, Colour> palette = createPalette(mtlFilename);
    std::vector<ModelTriangle> triangles = {};
    std::vector<glm::vec3> vertices = {};
    std::vector<std::vector<int> > faces = {};
    std::vector<std::string> faceColours = {};

    // open file
    std::ifstream file(objFilename);
    if (file.is_open()) {
        std::string currentColour;
        std::string line;
        while (std::getline(file, line)) {
            // add verticies
            if (line[0] == 'v') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    // get vetex coordinates and apply scale
                    float x = std::stof(nums.at(1)) * scale;
                    float y = std::stof(nums.at(2)) * scale;
                    float z = std::stof(nums.at(3)) * scale;

                    glm::vec3 vert(x, y, z);
                    vertices.push_back(vert);
                } catch (...) {
                    std::cerr << "Error parsing OBJ file: couldn't convert vertices from string to float" << std::endl;
                }
            }
            // add faces
            else if (line[0] == 'f') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    // remove backslashes
                    nums.at(1).pop_back(); // vertex 1
                    nums.at(2).pop_back(); // vertex 2
                    nums.at(3).pop_back(); // vertex 3

                    std::vector<int> face = {};

                    // add verticies to face
                    try {
                        face.push_back(std::stoi(nums[1]));
                        face.push_back(std::stoi(nums[2]));
                        face.push_back(std::stoi(nums[3]));
                    } catch (...) {
                        std::cerr << "Error parsing OBJ file: error converting face indices" << std::endl;
                    }

                    //add face
                    faces.push_back(face);
                    faceColours.push_back(currentColour);
                } catch (...) {
                    std::cerr << "Error parsing OBJ file: error parsing face" << std::endl;
                }
            } else if (line.substr(0, 6) == "usemtl") {
                // changing current colour
                std::vector<std::string> tokens = splitByDelimiter(line, ' ');
                try {
                    std::string colourName = tokens.at(1);
                    currentColour = colourName;
                } catch (...) {
                    std::cerr << "Error parsing OBJ file: error parsing usemtl line" << std::endl;
                }
            }
        }

        file.close();

        // add every face along with its colour
        for (int i = 0; i < faces.size(); i++) {
            std::vector<int> face = faces.at(i);
            std::string faceColour = faceColours.at(i);

            try {
                // get each vertex's index and correct for 1 based indexing in obj files then retrieve the vertex
                int v1index = face[0] - 1;
                glm::vec3 vert1(vertices.at(v1index));

                int v2index = face[1] - 1;
                glm::vec3 vert2(vertices.at(v2index));

                int v3index = face[2] - 1;
                glm::vec3 vert3(vertices.at(v3index));

                // create triangle
                ModelTriangle triangle = ModelTriangle(vert1, vert2, vert3, palette[faceColour]);

                // add triangle
                triangles.push_back(triangle);
            } catch (...) {
                std::cerr << "Error parsing OBJ file: error finding vertex for face" << std::endl;
            }
        }
    } else {
        std::cerr << "Could not open file " << objFilename << std::endl;
    }
    return triangles;
}
