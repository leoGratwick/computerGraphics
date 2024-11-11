#include <Colour.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "ModelTriangle.h"
#include "Helpers.h"
#include <tuple>


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

std::map<int, glm::vec3> getVertexNormals(
    std::map<int, std::vector<glm::vec3> > verticesMap) {
    std::map<int, glm::vec3> vertexNormals;

    // loop through each vertex
    for (const auto &[vertexIndex, triangleNormals]: verticesMap) {
        glm::vec3 triangleNormalSum = glm::vec3(0.0, 0.0, 0.0);

        // calculate the sum of triangle normals
        for (const auto &triangleNormal: triangleNormals) {
            triangleNormalSum += triangleNormal;
        }

        // calculate the average of the triangle normals
        glm::vec3 vertexNormal = triangleNormalSum / float(triangleNormals.size());

        // add to the vertax normals map
        vertexNormals[vertexIndex] = vertexNormal;
    }

    return vertexNormals;
}

std::tuple<std::vector<ModelTriangle>, std::map<int, glm::vec3> > parseObj(
    const std::string &objFilename, float scale, const std::string &mtlFilename = "") {
    // create colour palette
    std::string currentColour;
    std::map<std::string, Colour> palette;
    std::map<int, glm::vec3> vertexNormals;

    if (mtlFilename == "") {
        palette["grey"] = Colour(100, 100, 100);
        currentColour = "grey";
    } else {
        palette = createPalette(mtlFilename);
    }

    std::vector<ModelTriangle> triangles = {};
    std::vector<glm::vec3> vertices = {};
    std::vector<glm::vec3> normals = {};
    // map of verticies to each triangle normal it is associated with
    std::map<int, std::vector<glm::vec3> > verticesMap = {};
    std::vector<std::vector<int> > faces = {};
    std::vector<std::string> faceColours = {};

    // open file
    std::ifstream file(objFilename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // add verticies
            if (line[0] == 'v' and line[1] == ' ') {
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
                    std::cerr << line << std::endl;
                    for (const auto &num: nums) {
                        std::cout << "Parsed element: " << num << std::endl;
                    }
                    return {};
                }
            }
            // vertex normals given
            else if (line[0] == 'v' and line[1] == 'n') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    // get vetex normal values and apply scale
                    float x = std::stof(nums.at(1)) * scale;
                    float y = std::stof(nums.at(2)) * scale;
                    float z = std::stof(nums.at(3)) * scale;

                    glm::vec3 norm(x, y, z);
                    normals.push_back(normalize(norm));
                } catch (...) {
                    std::cerr << "Error parsing OBJ file: couldn't convert vertices from string to float" << std::endl;
                    std::cerr << line << std::endl;
                    for (const auto &num: nums) {
                        std::cout << "Parsed element: " << num << std::endl;
                    }
                    return {};
                }
            }
            // add faces
            else if (line[0] == 'f') {
                std::vector<std::string> nums = splitByDelimiter(line, ' ');
                try {
                    // remove backslashes if they exist
                    // if (nums.at(1).back() == '/') {
                    //     nums.at(1).pop_back(); // vertex 1
                    //     nums.at(2).pop_back(); // vertex 2
                    //     nums.at(3).pop_back(); // vertex 3
                    // }

                    // split the verticies into vertex index / texture point index / normal index
                    std::vector<std::string> numsSplit1 = splitByDelimiter(nums.at(1), '/');
                    std::vector<std::string> numsSplit2 = splitByDelimiter(nums.at(2), '/');
                    std::vector<std::string> numsSplit3 = splitByDelimiter(nums.at(3), '/');


                    std::vector<int> face = {};

                    // add verticies to face
                    try {
                        face.push_back(std::stoi(numsSplit1[0]));
                        face.push_back(std::stoi(numsSplit2[0]));
                        face.push_back(std::stoi(numsSplit3[0]));
                    } catch (...) {
                        std::cerr << "Error parsing OBJ file: error converting face indices" << std::endl;
                    }

                    // try to get normals
                    if (!normals.empty()) {
                        try {
                            vertexNormals[std::stoi(numsSplit1[0])] = normals.at(std::stoi(numsSplit1[2]));
                            vertexNormals[std::stoi(numsSplit2[0])] = normals.at(std::stoi(numsSplit2[2]));
                            vertexNormals[std::stoi(numsSplit3[0])] = normals.at(std::stoi(numsSplit3[2]));
                        } catch (...) {
                            std::cerr << "Error parsing OBJ file: error converting normals" << std::endl;
                        }
                    }


                    //add face
                    faces.push_back(face);
                    faceColours.push_back(currentColour);
                } catch (...) {
                    std::cerr << "Error parsing OBJ file: error parsing face" << std::endl;
                }
            } else if (line.substr(0, 6) == "usemtl" && mtlFilename != "") {
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
                triangle.verticeIndecies = {v1index, v2index, v3index};

                // add triangle
                triangles.push_back(triangle);

                // add to normal of the triangle to each vertex in vertex map
                glm::vec3 normal = triangleNormal(vert1, vert2, vert3);
                verticesMap[v1index].push_back(normal);
                verticesMap[v2index].push_back(normal);
                verticesMap[v3index].push_back(normal);
            } catch (...) {
                std::cerr << "Error parsing OBJ file: error finding vertex for face" << std::endl;
            }
        }
    } else {
        std::cerr << "Could not open file " << objFilename << std::endl;
    }

    // create vertex normals map
    if (vertexNormals.empty()) {
        vertexNormals = getVertexNormals(verticesMap);
    }


    std::cout << "finished loading model from " << objFilename << std::endl;
    return std::make_tuple(triangles, vertexNormals);
}
