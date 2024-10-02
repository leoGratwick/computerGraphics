#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Utils.h>
#include <fstream>
#include <iostream>
#include <TextureMap.h>
#include <vector>
#include <bits/stdc++.h>
#include <glm/detail/type_vec.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240

bool windowOpen = true;

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {

	float gap = (to - from) / (numberOfValues - 1 ) ;
	std::vector<float> vec = {};
	for (int i = 0; i < numberOfValues; i++) {
		vec.push_back(from + gap * i);
	}

	return vec;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> lerpd = {};
	std::vector<float>  col1 = interpolateSingleFloats(from[0], to[0], numberOfValues);
	std::vector<float>  col2 = interpolateSingleFloats(from[1], to[1], numberOfValues);
	std::vector<float>  col3 = interpolateSingleFloats(from[2], to[2], numberOfValues);

	for (int i = 0; i < numberOfValues; i++) {
		glm::vec3 add(col1.at(i), col2.at(i), col3.at(i))	;
		lerpd.push_back(add);
		// std::cout << glm::to_string(add) << std::endl;
	}
	return lerpd;

}

uint32_t vec3toColour(glm::vec3 vec) {
	float red = vec[0];
	float green = vec[1];
	float blue = vec[2];
	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

	return colour;
}

void drawCanvasPoint(CanvasPoint point, Colour col, DrawingWindow& window) {
	int x = std::round(point.x);
	int y = std::round(point.y);
	uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue);
	window.setPixelColour(x, y, colour);

}

void redNoise(DrawingWindow &window) {
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void grayScale1D(DrawingWindow &window) {
	// get horizonal line of values
	std::vector<float> lineCols = interpolateSingleFloats(255, 0, window.width);

	for(size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			int col = int(lineCols.at(x));
			uint32_t colour = (225 << 24) + (col << 16) + (col << 8) + col;
			window.setPixelColour(x, y, colour);
		}
	}

}

void rainbow(DrawingWindow &window) {

	glm::vec3 topLeft(255, 0, 0);        // red
	glm::vec3 topRight(0, 0, 255);       // blue
	glm::vec3 bottomRight(0, 255, 0);    // green
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	// left
	std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, window.height);

	// right
	std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, window.height);


	// fill
	for(size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> lineCols = interpolateThreeElementValues( left.at(y),right.at(y), window.width);
		for (size_t x = 0; x < window.width; x++) {
			int col = vec3toColour(lineCols.at(x));
			window.setPixelColour(x, y, col);
		}
	}


}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
	// std::cout << "drawing line from: " << from.x << ", " << from.y << " to: " << to.x << ", " << to.y << std::endl;
	int range = std::max(std::abs(to.x - from.x), std::abs(to.y - from.y))*5;
	std::vector<float> lerpX = interpolateSingleFloats(from.x, to.x, range);
	std::vector<float> lerpY = interpolateSingleFloats(from.y, to.y, range);

	for (int i = 0; i < range; i++) {
		CanvasPoint point = CanvasPoint(lerpX.at(i), lerpY.at(i));
		drawCanvasPoint(point, colour, window);

	}

}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle.v0(), triangle.v1(), colour);
	drawLine(window, triangle.v0(), triangle.v2(), colour);
	drawLine(window, triangle.v1(), triangle.v2(), colour);
}

void fillFlatBottomTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	std::cout << "drawing flat bottom triangle" << std::endl;
	// the first vertex is at the top of the triangle

	CanvasPoint v0 = triangle.v0();
	CanvasPoint v1 = triangle.v1();
	CanvasPoint v2 = triangle.v2();

	float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
	float slope2 = (v2.x - v0.x) / (v2.y - v0.y);

	float currentX1 = v0.x;
	float currentX2 = v0.x;

	for (int y = v0.y; y <= v1.y; y++) {
		CanvasPoint from = CanvasPoint(currentX1, y);
		CanvasPoint to = CanvasPoint(currentX2, y);
		drawLine(window, from, to, colour);
		currentX1 += slope1;
		currentX2 += slope2;
	}

}

void fillFlatTopTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	std::cout << "drawing flat top triangle" << std::endl;
	// v2 id the bottom vertex
	CanvasPoint v0 = triangle.v0();
	CanvasPoint v1 = triangle.v1();
	CanvasPoint v2 = triangle.v2();

	float slope1 = (v0.x - v2.x) / (v0.y - v2.y);
	float slope2 = (v1.x - v2.x) / (v1.y - v2.y);

	float currentX1 = v2.x;
	float currentX2 = v2.x;

	for (int y = v2.y; y >= v1.y; y--) {
		CanvasPoint from = CanvasPoint(currentX1, y);
		CanvasPoint to = CanvasPoint(currentX2, y);
		drawLine(window, from, to, colour);
		currentX1 -= slope1;
		currentX2 -= slope2;
	}
}

void fillTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	std::vector<CanvasPoint> verts = {triangle.v0(), triangle.v1(), triangle.v2()};

	// sort verticies by y coordinate
	std::sort(verts.begin(), verts.end(), [](const CanvasPoint& a, const CanvasPoint& b) {
		return a.y < b.y;  // Sort by y
	});

	CanvasPoint v0 = verts[0];;
	std::cout << v0 << std::endl;
	CanvasPoint v1 = verts[1];
	std::cout << v1 << std::endl;
	CanvasPoint v2 = verts[2];
	std::cout << v2 << std::endl;

	// flat bottom triangle
	if (v2.y == v1.y) {
		CanvasTriangle tri = CanvasTriangle(v0, v1, v2);
		fillFlatBottomTriangle(window, tri, colour);
	}
	// flat top triangle
	else if(v0.y == v1.y) {
		CanvasTriangle tri = CanvasTriangle(v0, v1, v2);
		fillFlatBottomTriangle(window, tri, colour);
	}
	// all other triangles
	else {
		float m = (v2.y-v0.y)/(v2.x-v0.x);
		float x = (v1.y - v0.y)/m + v0.x;
		CanvasPoint v3 = CanvasPoint(x, v1.y);
		std::cout << v3 << std::endl;
		CanvasTriangle flatBottomTri = CanvasTriangle(v0, v1, v3);
		fillFlatBottomTriangle(window, flatBottomTri, colour);

		CanvasTriangle flatTopTri = CanvasTriangle(v1, v3, v2);
		fillFlatTopTriangle(window, flatTopTri, colour);

	}
}

void randomTriangle(DrawingWindow &window) {
	CanvasPoint v0 = CanvasPoint(rand()%(window.width-1), rand()%(window.height-1));
	CanvasPoint v1 = CanvasPoint(rand()%(window.width-1), rand()%(window.height-1));
	CanvasPoint v2 = CanvasPoint(rand()%(window.width-1), rand()%(window.height-1));
	Colour randomCol = Colour(rand()%255, rand()%255, rand()%255);
	CanvasTriangle randomTriangle = CanvasTriangle( v0, v1, v2);
	fillTriangle(window, randomTriangle, randomCol);
	drawTriangle(window, randomTriangle, Colour(255,255,255));
	window.renderFrame();
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_c) {
			std::cout << "closing" << std::endl;
			windowOpen = false;
		}
		else if (event.key.keysym.sym == SDLK_u) {
			std::cout << "generate random triangle" << std::endl;
			randomTriangle(window);

		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	else if (event.type == SDL_WINDOWEVENT) {
		if (event.window.event == SDL_WINDOWEVENT_CLOSE) {
			std::cout << "Quiting" << std::endl;
			windowOpen = false;
		}

	}
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	rainbow(window);



}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	// std::string myText;
	//
	// // Read from the text file
	// std::ifstream MyReadFile("test.txt");
	//
	// // Use a while loop together with the getline() function to read the file line by line
	// while (getline (MyReadFile, myText)) {
	// 	// Output the text from the file
	// 	std::cout << myText;
	// }
	//
	// // Close the file
	// MyReadFile.close();

	const std::string textureFile = "texture.ppm";
	TextureMap brickTexture = TextureMap(textureFile);


	while (windowOpen) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

	return 0;
}
