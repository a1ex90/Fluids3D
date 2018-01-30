#include "FluidRenderer3d.h"
#include <iostream>
#include <string>
#include <fstream>
#include <SDL2/SDL.h>
#include "timing.h"

#define WIDTH 800
#define HEIGHT 600

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------
FluidRenderer3D::FluidRenderer3D(std::string geoFileName, int mode) {
	m_visualizationMode = mode;

	m_display = new Display{ WIDTH, HEIGHT, "2D Fluid Simulation" };

	m_pointShader = new Shader{ "./pointShader" };
	m_colorShader = new Shader{ "./basicShader" };
	m_opacityShader = new Shader{ "./opacityShader" };
	m_solidsTexShader = new Shader{ "./textureShader" };

	m_texSolids = new Texture("./sand-texture.jpg");

	m_geoFileName = geoFileName;

	initGeom(geoFileName);
	initBackground("./background-dune.png");

	m_currentFrame = 0;
}

FluidRenderer3D::FluidRenderer3D(std::string FileName, std::string geoFileName, int mode) {
	m_visualizationMode = mode;

	m_display = new Display{ WIDTH, HEIGHT, "2D Fluid Simulation" };

	m_pointShader = new Shader{"./pointShader"};
	m_colorShader = new Shader{ "./basicShader" };
	m_opacityShader = new Shader{ "./opacityShader" };
	m_solidsTexShader = new Shader{ "./textureShader" };

	m_texSolids = new Texture("./sand-texture.jpg");

	m_geoFileName = geoFileName;

	initGeom(geoFileName);
	initBackground("./background-dune.png");

	m_currentFrame = 0;

	if (mode == 0) {
		std::vector<std::string> linesFluid;
		readLines(FileName + "-lines.csv", linesFluid);
		m_numberOfFrames = linesFluid.size();
		for (int i = 0; i < linesFluid.size(); i++) {
			m_vertFluid.push_back(returnVecVertices(linesFluid[i]));
		}
	}
	else if (mode == 1) {
		std::vector<std::string> vertFluid;
		std::vector<std::string> indFluid;
		readLines(FileName + "-vert.csv", vertFluid);
		readLines(FileName + "-ind.csv", indFluid);
		m_numberOfFrames = vertFluid.size();
		for (int i = 0; i < vertFluid.size(); i++) {
			m_vertFluid.push_back(returnVecVertices(vertFluid[i]));
		}
		for (int i = 0; i < indFluid.size(); i++) {
			m_indFluid.push_back(returnIndices(indFluid[i]));
		}
	}
	else if (mode == 2) {
		std::vector<std::string> vertFluid;
		std::vector<std::string> indFluid;
		std::vector<std::string> opaFluid;
		readLines(FileName + "-vert.csv", vertFluid);
		readLines(FileName + "-ind.csv", indFluid);
		readLines(FileName + "-opacity.csv", opaFluid);
		m_numberOfFrames = vertFluid.size();
		for (int i = 0; i < vertFluid.size(); i++) {
			m_vertFluid.push_back(returnVecVertices(vertFluid[i]));
			m_opaFluid.push_back(returnOpacity(opaFluid[i]));
		}
		for (int i = 0; i < indFluid.size(); i++) {
			m_indFluid.push_back(returnIndices(indFluid[i]));
		}
	}
	else if (mode == 3) {
		std::vector<std::string> partFluid;
		readLines(FileName + "-part.csv", partFluid);
		m_numberOfFrames = partFluid.size();
		for (int i = 0; i < partFluid.size(); i++) {
			m_partFluid.push_back(returnVecVertices(partFluid[i]));
		}
	}
	else if (mode == 4) {
		std::vector<std::string> partFluid;
		readLines(FileName + "-part.csv", partFluid);
		m_numberOfFrames = partFluid.size();
		for (int i = 0; i < partFluid.size(); i++) {
			m_partFluid.push_back(returnVecVertices(partFluid[i]));
		}
		std::vector<std::string> linesFluid;
		readLines(FileName + "-lines.csv", linesFluid);
		for (int i = 0; i < linesFluid.size(); i++) {
			m_vertFluid.push_back(returnVecVertices(linesFluid[i]));
		}
	}
	else {
		std::cout << "Renderer invalid mode selected \n";
	}
}

FluidRenderer3D::~FluidRenderer3D() {
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidRenderer3D::run() {
	while (!m_display->isClosed()) {
		//m_display->clear(0.686f, 0.933f, 0.933f, 1.0f);
		m_solidsTexShader->bind();
		m_texImgBackground->bind(1);
		m_solidsTexShader->setTexture(1);
		m_meshBackground->draw();

		if (m_visualizationMode == 0) {
			drawLines(m_vertFluid[m_currentFrame]);
		}
		else if (m_visualizationMode == 1) {
			drawTriangles(m_vertFluid[m_currentFrame], m_indFluid[m_currentFrame]);
		}
		else if (m_visualizationMode == 2) {
			drawOpacityTriangles(m_vertFluid[m_currentFrame], m_indFluid[m_currentFrame], m_opaFluid[m_currentFrame]);
		}
		else if (m_visualizationMode == 3) {
			drawPoints(m_partFluid[m_currentFrame]);
		}
		else if (m_visualizationMode == 4) {
			drawPoints(m_partFluid[m_currentFrame]);
			drawLines(m_vertFluid[m_currentFrame]);
		}

		//m_colorShader->setColor(0.741f, 0.718f, 0.420f);
		m_solidsTexShader->bind();
		m_texSolids->bind(0);		
		m_solidsTexShader->setTexture(0);
		m_meshSolid->draw();

		m_display->update();

		//Picture Capturing
		//capturePicture(m_currentFrame);
		//End Picture Capturing

		m_currentFrame++;
		if (m_currentFrame >= m_numberOfFrames) {
			m_currentFrame = 0;
		}

		SDL_Delay(1000 / FRAME_RATE);
	}
}

void FluidRenderer3D::drawP(std::vector<glm::vec2> particles) {
	m_solidsTexShader->bind();
	m_texImgBackground->bind(1);
	m_solidsTexShader->setTexture(1);
	m_meshBackground->draw();

	drawPoints(particles);

	m_solidsTexShader->bind();
	m_texSolids->bind(0);
	m_solidsTexShader->setTexture(0);
	m_meshSolid->draw();

	m_display->update();
}

void FluidRenderer3D::draw(std::vector<glm::vec2> vertices) {
	m_solidsTexShader->bind();
	m_texImgBackground->bind(1);
	m_solidsTexShader->setTexture(1);
	m_meshBackground->draw();

	drawLines(vertices);

	m_solidsTexShader->bind();
	m_texSolids->bind(0);
	m_solidsTexShader->setTexture(0);
	m_meshSolid->draw();

	m_display->update();
}

void FluidRenderer3D::draw(std::vector<glm::vec2> vertices, std::vector<int> indicies) {
	m_solidsTexShader->bind();
	m_texImgBackground->bind(1);
	m_solidsTexShader->setTexture(1);
	m_meshBackground->draw();

	drawTriangles(vertices, indicies);

	m_solidsTexShader->bind();
	m_texSolids->bind(0);
	m_solidsTexShader->setTexture(0);
	m_meshSolid->draw();

	m_display->update();
}

void FluidRenderer3D::draw(std::vector<glm::vec2> vertices, std::vector<int> indicies, std::vector<float> opacities) {
	m_solidsTexShader->bind();
	m_texImgBackground->bind(1);
	m_solidsTexShader->setTexture(1);
	m_meshBackground->draw();

	drawOpacityTriangles(vertices, indicies, opacities);

	m_solidsTexShader->bind();
	m_texSolids->bind(0);
	m_solidsTexShader->setTexture(0);
	m_meshSolid->draw();

	m_display->update();
}

//----------------------------------------------------------------------
// Private Functions
//----------------------------------------------------------------------

/*
Draws the given points as dots in the scene
Args:
points - array of 2d points
*/
void FluidRenderer3D::drawPoints(std::vector<glm::vec2> points) {
	m_pointShader->bind();
	m_pointShader->setColor(0.000f, 0.000f, 0.804f);

	Point point{ points };
	point.draw();
}

/*
Draws the given vertices as lines in the scene
Args:
vertices - array with 2d points. each point is the beginning of one line and the end of another line
*/
void FluidRenderer3D::drawLines(std::vector<glm::vec2> vertices) {
	m_colorShader->bind();
	m_colorShader->setColor(0.000f, 0.000f, 0.804f);

	Line line{ vertices };
	line.draw();
}

/*
Draws the given vertices as a triangle mesh
Args:
vertices - vertices of the triangle mesh
indicies - indicies of the triangle mesh
*/
void FluidRenderer3D::drawTriangles(std::vector<glm::vec2> vertices, std::vector<int> indicies) {
	m_colorShader->bind();
	m_colorShader->setColor(0.000f, 0.000f, 0.804f);

	Mesh meshFluid{ vertices, indicies };
	meshFluid.draw();
}

/*
Draws the given vertices as a opaque triangle mesh
Args:
vertices - vertices of the triangle mesh
indicies - indicies of the triangle mesh
opacities - opacity values of each triangle
*/
void FluidRenderer3D::drawOpacityTriangles(std::vector<glm::vec2> vertices, std::vector<int> indicies, std::vector<float> opacities) {
	m_opacityShader->bind();
	m_opacityShader->setColor(0.000f, 0.000f, 0.804f);

	Mesh meshFluid{ vertices, indicies, opacities };
	meshFluid.draw();

	m_colorShader->bind();
}

/*
Initializes the lines array with the lines from a given file
Args:
file - filename
*/
void FluidRenderer3D::readLines(std::string file, std::vector<std::string>& lines) {
	std::ifstream linesFile(file);
	if (linesFile.is_open()) {
		while (linesFile.good()) {
			std::string lineStr;
			std::getline(linesFile, lineStr);
			lines.push_back(lineStr);
		}
		linesFile.close();
	}
}

/*
Reads in the inital geometry and generates the vertices of the solids
in a vector
Args:
file - name of the .txt file containing initial geometry
*/
void FluidRenderer3D::initGeom(std::string file) {
	//defines how often the texture should get repeated in x
	int xRepitions = 3;
	//defines how often the texture should get repeated in y
	int yRepitions = 3;

	std::vector<glm::vec2> vertSolid;
	std::vector<glm::vec2> texCoordsSolid;
	std::vector<int> indSolid;
	int currentInd = 0;

	std::vector <std::string> lines;
	std::ifstream geoFile(file);
	if (geoFile.is_open()) {
		while (geoFile.good()) {
			std::string lineStr;
			std::getline(geoFile, lineStr);
			lines.push_back(lineStr);
		}
		geoFile.close();
	}
	int height = lines.size();
	for (int j = 0; j < height; j++) {
		//iterate over in reverse order since i=0,j=0 is bottom-left
		std::string line = lines[height - 1 - j];
		int width = line.size();
		for (int i = 0; i < width; i++) {
			if (line[i] == 's') {
				float x1 = 2.0f * i / (width)-1.0f;
				float u1 = 1.0f * i / width * xRepitions - uvFloor(1.0f * i / width * xRepitions);
				float x2 = 2.0f * (i + 1) / (width)-1.0f;
				float u2 = 1.0f * (i + 1) / width * xRepitions - uvFloor(1.0f * (i + 1) / width * xRepitions);
				float y1 = 2.0f * j / (height)-1.0f;
				float v1 = 1.0 - 1.0f * j / height * yRepitions - uvFloor(1.0f * j / height * yRepitions);
				float y2 = 2.0f * (j + 1) / (height)-1.0f;
				float v2 = 1.0 - 1.0f * (j + 1) / height * yRepitions - uvFloor(1.0f * (j + 1) / height * yRepitions);
				//bottom-left
				vertSolid.push_back(glm::vec2(x1, y1));
				texCoordsSolid.push_back(glm::vec2(u1, v1));
				indSolid.push_back(currentInd);
				//top-left
				vertSolid.push_back(glm::vec2(x1, y2));
				texCoordsSolid.push_back(glm::vec2(u1, v2));
				indSolid.push_back(currentInd + 1);
				//top-right
				vertSolid.push_back(glm::vec2(x2, y2));
				texCoordsSolid.push_back(glm::vec2(u2, v2));
				indSolid.push_back(currentInd + 2);
				//bottom-left
				indSolid.push_back(currentInd);
				//top-right
				indSolid.push_back(currentInd + 2);
				//bottom-right
				vertSolid.push_back(glm::vec2(x2, y1));
				texCoordsSolid.push_back(glm::vec2(u2, v1));
				indSolid.push_back(currentInd + 3);
				
				currentInd += 4;
			}
		}
	}
	m_meshSolid = new Mesh{ vertSolid , texCoordsSolid, indSolid };
}

/*
sets a given image as a background texture for the whole frame
Args:
file - image file
*/
void FluidRenderer3D::initBackground(std::string file) {
	std::vector<glm::vec2> vertSolid;
	std::vector<glm::vec2> texCoordsSolid;
	std::vector<int> indSolid;
	//bottom-left
	vertSolid.push_back(glm::vec2(-1.0f, -1.0f));
	texCoordsSolid.push_back(glm::vec2(0.0f, 1.0f));
	indSolid.push_back(0);
	//top-left
	vertSolid.push_back(glm::vec2(-1.0f, 1.0f));
	texCoordsSolid.push_back(glm::vec2(0.0f, 0.0f));
	indSolid.push_back(1);
	//top-right
	vertSolid.push_back(glm::vec2(1.0f, 1.0f));
	texCoordsSolid.push_back(glm::vec2(1.0f, 0.0f));
	indSolid.push_back(2);
	//bottom-left
	indSolid.push_back(0);
	//top-right
	indSolid.push_back(2);
	//bottom-right
	vertSolid.push_back(glm::vec2(1.0f, -1.0f));
	texCoordsSolid.push_back(glm::vec2(1.0f, 1.0f));
	indSolid.push_back(3);
	m_meshBackground = new Mesh{ vertSolid , texCoordsSolid, indSolid };

	m_texImgBackground = new Texture(file);
}

/*
reads in a string with data and returns a vector with each seperate data point
Args:
lines - data string, format "a1 a2 a3 a4 ...."
*/
std::vector<glm::vec2> FluidRenderer3D::returnVecVertices(std::string lines) {
	std::vector<std::string> numbers = split(lines, " ");
	std::vector<glm::vec2> vertices;
	for (int i = 0; i < numbers.size() / 2; i++) {
		vertices.push_back(glm::vec2(std::stof(numbers[2 * i]), std::stof(numbers[2 * i + 1])));
	}
	return vertices;
}

std::vector<int> FluidRenderer3D::returnIndices(std::string lines) {
	std::vector<std::string> numbers = split(lines, " ");
	std::vector<int> indices;
	for (int i = 0; i < numbers.size(); i++) {
		indices.push_back(std::stoi(numbers[i]));
	}
	return indices;
}

std::vector<float> FluidRenderer3D::returnOpacity(std::string lines) {
	std::vector<std::string> numbers = split(lines, " ");
	std::vector<float> opacities;
	for (int i = 0; i < numbers.size(); i++) {
		opacities.push_back(std::stof(numbers[i]));
	}
	return opacities;
}


/*
	renders the current frame buffer to a bitmap picture. stores the picture
	as "'geoFileName'-frame-'frameNumber'.bmp"
	Args:
	frame - current frame number
*/
void FluidRenderer3D::capturePicture(int frame) {

	SDL_Surface * temp = SDL_CreateRGBSurface(SDL_SWSURFACE, WIDTH, HEIGHT, 24, 0x000000FF, 0x0000FF00, 0x00FF0000, 0);

	char * pixels = new char[3 * WIDTH * HEIGHT];

	glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	for (int i = 0; i < HEIGHT; i++)
		std::memcpy(((char *)temp->pixels) + temp->pitch * i, pixels + 3 * WIDTH * (HEIGHT - i - 1), WIDTH * 3);

	delete[] pixels;

	std::string filename = m_geoFileName.substr(0, m_geoFileName.length() - 4);
	filename.append("-frame-");
	filename.append(std::to_string(frame));
	filename.append(".bmp");


	SDL_SaveBMP(temp, filename.c_str());
}

//----------------------------------------------------------------------
// Private Helper Functions
//----------------------------------------------------------------------

/*
Splits a given string at the given token and returns a vector of the splitted strings
Args:
str - String to split
token - token at which the string should be split
*/
std::vector<std::string> FluidRenderer3D::split(std::string str, std::string token) {
	std::vector<std::string>result;
	while (str.size()) {
		int index = str.find(token);
		if (index != std::string::npos) {
			result.push_back(str.substr(0, index));
			str = str.substr(index + token.size());
			if (str.size() == 0)result.push_back(str);
		}
		else {
			result.push_back(str);
			str = "";
		}
	}
	return result;
}

float FluidRenderer3D::uvFloor(float number) {
	float a = floorf(number);
	if (a == number && a != 0.0f) {
		return a - 1.0f;
	}
	return a;
}