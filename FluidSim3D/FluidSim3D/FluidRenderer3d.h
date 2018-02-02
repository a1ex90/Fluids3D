#pragma once
#include <vector>
#include "RenderUtil.h"
#include "SimUtil.h"

class FluidRenderer3D
{
public:
	/*
	Constructor
	Args:
	lineFileName - name of the .csv file with the line data of the watersurface
	geoFileName - name of the .txt file containing initial geometry				
	*/
	FluidRenderer3D(std::string geoFileName, int mode);
	FluidRenderer3D::FluidRenderer3D(SimUtil::Mat3Di *labels, int gridWidth, int gridHeight, int gridDepth, int mode);
	~FluidRenderer3D();
	/*
	Draws the fluid for at each frame automatically on repeat. 
	applies the chosen visualisation mode
	*/
	void run();
	/*
	Draws the given particles as dots
	*/
	void drawP(std::vector<glm::vec2> particles);
	void drawP(std::vector<glm::vec3> particles);
	/*
	Draws the given vertices as a line
	*/
	void draw(std::vector<glm::vec2> vertices);
	/*
	Draws the given vertices as a mesh
	*/
	void draw(std::vector<glm::vec2> vertices, std::vector<int> indicies);
	void draw(std::vector<glm::vec3> vertices, std::vector<int> indicies);
	/*
	Draws the given vertices as a mesh with given opacities
	*/
	void draw(std::vector<glm::vec2> vertices, std::vector<int> indicies, std::vector<float> opacities);

private:
	//output window
	Display *m_display;
	//constant color shader
	Shader *m_pointShader;
	//constant color shader
	Shader *m_colorShader;
	//textured solids shader
	Shader *m_solidsTexShader;
	//varying opacity shader for fluid
	Shader *m_opacityShader;
	//texture for the solids
	Texture *m_texSolids;
	//texture for the background image
	Texture *m_texImgBackground;

	//visualisation mode switch
	int m_visualizationMode;
	//stores the currentFrame to be displayed
	int m_currentFrame;
	//stores the number of frames
	int m_numberOfFrames;
	/*each element of a vector containing the vertex, index and opacity
	  values for a specific frame*/
	std::vector<std::vector<glm::vec2>> m_partFluid;
	std::vector<std::vector<glm::vec2>> m_vertFluid;
	std::vector<std::vector<int>> m_indFluid;
	std::vector<std::vector<float>> m_opaFluid;

	//Mesh with the geometry of the solids
	Mesh *m_meshSolid;
	//Line for the edges of the solids
	Line *m_borderSolid;
	//Mesh for background display
	Mesh *m_meshBackground;

	/*framerate at which the simulation should be played.
	should match the framerate at which the simulation was calculated*/
	const int FRAME_RATE = 25;

	//filename of the geometry
	std::string m_geoFileName;

	void drawPoints(std::vector<glm::vec2> points);
	void drawPoints(std::vector<glm::vec3> points);
	void drawLines(std::vector<glm::vec2> vertices);
	void drawTriangles(std::vector<glm::vec2> vertices, std::vector<int> indicies);
	void drawOpacityTriangles(std::vector<glm::vec2> vertices, std::vector<int> indicies, std::vector<float> opacities);
	

	void readLines(std::string file, std::vector<std::string> &lines);
	void initGeom(std::string geoFileName);
	void initGeom(SimUtil::Mat3Di *label, int x, int y, int z);
	void initBorderLines(int x, int y, int z, int maxGridSize);
	void initBackground(std::string backgroundFileName);
	void capturePicture(int frame);
	std::vector<glm::vec2> returnVecVertices(std::string lines);
	std::vector<int> returnIndices(std::string lines);
	std::vector<float> returnOpacity(std::string lines);
	std::vector<std::string> split(std::string str, std::string token);
	float uvFloor(float number);

	void bottomLineAt(std::vector<glm::vec3> &lines, int yLoc, int crop, int x, int y, int z, int maxGridSize);
	void sideLineAt(std::vector<glm::vec3> &lines, int xLoc, int crop, int bottomCrop, int topCrop, int x, int y, int z, int maxGridSize);
	void frontLineAt(std::vector<glm::vec3> &lines, int zLoc, int crop, int bottomCrop, int topCrop, int x, int y, int z, int maxGridSize);


};

