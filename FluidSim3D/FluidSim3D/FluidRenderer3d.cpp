#include "FluidRenderer3d.h"


#define WIDTH 800
#define HEIGHT 600

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------

FluidRenderer3D::FluidRenderer3D(SimUtil::Mat3Di *labels, int gridWidth, int gridHeight, int gridDepth) {
	m_transform = new Transform();
	m_display = new Display{ WIDTH, HEIGHT, "2D Fluid Simulation", m_transform };
	m_camera = new Camera(glm::vec3(0, 0, -4), 70.0f, (float)WIDTH / (float)HEIGHT, 0.01f, 1000.0f);	
	m_colorShader = new Shader{ "./basicShader" };
	initGeom(labels, gridWidth, gridHeight, gridDepth);
	//initially pause the simulation
	m_isPaused = true;
	m_forwardPressed = false;
	m_orientation = glm::vec3(0.0f, -1.0f, 0.0f);
}

FluidRenderer3D::~FluidRenderer3D() {
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidRenderer3D::drawP(std::vector<glm::vec3> particles) {
	m_display->clear(0.686f, 0.933f, 0.933f, 1.0f);

	m_colorShader->bind();
	m_colorShader->update(m_transform, m_camera);
	m_colorShader->setColor(0.0f, 0.0f, 0.0f, 1.0f);

	m_borderSolid->draw();

	Point point{ particles };
	point.draw();

	m_colorShader->setColor(0.392f, 0.584f, 0.929f, 0.1f);

	m_meshSolid->draw();

	m_display->update(m_orientation, m_isPaused, m_forwardPressed);
}

void FluidRenderer3D::draw(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<int> indicies) {
	m_display->clear(0.686f, 0.933f, 0.933f, 1.0f);

	m_colorShader->bind();
	m_colorShader->update(m_transform, m_camera);
	m_colorShader->setColor(0.0f, 0.0f, 0.0f, 1.0f);

	m_borderSolid->draw();

	Mesh mesh{ vertices, normals, indicies };

	mesh.draw();

	m_colorShader->setColor(0.392f, 0.584f, 0.929f, 0.1f);

	m_meshSolid->draw();

	m_display->update(m_orientation, m_isPaused, m_forwardPressed);
}

void FluidRenderer3D::drawCubes(std::vector<glm::vec3> vertices, std::vector<int> indices, std::vector<glm::vec3> darkDots, std::vector<glm::vec3> brightDots) {
	m_display->clear(0.686f, 0.933f, 0.933f, 1.0f);

	m_transform->SetScale(glm::vec3(-1, 1, 1));

	m_colorShader->bind();
	m_colorShader->update(m_transform, m_camera);
	m_colorShader->setColor(0.0f, 0.0f, 0.0f, 1.0f);

	Point activeDots{ darkDots };
	activeDots.draw();

	m_colorShader->setColor(1.0f, 1.0f, 1.0f, 1.0f);

	Point inactiveDots{ brightDots };
	inactiveDots.draw();

	Mesh mesh{ vertices, std::vector<glm::vec3> {}, indices };
	m_colorShader->setColor(1.0f, 0.0f, 0.0f, 0.5f);
	mesh.draw();
	m_colorShader->setColor(1.0f, 0.9f, 0.9f, 1.0f);
	mesh.drawOutline();
	
	m_display->update(m_orientation, m_isPaused, m_forwardPressed);
}

//----------------------------------------------------------------------
// Private Functions
//----------------------------------------------------------------------

/*
Reads in the inital geometry and generates the vertices of the solids
in a vector
Args:
label - pointer reference to the grid containing the geometry info
x - gridwidth
y - gridheight
z - griddepth
*/
void FluidRenderer3D::initGeom(SimUtil::Mat3Di *label, int x, int y, int z) {
	std::vector<glm::vec3> vertSolid;
	std::vector<glm::vec3> normalsSolid;
	std::vector<int> indSolid;
	int maxGridSize;
	int currentInd = 0;
	//get maximum grid size for scaling
	if (x > y) {
		if (x > z)
			maxGridSize = x;
		else
			maxGridSize = z;
	}
	else {
		if (y > z)
			maxGridSize = y;
		else
			maxGridSize = z;
	}
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				//only draw solids
				if (label->get(i, j, k) == SimUtil::SOLID) {
					float x1 = 2.0f * i / (maxGridSize)-1.0f;
					float x2 = 2.0f * (i + 1) / (maxGridSize)-1.0f;
					float y1 = 2.0f * j / (maxGridSize)-1.0f;
					float y2 = 2.0f * (j + 1) / (maxGridSize)-1.0f;
					float z1 = 2.0f * k / (maxGridSize)-1.0f;
					float z2 = 2.0f * (k + 1) / (maxGridSize)-1.0f;
					//left-side
					if (i == 0 || label->get(i - 1, j, k) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x1, y1, z2));
						normalsSolid.push_back(glm::vec3(-1, 0, 0));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x1, y2, z2));
						normalsSolid.push_back(glm::vec3(-1, 0, 0));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x1, y2, z1));
						normalsSolid.push_back(glm::vec3(-1, 0, 0));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x1, y1, z1));
						normalsSolid.push_back(glm::vec3(-1, 0, 0));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
					//right-side
					if (i == x - 1 || label->get(i + 1, j, k) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x2, y1, z2));
						normalsSolid.push_back(glm::vec3(1, 0, 0));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x2, y2, z2));
						normalsSolid.push_back(glm::vec3(1, 0, 0));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x2, y2, z1));
						normalsSolid.push_back(glm::vec3(1, 0, 0));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x2, y1, z1));
						normalsSolid.push_back(glm::vec3(1, 0, 0));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
					//front-side
					if (k == 0 || label->get(i, j, k - 1) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x1, y1, z1));
						normalsSolid.push_back(glm::vec3(0, 0, -1));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x1, y2, z1));
						normalsSolid.push_back(glm::vec3(0, 0, -1));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x2, y2, z1));
						normalsSolid.push_back(glm::vec3(0, 0, -1));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x2, y1, z1));
						normalsSolid.push_back(glm::vec3(0, 0, -1));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
					//back-side
					if (k == z - 1 || label->get(i, j, k + 1) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x1, y1, z2));
						normalsSolid.push_back(glm::vec3(0, 0, 1));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x1, y2, z2));
						normalsSolid.push_back(glm::vec3(0, 0, 1));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x2, y2, z2));
						normalsSolid.push_back(glm::vec3(0, 0, 1));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x2, y1, z2));
						normalsSolid.push_back(glm::vec3(0, 0, 1));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
					//bottom-side
					if (j == 0 || label->get(i, j - 1, k) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x1, y1, z2));
						normalsSolid.push_back(glm::vec3(0, -1, 0));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x1, y1, z1));
						normalsSolid.push_back(glm::vec3(0, -1, 0));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x2, y1, z1));
						normalsSolid.push_back(glm::vec3(0, -1, 0));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x2, y1, z2));
						normalsSolid.push_back(glm::vec3(0, -1, 0));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
					//top-side
					if (j == y - 1 || label->get(i, j + 1, k) != SimUtil::SOLID) {
						//bottom-left
						vertSolid.push_back(glm::vec3(x1, y2, z2));
						normalsSolid.push_back(glm::vec3(0, 1, 0));
						indSolid.push_back(currentInd);
						//top-left
						vertSolid.push_back(glm::vec3(x1, y2, z1));
						normalsSolid.push_back(glm::vec3(0, 1, 0));
						indSolid.push_back(currentInd + 1);
						//top-right
						vertSolid.push_back(glm::vec3(x2, y2, z1));
						normalsSolid.push_back(glm::vec3(0, 1, 0));
						indSolid.push_back(currentInd + 2);
						//bottom-left
						indSolid.push_back(currentInd);
						//top-right
						indSolid.push_back(currentInd + 2);
						//bottom-right
						vertSolid.push_back(glm::vec3(x2, y2, z2));
						normalsSolid.push_back(glm::vec3(0, 1, 0));
						indSolid.push_back(currentInd + 3);

						currentInd += 4;
					}
				}
			}
		}
	}

	initBorderLines(x, y, z, maxGridSize);
	m_meshSolid = new Mesh{ vertSolid , normalsSolid, indSolid };
}

/*
Draws the boundary lines of the geometry
Args:
x - gridwidth
y - gridheight
z - griddepth
maxGridSize - maximum of x,y,z for scaling purposes
*/
void FluidRenderer3D::initBorderLines(int x, int y, int z, int maxGridSize) {
	std::vector<glm::vec3> lines;
	//inner
	bottomLineAt(lines, 3, 3, x, y, z, maxGridSize);
	sideLineAt(lines, 3, 3, 3, 3, x, y, z, maxGridSize);
	sideLineAt(lines, x - 3, 3, 3, 3, x, y, z, maxGridSize);
	frontLineAt(lines, 3, 3, 3, 3, x, y, z, maxGridSize);
	frontLineAt(lines, z - 3, 3, 3, 3, x, y, z, maxGridSize);
	//outer
	/*bottomLineAt(lines, 0, 0, x, y, z, maxGridSize);
	sideLineAt(lines, 0, 0, 0, 0, x, y, z, maxGridSize);
	sideLineAt(lines, x, 0, 0, 0, x, y, z, maxGridSize);
	frontLineAt(lines, 0, 0, 0, 0, x, y, z, maxGridSize);
	frontLineAt(lines, z, 0, 0, 0, x, y, z, maxGridSize);*/

	m_borderSolid = new Line{ lines };
}

/*
Draws the bottom line
Args:
lines - reference where line vertices should be stored in
yLoc - y coordinate for the line (line distance from the bottom)
crop - left and right cropping of the line
x - gridwidth
y - gridheight
z - griddepth
maxGridSize - maximum of x,y,z for scaling purposes
*/
void FluidRenderer3D::bottomLineAt(std::vector<glm::vec3> &lines, int yLoc, int crop, int x, int y, int z, int maxGridSize) {
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * yLoc / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
}

/*
Draws a side line
Args:
lines - reference where line vertices should be stored in
xLoc - x coordinate for the line
crop - front and back cropping of the line
bottomCrop - bottom crop of the line
topCrop - top crop of the line
x - gridwidth
y - gridheight
z - griddepth
maxGridSize - maximum of x,y,z for scaling purposes
*/
void FluidRenderer3D::sideLineAt(std::vector<glm::vec3> &lines, int xLoc, int crop, int bottomCrop, int topCrop, int x, int y, int z, int maxGridSize) {
	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * (z - crop) / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * xLoc / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * crop / maxGridSize - 1.0));
}

/*
Draws a side line
Args:
lines - reference where line vertices should be stored in
zLoc - z coordinate for the line
crop - front and back cropping of the line
bottomCrop - bottom crop of the line
topCrop - top crop of the line
x - gridwidth
y - gridheight
z - griddepth
maxGridSize - maximum of x,y,z for scaling purposes
*/
void FluidRenderer3D::frontLineAt(std::vector<glm::vec3> &lines, int zLoc, int crop, int bottomCrop, int topCrop, int x, int y, int z, int maxGridSize) {
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * (y - topCrop) / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));

	lines.push_back(glm::vec3(2.0 * (x - crop) / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));
	lines.push_back(glm::vec3(2.0 * crop / maxGridSize - 1.0, 2.0 * bottomCrop / maxGridSize - 1.0, 2.0 * zLoc / maxGridSize - 1.0));
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

	std::string filename = "pic";
	filename.append("-frame-");
	filename.append(std::to_string(frame));
	filename.append(".bmp");


	SDL_SaveBMP(temp, filename.c_str());
}