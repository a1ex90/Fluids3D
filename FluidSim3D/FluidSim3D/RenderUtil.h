#include <string>
#include <vector>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <glm/glm.hpp>

//----------------------------------------------------------------------
// Display Class
//----------------------------------------------------------------------
class Display
{
public:
	Display(int width, int height, const std::string& title);

	void update();
	void clear(float r, float g, float b, float a);
	bool isClosed();

	~Display();
private:
	Display(const Display& other) {}
	void operator=(const Display& other) {}

	SDL_Window* m_window;
	SDL_GLContext m_glContext;
	bool m_isClosed;
};

//----------------------------------------------------------------------
// Mesh Class
//----------------------------------------------------------------------
class Mesh {
public:
	Mesh(std::vector<glm::vec2> vertices, std::vector<int> indices);
	Mesh(std::vector<glm::vec2> vertices, std::vector<glm::vec2> textureCoords, std::vector<int> indices);
	Mesh(std::vector<glm::vec2> vertices, std::vector<int> indices, std::vector<float> opacites);

	void draw();

	~Mesh();
private:
	enum {
		POSITION_VB,
		TEXCOORD_VB,
		OPACITY_VB,	
		INDEX_VB,
		NUM_BUFFERS
	};

	GLuint m_vertexArrayObject;
	GLuint m_vertexArrayBuffers[NUM_BUFFERS];
	unsigned int m_drawCount;
};

//----------------------------------------------------------------------
// Line Class
//----------------------------------------------------------------------
class Line {
public:
	Line(std::vector<glm::vec2> vertices);
	
	void draw();

	~Line();
private:
	enum {
		POSITION_VB,
		NUM_BUFFERS
	};
	GLuint m_vertexArrayObject;
	GLuint m_vertexArrayBuffers[NUM_BUFFERS];
	unsigned int m_verticesCount;
};

//----------------------------------------------------------------------
// Point Class
//----------------------------------------------------------------------
class Point {
public:
	Point(std::vector<glm::vec2> points);

	void draw();

	~Point();
private:
	enum {
		POSITION_VB,
		NUM_BUFFERS
	};
	GLuint m_vertexArrayObject;
	GLuint m_vertexArrayBuffers[NUM_BUFFERS];
	unsigned int m_pointsCount;
};

//----------------------------------------------------------------------
// Shader Class
//----------------------------------------------------------------------
class Shader {
public:
	Shader(const std::string& fileName);

	void bind();
	void setColor(float r, float g, float b);
	void setTexture(int unit);
	
	~Shader();
private:
	static const unsigned int NUM_SHADERS = 2;
	static const unsigned int NUM_UNIFORMS = 2;

	GLuint m_program;
	GLuint m_shaders[NUM_SHADERS];
	GLuint m_uniforms[NUM_UNIFORMS];	
};

//----------------------------------------------------------------------
// Texture Class
//----------------------------------------------------------------------
class Texture
{
public:
	Texture(const std::string& fileName);

	void bind(unsigned int unit);

	~Texture();
private:
	Texture(const Texture& other) {};
	void operator=(const Texture& other) {};

	GLuint m_texture;
};