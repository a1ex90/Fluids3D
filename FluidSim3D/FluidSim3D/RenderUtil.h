#include <string>
#include <vector>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <iostream>
#include <fstream>
#include <cassert>
#define GLM_ENABLE_EXPERIMENTAL 
#include <glm/gtx/transform.hpp>
#include <glm/gtc/quaternion.hpp>


//----------------------------------------------------------------------
//  Light Structure
//----------------------------------------------------------------------
struct Light {
	glm::vec3 position;
	glm::vec3 intensities; //a.k.a. the color of the light
	float attenuation;
	float ambientCoefficient;
};

//----------------------------------------------------------------------
// Transform Class
//----------------------------------------------------------------------
class Transform
{
public:
	Transform(const glm::vec3& pos = glm::vec3(), const glm::vec3& rot = glm::vec3(), const glm::vec3& scale = glm::vec3(1.0f, 1.0f, 1.0f)) :
		m_pos(pos),
		m_rot(rot),
		m_scale(scale) {}

	inline glm::mat4 GetModel() const {
		glm::mat4 posMatrix = glm::translate(m_pos);
		glm::mat4 scaleMatrix = glm::scale(m_scale);
		glm::mat4 rotMatrix = glm::mat4_cast(m_rot);

		return posMatrix * rotMatrix * scaleMatrix;
	}

	inline glm::vec3& GetPos() { return m_pos; }
	inline glm::quat& GetRot() { return m_rot; }
	inline glm::vec3& GetScale() { return m_scale; }

	inline void SetPos(const glm::vec3& pos) { m_pos = pos; }
	inline void SetRot(const glm::quat& rot) { m_rot = rot; }
	inline void SetScale(const glm::vec3& scale) { m_scale = scale; }
private:
	glm::vec3 m_pos;
	glm::quat m_rot;
	glm::vec3 m_scale;
};

//----------------------------------------------------------------------
// Display Class
//----------------------------------------------------------------------
class Display
{
public:
	Display(int width, int height, const std::string& title, Transform* transform);

	void update(glm::vec3 &orientation, bool &pausePressed, bool &forwardPressed, int &visualMode, bool &manipulation);
	void clear(float r, float g, float b, float a);
	bool isClosed();

	~Display();
private:
	Display(const Display& other) {}
	void operator=(const Display& other) {}

	const float PI = 3.1415927410125732421875f;
	//Defines the size of the arcball for rotation
	float m_r;
	//Stores the startpoint of an arcball rotation
	glm::vec3 m_startP;
	//Stores the axis for the arcball rotation
	glm::vec3 m_axis;
	//stores the angle of the arcball rotation around axis
	float m_angle;
	//Stores if mouse is clicked
	bool m_doRotation;

	glm::vec3 projectOnSphere(float x, float y, float r);
	glm::vec3 updateOrientation(glm::quat rotation);

	Transform* m_transform;

	SDL_Window* m_window;
	SDL_GLContext m_glContext;
	bool m_isClosed;

	int m_width;
	int m_height;
};

//----------------------------------------------------------------------
// Mesh Class
//----------------------------------------------------------------------
class Mesh {
public:
	Mesh::Mesh(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<int> indices);

	void draw();
	void drawOutline();

	~Mesh();
private:
	enum {
		POSITION_VB,
		TEXCOORD_VB,
		NORMAL_VB,	
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
	Line(std::vector<glm::vec3> vertices);
	
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
	Point(std::vector<glm::vec3> points);

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
// Camera Class
//----------------------------------------------------------------------
class Camera
{
public:
	Camera(const glm::vec3& pos, float fov, float aspect, float zNear, float zFar) {
		m_perspective = glm::perspective(fov, aspect, zNear, zFar);
		m_position = pos;
		m_forward = glm::vec3(0, 0, 1);
		m_up = glm::vec3(0, 1, 0);

	}

	inline glm::mat4 getViewProjection() const {
		return m_perspective * glm::lookAt(m_position, m_position + m_forward, m_up);
	}

	inline glm::vec3 getPos() const {
		return m_position;
	}
private:
	glm::mat4 m_perspective;
	glm::vec3 m_position;
	//for camera roation
	glm::vec3 m_forward;
	glm::vec3 m_up;
};

//----------------------------------------------------------------------
// Shader Class
//----------------------------------------------------------------------
class Shader {
public:
	Shader(const std::string& fileName);

	void bind();
	void update(const Transform* transform, const Camera* camera);
	void setLight(const Light& light);
	void setMaterialSettings(float shininess, glm::vec3 specularColor);
	void setColor(float r, float g, float b, float a);
	
	~Shader();
private:
	static const unsigned int NUM_SHADERS = 2;

	enum {
		MODEL_U,
		CAMERA_U,
		CAMERA_POS_U,
		COLOR,
		SHININESS,
		SPECULARCOLOR,

		NUM_UNIFORMS
	};

	GLuint m_program;
	GLuint m_shaders[NUM_SHADERS];
	GLuint m_uniforms[NUM_UNIFORMS];	
};

