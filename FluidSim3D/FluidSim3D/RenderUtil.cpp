#include "RenderUtil.h"


//----------------------------------------------------------------------
// Display Class
//----------------------------------------------------------------------
/*
creates an OpenGL Window with given attributes using cross-platform OpenGL Loader SDL
---
Arguments:
width - Window width
height - Window height
title - Window title

*/
Display::Display(int width, int height, const std::string& title, Transform* transform)
{
	SDL_Init(SDL_INIT_EVERYTHING);

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE, 32);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	m_window = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_OPENGL);
	m_glContext = SDL_GL_CreateContext(m_window);

	GLenum status = glewInit();

	if (status != GLEW_OK) {
		std::cerr << "Glew failed to initialize" << std::endl;
	}

	m_transform = transform;

	m_isClosed = false;
}


Display::~Display()
{
	SDL_GL_DeleteContext(m_glContext);
	SDL_DestroyWindow(m_window);
	SDL_Quit();
}

/*
Fills the window with a given solid color
---
Arguments:
r,g,b,a - rot, green, blue and alpha values that define the color
values between 0 and 1
*/
void Display::clear(float r, float g, float b, float a) {
	glClearColor(r, g, b, a);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

/*
Returns true if the window was closed false else
*/
bool Display::isClosed() {
	return m_isClosed;
}

/*
Updates the window
*/
void Display::update(glm::vec3 &orientation, bool &pausePressed, bool &forwardPressed, int &visualMode) {
	SDL_GL_SwapWindow(m_window);
	SDL_Event e;
	int handled;

	float PI = 3.1415927410125732421875f;
	//Defines how much the scene should be rotated by one keypress
	int increments = 16;
	//Defines how much the scene should be moved by one keypress
	float steps = 0.2f;

	while (SDL_PollEvent(&e)) {
		if (e.type == SDL_QUIT) {
			m_isClosed = true;
		}
		switch (e.type)
		{
		case SDL_KEYDOWN:
			switch (e.key.keysym.sym)
			{
			case SDLK_1:
				visualMode = 1;
				break;
			case SDLK_2:
				visualMode = 2;
				break;
			case SDLK_3:
				visualMode = 3;
				break;
			case SDLK_4:
				visualMode = 4;
				break;
			case SDLK_w:
				//move screen up
				m_transform->SetPos(m_transform->GetPos() + glm::vec3(0, steps, 0));
				break;
			case SDLK_s:
				//move screen down
				m_transform->SetPos(m_transform->GetPos() - glm::vec3(0, steps, 0));
				break;
			case SDLK_d:
				//move screen right
				m_transform->SetPos(m_transform->GetPos() - glm::vec3(steps, 0, 0));
				break;
			case SDLK_a:
				//move screen left
				m_transform->SetPos(m_transform->GetPos() + glm::vec3(steps, 0, 0));
				break;
			case SDLK_q:
				//move closer
				m_transform->SetPos(m_transform->GetPos() + glm::vec3(0, 0, steps));
				break;
			case SDLK_e:
				//move screen left
				m_transform->SetPos(m_transform->GetPos() - glm::vec3(0, 0, steps));
				break;
			case SDLK_r:
				//move screen left
				m_transform->SetPos(glm::vec3(0, 0, 0));
				m_transform->SetRot(glm::vec3(0, 0, 0));
				orientation.x = 0.0f;
				orientation.y = -1.0f;
				orientation.z = 0.0f;
				break;
			case SDLK_ESCAPE:
				exit(0);
				break;
			case SDLK_RIGHT: {
				m_transform->SetRot((m_transform->GetRot() + glm::vec3(0, 2 * PI / increments, 0)));
				float z = orientation.z * cos(2 * PI / increments) + orientation.x * sin(2 * PI / increments);
				float x = -orientation.z * sin(2 * PI / increments) + orientation.x * cos(2 * PI / increments);
				orientation.z = z;
				orientation.x = x;
				break; }
			case SDLK_LEFT: {
				m_transform->SetRot((m_transform->GetRot() - glm::vec3(0, 2 * PI / increments, 0)));
				float z = orientation.z * cos(2 * PI / increments) - orientation.x * sin(2 * PI / increments);
				float x = orientation.z * sin(2 * PI / increments) + orientation.x * cos(2 * PI / increments);
				orientation.z = z;
				orientation.x = x;
				break; }
			case SDLK_DOWN: {
				m_transform->SetRot((m_transform->GetRot() - glm::vec3(2 * PI / increments, 0, 0)));
				float y = orientation.y * cos(2 * PI / increments) - orientation.z * sin(2 * PI / increments);
				float z = orientation.y * sin(2 * PI / increments) + orientation.z * cos(2 * PI / increments);
				orientation.y = y;
				orientation.z = z;
				break; }
			case SDLK_UP: {
				m_transform->SetRot((m_transform->GetRot() + glm::vec3(2 * PI / increments, 0, 0)));			
				float y = orientation.y * cos(2 * PI / increments) + orientation.z * sin(2 * PI / increments);
				float z = -orientation.y * sin(2 * PI / increments) + orientation.z * cos(2 * PI / increments);
				orientation.y = y;
				orientation.z = z;
				break; }
			case SDLK_o: {
				m_transform->SetRot((m_transform->GetRot() + glm::vec3(0, 0, 2 * PI / increments)));
				float x = orientation.x * cos(2 * PI / increments) + orientation.y * sin(2 * PI / increments);
				float y = -orientation.x * sin(2 * PI / increments) + orientation.y * cos(2 * PI / increments);
				orientation.x = x;
				orientation.y = y;
				break; }
			case SDLK_l: {
				m_transform->SetRot((m_transform->GetRot() - glm::vec3(0, 0, 2 * PI / increments)));
				float x = orientation.x * cos(2 * PI / increments) - orientation.y * sin(2 * PI / increments);
				float y = orientation.x * sin(2 * PI / increments) + orientation.y * cos(2 * PI / increments);
				orientation.x = x;
				orientation.y = y;
				break; }
			case SDLK_SPACE:
				pausePressed = !pausePressed;
				break;
			case SDLK_m:
				forwardPressed = true;
				pausePressed = true;
				break;
			case SDLK_n:
				orientation.x = 0.0f;
				orientation.y = 0.0f;
				orientation.z = 1.0f;
				break;
			}
		}
	}
}

//----------------------------------------------------------------------
// Mesh Class
//----------------------------------------------------------------------

Mesh::Mesh(std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<int> indices) {
	m_drawCount = indices.size();

	glGenVertexArrays(1, &m_vertexArrayObject);
	glBindVertexArray(m_vertexArrayObject);

	glGenBuffers(NUM_BUFFERS, m_vertexArrayBuffers);

	glBindBuffer(GL_ARRAY_BUFFER, m_vertexArrayBuffers[POSITION_VB]);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0]), vertices.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_vertexArrayBuffers[NORMAL_VB]);
	glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(normals[0]), normals.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_vertexArrayBuffers[INDEX_VB]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(indices[0]), indices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);

	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	//needs to be enabled to active opacity rendering
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

Mesh::~Mesh() {
	glDeleteBuffers(NUM_BUFFERS, m_vertexArrayBuffers);
	glDeleteVertexArrays(1, &m_vertexArrayObject);
}

void Mesh::draw() {
	glBindVertexArray(m_vertexArrayObject);
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElementsBaseVertex(GL_TRIANGLES, m_drawCount, GL_UNSIGNED_INT, 0, 0);

	glBindVertexArray(0);
}

void Mesh::drawOutline() {
	glBindVertexArray(m_vertexArrayObject);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElementsBaseVertex(GL_TRIANGLES, m_drawCount, GL_UNSIGNED_INT, 0, 0);

	glBindVertexArray(0);
}

//----------------------------------------------------------------------
// Line Class
//----------------------------------------------------------------------
/*
	Constructor
	---
	Arguments:
	vertices - array with 2D vertice data. for vertex 1 per example vertices[0] = x-coord and vertices[1] = y-coord
			   so vertices is twice as long as vertices exist
	verticesCount - number of vertices

*/
Line::Line(std::vector<glm::vec3> vertices) {
	m_verticesCount = vertices.size();

	glGenVertexArrays(1, &m_vertexArrayObject);
	glBindVertexArray(m_vertexArrayObject);

	glGenBuffers(NUM_BUFFERS, m_vertexArrayBuffers);

	glBindBuffer(GL_ARRAY_BUFFER, m_vertexArrayBuffers[POSITION_VB]);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0]), vertices.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindVertexArray(0);
}

Line::~Line() {
	glDeleteBuffers(NUM_BUFFERS, m_vertexArrayBuffers);
	glDeleteVertexArrays(1, &m_vertexArrayObject);
}

void Line::draw() {
	glBindVertexArray(m_vertexArrayObject);

	glEnable(GL_LINE_SMOOTH);
	glLineWidth(1);
	glDrawArrays(GL_LINES, 0, m_verticesCount);

	glBindVertexArray(0);

}

//----------------------------------------------------------------------
// Point Class
//----------------------------------------------------------------------
/*
Constructor
---
Arguments:
vertices - array with 2D vertice data. for vertex 1 per example vertices[0] = x-coord and vertices[1] = y-coord
so vertices is twice as long as vertices exist
verticesCount - number of vertices

*/
Point::Point(std::vector<glm::vec3> points) {
	m_pointsCount = points.size();

	glGenVertexArrays(1, &m_vertexArrayObject);
	glBindVertexArray(m_vertexArrayObject);

	glGenBuffers(NUM_BUFFERS, m_vertexArrayBuffers);

	glBindBuffer(GL_ARRAY_BUFFER, m_vertexArrayBuffers[POSITION_VB]);
	glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(points[0]), points.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindVertexArray(0);
}

Point::~Point() {
	glDeleteBuffers(NUM_BUFFERS, m_vertexArrayBuffers);
	glDeleteVertexArrays(1, &m_vertexArrayObject);
}

void Point::draw() {
	glBindVertexArray(m_vertexArrayObject);

	glEnable(GL_PROGRAM_POINT_SIZE);
	//DISABLE MEE
	//glPointSize(5.0f);
	glEnable(GL_POINT_SMOOTH);
	glDrawArrays(GL_POINTS, 0, m_pointsCount);

	glBindVertexArray(0);

}

//----------------------------------------------------------------------
// Shader Class
//----------------------------------------------------------------------
void checkShaderError(GLuint shader, GLuint flag, bool isProgram, const std::string& errorMessage);
std::string loadShader(const std::string& fileName);
GLuint createShader(const std::string& text, GLenum shaderType);

Shader::Shader(const std::string& fileName) {
	m_program = glCreateProgram();
	m_shaders[0] = createShader(loadShader(fileName + ".vs"), GL_VERTEX_SHADER);
	m_shaders[1] = createShader(loadShader(fileName + ".fs"), GL_FRAGMENT_SHADER);

	for (int i = 0; i < NUM_SHADERS; i++) {
		glAttachShader(m_program, m_shaders[i]);
	}

	glBindAttribLocation(m_program, 0, "position");
	glBindAttribLocation(m_program, 1, "texCoord");
	glBindAttribLocation(m_program, 2, "normal");

	glLinkProgram(m_program);
	glValidateProgram(m_program);

	//m_uniforms[0] = glGetUniformLocation(m_program, "Color");
	//needed for texture
	//m_uniforms[1] = glGetUniformLocation(m_program, "Diffuse");
	m_uniforms[TRANSFORM_U] = glGetUniformLocation(m_program, "transform");
	m_uniforms[COLOR] = glGetUniformLocation(m_program, "color");
}

Shader::~Shader() {
	for (int i = 0; i < NUM_SHADERS; i++) {
		glDetachShader(m_program, m_shaders[i]);
		glDeleteShader(m_shaders[i]);
	}
	glDeleteProgram(m_program);
}

void Shader::bind() {
	glUseProgram(m_program);
}

void Shader::update(const Transform* transform, const Camera* camera) {
	glm::mat4 model = camera->getViewProjection() * transform->GetModel();

	glUniformMatrix4fv(m_uniforms[TRANSFORM_U], 1, GL_FALSE, &model[0][0]);
}

void Shader::setColor(float r, float g, float b, float a) {
	glUniform4f(m_uniforms[COLOR], r, g, b, a);
}

std::string loadShader(const std::string& fileName)
{
	std::ifstream file;
	file.open((fileName).c_str());

	std::string output;
	std::string line;

	if (file.is_open())
	{
		while (file.good())
		{
			getline(file, line);
			output.append(line + "\n");
		}
	}
	else
	{
		std::cerr << "Unable to load shader: " << fileName << std::endl;
	}

	return output;
}

GLuint createShader(const std::string& text, GLenum shaderType) {
	GLuint shader = glCreateShader(shaderType);

	if (shader == 0) {
		std::cerr << "Error: Shader creation failed" << std::endl;
	}

	const GLchar* shaderSourceStrings[1];
	GLint shaderSourceStringLengths[1];
	shaderSourceStrings[0] = text.c_str();
	shaderSourceStringLengths[0] = text.length();

	glShaderSource(shader, 1, shaderSourceStrings, shaderSourceStringLengths);
	glCompileShader(shader);

	checkShaderError(shader, GL_COMPILE_STATUS, true, "Error: Program is compilation failed: ");

	return shader;
}

void checkShaderError(GLuint shader, GLuint flag, bool isProgram, const std::string& errorMessage)
{
	GLint success = 0;
	GLchar error[1024] = { 0 };

	if (isProgram)
		glGetProgramiv(shader, flag, &success);
	else
		glGetShaderiv(shader, flag, &success);

	if (success == GL_FALSE)
	{
		if (isProgram)
			glGetProgramInfoLog(shader, sizeof(error), NULL, error);
		else
			glGetShaderInfoLog(shader, sizeof(error), NULL, error);

		std::cerr << errorMessage << ": '" << error << "'" << std::endl;
	}
}