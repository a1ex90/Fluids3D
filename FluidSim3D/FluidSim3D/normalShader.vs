#version 150

attribute vec3 position;
attribute vec3 normal;

out vec3 fragNormal;
out vec3 fragVert;

uniform mat4 model;
uniform mat4 camera; 

void main() {
	fragNormal = normal;
	fragVert = position;
	gl_Position = camera * model * vec4(position, 1.0);
}