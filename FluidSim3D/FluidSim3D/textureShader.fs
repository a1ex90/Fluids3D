#version 120

varying vec2 texCoord0;

uniform sampler2D Diffuse;

void main() {
	gl_FragColor = texture2D(Diffuse, texCoord0);
}