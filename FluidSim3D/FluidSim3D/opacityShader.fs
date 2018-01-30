#version 120

varying float opacity0;

uniform vec3 Color;

void main() {
	float opaMax = 0.95;
	float opaMin = 0.5;

	//gl_FragColor = vec4(Color, mix(opaMin, opaMax, opacity0));

	vec3 color2 = vec3(0.392, 0.584, 0.929);
	gl_FragColor = vec4(mix(color2, Color, opacity0), mix(opaMin, opaMax, opacity0));
}