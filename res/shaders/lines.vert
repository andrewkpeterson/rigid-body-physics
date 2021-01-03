#version 330 core

uniform mat4 m;
uniform mat4 v;
uniform mat4 p;

layout(location = 0) in vec3 position;

void main(void)
{
    gl_Position = p * v * m * vec4(position,1);
}
