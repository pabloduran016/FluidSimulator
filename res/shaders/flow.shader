#shader vertex

#version 330 core

uniform vec2 window_resolution;
uniform float max_arrow_size;
uniform float arrow_width;

layout(location=0) in vec2 pos;
layout(location=1) in vec2 vel;


float sigmoid(float x) {return 1 / (1 + exp(-x));}
vec2 sigmoid(vec2 a) {return 1 / (1 + exp(-a));}


const float MAX_VEL = 100;


void main() {

    
    vec2 u_displacement = sigmoid(vel/MAX_VEL)*2 - 1;
    u_displacement = (((gl_VertexID >> 1) & 1) * 2 - 1) * u_displacement / 2;
    vec2 world_pos = pos + ((gl_VertexID >> 1) & 1) * 2 * u_displacement * max_arrow_size / 2 + ((gl_VertexID & 1) * 2 - 1) * (((gl_VertexID >> 1) & 1) * 2 - 1) * vec2(- u_displacement.y, u_displacement.x) * arrow_width;

    vec2 screen_pos = (world_pos) / window_resolution * 2;
   
    gl_Position = vec4(screen_pos, 0, 1);
    //vec2 screen_pos = (pos + arrow_width/2*(2*uv - vec2(1, 1))) / window_resolution * 2;
    //gl_Position = vec4(screen_pos, 0, 1);
   
}

#shader fragment

#version 330 core

uniform vec4 flow_color;

out vec4 color;

void main()
{
    color = flow_color;
}
