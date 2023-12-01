#shader vertex

#version 330 core

uniform vec2 window_resolution;
uniform float max_arrow_size;
uniform float arrow_width;

layout(location=0) in vec2 pos;
layout(location=1) in vec2 vel;

out vec3 flow_color;

float sigmoid(float x) {return 1 / (1 + exp(-x));}
vec2 sigmoid(vec2 a) {return 1 / (1 + exp(-a));}

vec3 lerp(vec3 c0, vec3 cf, float t) {return c0 * (1 - t) + cf * t;}



const float MAX_VEL = 100;

#define RED_PASTEL vec3(0.9882352941176471, 0.43529411764705883, 0.43529411764705883)
#define GREEN_PASTEL vec3(0.7450980392156863, 1.0, 0.6078431372549019)
#define BLUE_PASTEL vec3(0.3333333333333333, 0.40784313725490196, 0.9686274509803922)

#define COLORING_VELOCITY_COEFICIENT 40

void main() {
    vec2 u_displacement = sigmoid(vel/MAX_VEL)*2 - 1;
    u_displacement = (((gl_VertexID >> 1) & 1) * 2 - 1) * u_displacement / 2;
    vec2 world_pos = pos + ((gl_VertexID >> 1) & 1) * 2 * u_displacement * max_arrow_size / 2 + ((gl_VertexID & 1) * 2 - 1) * (((gl_VertexID >> 1) & 1) * 2 - 1) * vec2(- u_displacement.y, u_displacement.x) * arrow_width;

    vec2 screen_pos = (world_pos) / window_resolution * 2;
   
    gl_Position = vec4(screen_pos, 0, 1);

    float t = sigmoid(length(vel) / COLORING_VELOCITY_COEFICIENT - 2.5);
    if (t <= 0.5) {
        flow_color = lerp(BLUE_PASTEL, GREEN_PASTEL, 2*t);
    } else {
        flow_color = lerp(GREEN_PASTEL, RED_PASTEL, 2*(t - 0.5));
    }
    //vec2 screen_pos = (pos + arrow_width/2*(2*uv - vec2(1, 1))) / window_resolution * 2;
    //gl_Position = vec4(screen_pos, 0, 1);
   
}

#shader fragment

#version 330 core

// uniform vec4 flow_color;

in vec3 flow_color;
out vec4 color;

void main()
{
    color = vec4(flow_color, 1);
}
