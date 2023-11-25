#shader vertex

#version 330 core

uniform float particle_radius;
uniform vec2 window_resolution;
uniform bool draw_density;

layout(location=0) in vec2 pos;
layout(location=1) in vec2 vel;

out vec3 vertex_color;
out vec2 uv;


vec3 lerp(vec3 c0, vec3 cf, float t)
{
    return c0 * (1 - t) + cf * t;
}


float sigmoid(float x)
{
    return 1 / (1 + exp(-x));
}


#define RED_PASTEL vec3(0.9882352941176471, 0.43529411764705883, 0.43529411764705883)
#define GREEN_PASTEL vec3(0.7450980392156863, 1.0, 0.6078431372549019)
#define BLUE_PASTEL vec3(0.3333333333333333, 0.40784313725490196, 0.9686274509803922)

#define COLORING_VELOCITY_COEFICIENT 40
#define DRAW_DENSITY_COLOR 0, 0, 0

void main() {
    uv = vec2(float(gl_VertexID & 1), float((gl_VertexID >> 1) & 1));
    vec2 screen_pos = (pos + particle_radius*(2*uv - vec2(1, 1))) / window_resolution * 2;
    gl_Position = vec4(screen_pos, 0, 1);
    float t = sigmoid(length(vel) / COLORING_VELOCITY_COEFICIENT - 2.5);
    if (draw_density) {
        vertex_color = vec3(DRAW_DENSITY_COLOR);  
    } else {
        if (t <= 0.5) {
            vertex_color = lerp(BLUE_PASTEL, GREEN_PASTEL, 2*t);
        } else {
            vertex_color = lerp(GREEN_PASTEL, RED_PASTEL, 2*(t - 0.5));
        }
    }
}

#shader fragment

#version 330 core

uniform bool draw_density;

in vec3 vertex_color;
in vec2 uv;

out vec4 color;

void main()
{
    if (length(uv - vec2(0.5, 0.5)) > 0.5) {    
        color = vec4(0, 0, 0, 0);
        return;
    } 
    color = vec4(vertex_color.rgb, 0.5);
}
