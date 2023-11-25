#shader vertex

#version 330 core

uniform float interaction_radius;
uniform vec2 window_resolution;

layout(location=0) in vec2 pos;

out vec3 vertex_color;
out vec2 uv;

void main() {
    uv = vec2(float(gl_VertexID & 1), float((gl_VertexID >> 1) & 1));
    vec2 screen_pos = (pos + interaction_radius*(2*uv - vec2(1, 1))) / window_resolution * 2;
    gl_Position = vec4(screen_pos, 0, 1);
    vertex_color = vec3(1, 1, 1);
}

#shader fragment

#version 330 core

in vec3 vertex_color;
in vec2 uv;

out vec4 color;

void main()
{
    vec2 r = uv - vec2(0.5, 0.5);
    float radius = length(r);
    if (radius > 0.5) {    
        color = vec4(0, 0, 0, 0);
        return;
    } 
    color = vec4(vertex_color.rgb, pow(2*radius - 0.5, 2));
}
