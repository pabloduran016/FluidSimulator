#shader vertex

#version 330 core

layout (location = 0) in vec2 pos; 
layout (location = 1) in vec2 uv; 

out vec2 frag_uv;

uniform vec2 window_resolution;

void main()
{
    frag_uv = uv;
    vec2 screen_pos = pos / window_resolution * 2;
    gl_Position = vec4(screen_pos, 0, 1);
}

#shader fragment

#version 330 core

in vec2 frag_uv;
out vec4 color;

uniform sampler2D text;
uniform vec4 text_color;

void main()
{    
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, frag_uv).r);
    color = text_color * sampled;
}