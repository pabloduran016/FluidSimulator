#shader vertex

#version 330 core

layout (location = 0) in vec2 pos; 

uniform vec2 window_resolution;
uniform float particle_mass;
uniform float smoothing_radius;
uniform int n_particles;

uniform sampler1D positions_encoded;

const int MAX_NUMBER_OF_PARTICLES = 3000;

// uniform vec2 positions[1000];

out float density_value;

const float PI = 3.141593;

float calculate_density(vec2 pos)
{
    float density = 0;
    for (int i = 0; i < MAX_NUMBER_OF_PARTICLES; ++i) {
        if (i >= n_particles) break;
        // vec2 pos_i = positions[i];
        vec2 pos_i = texelFetch(positions_encoded, i, 0).rg;
        float length_squared = dot(pos_i - pos, pos_i - pos);
        if (length_squared >= smoothing_radius * smoothing_radius) continue;
        float length = sqrt(length_squared);
        density += pow((smoothing_radius - length), 2);
    }
    return particle_mass * density * 6 / (PI * pow(smoothing_radius, 4));  // SpikyKernelPow2
}

void main()
{
    density_value = calculate_density(pos);
    vec2 screen_pos = pos / window_resolution * 2;
    gl_Position = vec4(screen_pos, 0, 1);
}

#shader fragment

#version 330 core

in float density_value;

uniform float target_density;

out vec4 color;

#define DENSITY_COLOR0 0, 0.2, 1, 0.8
#define DENSITY_COLOR1 1, 1, 1, 0.8
#define DENSITY_COLOR2 0, 0, 0, 0.8

#define DRAW_DENSITY_EXPONENT 0.5

void main()
{    
    float t = (2/(1 + exp(1 - pow(density_value/target_density, DRAW_DENSITY_EXPONENT))));  // in (0, 2)
    if (t > 1) {
        color = vec4(DENSITY_COLOR0) * (t - 1) + vec4(DENSITY_COLOR1) * (2 - t);
    } else {
        color = vec4(DENSITY_COLOR1) * t + vec4(DENSITY_COLOR2) * (1 - t);
    }
}