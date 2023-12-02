#define GLAD_GL_IMPLEMENTATION
#include <GL/glew.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <ft2build.h>
#include FT_FREETYPE_H
 
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>


// BEGIN TIMING
#define SHOULD_TIME 0
#if SHOULD_TIME 
    #define TIME_THIS_INIT(name) \
        long int _total_time##name = 0;
    #define TIME_THIS_BEGIN(name) \
        struct timespec _before_##name, _after_##name; \
        timespec_get(&_before_##name, TIME_UTC);
    #define TIME_THIS_END(name) \
        timespec_get(&_after_##name, TIME_UTC); \
        long int _time_dif##name = time_difference_ns(_before_##name, _after_##name); \
        _total_time##name += _time_dif##name; \
        printf("took %ld ns (%.04f s) to run code "#name"\n", _time_dif##name, (float) _time_dif##name*1e-9);
    #define TIME_THIS_RECAP(name) \
        printf("total time: took %ld ns (%.04f s) to run code "#name"\n", _total_time##name, (float) _total_time##name*1e-9);
#else
    #define TIME_THIS_INIT(name) ;
    #define TIME_THIS_BEGIN(name) ;
    #define TIME_THIS_END(name) ;
    #define TIME_THIS_RECAP(name) ;
#endif
// END TIMING


// BEGIN GL ERRORS
#ifdef NOT_GL_LOG_CALL
    #define GLCall(x) x
#else
    #define GLCall(x) GLClearError();x;GLLogCall(__LINE__, __FILE__)
#endif

static void GLClearError(void) {
    while (glGetError() != GL_NO_ERROR) {};    
}

static int GLLogCall(unsigned int line, const char* file) {
    unsigned int error;
    while ((error = glGetError())) {
        fprintf(stderr, "Open GL Error: %s:%u: %x\n", file, line, error);
        return 0;
    }
    return 1;
} 

void debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam)
{
    (void) source;
    (void) id;
    (void) length;
    (void) userParam;
    fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
            (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
            type, severity, message);
}
 
static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}
// END GL ERRORS


#define randnum(min, max) (float)rand()/(float)(RAND_MAX) * (max - min) + min


// BEGIN LINEAR ALGEBRA
typedef struct {
    int x, y;
} iVec2;

static inline iVec2 ivec2(int x, int y) {
    return (iVec2) {.x = x, .y = y};
}

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    float x, y, z, w;
} Vec4;

typedef struct {
    float x, y;
} Vec2;

static inline Vec2 vec2(float x, float y) {
    return (Vec2) {.x = x, .y = y};
}

static inline Vec2 vec2_add(Vec2 a, Vec2 b) {
    return (Vec2) {.x = a.x + b.x, .y = a.y + b.y};
}
static inline Vec2 vec2_sub(Vec2 a, Vec2 b) {
    return (Vec2) {.x = a.x - b.x, .y = a.y - b.y};
}
static inline Vec2 vec2_mul(Vec2 a, Vec2 b) {
    return (Vec2) {.x = a.x * b.x, .y = a.y * b.y};
}
static inline Vec2 vec2_scale(Vec2 a, float b) {
    return (Vec2) {.x = a.x * b, .y = a.y * b};
}

static inline float vec2_length(Vec2 a) {
    return sqrt(a.x * a.x + a.y * a.y);
}

static inline float vec2_squared_length(Vec2 a) {
    return a.x * a.x + a.y * a.y;
}

static inline Vec2 vec2_normalize(Vec2 a) {
    return vec2_scale(a, 1/vec2_length(a));
}
// END LINEAR ALGEBRA
 

// BEGIN OpenGL
#define MAX_SHADER_LENGTH 10000

static unsigned int compile_shader(const char* content, unsigned int type) {
    unsigned int shader =  glCreateShader(type);
    glShaderSource(shader, 1, &content, NULL);
    glCompileShader(shader);
    int result;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE) {
        int length;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        char* message = (char*) alloca(length * sizeof (char));
        glGetShaderInfoLog(shader, length, &length, message);
        fprintf(stderr, "Failed to compile %s shader!\n", (type == GL_VERTEX_SHADER ? "vertex" : "fragment"));
        fprintf(stderr, "%s", message);
        return 0;
    }
    return shader;
}

static unsigned int create_program(const char* shader_filepath) 
{
    char vertex_shader_content[MAX_SHADER_LENGTH + 1] = {0};
    size_t vertex_shader_length = 0;
    char fragment_shader_content[MAX_SHADER_LENGTH + 1] = {0};
    size_t fragment_shader_length = 0;

    FILE* file = fopen(shader_filepath, "r");
    if (file == NULL) {
        fprintf(stderr, "Failed to open shader file: %s\n", shader_filepath);
        fclose(file);
        return 0;
    }
    int parsing = -1;
    char* line = NULL;
    size_t line_buf_size = 0;
    int line_size;
    while ((line_size = getline(&line, &line_buf_size, file)) >= 0)
    {
        if (strstr(line, "#shader vertex") != NULL) {
            parsing = 1;
            line_buf_size = 0;
            line = NULL;
            continue;
        } else if (strstr(line, "#shader fragment") != NULL) {
            parsing = 2;
            line_buf_size = 0;
            line = NULL;
            continue;
        }
        if (parsing < 0) {}
        else {
            char* shader_start = (char*) (parsing == 1 ? &vertex_shader_content[vertex_shader_length] : &fragment_shader_content[fragment_shader_length]);
            size_t* length = parsing == 1 ? &vertex_shader_length : &fragment_shader_length;
            if (*length + line_size > MAX_SHADER_LENGTH) {
                fprintf(stderr, "Could not load shader: exceeded maximum length on line `%s` (%zu > %d)!\n", line, line_size + *length, MAX_SHADER_LENGTH);
                return 0;
            }
            memcpy(shader_start, line, line_size);
            *length += line_size;
        }
        line_buf_size = 0;
        line = NULL;
    }
    fclose(file);
    // printf("VERTEX\n");
    // printf("%s", vertex_shader_content);
    // printf("FRAGMENT\n");
    // printf("%s", fragment_shader_content);
    unsigned int program = glCreateProgram();

    unsigned int vertex_shader = compile_shader(&vertex_shader_content[0], GL_VERTEX_SHADER);
    if (vertex_shader == 0) return 0;
    unsigned int fragment_shader = compile_shader(&fragment_shader_content[0], GL_FRAGMENT_SHADER);
    if (fragment_shader == 0) return 0;

    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    int program_linked;
    glGetProgramiv(program, GL_LINK_STATUS, &program_linked);
    if (program_linked != GL_TRUE)
    {
        int length;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
        char* message = (char*) alloca(length * sizeof (char));
        glGetProgramInfoLog(program, length, &length, message);
        fprintf(stderr, "Failed to compile program shader!\n");
        fprintf(stderr, "%s", message);
        return 0;
    }
    glValidateProgram(program);
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
	return program;
}

typedef struct {
    float x, y, vx, vy;
} ParticleVertex;


typedef struct {
    unsigned short first, second, third;
} Triangle;

typedef struct {
    Triangle top, bottom;
} Quad;


/// Holds all state information relevant to a character as loaded using FreeType
typedef struct {
    unsigned int texture_id;  // ID handle of the glyph texture
    iVec2 size;  // Size of glyph
    iVec2 bearing;  // Offset from baseline to left/top of glyph
    unsigned int advance;  // Horizontal offset to advance to next glyph
} Character;

static Character characters_map[128];

void render_text(char* text, float left, float top, float scale, int program, unsigned int vao, unsigned int vbo, int text_color_location, Vec4 text_color) 
{
    // printf("text: %s, size: %zu, program: %d, vao: %u, vbo: %u, text loc: %d, r: %f, g: %f, b: %f\n", text, text_size, program, vao, vbo, text_color_location, r, g, b);
    GLCall(glUseProgram(program));
    GLCall(glUniform4f(text_color_location, text_color.x, text_color.y, text_color.z, text_color.w));
    GLCall(glActiveTexture(GL_TEXTURE0));
    GLCall(glBindVertexArray(vao));
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, vbo));

    // iterate through all characters
    float x = left;
    size_t text_size = strlen(text);
    for (size_t i = 0; i < text_size; i++) 
    {
        unsigned char c = text[i];
        if (c >= 128) continue;
        Character ch = characters_map[(unsigned char) c];

        float xpos = x + ch.bearing.x * scale;
        float ypos = top - (ch.size.y - ch.bearing.y) * scale;

        float w = ch.size.x * scale;
        float h = ch.size.y * scale;
        // update VBO for each character
        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },            
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },

            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }           
        };
        // render glyph texture over quad
        GLCall(glBindTexture(GL_TEXTURE_2D, ch.texture_id));
        // update content of VBO memory
        GLCall(glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices)); // be sure to use glBufferSubData and not glBufferData
        // render quad
        GLCall(glDrawArrays(GL_TRIANGLES, 0, 6));
        // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
    GLCall(glBindVertexArray(0));
    GLCall(glBindTexture(GL_TEXTURE_2D, 0));
}
// END Open GL


// BEGIN PHYSICS
typedef struct {
    Vec2 pos;
    Vec2 vel;
} Particle;

typedef struct {
    float density;
    float near_density;
} Density;

Density calculate_density(
    Vec2 pos,
    int* cells,
    int* cells_start,
    size_t n_cells_rows,
    size_t n_cells_cols,
    float window_width,
    float window_height,
    Vec2* positions,
    float particle_mass,
    float smoothing_radius,
    size_t n_particles)
{
    float density = 0;
    float near_density = 0;
    float squared_radius = smoothing_radius * smoothing_radius;
    size_t cell_x = floor((pos.x / window_width + 0.5) * n_cells_cols);
    size_t cell_y = floor((pos.y / window_height + 0.5) * n_cells_rows);

    // Iterate through neighbourghood only
    //printf("pos (%.02f, %.02f)\n", pos.x, pos.y);
    for (int offset_x = -1; offset_x < 2; offset_x++) {
        for (int offset_y = -1; offset_y < 2; offset_y++) {
            size_t cell_x_i = cell_x + offset_x;
            size_t cell_y_i = cell_y + offset_y;
            if ((cell_x_i < 0) || (cell_y_i < 0) || (cell_x_i >= n_cells_cols) || (cell_y_i >= n_cells_rows)) continue;
            //printf("offset(%zu, %zu)\n", offset_x, offset_y);
            int start = cells_start[cell_y_i*n_cells_cols + cell_x_i];
            while (start >= 0)
            {
                Vec2 pos_i = positions[start];
                //printf("\t (%d) offset(%zu, %zu) pos_i (%.02f, %.02f)\n", start, offset_x, offset_y, pos_i.x, pos_i.y);
                start = cells[start];

                float length_squared = vec2_squared_length(vec2_sub(pos_i, pos));
                if (length_squared >= squared_radius) continue;
                float length = sqrt(length_squared);
                density += pow((smoothing_radius - length), 2);
                near_density += pow((smoothing_radius - length), 3);
            }
            
        }
    }
    return (Density) {
        .density = particle_mass * density * 6 / (M_PI * pow(smoothing_radius, 4)),  // SpikyKernelPow2
        .near_density = particle_mass * near_density * 10 / (M_PI * pow(smoothing_radius, 5)),  // SpikyKernelPow3
    };
}


typedef struct {
    Vec2 viscosity;
    Vec2 pressure;
} Forces;


Forces calculate_forces(
    size_t particle_index,
    int* cells,
    int* cells_start,
    size_t n_cells_rows,
    size_t n_cells_cols,
    float window_width,
    float window_height,
    Vec2* positions,
    Vec2* velocities,
    float particle_mass,
    float smoothing_radius,
    float pressure_multiplier,
    float target_density,
    Density* densities,
    size_t n_particles,
    float viscosity_multiplier,
    float near_pressure_multiplier)
{
    Forces forces = {0};

	float density = densities[particle_index].density;
	float near_density = densities[particle_index].near_density;
	float pressure = (density - target_density) * pressure_multiplier;  // pressure from density
	float near_pressure = near_pressure_multiplier * near_density;  // near pressure from density
	Vec2 pressure_force = {0};
    float squared_radius = smoothing_radius * smoothing_radius;
	
	Vec2 pos = positions[particle_index];
	Vec2 velocity = velocities[particle_index];
    size_t cell_x = floor((pos.x / window_width + 0.5) * n_cells_cols);
    size_t cell_y = floor((pos.y / window_height + 0.5) * n_cells_rows);

    // Iterate through neighbourghood only
    for (int offset_x = -1; offset_x < 2; offset_x++) {
        for (int offset_y = -1; offset_y < 2; offset_y++) {
            size_t cell_x_i = cell_x + offset_x;
            size_t cell_y_i = cell_y + offset_y;
            if ((cell_x_i < 0) || (cell_y_i < 0) || (cell_x_i >= n_cells_cols) || (cell_y_i >= n_cells_rows)) continue;
            int start = cells_start[cell_y_i*n_cells_cols + cell_x_i];
            while (start >= 0)
            {
                size_t i = start;
                start = cells[start];
                if (i == particle_index) continue;

                Vec2 pos_i = positions[i];
                Vec2 velocity_i = velocities[i];
                float density_i = densities[i].density;
                float near_density_i = densities[i].near_density;

                Vec2 d_pos = vec2_sub(pos_i, pos);

                float length_squared = vec2_squared_length(d_pos);
                if (length_squared > squared_radius) continue;
                float length = sqrt(length_squared);

                Vec2 dir;
                if (length_squared == 0) {
                    dir = vec2_normalize((Vec2) {randnum(-1, 1), randnum(-1, 1)});
                } else {
                    dir = vec2_normalize(d_pos);
                }
                
                float pressure_i = (density_i - target_density) * pressure_multiplier;  // pressure from density
                float near_pressure_i = near_pressure_multiplier * near_density_i;  // near pressure from density

                float shared_pressure = (pressure + pressure_i) * 0.5;
                float shared_near_pressure = (near_pressure + near_pressure_i) * 0.5;

                forces.viscosity = vec2_add(forces.viscosity, 
                    vec2_scale(
                        vec2_sub(velocity_i, velocity), 
                        pow(smoothing_radius * smoothing_radius - length * length, 3) * viscosity_multiplier * 4 / (M_PI * pow(smoothing_radius, 8))
                    )
                );  // SmoothingKernelPoly6

                forces.pressure = vec2_add(
                    forces.pressure, 
                    vec2_scale(
                        dir, 
                        (
                            - (smoothing_radius - length) * 12 / (pow(smoothing_radius, 4) * M_PI) * shared_pressure / density_i +  // DensityDerivative
                            - pow(smoothing_radius - length, 2) * 30 / (pow(smoothing_radius, 5) * M_PI) * shared_near_pressure / near_density_i  // NearDensityDerivative
                )));
                // printf("(%zu) -> pressure_force (%f, %f); length %f; smoothing_radius %f; position (%f, %f); d_pos (%f, %f)\n", particle_index, pressure_force.x, pressure_force.y, length, smoothing_radius, pos.x, pos.y, d_pos.x, d_pos.y);
            }
        }
    }

    return forces;
}


Vec2 calculate_wall_force(Vec2 pos, float window_width, float window_height, float wall_force_radius, float wall_force_multiplier) 
{
    Vec2 wall_force = {0};
    float distance_to_left = (pos.x + window_width / 2);
    wall_force.x += distance_to_left < wall_force_radius ? wall_force_multiplier : 0;
    float distance_to_right = (window_width / 2 - pos.x);
    wall_force.x -= distance_to_right < wall_force_radius ? wall_force_multiplier: 0;

    float distance_to_bottom = (pos.y + window_height / 2);
    wall_force.y += distance_to_bottom < wall_force_radius ? wall_force_multiplier: 0;
    float distance_to_top = (window_height / 2 - pos.y);
    wall_force.y -= distance_to_top < wall_force_radius ? wall_force_multiplier: 0;
    // printf("(%f, %f)\n", wall_force.x, wall_force.y);
    return wall_force;
}

Vec2 calculate_interaction_force(size_t particle_index, Vec2 *positions, Vec2 *velocities, Vec2 interaction_position, float interaction_radius, float interaction_multiplier)
{
    Vec2 pos = positions[particle_index];
    Vec2 vel = velocities[particle_index];
    Vec2 d_pos = vec2_sub(interaction_position, pos);
    float length_squared = vec2_squared_length(d_pos);
    if (length_squared > interaction_radius*interaction_radius) {
        return (Vec2) {0};
    }
    float length = sqrt(length_squared);
    Vec2 interaction_force = vec2_add(
        vec2_scale(d_pos, (1 - length/interaction_radius)*interaction_multiplier),
        vec2_scale(vel, -(1 - length/interaction_radius))
    );
    return interaction_force;
}
// END PHYSICS

// BEGIN CELLS
static void fill_cells(int* cells_start, size_t number_of_items, int fill_val)
{
    memset(cells_start, -1, sizeof(int) * number_of_items);
}
// END CELLS


// BEGIN CLOCK
typedef struct {
    struct timespec last_call_ts;
    double dt;
} Clock;

static inline long int time_difference_ns(const struct timespec start, const struct timespec end)
{
    return 1e9*(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
}

static void clock_init(Clock *clock) 
{
    timespec_get(&clock->last_call_ts, TIME_UTC);
    clock->dt = 0;
}

// Ensure that between two consecutive calls to clock_tick it passes at least 1/FPS seconds
static inline void clock_tick(Clock *clock, const long int time_per_frame_ns) 
{
    struct timespec last_call_ts = clock->last_call_ts;
    struct timespec intermediate_ts;
    timespec_get(&intermediate_ts, TIME_UTC);
    long int elapsed_ns = time_difference_ns(last_call_ts, intermediate_ts);
    // printf("elpased_ns=%ld, time_per_frame_ns=%ld\n", elapsed_ns, time_per_frame_ns);
    long int waiting_ns = time_per_frame_ns - elapsed_ns;
    // printf("waiting_ns = %ld; elapsed_ns = %ld\n", waiting_ns, elapsed_ns);
    if (waiting_ns < 0) {
        timespec_get(&clock->last_call_ts, TIME_UTC);
        clock->dt = 1.0e-9 * time_difference_ns(last_call_ts, clock->last_call_ts);
        return;
    }
    struct timespec waiting_ts;
    waiting_ts.tv_sec = (waiting_ns / 1e9);
    waiting_ts.tv_nsec = waiting_ns - waiting_ts.tv_sec * 1e9;
    // printf("waiting_ts = { .tv_sec = %ld, .tv_nsec = %ld; }\n", waiting_ts.tv_sec, waiting_ts.tv_nsec);
    nanosleep(&waiting_ts, NULL);

    timespec_get(&clock->last_call_ts, TIME_UTC);
    clock->dt = 1.0e-9 * time_difference_ns(last_call_ts, clock->last_call_ts);
    // printf("dt: %f; fps: %f\n", clock->dt, 1.f / clock->dt);
    return;
}
// END CLOCK


// BEGIN GLOBAL SETTINGS
typedef struct {
    bool draw_density;
    bool draw_flow;
    bool draw_particles;
    bool draw_help;
    bool paused;
    bool gravity_on;
} GlobalSettings;

GlobalSettings global_settings;
bool need_reload_global_settings = false;
bool should_reset = false;
GlobalSettings updated_global_settings;

void reload_global_settings(unsigned int particles_program, int particle_draw_density_location, Vec2* gravity, const Vec2 gravity0) 
{
    if (updated_global_settings.draw_density != global_settings.draw_density) {
        GLCall(glUseProgram(particles_program));
        GLCall(glUniform1i(particle_draw_density_location, updated_global_settings.draw_density));
        GLCall(glUseProgram(0));
    } 
    if (updated_global_settings.gravity_on != global_settings.gravity_on) {
        *gravity = updated_global_settings.gravity_on ? gravity0 : vec2(0, 0);
    } 
    global_settings = updated_global_settings;
    return;
}
// END GLOBAL SETTINGS


// BEGIN USER INTERACTION
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(window, GLFW_TRUE);
    if (key == GLFW_KEY_D && action == GLFW_PRESS) {
        updated_global_settings.draw_density ^= 1; 
        need_reload_global_settings = true;
    }
    if (key == GLFW_KEY_P && action == GLFW_PRESS) {
        updated_global_settings.draw_particles ^= 1; 
        need_reload_global_settings = true;
    }
    if (key == GLFW_KEY_F && action == GLFW_PRESS) {
        updated_global_settings.draw_flow ^= 1; 
        need_reload_global_settings = true;
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS) {
        should_reset = true;
    }
    if (key == GLFW_KEY_H && action == GLFW_PRESS) {
        updated_global_settings.draw_help ^= 1; 
        need_reload_global_settings = true;
    }
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        updated_global_settings.paused ^= 1; 
        need_reload_global_settings = true;
    }
    if (key == GLFW_KEY_G && action == GLFW_PRESS) {
        updated_global_settings.gravity_on ^= 1; 
        need_reload_global_settings = true;
    }
}
// END USER INTERACTION


// BEGIN INITIALIZATION
static void init_particles(Particle* particles, size_t n_particles, float max_velocity, float spacing)
{
    size_t ncols = (int) sqrt(n_particles);
    size_t nrows = n_particles / ncols + 1;
    for (size_t i = 0; i < ncols; i++) {
        for (size_t j = 0; (j < nrows) && (i*nrows + j) < n_particles; j++) {
            size_t p_indx = i*nrows + j;
            Vec2 pos = {(i - ncols/2.f) * spacing, (nrows/2.f - j) * spacing};
            particles[p_indx] = (Particle) {.pos = pos, .vel = vec2(randnum(-max_velocity, max_velocity), randnum(-max_velocity, max_velocity))};
        }
    }
}

bool initialized_glfw = false;

static void terminate()
{
    if (initialized_glfw) {
        glfwTerminate();
    }
    return;
}

static GLFWwindow* initialize_glfw(const float window_width, const float window_height)
{
    glfwSetErrorCallback(glfw_error_callback);
 
    if (!glfwInit())
    {
        fprintf(stderr, "Failed to initialize GLFW!\n");
        exit(EXIT_FAILURE);
    }
    initialized_glfw = true;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
 
    GLFWwindow* window = glfwCreateWindow(window_width, window_height, "Fluid Simulator", NULL, NULL);
    if (!window)
    {
        fprintf(stderr, "Failed to create window!\n");
        terminate();
        exit(EXIT_FAILURE);
    }
    
    glfwSetKeyCallback(window, key_callback);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    printf("Using OpenGL version %s\n", glGetString(GL_VERSION));    

    return window;
}

static void initialize_glew(void) 
{
    unsigned int err = glewInit();
    if (err != GLEW_OK)
    {
        fprintf(stderr, "Failed to initialize GLEW: %s\n", glewGetErrorString(err));
        terminate();
        exit(EXIT_FAILURE);
    }
    printf("Using GLEW version %s\n", glewGetString(GLEW_VERSION));
}

static void initialize_fonts(const char* font_path)
{
    // FreeType
    // --------
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        fprintf(stderr, "ERROR::FREETYPE: Could not init FreeType Library");
        terminate();
        exit(EXIT_FAILURE);
    }
	// find path to font
    FT_Face face;
    if (FT_New_Face(ft, font_path, 0, &face)) {
        fprintf(stderr, "ERROR::FREETYPE: Failed to load font");
        terminate();
        exit(EXIT_FAILURE);
    }
    // set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 24);

    // disable byte-alignment restriction
    GLCall(glPixelStorei(GL_UNPACK_ALIGNMENT, 1));

    // load first 128 characters of ASCII set
    for (unsigned char c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            fprintf(stderr, "ERROR::FREETYTPE: Failed to load Glyph");
            continue;
        }
        // generate texture
        unsigned int texture;
        GLCall(glGenTextures(1, &texture));
        GLCall(glActiveTexture(GL_TEXTURE0));
        GLCall(glBindTexture(GL_TEXTURE_2D, texture));
        GLCall(glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer
        ));
        // set texture options
        GLCall(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
        GLCall(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
        GLCall(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
        GLCall(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
        // now store character for later use
        Character character = {
            texture,
            ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            (unsigned int) (face->glyph->advance.x)
        };
        characters_map[c] = character;
    }
    GLCall(glBindTexture(GL_TEXTURE_2D, 0));

    // destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);
}

static void initialize_density_vertices(size_t number_of_columns, size_t number_of_rows, Vec2 density_vertices[number_of_rows + 1][number_of_columns + 1], Quad density_indices[number_of_rows][number_of_columns], float window_width, float window_height)
{
    for (size_t j = 0; j < number_of_rows + 1; j++) {
        for (size_t i = 0; i < number_of_columns + 1; i++) {
            float v_pos_x = i * window_width / (float) number_of_columns - window_width / 2;
            float v_pos_y = window_height / 2 - j * window_height / (float) number_of_rows;
            density_vertices[j][i].x = v_pos_x;
            density_vertices[j][i].y = v_pos_y;
        }
    }
    for (size_t j = 0; (j < number_of_rows); j++) {
        for (size_t i = 0; i < number_of_columns; i++) {
            unsigned short top_left = i + j*(number_of_columns + 1);
            unsigned short top_right = i + 1 + j*(number_of_columns + 1);
            unsigned short bottom_left = i + (j + 1)*(number_of_columns + 1);
            unsigned short bottom_right = i + 1 + (j + 1)*(number_of_columns + 1);
            // printf("top_left %u (%f, %f)\n", top_left, ((Vec2*) density_vertices)[top_left].x, ((Vec2*) density_vertices)[top_left].y);
            density_indices[j][i] = (Quad) {
                .top = {
                    .first = top_left, .second = top_right, .third = bottom_left,
                }, .bottom = {
                    .first = bottom_left, .second = top_right, .third = bottom_right,
                },
            };
        }
    }
    return;
}
// END INITIALIZATION


// BEGIN PROGRAMS
typedef enum {
    UNIFORM_TYPE_INT,
    UNIFORM_TYPE_FLOAT,
    UNIFORM_TYPE_VEC2,
    UNIFORM_TYPE_VEC3,
    UNIFORM_TYPE_VEC4,
    UNIFORM_TYPE_TEXTURE1D
} UNIFORM_TYPE;

typedef struct {
    unsigned int id;
    unsigned int unit;
    GLenum internal_format;
    size_t width;
    GLenum format;
    GLenum type;
} Texture1D;

typedef struct {
    const char* name; 
    int location; 
    UNIFORM_TYPE type; 
    union {
        int initial_int;
        float initial_float;
        Vec2 initial_vec2;
        Vec3 initial_vec3;
        Vec4 initial_vec4;
        Texture1D texture1D;
    };
} Uniform;

typedef struct {
    const char* name; 
    size_t size_items;
    GLenum type;
    GLboolean normalized;
    size_t size_bytes;
    unsigned int divisor;
} Attribute;

typedef struct {
    unsigned int id; 
} VAO;

typedef struct {
    unsigned int id; 
    size_t size;
    void* data;
    GLenum usage;
} VBO;

typedef struct {
    unsigned int id; 
    size_t size;
    void* data;
    GLenum usage;
} IBO;

typedef struct {
    const char* program_name;
    const char* shader_path;
    Attribute *attributes;
    size_t n_attributes;
    Uniform *uniforms;
    size_t n_uniforms;
    VAO vao;
    VBO vbo;
    IBO ibo;
    unsigned int program_id;
} Program;


// May seg-fault but it is quick way in deed. 
// NOT MENT TO BE CALLED DURING LOOP, JUST AT THE BEGINNING OF THE CODE SO THAT IT IS ALMOST NEVER CHANGED.
static Uniform* get_uniform(Uniform* uniforms, char* name) 
{   
    const size_t max_iter = 15;
    for (size_t i = 0; i < max_iter; i++){
        if (strcmp(uniforms[i].name, name) == 0) {
            return &uniforms[i];
        } 
    }
    fprintf(stderr, "ERROR: could not find uniform %s\n", name);
    exit(EXIT_FAILURE);
}

static void init_program(Program* program)
{
    if (program->vao.id == 0) {
        GLCall(glGenVertexArrays(1, &program->vao.id));
    };
    GLCall(glBindVertexArray(program->vao.id));
    if (program->vbo.id == 0) {
        GLCall(glGenBuffers(1, &program->vbo.id));
    };
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, program->vbo.id));
    GLCall(glBufferData(GL_ARRAY_BUFFER, program->vbo.size, program->vbo.data, program->vbo.usage));
    if (program->ibo.size > 0) {
        GLCall(glGenBuffers(1, &program->ibo.id));
        GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, program->ibo.id));
        GLCall(glBufferData(GL_ELEMENT_ARRAY_BUFFER, program->ibo.size, program->ibo.data, program->ibo.usage));
    }

    size_t vertex_size = 0;
    for (size_t i = 0; i < program->n_attributes; i++) vertex_size += program->attributes[i].size_bytes;

    size_t offset = 0;
    for (size_t i = 0; i < program->n_attributes; i++) {
        GLCall(glEnableVertexAttribArray(i));
        GLCall(glVertexAttribPointer(
            i, 
            program->attributes[i].size_items, 
            program->attributes[i].type, 
            program->attributes[i].normalized, 
            vertex_size, 
            (const void *) offset));
        offset += program->attributes[i].size_bytes; 
        GLCall(glVertexAttribDivisor(i, program->attributes[i].divisor));
    }
    GLCall(glBindVertexArray(0));
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
    GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

    program->program_id = create_program(program->shader_path);
    if (program->program_id == 0) {
        fprintf(stderr, "Failed to create program %s!\n", program->program_name);
        terminate();
        exit(EXIT_FAILURE);
    }
    GLCall(glUseProgram(program->program_id));
    for (size_t i = 0; i < program->n_uniforms; i++)
    {   
        GLCall(int loc = glGetUniformLocation(program->program_id, program->uniforms[i].name));
        if (loc < 0) {
            fprintf(stderr, "WARNING: Could not get uniform location for `%s` in program %s: maybe it is not in use!\n",
            program->uniforms[i].name, program->program_name);
        }
        program->uniforms[i].location = loc;
        switch (program->uniforms[i].type) {
            case UNIFORM_TYPE_INT:
                GLCall(glUniform1i(loc, program->uniforms[i].initial_int));
                break;
            case UNIFORM_TYPE_FLOAT:
                GLCall(glUniform1f(loc, program->uniforms[i].initial_float));
                break;
            case UNIFORM_TYPE_VEC2:
                GLCall(glUniform2f(loc, program->uniforms[i].initial_vec2.x, program->uniforms[i].initial_vec2.y));
                break;
            case UNIFORM_TYPE_VEC3:
                GLCall(glUniform3f(loc, program->uniforms[i].initial_vec3.x, program->uniforms[i].initial_vec3.y, program->uniforms[i].initial_vec3.z));
                break;
            case UNIFORM_TYPE_VEC4:
                GLCall(glUniform4f(loc, program->uniforms[i].initial_vec4.x, program->uniforms[i].initial_vec4.y, program->uniforms[i].initial_vec4.z, program->uniforms[i].initial_vec4.w));
                break;
            case UNIFORM_TYPE_TEXTURE1D:
                GLCall(glUniform1i(loc, program->uniforms[i].texture1D.unit));
                GLCall(glGenTextures(1, &program->uniforms[i].texture1D.id));
                GLCall(glActiveTexture(GL_TEXTURE0 + program->uniforms[i].texture1D.unit));
                GLCall(glBindTexture(GL_TEXTURE_1D, program->uniforms[i].texture1D.id));
                GLCall(glTexImage1D(GL_TEXTURE_1D, 0, program->uniforms[i].texture1D.internal_format, program->uniforms[i].texture1D.width, 0, program->uniforms[i].texture1D.format, program->uniforms[i].texture1D.type, NULL));
                GLCall(glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
                GLCall(glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
                GLCall(glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
                break;
            default:
                assert(0);
        }
    }
    GLCall(glUseProgram(0));
    return;
}
// END PROGRAMS

// BEGIN MENUS
static void draw_help(const Program* text_program, const Uniform* text_color_uniform, float x, float y, int offset, Vec4 text_color)
{
    render_text("`ESC` > CLOSE", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text(global_settings.draw_particles ? "`P`   > HIDE PARTICLES": "`P`   > SHOW PARTICLES", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text(global_settings.draw_density ? "`D`   > HIDE DENSITY": "`D`   > DRAW DENSITY", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text(global_settings.draw_flow ? "`F`   > HIDE FLOW": "`F`   > DRAW FLOW", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text(global_settings.paused ? "`S`   > RESUME": "`S`   > STOP", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text(global_settings.gravity_on ? "`G`   > DISABLE GRAVITY": "`G`   > ENABLE_GRAVITY", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
    offset += 20;
    render_text("`R`   > RESET", x, y - offset, 0.5f, text_program->program_id, text_program->vao.id, text_program->vbo.id, text_color_uniform->location, text_color);
} 
// END MENUS


int main(void)
{
    float window_height = 800;
    float window_width = 1200;

    GLFWwindow* window = initialize_glfw(window_width, window_height);
    initialize_glew();
    
    // glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEBUG_OUTPUT);
    // glDebugMessageCallback(debug_callback, 0);  

    char* monospaced_font = "./res/fonts/AnonymousMonospaced.ttf";
    initialize_fonts(monospaced_font);
    
    // TEXT PROGRAM
    Uniform text_uniforms[] = {
        {.name = "text_color", .type = UNIFORM_TYPE_VEC4, .initial_vec4 = {1, 1, 1, 1}},
        {.name = "window_resolution", .type = UNIFORM_TYPE_VEC2, .initial_vec2 = {window_width, window_height}}
    };
    Uniform *text_color_uniform = get_uniform(text_uniforms, "text_color");

    Attribute text_attributes[] = {
        {.name = "position", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float)},
        {.name = "uv", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float)},
    };
    Program text_program = {
        .program_name = "text",
        .shader_path = "res/shaders/text.shader",
        .vao = {0},
        .ibo = {0},
        .vbo = {
            .size = sizeof(float) * 6 * 4,
            .data = NULL,
            .usage = GL_DYNAMIC_DRAW
        },
        .attributes = text_attributes,
        .n_attributes = sizeof(text_attributes) / sizeof(text_attributes[0]),
        .uniforms = text_uniforms,
        .n_uniforms = sizeof(text_uniforms) / sizeof(text_uniforms[0])
    };
    init_program(&text_program);

    // INTERACTION_FORCE PROGRAM
    Uniform interaction_force_uniforms[] = {
        {.name = "interaction_radius", .type = UNIFORM_TYPE_FLOAT, .initial_float = 0},
        {.name = "window_resolution", .type = UNIFORM_TYPE_VEC2, .initial_vec2 = {window_width, window_height}}
    };
    Uniform *interaction_radius_uniform = get_uniform(interaction_force_uniforms, "interaction_radius");
    Attribute interaction_force_attributes[] = {
        {.name = "position", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float), .divisor = 1},
    };
    Program interaction_force_program = {
        .program_name = "interaction_force",
        .shader_path = "res/shaders/interaction_force.shader",
        .vao = {0},
        .ibo = {0},
        .vbo = {
            .size = sizeof(Vec2) * 1,
            .data = NULL,
            .usage = GL_DYNAMIC_DRAW
        },
        .attributes = interaction_force_attributes,
        .n_attributes = sizeof(interaction_force_attributes) / sizeof(interaction_force_attributes[0]),
        .uniforms = interaction_force_uniforms,
        .n_uniforms = sizeof(interaction_force_uniforms) / sizeof(interaction_force_uniforms[0])
    };
    init_program(&interaction_force_program);

    // PARTICLES PROGRAM
#define MAX_NUMBER_OF_PARTICLES 3000
    ParticleVertex particle_vertices[MAX_NUMBER_OF_PARTICLES] = {0};
    size_t number_of_vertices;
    Uniform particles_uniforms[] = {
        {.name = "particle_radius", .type = UNIFORM_TYPE_FLOAT, .initial_float = 0},
        {.name = "draw_density", .type = UNIFORM_TYPE_INT, .initial_float = 0},
        {.name = "window_resolution", .type = UNIFORM_TYPE_VEC2, .initial_vec2 = {window_width, window_height}}
    };
    Uniform* particle_radius_uniform = get_uniform(particles_uniforms, "particle_radius");
    Uniform* draw_density_uniform = get_uniform(particles_uniforms, "draw_density");
    Attribute particles_attributes[] = {
        {.name = "position", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float), .divisor = 1},
        {.name = "velocity", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float), .divisor = 1},
    };
    Program particles_program = {
        .program_name = "particles",
        .shader_path = "res/shaders/particles.shader",
        .vao = {0},
        .ibo = {0},
        .vbo = {
            .size = sizeof(particle_vertices[0]) * MAX_NUMBER_OF_PARTICLES,
            .data = NULL,
            .usage = GL_DYNAMIC_DRAW
        },
        .attributes = particles_attributes,
        .n_attributes = sizeof(particles_attributes) / sizeof(particles_attributes[0]),
        .uniforms = particles_uniforms,
        .n_uniforms = sizeof(particles_uniforms) / sizeof(particles_uniforms[0])
    };
    init_program(&particles_program);

    Uniform flow_uniforms[] = {
        {.name = "max_arrow_size", .type = UNIFORM_TYPE_FLOAT, .initial_float = 50},
        {.name = "arrow_width", .type = UNIFORM_TYPE_FLOAT, .initial_float = 2},
        //{.name = "flow_color", .type = UNIFORM_TYPE_VEC4, .initial_vec4 = {0.5, 0.4, 0.7, 1}},
        {.name = "window_resolution", .type = UNIFORM_TYPE_VEC2, .initial_vec2 = {window_width, window_height}}
    };
    Uniform* max_arrow_size_uniform = get_uniform(particles_uniforms, "particle_radius");
    Program flow_program = {
        .program_name = "flow",
        .shader_path = "res/shaders/flow.shader",
        .vao = particles_program.vao,
        .ibo = {0},
        .vbo = particles_program.vbo,
        .attributes = particles_program.attributes,
        .n_attributes = particles_program.n_attributes,
        .uniforms = flow_uniforms,
        .n_uniforms = sizeof(flow_uniforms) / sizeof(flow_uniforms[0])
    };
    init_program(&flow_program);
    
#define NUMBER_OF_ROWS 200
#define NUMBER_OF_COLUMNS 200
    // TODO: This could be done with draw instanced...just save veertices for oone cell
    Vec2 density_vertices[NUMBER_OF_ROWS + 1][NUMBER_OF_COLUMNS + 1];
    // row x column x (top or bottom) x (triangle vertex)
    Quad density_indices[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
    initialize_density_vertices(NUMBER_OF_COLUMNS, NUMBER_OF_ROWS, density_vertices, density_indices, window_width, window_height);

#define MIN_SMOOTHING_RADIUS 20
#define MAX_WINDOW_WIDTH 2000
#define MAX_WINDOW_HEIGHT 2000
#define MAX_NUMBER_OF_CELLS MAX_WINDOW_HEIGHT / MIN_SMOOTHING_RADIUS * MAX_WINDOW_WIDTH / MIN_SMOOTHING_RADIUS

    Uniform density_uniforms[] = {
        {.name = "particle_mass", .type = UNIFORM_TYPE_FLOAT, .initial_float = 0},
        {.name = "smoothing_radius", .type = UNIFORM_TYPE_FLOAT, .initial_float = 0},
        {.name = "target_density", .type = UNIFORM_TYPE_FLOAT, .initial_float = 1},
        {.name = "window_resolution", .type = UNIFORM_TYPE_VEC2, .initial_vec2 = {window_width, window_height}},
        {.name = "n_cells_cols", .type = UNIFORM_TYPE_INT, .initial_int = 0},
        {.name = "n_cells_rows", .type = UNIFORM_TYPE_INT, .initial_int = 0},

        {.name = "positions_encoded", .type = UNIFORM_TYPE_TEXTURE1D,   .texture1D = {.unit = 0, .internal_format = GL_RG16F, .width = MAX_NUMBER_OF_PARTICLES, .format = GL_RG,          .type = GL_FLOAT}},
        {.name = "cells_encoded", .type = UNIFORM_TYPE_TEXTURE1D,       .texture1D = {.unit = 1, .internal_format = GL_R32I,  .width = MAX_NUMBER_OF_PARTICLES, .format = GL_RED_INTEGER, .type = GL_INT}},
        {.name = "cells_start_encoded", .type = UNIFORM_TYPE_TEXTURE1D, .texture1D = {.unit = 2, .internal_format = GL_R32I,  .width = MAX_NUMBER_OF_CELLS,     .format = GL_RED_INTEGER, .type = GL_INT}},
    };
    Uniform *particle_mass_uniform = get_uniform(density_uniforms, "particle_mass");
    Uniform *smoothing_radius_uniform = get_uniform(density_uniforms, "smoothing_radius");
    Uniform *target_density_uniform = get_uniform(density_uniforms, "target_density");
    Uniform *n_cells_cols_uniform = get_uniform(density_uniforms, "n_cells_cols");
    Uniform *n_cells_rows_uniform = get_uniform(density_uniforms, "n_cells_rows");

    Uniform *positions_encoded_uniform = get_uniform(density_uniforms, "positions_encoded");
    Uniform *cells_encoded_uniform = get_uniform(density_uniforms, "cells_encoded");
    Uniform *cells_start_encoded_uniform = get_uniform(density_uniforms, "cells_start_encoded");
    Attribute density_attributes[] = {
        {.name = "position", .size_items = 2, .type = GL_FLOAT, .normalized = GL_FALSE, .size_bytes = 2 * sizeof(float)},
    };
    Program density_program = {
        .program_name = "density",
        .shader_path = "res/shaders/density.shader",
        .vao = {0},
        .vbo = {
            .size = sizeof(density_vertices[0][0]) * (NUMBER_OF_ROWS + 1) * (NUMBER_OF_COLUMNS + 1),
            .data = density_vertices,
            .usage = GL_DYNAMIC_DRAW
        },
        .ibo = {
            .size = sizeof(density_indices[0][0]) * NUMBER_OF_COLUMNS * NUMBER_OF_ROWS,
            .data = (unsigned short*) density_indices,
            .usage = GL_STATIC_DRAW
        },
        .attributes = density_attributes,
        .n_attributes = sizeof(density_attributes) / sizeof(density_attributes[0]),
        .uniforms = density_uniforms,
        .n_uniforms = sizeof(density_uniforms) / sizeof(density_uniforms[0])
    };
    init_program(&density_program);
   
    // SETTINGS
    size_t number_of_particles = 2000;
    float particle_radius = 2.0f;
    float particle_mass = 0.5f;
    Vec2 gravity0 = (Vec2) {.x = 0, .y = -300};
    Vec2 gravity = gravity0;
    float dampig_coefficient = 0.7;
    float smoothing_radius = 50;
    float pressure_multiplier      = 35000;
    float near_pressure_multiplier = 5000;
    float viscosity_multiplier = 700;
    float target_density = 0.0003;
    float wall_force_radius = particle_radius*1.1;
    float wall_force_multiplier = 200;
    
    Vec2 interaction_position;
    float interaction_radius = 90;
    float interaction_multiplier_attraction = 80;
    float interaction_multiplier_repulsion = -80;
    float interaction_multiplier = 0;

    Particle particles[MAX_NUMBER_OF_PARTICLES] = {0};

    // srand(time(NULL));
    // for (size_t i = 0; i < number_of_particles; ++i) {
    //     particles[i] = (Particle) {.pos = {randnum(- window_width / 2, window_width / 2), randnum(- window_height / 2, window_height / 2)}, .vel = {randnum(-MAX_INITIAL_VELOCITY, MAX_INITIAL_VELOCITY), randnum(-MAX_INITIAL_VELOCITY, MAX_INITIAL_VELOCITY)}};
    // }
    int spacing = smoothing_radius*0.12;
    float max_velocity = 0;
    init_particles(particles, number_of_particles, max_velocity, spacing);

#define FPS 200
#define TIME_PER_FRAME_NS 1e9 / FPS

    Density densities[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 pressure_forces[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 viscosity_forces[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 interaction_forces[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 predicted_positions[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 velocities[MAX_NUMBER_OF_PARTICLES] = {0};
    Vec2 positions_uniform[MAX_NUMBER_OF_PARTICLES] = {0};
    
    // Number of cells depends on smoothing radius. Lets say that smoothing radius is never less than 10 and window size is less that 2000 x 2000
    assert(smoothing_radius >= MIN_SMOOTHING_RADIUS);
    assert(window_width <= MAX_WINDOW_WIDTH);
    assert(window_height <= MAX_WINDOW_HEIGHT);
    size_t n_cells_cols = ceil(window_width / smoothing_radius);
    size_t n_cells_rows = ceil(window_height / smoothing_radius);
    // This is going to store a linked list in the following way:
    // To get the index of the first particle at cell (i, j) you do: start = cells_start[j * n_cells_cols + i]
    // The first particle at that cell will be in in particles[start];
    // To get the next element you reach next = cells[start];
    // Until next is < 0 that you break.
    int cells_start[MAX_WINDOW_HEIGHT / MIN_SMOOTHING_RADIUS * MAX_WINDOW_WIDTH / MIN_SMOOTHING_RADIUS] = {-1};
    fill_cells(cells_start, n_cells_rows * n_cells_cols, -1);
    int cells[MAX_NUMBER_OF_PARTICLES] = {-1};

    Clock clock;
    clock_init(&clock);

    global_settings.draw_density = false;
    global_settings.draw_flow = true;
    global_settings.draw_particles = true;
    global_settings.gravity_on = true;
    updated_global_settings = global_settings;
    reload_global_settings(particles_program.program_id, draw_density_uniform->location, &gravity, gravity0);

    GLCall(glUseProgram(particles_program.program_id));
    GLCall(glUniform1f(particle_radius_uniform->location, particle_radius));
    GLCall(glUniform1i(draw_density_uniform->location, global_settings.draw_density));
  
    GLCall(glUseProgram(density_program.program_id));
    GLCall(glUniform1f(particle_mass_uniform->location, particle_mass));
    GLCall(glUniform1f(smoothing_radius_uniform->location, smoothing_radius));
    GLCall(glUniform1f(target_density_uniform->location, target_density));
    GLCall(glUniform1i(n_cells_cols_uniform->location, n_cells_cols));
    GLCall(glUniform1i(n_cells_rows_uniform->location, n_cells_rows));
    glActiveTexture(GL_TEXTURE0 + positions_encoded_uniform->texture1D.unit);
    glBindTexture(GL_TEXTURE_1D, positions_encoded_uniform->texture1D.id);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RG16F, MAX_NUMBER_OF_PARTICLES, 0, GL_RG,
  GL_FLOAT, NULL);


    GLCall(glUseProgram(interaction_force_program.program_id));
    GLCall(glUniform1f(interaction_radius_uniform->location, interaction_radius));

    GLCall(glUseProgram(0));

    TIME_THIS_INIT(total_execution)
    TIME_THIS_INIT(physics)
    TIME_THIS_INIT(physics_loop_calculate_densities)
    TIME_THIS_INIT(physics_loop_calculate_interaction_forces)
    TIME_THIS_INIT(physics_loop_calculate_viscosity)
    TIME_THIS_INIT(physics_loop_calculate_pressure)
    TIME_THIS_INIT(physics_loop_calculate_forces)

    TIME_THIS_BEGIN(total_execution)

    while (!glfwWindowShouldClose(window))
    {
        // Clock tick FPS
        clock_tick(&clock, TIME_PER_FRAME_NS);
            
        if (need_reload_global_settings) {
            reload_global_settings(particles_program.program_id, draw_density_uniform->location, &gravity, gravity0);
            need_reload_global_settings = false;
        }
        if (should_reset)
        {
            init_particles(particles, number_of_particles, max_velocity, spacing);
            should_reset = false;
        }
        
        int pressed_left = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);	
        int pressed_right = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);	
        interaction_multiplier = 0;
        if (pressed_right == GLFW_PRESS) {
            interaction_multiplier = interaction_multiplier_repulsion;
        } else if (pressed_left == GLFW_PRESS) {
            interaction_multiplier = interaction_multiplier_attraction;
        }
        double mouse_xpos, mouse_ypos;
        glfwGetCursorPos(window, &mouse_xpos, &mouse_ypos);
        interaction_position = vec2(mouse_xpos - window_width / 2, window_height / 2 - mouse_ypos);	
        
        if (!global_settings.paused)
        {
            TIME_THIS_BEGIN(physics)
            // Physics -> slow (6.4 / 7.2)

            // Physics-Predicted Positions -> fast (0.0021 / 7.2)
            for (size_t i = 0; i < number_of_particles; ++i) {
                predicted_positions[i] = vec2_add(particles[i].pos, vec2_scale(particles[i].vel, 1/120));
                velocities[i] = particles[i].vel;
            }

            // Organise in cells
            fill_cells(cells_start, n_cells_rows * n_cells_cols, -1);
            for (size_t i = 0; i < number_of_particles; ++i) {
                Vec2 pos = predicted_positions[i];
                size_t cell_x = floor((pos.x / window_width + 0.5) * n_cells_cols);
                size_t cell_y = floor((pos.y / window_height + 0.5) * n_cells_rows);
                int start = cells_start[cell_y * n_cells_cols + cell_x];
                // printf("n_pos_x %f, n_pos_y %f ;cell_x %zu, cell_y %zu -> %d\n", normalised_pos.x, normalised_pos.y, cell_x, cell_y, start);
                if (start < 0) {
                    cells_start[cell_y * n_cells_cols + cell_x] = i;
                    cells[i] = -1;
                } else {
                    int next = start;
                    while (next >= 0) {
                        start = next;
                        next = cells[start];
                    }
                    cells[start] = i;
                    cells[i] = -1;
                }
            }
            
            // Physics-Densities -> slow (2.0 / 7.2)
            TIME_THIS_BEGIN(physics_loop_calculate_densities)
            for (size_t i = 0; i < number_of_particles; ++i) {
                densities[i] = calculate_density(predicted_positions[i], cells, cells_start, n_cells_rows, n_cells_cols, window_width, window_height, predicted_positions, particle_mass, smoothing_radius, number_of_particles);
            }      
            TIME_THIS_END(physics_loop_calculate_densities)

            // Physics-Forces -> slow (4.4 / 7.2)

            TIME_THIS_BEGIN(physics_loop_calculate_forces)
            for (size_t i = 0; i < number_of_particles; ++i) {
                Forces forces = calculate_forces(i, cells, cells_start, n_cells_rows, n_cells_cols, window_width, window_height, predicted_positions, velocities, particle_mass, smoothing_radius, pressure_multiplier, target_density, densities, number_of_particles, viscosity_multiplier, near_pressure_multiplier);
                pressure_forces[i] = forces.pressure;
                viscosity_forces[i] = forces.viscosity;
            }
            TIME_THIS_END(physics_loop_calculate_forces)

            for (size_t i = 0; i < number_of_particles; ++i) {
                interaction_forces[i] = interaction_multiplier != 0 ? calculate_interaction_force(i, predicted_positions, velocities, interaction_position, interaction_radius, interaction_multiplier) : vec2(0, 0);
            }
            
            // Physics-Update Particles -> normal-fast (0.02 / 7.2)
            for (size_t i = 0; i < number_of_particles; ++i) {
                Vec2 pos = particles[i].pos;
                Vec2 vel = particles[i].vel;
                Vec2 interaction_force = interaction_forces[i];
                Vec2 pressure_force = pressure_forces[i];
                Vec2 viscosity_force = viscosity_forces[i];
                Vec2 wall_force = calculate_wall_force(pos, window_width, window_height, wall_force_radius, wall_force_multiplier);
                Vec2 acceleration = vec2_add(vec2_add(vec2_add(vec2_add(
                    gravity, 
                    vec2_scale(pressure_force, 1/densities[i].density)), 
                    vec2_scale(wall_force, 1/particle_mass)),
                    vec2_scale(viscosity_force, 1/particle_mass)),
                    interaction_force
                );
                Vec2 new_vel = vec2_add(vel, vec2_scale(acceleration, clock.dt));
                Vec2 new_pos = vec2_add(pos, vec2_scale(vec2_add(new_vel, vel), clock.dt/2));
                if ((new_pos.x + particle_radius) >= window_width / 2) {
                    new_pos.x = window_width / 2 - particle_radius;
                    new_vel.x *= -dampig_coefficient;
                } else if ((new_pos.x - particle_radius) <= - window_width / 2) {
                    new_pos.x = - window_width / 2 + particle_radius;
                    new_vel.x *= -dampig_coefficient;
                }
                if ((new_pos.y + particle_radius) >= window_height / 2) {
                    new_pos.y = window_height / 2 - particle_radius;
                    new_vel.y *= -dampig_coefficient;
                } else if ((new_pos.y - particle_radius) <= - window_height / 2) {
                    new_pos.y = - window_height / 2 + particle_radius;
                    new_vel.y *= -dampig_coefficient;
                }
                particles[i].pos = new_pos;
                particles[i].vel = new_vel;
            }
            TIME_THIS_END(physics)
        }

        // FILL Buffers -> fast (0.0008 / 7.2)
        for (size_t i = 0; i < number_of_particles; ++i) {
            particle_vertices[i] = (ParticleVertex) {particles[i].pos.x, particles[i].pos.y, particles[i].vel.x, particles[i].vel.y};
            positions_uniform[i] = (Vec2) {particles[i].pos.x, particles[i].pos.y};
        }
        number_of_vertices = number_of_particles;


        // Buffer Sub Data -> fast (0.0020 / 7.2)
        GLCall(glBindBuffer(GL_ARRAY_BUFFER, particles_program.vbo.id));
            GLCall(glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(particle_vertices[0]) * number_of_vertices, particle_vertices));
        GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));

        // Draw -> normal-slow (0.58 / 7.2)
        GLCall(glClear(GL_COLOR_BUFFER_BIT));    

        Vec4 text_color = {1, 1, 1, 0.6};
        if (global_settings.draw_density) {
            text_color = (Vec4) {0.8, 0.8, 0.8, 0.6};
            GLCall(glUseProgram(density_program.program_id));

            GLCall(glActiveTexture(GL_TEXTURE0 + positions_encoded_uniform->texture1D.unit));
            GLCall(glBindTexture(GL_TEXTURE_1D, positions_encoded_uniform->texture1D.id));
            GLCall(glTexSubImage1D(
                GL_TEXTURE_1D,
                0,
                0,
                MAX_NUMBER_OF_PARTICLES,
                GL_RG,
                GL_FLOAT,
                positions_uniform
            ));

            //int c[] = {1, 2, 3};
            // float c[] = {{1, 2}, {3, 4}, {5, 6}};
            GLCall(glActiveTexture(GL_TEXTURE0 + cells_encoded_uniform->texture1D.unit));
            GLCall(glBindTexture(GL_TEXTURE_1D, cells_encoded_uniform->texture1D.id));
            GLCall(glTexSubImage1D(    
                GL_TEXTURE_1D,
                0,
                0,
                MAX_NUMBER_OF_PARTICLES,
                GL_RED_INTEGER,
                GL_INT,
                cells
            ));

            GLCall(glActiveTexture(GL_TEXTURE0 + cells_start_encoded_uniform->texture1D.unit));
            GLCall(glBindTexture(GL_TEXTURE_1D, cells_start_encoded_uniform->texture1D.id));
            GLCall(glTexSubImage1D(
                GL_TEXTURE_1D,
                0,
                0,
                n_cells_cols * n_cells_rows,
                GL_RED_INTEGER,
                GL_INT,
                cells_start
            ));
            
            GLCall(glBindVertexArray(density_program.vao.id));
            //GLCall(glBindBuffer(GL_ARRAY_BUFFER, density_vbo));
            GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, density_program.ibo.id));
                GLCall(glDrawElements(GL_TRIANGLES, NUMBER_OF_COLUMNS * NUMBER_OF_ROWS * 3 * 2, GL_UNSIGNED_SHORT, (const void*) (0 * sizeof(unsigned short))));
            //GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
            GLCall(glBindTexture(GL_TEXTURE_1D, 0));
            GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
            GLCall(glBindVertexArray(0));
        }

        if (global_settings.draw_particles) {
            GLCall(glUseProgram(particles_program.program_id));
            GLCall(glBindVertexArray(particles_program.vao.id));
                GLCall(glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, number_of_vertices));
            GLCall(glBindVertexArray(0));
        }

        if (global_settings.draw_flow) {
            GLCall(glUseProgram(flow_program.program_id));
            GLCall(glBindVertexArray(flow_program.vao.id));
                GLCall(glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, number_of_vertices));
            GLCall(glBindVertexArray(0));
        }


        if (interaction_multiplier != 0) {
            GLCall(glUseProgram(interaction_force_program.program_id));
            GLCall(glBindBuffer(GL_ARRAY_BUFFER, interaction_force_program.vbo.id));
                GLCall(glBufferSubData(GL_ARRAY_BUFFER, 0, 1 * sizeof(Vec2), &interaction_position));
            GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
            GLCall(glBindVertexArray(interaction_force_program.vao.id));
                GLCall(glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1));
            GLCall(glBindVertexArray(0));
        }
        
        char fps_text[50];
        sprintf(fps_text, "FPS: %.02f", 1.f/clock.dt);
        render_text(fps_text, -window_width / 2 + 10, window_height/2 - 18, 0.6f, text_program.program_id, text_program.vao.id, text_program.vbo.id, text_color_uniform->location, text_color);

        if (global_settings.draw_help) {
            float x = window_width / 2 - 190;
            float y = window_height/2 - 18;
            render_text("PRESS `H` TO HIDE HELP", x, y, 0.5f, text_program.program_id, text_program.vao.id, text_program.vbo.id, text_color_uniform->location, text_color);
            draw_help(&text_program, text_color_uniform, window_width / 2 - 190, window_height/2 - 18, 20, text_color);
        } else {
            render_text("PRESS `H` FOR HELP", window_width / 2 - 160, window_height/2 - 18, 0.5f, text_program.program_id, text_program.vao.id, text_program.vbo.id, text_color_uniform->location, text_color);
        }

        GLCall(glUseProgram(0));
        // Swap Buffers -> normal-slow (0.19 / 7.2)
        glfwSwapBuffers(window);
        glfwPollEvents();

    }
    TIME_THIS_END(total_execution)

    TIME_THIS_RECAP(total_execution)
    TIME_THIS_RECAP(physics)
    TIME_THIS_RECAP(physics_loop_calculate_densities)
    TIME_THIS_RECAP(physics_loop_calculate_interaction_forces)
    TIME_THIS_RECAP(physics_loop_calculate_viscosity)
    TIME_THIS_RECAP(physics_loop_calculate_pressure)
    TIME_THIS_RECAP(physics_loop_calculate_forces)
    GLCall(glDeleteProgram(particles_program.program_id));
    GLCall(glDeleteProgram(density_program.program_id));
    GLCall(glDeleteProgram(text_program.program_id));
    glfwDestroyWindow(window);
    initialized_glfw = false;
    terminate(); 
    exit(EXIT_SUCCESS);
}
