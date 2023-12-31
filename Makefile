APP_NAME = fluid_simulator
BUILD_DIR = ./bin
C_FILES = ./src/*.c
CC = clang
FLAGS = -g -fcoverage-mapping -fprofile-instr-generate

APP_DEFINES:= -DGL_SILENCE_DEPRECATION
APP_INCLUDES:= -I./src/vendors/ -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -framework Carbon 
APP_LINKERS:= -L./src/vendors/GLFW/lib -lglfw3 -L./src/vendors/GLEW/lib -lGLEW -L./src/vendors/freetype/lib -lfreetype


build:
	$(CC) $(C_FILES) -o $(BUILD_DIR)/$(APP_NAME) $(APP_INCLUDES) $(APP_LINKERS) $(APP_DEFINES) $(FLAGS)

run: build
	$(BUILD_DIR)/$(APP_NAME)
