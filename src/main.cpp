#include <SDL2/SDL.h>
#include <iostream>
#include "window_manager.h"

// Screen dimension constants
const int SCREEN_WIDTH = 1080;
const int SCREEN_HEIGHT = 960;

int main(int argc, char* args[]) {
    WindowManager windowmanager(SCREEN_WIDTH, SCREEN_HEIGHT);

    // Initialize SDL and create window and renderer
    if (!windowmanager.init()) {
        std::cerr << "Failed to initialize!" << std::endl;
        return -1;
    }

    SDL_Window* window = windowmanager.getWindow();
    SDL_Renderer* renderer = windowmanager.getRenderer();

    // Set draw color to blue
    windowmanager.setColor(0, 0, 255, 255);

    // Draw a rectangle (cuboid)
    SDL_Rect rect;
    rect.x = SCREEN_WIDTH / 4;
    rect.y = SCREEN_HEIGHT / 4;
    rect.w = SCREEN_WIDTH / 2;
    rect.h = SCREEN_HEIGHT / 2;
    SDL_RenderFillRect(renderer, &rect);

    // Present the renderer
    SDL_RenderPresent(renderer);

    // Wait for user to close the window
    bool quit = false;
    SDL_Event e;
    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
        }
    }

    // Free resources and close SDL
    windowmanager.close();

    return 0;
}