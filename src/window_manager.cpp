#include "window_manager.h"

WindowManager::WindowManager(int width, int height, double offset_x, double offset_y, double scale) 
    : SCREEN_WIDTH_(width), SCREEN_HEIGHT_(height), x_offset_(offset_x), y_offset_(offset_y), 
        window(nullptr), renderer(nullptr), scale_(scale) {
            x_offset_ += SCREEN_WIDTH_ >> 1;
            y_offset_ += SCREEN_HEIGHT_ >> 1;
            scale_ *= 100.0;
        }

WindowManager::~WindowManager() {
    close();
}

bool WindowManager::init() {
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        return false;
    }

    // Create window
    window = SDL_CreateWindow("SDL2 Drawing Example", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH_, SCREEN_HEIGHT_, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        std::cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return false;
    }

    // Create renderer
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return false;
    }

    // Set renderer color
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // White background

    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    return true;
}

void WindowManager::close() {
    if (renderer) {
        SDL_DestroyRenderer(renderer);
        renderer = nullptr;
    }
    if (window) {
        SDL_DestroyWindow(window);
        window = nullptr;
    }
    SDL_Quit();
}

SDL_Window* WindowManager::getWindow() const {
    return window;
}

SDL_Renderer* WindowManager::getRenderer() const {
    return renderer;
}

void WindowManager::setColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
    // 0 is fully transparent, 255 is fully opaque
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
}

void WindowManager::setOffset(int x, int y) {
    x_offset_ = x + SCREEN_WIDTH_ / 2;
    y_offset_ = y + SCREEN_HEIGHT_ / 2;
}

void WindowManager::setScale(double scale) {
    scale_ = scale;
}

int WindowManager::getOffsetX() const {
    return x_offset_;
}

int WindowManager::getOffsetY() const {
    return y_offset_; 
}

double WindowManager::getScale() const {
    return scale_;
}