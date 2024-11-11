#ifndef WINDOW_MANAGER_H
#define WINDOW_MANAGER_H

#include <SDL2/SDL.h>
#include <iostream>

class WindowManager {
public:
    WindowManager(int, int, double offset_x = 0, double offset_y = 0, double scale = 1);
    ~WindowManager();
    bool init();
    void close();
    void setOffset(int x, int y);
    void setScale(double scale);
    int getOffsetX() const ;
    int getOffsetY() const ;
    double getScale() const ;
    SDL_Window* getWindow() const ;
    SDL_Renderer* getRenderer() const ;
    void setColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a);

private:
    SDL_Window* window;
    SDL_Renderer* renderer;
    int SCREEN_WIDTH_, SCREEN_HEIGHT_;
    int x_offset_, y_offset_; // (0,0) actual screen coordinates
    double scale_;
};

#endif // WINDOW_MANAGER_H