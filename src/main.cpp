#include "render.h"
#include <ctime>

void InitCamera(Camera &camera) {
    camera.position = glm::vec<3, float>(-50, 0, 0);
    camera.accuracy = 0.01;
    camera.direction = glm::vec<3, float>(1, 0, 0);
    camera.up = glm::vec<3, float>(0, 1, 0);
    camera.right = glm::vec<3, float>(0, 0, 1);
}
int main() {
#if defined(__OPTIMIZE__) && !defined(__OPTIMIZE_SIZE__)
    std::cout << "O2 optimization is enabled." << std::endl;
#else
    std::cout << "O2 optimization is not enabled." << std::endl;
#endif
    int sav = clock();
    Renderer renderer(256, 256, "fbx/chair.fbx");
    InitCamera(renderer.camera);
    renderer.render();
    renderer.saveImage("test.png");
    printf("%d", clock() - sav);
    return 0;
}
