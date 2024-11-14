#include "render.h"
#include <ctime>

void InitCamera(Camera &camera) {
    camera.accuracy = 0.001;
    camera.direction = glm::vec<3, float>(1/sqrt(2), 1/sqrt(2), 0);
    camera.up = glm::vec<3, float>(1/sqrt(2), -1/sqrt(2), 0);
    camera.right = glm::vec<3, float>(0, 0, 1);
    camera.position = -(float)64 * camera.direction;
}
int main() {
#if defined(__OPTIMIZE__) && !defined(__OPTIMIZE_SIZE__)
    std::cout << "O2 optimization is enabled." << std::endl;
#else
    std::cout << "O2 optimization is not enabled." << std::endl;
#endif
    int sav = clock();
    Renderer renderer(1024, 1024, "fbx/chair.fbx");
    InitCamera(renderer.camera);
    renderer.render();
    renderer.saveImage("test.png");
    printf("%d", clock() - sav);
    return 0;
}
