#include "render.h"
#include <ctime>

void InitCamera(Camera &camera) {
    camera.Oversampling = 1;
    camera.accuracy = 0.0005;
    camera.direction = glm::vec<3, float>(0, 0, -1);
    camera.up = glm::vec<3, float>(0, -1, 0);
    camera.right  = glm::vec<3, float>(1, 0, 0);
}
void check(Renderer &renderer) {
    float D,U,R;
    scanf("%f%f%f", &D,&U,&R);
    int sav = clock();
    InitCamera(renderer.camera);
    renderer.camera.position = D * renderer.camera.direction + U * renderer.camera.up + R * renderer.camera.right;
    renderer.render();
    printf("%d\n", clock() - sav);
    renderer.saveImage("test.png");
    renderer.clearImage();
}

int main() {
#if defined(__OPTIMIZE__) && !defined(__OPTIMIZE_SIZE__)
    std::cout << "O2 optimization is enabled." << std::endl;
#else
    std::cout << "O2 optimization is not enabled." << std::endl;
#endif
    int sav = clock();
    Renderer renderer(3200, 2400, "fbx/model.fbx");
    printf("%d\n", clock() - sav);
    while(true) check(renderer);
    return 0;
}
