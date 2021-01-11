#include "view.h"
#include "src/mainwindow.h"
#include "viewformat.h"
#include "engine/graphics/GraphicsDebug.h"
#include "engine/graphics/Graphics.h"
#include "engine/graphics/Camera.h"
#include "engine/graphics/Material.h"
#include "engine/graphics/Shape.h"
#include "engine/graphics/GraphicsDebug.h"

#include <QApplication>
#include <QKeyEvent>
#include <QWindow>

#include <algorithm>

using namespace std;
using namespace glm;

View::View(QWidget *parent) : QGLWidget(ViewFormat(), parent),
    m_window(parent->parentWidget()),
    m_time(), m_timer(),
    m_captureMouse(true),
    m_right_mouse_pressed(false),
    m_mouse_pos(glm::vec2(0,0)),
    m_distance_along_look(5),
    m_fps(0), m_frameIndex(0),
    m_graphics(nullptr),
    m_camera(nullptr)
{
    /** SUPPORT CODE START **/

    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    if(m_captureMouse) {
        QApplication::setOverrideCursor(Qt::ArrowCursor);
    }
    else {
        QApplication::setOverrideCursor(Qt::ArrowCursor);
    }

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));

    // Initialize frame times for last FRAMES_TO_AVERAGE frames
    for (int i = 0; i < FRAMES_TO_AVERAGE; i++) {
        m_frameTimes[i] = 0;
    }

    m_frameIndex = 0;

    /** SUPPORT CODE END **/

}

View::~View()
{
}

void View::initializeGL()
{
    /** SUPPORT CODE START **/

    // Initialize graphics object
    m_graphics = Graphics::getGlobalInstance();

    // Enable depth testing, so that objects are occluded based on depth instead of drawing order.
    m_graphics->enableDepthTest();

    // Enable back-face culling, meaning only the front side of every face is rendered.
    // Also specify that the front face is represented by vertices in counterclockwise order (this is
    // the default).
    m_graphics->enableBackfaceCulling();

    // Enable alpha blending, so that materials with an alpha value < 1 are not totally opaque.
    m_graphics->enableBlendTest();

    // Disable stencil test for now (students may change this for final projects)
    m_graphics->disableStencilTest();

    // Start a timer that will try to get 60 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
    m_time.start();
    m_timer.start(1000 / 60);

    m_camera = std::make_shared<Camera>();
    m_camera->setEye(glm::vec3(0,0,0));
    m_graphics->setCamera(m_camera);

    m_graphics->addShader("lines", ":/shaders/lines.vert", ":/shaders/lines.frag");

    // camera for rendering axes in bottom left corner of screen
    m_axes_capture_camera = std::make_shared<Camera>(glm::vec2(150,150));

    std::vector<float> grid_positions;
    int s = 20;
    float i = -s;
    while (i <= s) {
        grid_positions.push_back(i); grid_positions.push_back(0); grid_positions.push_back(-s);
        grid_positions.push_back(i); grid_positions.push_back(0); grid_positions.push_back(s);
        grid_positions.push_back(-s); grid_positions.push_back(0); grid_positions.push_back(i);
        grid_positions.push_back(+s); grid_positions.push_back(0); grid_positions.push_back(i);
        i += 1;
    }
    std::vector<VBOAttribMarker> posMarkers;
    posMarkers.push_back(VBOAttribMarker(ShaderAttrib::POSITION, 3, 0));
    VBO posVBO(grid_positions.data(), grid_positions.size(), posMarkers, VBO::GEOMETRY_LAYOUT::LAYOUT_LINES);
    m_grid_vao = std::make_unique<VAO>(posVBO, grid_positions.size() / 3);

    // add light
    Light light1(Light::LIGHT_TYPE::POINT, glm::vec3(1,1,1), glm::vec3(0), glm::vec3(2,2,2), glm::vec2(0,0));
    m_graphics->addLight(light1);
    Light light2(Light::LIGHT_TYPE::POINT, glm::vec3(1,1,1), glm::vec3(0), glm::vec3(-2,2,-2), glm::vec2(0,0));
    m_graphics->addLight(light2);
    /** SUPPORT CODE END **/

    // create your physics system
    m_physics_system = std::make_unique<PhysicsSystem>();
}

void View::paintGL()
{
    /** SUPPORT CODE START **/

    checkError();
    glViewport(0,0,devicePixelRatio()*width(),devicePixelRatio()*height());
    m_graphics->setClearColor(glm::vec3(1,1,1));
    m_graphics->clearScreen(Graphics::CLEAR_FLAG::ALL);
    m_graphics->clearShader();
    m_graphics->setDefaultMaterial();

    checkError();

    m_graphics->setShader(m_graphics->getShader("lines"));
    checkError();
    m_graphics->setCamera(m_camera);
    m_graphics->clearTransform();
    m_grid_vao->draw();

    checkError();

    /** SUPPORT CODE END **/

    m_physics_system->draw();
    checkError();


    /** SUPPORT CODE START **/

    // draw some axes in the bottom left corner
    glViewport(50,50,200,200);
    m_graphics->setShader("default");
    m_graphics->setUseLighting(0);
    m_graphics->setCamera(m_axes_capture_camera);
    m_graphics->clearTransform();
    m_graphics->scale(glm::vec3(.2,2,.2));
    m_graphics->setColor(glm::vec3(0,1,0));
    m_graphics->drawShape("cylinder");
    m_graphics->clearTransform();
    m_graphics->rotate(M_PI / 2, glm::vec3(1,0,0));
    m_graphics->scale(glm::vec3(.2,2,.2));
    m_graphics->setColor(glm::vec3(0,0,1));
    m_graphics->drawShape("cylinder");
    m_graphics->clearTransform();
    m_graphics->rotate(3*M_PI / 2, glm::vec3(0,0,1));
    m_graphics->scale(glm::vec3(.2,2,.2));
    m_graphics->setColor(glm::vec3(1,0,0));
    m_graphics->drawShape("cylinder");
    m_graphics->setUseLighting(1);

    checkError();

#if GRAPHICS_DEBUG_LEVEL > 0
    m_graphics->printDebug();
    m_graphics->printShaderDebug();
    m_graphics->printFBODebug();
#endif

    /** SUPPORT CODE END **/
}

void View::resizeGL(int w, int h)
{
    /** SUPPORT CODE START **/

    m_graphics->setViewport(glm::vec2(0, 0), glm::vec2(w, h));
    m_camera->setScreenSize(glm::vec2(w,h));

    /** SUPPORT CODE END **/
}

void View::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        m_right_mouse_pressed = true;
    }
    MainWindow * win = (MainWindow *) QApplication::activeWindow();
    if (win->getMode() == PhysicsDebuggerMode::COLLISION_MODE) {

    }
}

void View::mouseMoveEvent(QMouseEvent *event)
{
    /** SUPPORT CODE START **/

    int deltaX = m_mouse_pos.x - event->x();
    int deltaY = m_mouse_pos.y - event->y();

    m_mouse_pos.x = event->x();
    m_mouse_pos.y = event->y();

    if (m_right_mouse_pressed) {
        m_camera->rotate(deltaX / 100.0, deltaY / 100.0);
        m_axes_capture_camera->rotate(deltaX / 100.0, deltaY / 100.0);
    }

    /** SUPPORT CODE END **/
}

void View::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        m_right_mouse_pressed = false;
    }
}

void View::wheelEvent(QWheelEvent *event)
{
    if (event->delta() > 0) {
        m_distance_along_look -= 1;
    } else {
        m_distance_along_look += 1;
    }
    m_distance_along_look = std::min(35.0f, std::max(5.0f, m_distance_along_look));
}

void View::keyPressEvent(QKeyEvent *event)
{
    /** SUPPORT CODE START **/

    // Don't remove this -- helper code for key repeat events
    if(event->isAutoRepeat()) {
        keyRepeatEvent(event);
        return;
    }

    // Feel free to remove this
    if (event->key() == Qt::Key_Escape) QApplication::quit();

    /** SUPPORT CODE END **/
}

void View::keyRepeatEvent(QKeyEvent *event)
{

}

void View::keyReleaseEvent(QKeyEvent *event)
{
    /** SUPPORT CODE START **/

    // Don't remove this -- helper code for key repeat events
    if(event->isAutoRepeat()) {
        return;
    }

    if (event->key() == Qt::RightButton) {
        m_right_mouse_pressed = false;
    }

    /** SUPPORT CODE END **/
}

void View::tick()
{
    /** SUPPORT CODE START **/

    // Get the number of seconds since the last tick (variable update rate)
    float seconds = m_time.restart() * 0.001f;

    m_frameTimes[m_frameIndex] = seconds;
    m_frameIndex = (m_frameIndex + 1) % FRAMES_TO_AVERAGE;
    m_fps = 0;
    for (int i = 0; i < FRAMES_TO_AVERAGE; i++) {
        m_fps += m_frameTimes[i];
    }
    m_fps /= FRAMES_TO_AVERAGE;
    m_fps = 1.f / m_fps;

    // Display fps
    QString title = "CS195U Engine";
    m_window->setWindowTitle(title + ", FPS: " + QString::number(m_fps, 'f', 3));

    // set camera
    glm::vec3 look = m_camera->getLook();
    m_camera->setEye(glm::vec3(0,2,0) - m_distance_along_look*look);
    m_axes_capture_camera->setEye(-5.0f*look);

    /** SUPPORT CODE END **/

    MainWindow * win = (MainWindow *) QApplication::activeWindow();
    if (win && win->getMode() == PhysicsDebuggerMode::PHYSICS_MODE) {
        m_physics_system->tick(seconds);
    } else if (win && win->getMode() == PhysicsDebuggerMode::COLLISION_MODE) {

    }

    /** SUPPORT CODE START **/

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();

    /** SUPPORT CODE END **/
}
