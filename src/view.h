#ifndef VIEW_H
#define VIEW_H

#include "engine/util/CommonIncludes.h"

#include <QGLWidget>
#include <QTime>
#include <QTimer>
#include <memory>
#include <glm/glm.hpp>
#include "src/engine/graphics/ShaderAttribLocations.h"
#include "src/engine/graphics/VAO.h"
#include "src/engine/graphics/VBOAttribMarker.h"
#include "src/engine/graphics/VBO.h"
#include "src/engine/physics/PhysicsSystem.h"

class Graphics;
class Camera;

/**
 * This is similar to your "CS1971FrontEnd" class. Here you will receive all of the input events
 * to forward to your game.
 *
 * @brief The View class
 */
class View : public QGLWidget
{
    Q_OBJECT

public:
    View(QWidget *parent);
    ~View();

private:
    static const int FRAMES_TO_AVERAGE = 30;

private:
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void keyRepeatEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

private:
    QWidget *m_window;
    QTime m_time;
    QTimer m_timer;
    bool m_captureMouse;
    bool m_right_mouse_pressed;
    glm::vec2 m_mouse_pos;
    float m_distance_along_look;

    float m_fps;
    int m_frameIndex;
    float m_frameTimes[FRAMES_TO_AVERAGE];

    Graphics* m_graphics;

    std::unique_ptr<VAO> m_grid_vao;

    std::shared_ptr<Camera> m_axes_capture_camera;

    std::shared_ptr<Camera> m_camera;

    std::unique_ptr<PhysicsSystem> m_physics_system;

private slots:
    void tick();
};

#endif // VIEW_H

