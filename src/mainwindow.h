#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "src/engine/util/CommonIncludes.h"

namespace Ui {
    class MainWindow;
}

enum class PhysicsDebuggerMode {
    COLLISION_MODE, PHYSICS_MODE
};

enum class PhysicsDebuggerObject {
    CYLINDER, CUBOID, ELLIPSOID
};

/**
 * QTs "frame" class, opens a window for you to draw in and interact with.
 *
 * @brief The MainWindow class
 */
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    PhysicsDebuggerMode getMode() { return mode; }

    glm::vec3 getPos1() { return pos1; }
    glm::vec3 getRot1() { return rot1; }
    glm::vec3 getScale1() { return scale1; }

    glm::vec3 getPos2() { return pos2; }
    glm::vec3 getRot2() { return rot2; }
    glm::vec3 getScale2() { return scale2; }

    PhysicsDebuggerObject getObj1() { return object1; }
    PhysicsDebuggerObject getObj2() { return object2; }

private slots:
    void on_collisionModeRadioButton_toggled(bool checked);
    void on_physicsModeRadioButton_toggled(bool checked);

    void on_xPos1_valueChanged(double arg1);

    void on_yPos1_valueChanged(double arg1);

    void on_zPos1_valueChanged(double arg1);

    void on_roll1_valueChanged(double arg1);

    void on_pitch1_valueChanged(double arg1);

    void on_yaw1_valueChanged(double arg1);

    void on_xScale1_valueChanged(double arg1);

    void on_yScale1_valueChanged(double arg1);

    void on_zScale1_valueChanged(double arg1);

    void on_xScale2_valueChanged(double arg1);

    void on_yPos2_valueChanged(double arg1);

    void on_zPos2_valueChanged(double arg1);

    void on_roll2_valueChanged(double arg1);

    void on_pitch2_valueChanged(double arg1);

    void on_yaw2_valueChanged(double arg1);

    void on_yScale2_valueChanged(double arg1);

    void on_zScale2_valueChanged(double arg1);

    void on_xPos2_valueChanged(double arg1);

    void on_radioButton_4_toggled(bool checked);

    void on_radioButton_5_toggled(bool checked);

    void on_radioButton_6_toggled(bool checked);

    void on_radioButton_toggled(bool checked);

    void on_radioButton_2_toggled(bool checked);

    void on_radioButton_3_toggled(bool checked);

private:
    Ui::MainWindow *ui;
    PhysicsDebuggerMode mode;

    PhysicsDebuggerObject object1 = PhysicsDebuggerObject::CUBOID;
    PhysicsDebuggerObject object2 = PhysicsDebuggerObject::CUBOID;

    glm::vec3 pos1 = glm::vec3(1,0,0);
    glm::vec3 rot1 = glm::vec3(0,0,0);
    glm::vec3 scale1 = glm::vec3(1,1,1);

    glm::vec3 pos2 = glm::vec3(-1,0,0);
    glm::vec3 rot2 = glm::vec3(0,0,0);
    glm::vec3 scale2 = glm::vec3(1,1,1);
};

#endif // MAINWINDOW_H

