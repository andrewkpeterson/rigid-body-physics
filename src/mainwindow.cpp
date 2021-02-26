#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_collisionModeRadioButton_toggled(bool checked)
{
    if (checked) {
        mode = PhysicsDebuggerMode::COLLISION_MODE;
    } else {
        mode = PhysicsDebuggerMode::PHYSICS_MODE;
    }
}

void MainWindow::on_physicsModeRadioButton_toggled(bool checked)
{
    if (checked) {
        mode = PhysicsDebuggerMode::PHYSICS_MODE;
    } else {
        mode = PhysicsDebuggerMode::COLLISION_MODE;
    }
}

void MainWindow::on_xPos1_valueChanged(double arg1)
{
    pos1 = glm::vec3(arg1, pos1.y, pos1.z);
}

void MainWindow::on_yPos1_valueChanged(double arg1)
{
    pos1 = glm::vec3(pos1.x, arg1, pos1.z);
}

void MainWindow::on_zPos1_valueChanged(double arg1)
{
    pos1 = glm::vec3(pos1.x, pos1.y, arg1);
}

void MainWindow::on_roll1_valueChanged(double arg1)
{
    rot1 = glm::vec3(arg1, rot1.y, rot1.z);
}

void MainWindow::on_pitch1_valueChanged(double arg1)
{
    rot1 = glm::vec3(rot1.x, rot1.y, arg1);
}

void MainWindow::on_yaw1_valueChanged(double arg1)
{
    rot1 = glm::vec3(rot1.x, arg1, rot1.z);
}

void MainWindow::on_xScale1_valueChanged(double arg1)
{
    scale1 = glm::vec3(arg1, scale1.y, scale1.z);
}

void MainWindow::on_yScale1_valueChanged(double arg1)
{
    scale1 = glm::vec3(scale1.x, arg1, scale1.z);
}

void MainWindow::on_zScale1_valueChanged(double arg1)
{
    scale1 = glm::vec3(scale1.x, scale1.y, arg1);
}

void MainWindow::on_xScale2_valueChanged(double arg1)
{
    scale2 = glm::vec3(arg1, scale2.y, scale2.z);
}

void MainWindow::on_yPos2_valueChanged(double arg1)
{
    pos2 = glm::vec3(pos2.x, arg1, pos2.z);
}

void MainWindow::on_zPos2_valueChanged(double arg1)
{
    pos2 = glm::vec3(pos2.x, pos2.y, arg1);
}

void MainWindow::on_roll2_valueChanged(double arg1)
{
    rot2 = glm::vec3(arg1, rot2.y, rot2.z);
}

void MainWindow::on_pitch2_valueChanged(double arg1)
{
    rot2 = glm::vec3(rot2.x, rot2.y, arg1);
}

void MainWindow::on_yaw2_valueChanged(double arg1)
{
    rot2 = glm::vec3(rot2.x, arg1, rot2.z);
}

void MainWindow::on_yScale2_valueChanged(double arg1)
{
    scale2 = glm::vec3(scale2.x, arg1, scale2.z);
}

void MainWindow::on_zScale2_valueChanged(double arg1)
{
    scale2 = glm::vec3(scale2.x, scale2.y, arg1);
}

void MainWindow::on_xPos2_valueChanged(double arg1)
{
    pos2 = glm::vec3(arg1, pos2.y, pos2.z);
}

void MainWindow::on_radioButton_4_toggled(bool checked)
{
    object1 = PhysicsDebuggerObject::CUBOID;
}

void MainWindow::on_radioButton_5_toggled(bool checked)
{
    object1 = PhysicsDebuggerObject::CYLINDER;
}

void MainWindow::on_radioButton_6_toggled(bool checked)
{
    object1 = PhysicsDebuggerObject::ELLIPSOID;
}

void MainWindow::on_radioButton_toggled(bool checked)
{
    object2 = PhysicsDebuggerObject::CUBOID;
}

void MainWindow::on_radioButton_2_toggled(bool checked)
{
    object2 = PhysicsDebuggerObject::CYLINDER;
}

void MainWindow::on_radioButton_3_toggled(bool checked)
{
    object2 = PhysicsDebuggerObject::ELLIPSOID;
}
