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
