#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

enum class PhysicsDebuggerMode {
    COLLISION_MODE, PHYSICS_MODE
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

private slots:
    void on_collisionModeRadioButton_toggled(bool checked);
    void on_physicsModeRadioButton_toggled(bool checked);

private:
    Ui::MainWindow *ui;
    PhysicsDebuggerMode mode;

};

#endif // MAINWINDOW_H

