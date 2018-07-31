#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "getxswindow.h"
#include "trainwindow.h"
#include "testwindow.h"
#include "vswindow.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

private slots:

    void on_actionGetXs_triggered();

    void on_actionTrain_triggered();

    void on_actionExit_triggered();

    void on_actionScore_triggered();

    void on_actionTest_triggered();

public slots:
    void onGetXsWindowClosed();
    void onGetTrainWindowClosed();
    void onGetTestWindowClosed();
    void onGetVSWindowClosed();

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    GetXsWindow *xsWindow;
    TrainWindow *train_window;
    TestWindow *test_window;
    VSWindow *vs_window;

    QString last_dir;

};

#endif // MAINWINDOW_H
