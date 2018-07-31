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
    QCoreApplication::exit();
}

void MainWindow::onGetXsWindowClosed()
{
    this->ui->actionGetXs->setEnabled(true);
}

void MainWindow::onGetTrainWindowClosed()
{
    this->ui->actionTrain->setEnabled(true);
}

void MainWindow::onGetTestWindowClosed()
{
    this->ui->actionTest->setEnabled(true);
}

void MainWindow::onGetVSWindowClosed()
{
    this->ui->actionScore->setEnabled(true);
}

void MainWindow::on_actionGetXs_triggered()
{
    xsWindow = new GetXsWindow(this);
    xsWindow->setWindowTitle("Get Xs Window");
    xsWindow->setWindowFlags(xsWindow->windowFlags() | Qt::WindowMinimizeButtonHint);
    xsWindow->setAttribute(Qt::WA_DeleteOnClose);
    xsWindow->show();

    this->ui->actionGetXs->setEnabled(false);
}

void MainWindow::on_actionTrain_triggered()
{
    train_window = new TrainWindow(this);
    train_window->setWindowTitle("Train Window");
    train_window->setWindowFlags(train_window->windowFlags() | Qt::WindowMinimizeButtonHint);
    train_window->setAttribute(Qt::WA_DeleteOnClose);
    train_window->show();

    this->ui->actionTrain->setEnabled(false);
}

void MainWindow::on_actionTest_triggered()
{
    test_window = new TestWindow(this);
    test_window->setWindowTitle("Test Window");
    test_window->setWindowFlags(test_window->windowFlags() | Qt::WindowMinimizeButtonHint);
    test_window->setAttribute(Qt::WA_DeleteOnClose);
    test_window->show();

    this->ui->actionTest->setEnabled(false);
}

void MainWindow::on_actionScore_triggered()
{
    vs_window = new VSWindow(this);
    vs_window->setWindowTitle("Score Window");
    vs_window->setWindowFlags(vs_window->windowFlags() | Qt::WindowMinimizeButtonHint);
    vs_window->setAttribute(Qt::WA_DeleteOnClose);
    vs_window->show();

    this->ui->actionScore->setEnabled(false);
}

void MainWindow::on_actionExit_triggered()
{
    delete ui;
    QCoreApplication::exit();
}
