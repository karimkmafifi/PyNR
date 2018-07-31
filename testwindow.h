#ifndef TESTWINDOW_H
#define TESTWINDOW_H

#include "testpdbbind.h"

namespace Ui {
class TestWindow;
}

class TestWindow : public QDialog
{
    Q_OBJECT

public:
    explicit TestWindow(QWidget *parent = 0);
    ~TestWindow();

private:
    Ui::TestWindow *ui;
    QString last_dir;
    TestPDBbind *test_pdb_bind;
    bool testIsInitiated;
    void changeWorkStyle();

signals:
    void getTestWindowClosed();

private slots:
    void on_expPushButton_pressed();
    void on_expPushButton_released();
    void on_testPushButton_pressed();
    void on_testPushButton_released();
    void on_logPushButton_pressed();
    void on_logPushButton_released();
    void on_testSubmitButton_pressed();
    void on_testSubmitButton_released();

public slots:
    void ontestStringEmission(QString);
    void ontestTaskFinished();
    void ontestErrorEmission(QString);
};

#endif // TESTWINDOW_H
