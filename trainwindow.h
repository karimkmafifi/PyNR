#ifndef TRAINWINDOW_H
#define TRAINWINDOW_H

#include "scorepdbbind.h"

namespace Ui {
class TrainWindow;
}

class TrainWindow : public QDialog
{
    Q_OBJECT

public:
    explicit TrainWindow(QWidget *parent = 0);
    ~TrainWindow();

public slots:
    void onlrStringEmission(QString);
    void onlrTaskFinished();
    void onlrErrorEmission(QString);
    void onrfStringEmission(QString);
    void onrfTaskFinished();
    void onrfErrorEmission(QString);

private:
    Ui::TrainWindow *ui;
    QString last_dir;
    LinearRegression *linear_regression;
    RandomForest *random_forest;
    bool lrIsInitiated;
    bool rfIsInitiated;
    void changeWorkStyle1();
    void changeWorkStyle2();

signals:
    void getTrainWindowClosed();

private slots:
    void on_lrTrainButton_pressed();
    void on_lrTrainButton_released();
    void on_lrOutButton_pressed();
    void on_lrOutButton_released();
    void on_lrSubmitButton_pressed();
    void on_lrSubmitButton_released();
    void on_rfTrainButton_pressed();
    void on_rfTrainButton_released();
    void on_rfOutButton_pressed();
    void on_rfOutButton_released();
    void on_rfDefaultCheck_clicked(bool checked);
    void on_rfSubmitButton_pressed();
    void on_rfSubmitButton_released();
};

#endif // TRAINWINDOW_H
