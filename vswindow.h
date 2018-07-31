#ifndef VSWINDOW_H
#define VSWINDOW_H

#include "scorevs.h"

namespace Ui {
class VSWindow;
}

class VSWindow : public QDialog
{
    Q_OBJECT

public:
    explicit VSWindow(QWidget *parent = 0);
    ~VSWindow();

public slots:
    void onvsStringEmission(QString);
    void onvsTaskFinished();
    void onvsErrorEmission(QString);

private:
    Ui::VSWindow *ui;
    QString last_dir;
    ScoreVS *score_vs;
    bool vsIsInitiated;
    void changeWorkStyle();

signals:
    void getVSWindowClosed();

private slots:
    void on_folderCheckBox_clicked(bool checked);
    void on_typeComboBox_currentIndexChanged(int index);
    void on_recPushButton_pressed();
    void on_recPushButton_released();
    void on_ligPushButton_pressed();
    void on_ligPushButton_released();
    void on_lrrfPushButton_pressed();
    void on_lrrfPushButton_released();
    void on_scoreSubmitButton_pressed();
    void on_scoreSubmitButton_released();
    void on_logPushButton_pressed();
    void on_logPushButton_released();
};

#endif // VSWINDOW_H
