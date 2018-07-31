#ifndef GETXSWINDOW_H
#define GETXSWINDOW_H

#include "scorepdbbind.h"

namespace Ui {
class GetXsWindow;
}

class GetXsWindow : public QDialog
{
    Q_OBJECT

public:
    explicit GetXsWindow(QWidget *parent = 0);
    ~GetXsWindow();

public slots:
    void onStringEmission(QString);
    void onTaskFinished();
    void onErrorEmission(QString);

private slots:
    void on_expDataButton_pressed();

    void on_expDataButton_released();

    void on_compParentButton_pressed();

    void on_compParentButton_released();

    void on_outXsButton_pressed();

    void on_outXsButton_released();

    void on_getXsStartButton_pressed();

    void on_getXsStartButton_released();

signals:
    void getXsWindowClosed();

private:
    Ui::GetXsWindow *ui;
    QString last_dir;
    ScorePDBbind *score_pdb_bind;
    bool isInitiated;
    void changeWorkStyle();
};

#endif // GETXSWINDOW_H
