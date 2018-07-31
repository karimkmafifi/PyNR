#include "getxswindow.h"
#include "ui_getxswindow.h"

GetXsWindow::GetXsWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::GetXsWindow), isInitiated(false)
{
    ui->setupUi(this);

    GetXsWindow::last_dir = QDir::currentPath();

    connect(this, SIGNAL(getXsWindowClosed()),
                this->parent(), SLOT(onGetXsWindowClosed()));

    GetXsWindow::changeWorkStyle();

}

GetXsWindow::~GetXsWindow()
{
    if(isInitiated)
    {
        score_pdb_bind->stopThread = true;
        score_pdb_bind->quit();

        if(!score_pdb_bind->wait(500))
        {
            score_pdb_bind->terminate();
            score_pdb_bind->wait();
        }

        delete score_pdb_bind;
    }

    isInitiated = false;

    emit GetXsWindow::getXsWindowClosed();

    delete ui;
}

void GetXsWindow::changeWorkStyle()
{
    this->ui->expDataButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->compParentButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->outXsButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->getXslineEdit1->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->getXslineEdit2->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->getXslineEdit3->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->label->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_2->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_3->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_4->setStyleSheet("font: bold 10.5pt; color: rgb(0, 0, 0);");
    this->ui->getXscomboBox->setStyleSheet("color: rgb(0, 0, 0); font: 9pt; selection-background-color: rgb(10, 102, 204); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
    this->ui->getXsStartButton->setStyleSheet("color: rgb(0, 0, 0); border: 1px solid rgb(0,0,0); background-color: rgb(50, 150, 0);");
}

void GetXsWindow::on_expDataButton_pressed()
{
    this->ui->expDataButton->setFont(QFont(this->ui->expDataButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("PDBbind Refined Data File"), GetXsWindow::last_dir, tr("All Files(*.*)"), 0, QFileDialog::DontUseNativeDialog);


    if(path != "")
    {
        GetXsWindow::last_dir = QFileInfo(path).path();
        this->ui->getXslineEdit1->setText(path);
    }
}

void GetXsWindow::onErrorEmission(QString err)
{
    QMessageBox::critical(this, "Get Xs Error: ", err, QMessageBox::Ok);
}

void GetXsWindow::onTaskFinished()
{
    if(isInitiated)
    {
        score_pdb_bind->stopThread = true;
        score_pdb_bind->quit();

        if(!score_pdb_bind->wait(500))
        {
            score_pdb_bind->terminate();
            score_pdb_bind->wait();
        }

        delete score_pdb_bind;
    }

    isInitiated = false;

    this->ui->frame->setEnabled(true);

    GetXsWindow::changeWorkStyle();
    this->ui->getXsoutScreen->setText("");
}

void GetXsWindow::on_expDataButton_released()
{
    this->ui->expDataButton->setFont(QFont(this->ui->expDataButton->font().family(), 8));
}

void GetXsWindow::on_compParentButton_pressed()
{
    this->ui->compParentButton->setFont(QFont(this->ui->compParentButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getExistingDirectory(this, tr("Complexes Parent Dir"), GetXsWindow::last_dir, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        GetXsWindow::last_dir = QFileInfo(path).path();
        this->ui->getXslineEdit2->setText(path);
    }
}

void GetXsWindow::on_compParentButton_released()
{
    this->ui->compParentButton->setFont(QFont(this->ui->compParentButton->font().family(), 8));
}

void GetXsWindow::on_outXsButton_pressed()
{
    this->ui->outXsButton->setFont(QFont(this->ui->outXsButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("Xs Output File"), GetXsWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        GetXsWindow::last_dir = QFileInfo(path).path();
        this->ui->getXslineEdit3->setText(path);
    }
}

void GetXsWindow::on_outXsButton_released()
{
    this->ui->outXsButton->setFont(QFont(this->ui->outXsButton->font().family(), 8));
}

void GetXsWindow::onStringEmission(QString log)
{
    this->ui->getXsoutScreen->setText(log);
}

void GetXsWindow::on_getXsStartButton_pressed()
{
    this->ui->getXsStartButton->setFont(QFont(this->ui->getXsStartButton->font().family(), 7, QFont::Bold));

    bool correctInput = true;

    std::string expDataFilePath = this->ui->getXslineEdit1->text().toStdString();
    std::string complParentDirPath = this->ui->getXslineEdit2->text().toStdString();
    std::string outXsFilePath = this->ui->getXslineEdit3->text().toStdString();

    Common::rmvBackslash(expDataFilePath);
    Common::rmvBackslash(complParentDirPath);
    Common::rmvBackslash(outXsFilePath);

    Common::bothTrim(expDataFilePath);
    Common::bothTrim(complParentDirPath);
    Common::bothTrim(outXsFilePath);

    if(expDataFilePath == "" || complParentDirPath == "" || outXsFilePath == "")
    {
        QMessageBox::critical(this, "Get Xs Window Error: ", "All entries should not be blank", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QDir(QString::fromStdString(complParentDirPath)).exists())
    {
        QMessageBox::critical(this, "Get Xs Window Error: ", "Complexes Parent Dir does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QFileInfo::exists(QString::fromStdString(expDataFilePath)))
    {
        QMessageBox::critical(this, "Get Xs Window Error: ", "PDBbind Refined data file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    std::vector<std::string> elems;
    Common::split(outXsFilePath, '/', elems);

    std::vector<std::string> elems2;

    Common::split(elems[elems.size() - 1], '.', elems2);

    if(elems2.size() > 1)
    {
        if(elems2[1] != "txt" && elems2[1] != "csv")
        {
            QMessageBox::critical(this, "Get Xs Window Error: ", "Xs Output File should only be a txt or a csv file", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        QMessageBox::critical(this, "Get Xs Window Error: ", "Error finding Xs Output File", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(correctInput)
    {
        int num_threads = QThread::idealThreadCount();

        if(num_threads <= 1)
        {
            QMessageBox::StandardButton reply = QMessageBox::question(this, "WARNING", "No other thread available or can be detected. This will cause the GUI to be unresponsive while performing this task. Continue?", QMessageBox::Yes | QMessageBox::No);

            if (reply == QMessageBox::No)
            {
                return;
            }
        }

        this->ui->frame->setEnabled(false);

        this->ui->expDataButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->compParentButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->outXsButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->getXslineEdit1->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->getXslineEdit2->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->getXslineEdit3->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->label->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_2->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_3->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_4->setStyleSheet("font: bold 10.5pt; color: rgb(150, 150, 150);");
        this->ui->getXscomboBox->setStyleSheet("color: rgb(150, 150, 150); font: 9pt; selection-background-color: rgb(20, 129, 255); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
        this->ui->getXsStartButton->setStyleSheet("color: rgb(150, 150, 150); border: 1px solid rgb(150,150,150); background-color: rgb(64, 193, 0);");

        unsigned int scoringFunctionType = this->ui->getXscomboBox->currentIndex();

        score_pdb_bind = new ScorePDBbind(this, expDataFilePath, complParentDirPath, outXsFilePath, scoringFunctionType);

        connect(score_pdb_bind, SIGNAL(stringEmission(QString)),
                    this, SLOT(onStringEmission(QString)));

        connect(score_pdb_bind, SIGNAL(taskFinished()),
                    this, SLOT(onTaskFinished()));

        connect(score_pdb_bind, SIGNAL(errorEmission(QString)),
                    this, SLOT(onErrorEmission(QString)));

        isInitiated = true;

        if(num_threads > 1)
        {
            score_pdb_bind->start();
        }
        else
        {
            score_pdb_bind->run();
        }
    }
}

void GetXsWindow::on_getXsStartButton_released()
{
    this->ui->getXsStartButton->setFont(QFont(this->ui->getXsStartButton->font().family(), 8, QFont::Bold));
}
