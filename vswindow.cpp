#include "vswindow.h"
#include "ui_vswindow.h"

VSWindow::VSWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VSWindow)
{
    ui->setupUi(this);

    VSWindow::last_dir = QDir::currentPath();

    connect(this, SIGNAL(getVSWindowClosed()),
                this->parent(), SLOT(onGetVSWindowClosed()));

    VSWindow::changeWorkStyle();
}

VSWindow::~VSWindow()
{
    emit VSWindow::getVSWindowClosed();

    delete ui;
}

void VSWindow::changeWorkStyle()
{
    this->ui->recPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->ligPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->lrrfPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->logPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->recLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->ligLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->lrrfLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->logLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->recLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->ligLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->lrrfLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->logLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->folderCheckBox->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->typeComboBox->setStyleSheet("color: rgb(0, 0, 0); font: 9pt; selection-background-color: rgb(10, 102, 204); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
    this->ui->sfTypeComboBox->setStyleSheet("color: rgb(0, 0, 0); font: 9pt; selection-background-color: rgb(10, 102, 204); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
    this->ui->scoreSubmitButton->setStyleSheet("color: rgb(0, 0, 0); border: 1px solid rgb(0,0,0); background-color: rgb(50, 150, 0);");
}

void VSWindow::onvsTaskFinished()
{
    if(vsIsInitiated)
    {
        score_vs->stopThread = true;
        score_vs->quit();

        if(!score_vs->wait(500))
        {
            score_vs->terminate();
            score_vs->wait();
        }

        delete score_vs;
    }

    vsIsInitiated = false;

    this->ui->frame_5->setEnabled(true);

    VSWindow::changeWorkStyle();
}

void VSWindow::onvsErrorEmission(QString err)
{
    QMessageBox::critical(this, "Score Error: ", err, QMessageBox::Ok);
}

void VSWindow::onvsStringEmission(QString log)
{
    this->ui->scoreLogViewer->moveCursor(QTextCursor::End);
    this->ui->scoreLogViewer->insertPlainText(log);
}

void VSWindow::on_folderCheckBox_clicked(bool checked)
{
    if(checked)
    {
        this->ui->recLabel->setText("Complexes Dir:");
        this->ui->recLineEdit->setText("");
        this->ui->frame_13->setEnabled(false);
        this->ui->ligLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->ligLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->ligPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");

        this->ui->ligLineEdit->setText("");
    }
    else
    {
        this->ui->recLabel->setText("Receptor Path:");
        this->ui->recLineEdit->setText("");
        this->ui->frame_13->setEnabled(true);
        this->ui->ligLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->ligLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
        this->ui->ligPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    }
}

void VSWindow::on_typeComboBox_currentIndexChanged(int index)
{
    if(index == 0)
    {
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->addItem("X-Score HC");
        this->ui->sfTypeComboBox->addItem("VINA");
        this->ui->sfTypeComboBox->addItem("VINA-Halogen");
        this->ui->sfTypeComboBox->addItem("AutoDock");
    }
    else if(index == 1)
    {
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->removeItem(0);
        this->ui->sfTypeComboBox->addItem("RF-SCORE");
        this->ui->sfTypeComboBox->addItem("RF-X-SCORE HC");
        this->ui->sfTypeComboBox->addItem("RF-VINA");
        this->ui->sfTypeComboBox->addItem("RF-VINA-Halogen");
        this->ui->sfTypeComboBox->addItem("RF-AutoDock");
    }
}

void VSWindow::on_recPushButton_pressed()
{
    this->ui->recPushButton->setFont(QFont(this->ui->recPushButton->font().family(), 7));

    QString path = "";

    if (this->ui->folderCheckBox->checkState() == Qt::Checked)
    {
        path = QFileDialog::getExistingDirectory(this, tr("Complexes Parent Dir"), VSWindow::last_dir, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);

        if(path != "")
        {
            VSWindow::last_dir = QFileInfo(path).path();
            this->ui->recLineEdit->setText(path);
        }
    }
    else
    {
        path = QFileDialog::getOpenFileName(this, tr("Receptor File"), VSWindow::last_dir, tr("pdbqt Files(*.pdbqt);;pdb Files(*.pdb)"), 0, QFileDialog::DontUseNativeDialog);

        if(path != "")
        {
            VSWindow::last_dir = QFileInfo(path).path();
            this->ui->recLineEdit->setText(path);
        }
    }
}

void VSWindow::on_recPushButton_released()
{
    this->ui->recPushButton->setFont(QFont(this->ui->recPushButton->font().family(), 8));
}

void VSWindow::on_ligPushButton_pressed()
{
    this->ui->ligPushButton->setFont(QFont(this->ui->ligPushButton->font().family(), 7));

    QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Ligands Files"), VSWindow::last_dir, tr("pdbqt Files(*.pdbqt);;sdf Files(*.sdf)"), 0, QFileDialog::DontUseNativeDialog);

    if(!filenames.isEmpty())
    {
        QString pathsToView = "";

        VSWindow::last_dir = QFileInfo(filenames[0]).path();

        for (int i = 0; i < filenames.count(); ++i)
        {
            pathsToView.append(filenames[i]);

            if(i < filenames.count() - 1)
            {
                pathsToView.append("|");
            }
        }

        this->ui->ligLineEdit->setText(pathsToView);
    }
}

void VSWindow::on_ligPushButton_released()
{
    this->ui->ligPushButton->setFont(QFont(this->ui->ligPushButton->font().family(), 8));
}

void VSWindow::on_lrrfPushButton_pressed()
{
    this->ui->lrrfPushButton->setFont(QFont(this->ui->lrrfPushButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("LR Weights or RF Trees File"), VSWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        VSWindow::last_dir = QFileInfo(path).path();
        this->ui->lrrfLineEdit->setText(path);
    }
}

void VSWindow::on_lrrfPushButton_released()
{
    this->ui->lrrfPushButton->setFont(QFont(this->ui->lrrfPushButton->font().family(), 8));
}

void VSWindow::on_scoreSubmitButton_pressed()
{
    this->ui->scoreSubmitButton->setFont(QFont(this->ui->scoreSubmitButton->font().family(), 7, QFont::Bold));
    this->ui->scoreLogViewer->setText("");

    bool correctInput = true;

    std::string recFilePath = this->ui->recLineEdit->text().toStdString();
    std::string ligFilesPaths = this->ui->ligLineEdit->text().toStdString();
    std::string lrrfFilePath = this->ui->lrrfLineEdit->text().toStdString();
    std::string outFilePath = this->ui->logLineEdit->text().toStdString();

    std::vector<std::string> ligandsPathsElems;

    Common::rmvBackslash(recFilePath);
    Common::rmvBackslash(ligFilesPaths);
    Common::rmvBackslash(lrrfFilePath);
    Common::rmvBackslash(outFilePath);

    Common::bothTrim(recFilePath);
    Common::bothTrim(ligFilesPaths);
    Common::bothTrim(lrrfFilePath);
    Common::bothTrim(outFilePath);

    if (this->ui->folderCheckBox->checkState() == Qt::Checked)
    {
        if(recFilePath == "" || lrrfFilePath == "" || outFilePath == "")
        {
            QMessageBox::critical(this, "Score Window Error: ", "All required entries should not be blank", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        if(recFilePath == "" || ligFilesPaths == "" || lrrfFilePath == "" || outFilePath == "")
        {
            QMessageBox::critical(this, "Score Window Error: ", "All required entries should not be blank", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }

    if (this->ui->folderCheckBox->checkState() == Qt::Checked)
    {
        if(!QDir(QString::fromStdString(recFilePath)).exists())
        {
            QMessageBox::critical(this, "Score Window Error: ", "Complexes Parent Dir does not exist", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        if(!QFileInfo::exists(QString::fromStdString(recFilePath)))
        {
            QMessageBox::critical(this, "Score Window Error: ", "Receptor file does not exist", QMessageBox::Ok);
            correctInput = false;
            return;
        }

        Common::split(ligFilesPaths, '|', ligandsPathsElems);

        for(unsigned int i = 0; i < ligandsPathsElems.size(); ++i)
        {
            if(!QFileInfo::exists(QString::fromStdString(ligandsPathsElems[i])))
            {
                QMessageBox::critical(this, "Score Window Error: ", "A Ligand file does not exist", QMessageBox::Ok);
                correctInput = false;
                return;
            }
        }
    }

    if(!QFileInfo::exists(QString::fromStdString(lrrfFilePath)))
    {
        QMessageBox::critical(this, "Score Window Error: ", "LR/RF file file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }


    std::vector<std::string> elems;
    Common::split(outFilePath, '/', elems);

    std::vector<std::string> elems2;

    Common::split(elems[elems.size() - 1], '.', elems2);

    if(elems2.size() > 1)
    {
        if(elems2[1] != "txt" && elems2[1] != "csv")
        {
            QMessageBox::critical(this, "Score Window Error: ", "Output File should only be a txt or a csv file", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        QMessageBox::critical(this, "Score Window Error: ", "Error finding Output File", QMessageBox::Ok);
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

        this->ui->frame_5->setEnabled(false);

        this->ui->recPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->ligPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->lrrfPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->logPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->recLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->ligLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->lrrfLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->logLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->recLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->ligLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->lrrfLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->logLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->folderCheckBox->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->typeComboBox->setStyleSheet("color: rgb(150, 150, 150); font: 9pt; selection-background-color: rgb(20, 129, 255); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
        this->ui->sfTypeComboBox->setStyleSheet("color: rgb(150, 150, 150); font: 9pt; selection-background-color: rgb(20, 129, 255); selection-color: rgb(255, 255, 255); background-color: rgb(255, 255, 255);");
        this->ui->scoreSubmitButton->setStyleSheet("color: rgb(150, 150, 150); border: 1px solid rgb(150,150,150); background-color: rgb(64, 193, 0);");

        unsigned int scoringFunctionType = 100; //Just for throwing an error later if the SF Type is not right somehow...

        if(this->ui->typeComboBox->currentIndex() == 0)
        {
            scoringFunctionType = this->ui->sfTypeComboBox->currentIndex();
        }
        else if(this->ui->typeComboBox->currentIndex() == 1)
        {
            scoringFunctionType = (this->ui->sfTypeComboBox->currentIndex()) + 4;
        }

        score_vs = new ScoreVS(this, recFilePath, ligandsPathsElems, lrrfFilePath, outFilePath, !(this->ui->folderCheckBox->checkState() == Qt::Checked), scoringFunctionType);

        connect(score_vs, SIGNAL(vsStringEmission(QString)),
                    this, SLOT(onvsStringEmission(QString)));

        connect(score_vs, SIGNAL(vsTaskFinished()),
                    this, SLOT(onvsTaskFinished()));

        connect(score_vs, SIGNAL(vsErrorEmission(QString)),
                    this, SLOT(onvsErrorEmission(QString)));

        vsIsInitiated = true;

        if(num_threads > 1)
        {
            score_vs->start();
        }
        else
        {
            score_vs->run();
        }
    }
}

void VSWindow::on_scoreSubmitButton_released()
{
    this->ui->scoreSubmitButton->setFont(QFont(this->ui->scoreSubmitButton->font().family(), 8, QFont::Bold));
}

void VSWindow::on_logPushButton_pressed()
{
    this->ui->logPushButton->setFont(QFont(this->ui->logPushButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("Output File"), VSWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        VSWindow::last_dir = QFileInfo(path).path();
        this->ui->logLineEdit->setText(path);
    }
}

void VSWindow::on_logPushButton_released()
{
    this->ui->logPushButton->setFont(QFont(this->ui->logPushButton->font().family(), 8));
}
