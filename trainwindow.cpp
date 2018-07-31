#include "trainwindow.h"
#include "ui_trainwindow.h"

TrainWindow::TrainWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TrainWindow), lrIsInitiated(false)
{
    ui->setupUi(this);

    TrainWindow::last_dir = QDir::currentPath();

    connect(this, SIGNAL(getTrainWindowClosed()),
                this->parent(), SLOT(onGetTrainWindowClosed()));

    TrainWindow::changeWorkStyle1();
    TrainWindow::changeWorkStyle2();
}


TrainWindow::~TrainWindow()
{
    if(lrIsInitiated)
    {
        linear_regression->saveWeights();
        linear_regression->stopThread = true;
        linear_regression->quit();

        if(!linear_regression->wait(500))
        {
            linear_regression->terminate();
            linear_regression->wait();
        }

        delete linear_regression;
    }

    lrIsInitiated = false;

    if(rfIsInitiated)
    {
        random_forest->save_forest();
        random_forest->stopThread = true;
        random_forest->quit();

        if(!random_forest->wait(500))
        {
            random_forest->terminate();
            random_forest->wait();
        }

        delete random_forest;
    }

    rfIsInitiated = false;

    emit TrainWindow::getTrainWindowClosed();

    delete ui;
}

void TrainWindow::changeWorkStyle1()
{
    this->ui->lrTrainButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->lrOutButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->lrTrainLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->lrOutLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->lrLearningRate->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->lrNoIteration->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->label->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_2->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_3->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_4->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->lrCheckMSE->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->lrSubmitButton->setStyleSheet("color: rgb(0, 0, 0); border: 1px solid rgb(0,0,0); background-color: rgb(50, 150, 0);");
}

void TrainWindow::changeWorkStyle2()
{
    this->ui->rfTrainButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->rfOutButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->rfTrainLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfOutLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfDefaultCheck->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->rfSubmitButton->setStyleSheet("color: rgb(0, 0, 0); border: 1px solid rgb(0,0,0); background-color: rgb(50, 150, 0);");
    this->ui->rfLogViewer->setText("");

    this->ui->frame_16->setEnabled(true);
    this->ui->label_7->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_8->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_9->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_10->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->label_11->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->rfNoTrees->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfMaxDepth->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfNoBootStrap->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfMtry->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->rfIntMinSplit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
}

void TrainWindow::on_lrTrainButton_pressed()
{
    this->ui->lrTrainButton->setFont(QFont(this->ui->lrTrainButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("LR Training File"), TrainWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TrainWindow::last_dir = QFileInfo(path).path();
        this->ui->lrTrainLineEdit->setText(path);
    }
}

void TrainWindow::on_lrTrainButton_released()
{
    this->ui->lrTrainButton->setFont(QFont(this->ui->lrTrainButton->font().family(), 8));
}

void TrainWindow::on_lrOutButton_pressed()
{
    this->ui->lrOutButton->setFont(QFont(this->ui->lrOutButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("LR Out File"), TrainWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TrainWindow::last_dir = QFileInfo(path).path();
        this->ui->lrOutLineEdit->setText(path);
    }
}

void TrainWindow::on_lrOutButton_released()
{
    this->ui->lrOutButton->setFont(QFont(this->ui->lrOutButton->font().family(), 8));
}

void TrainWindow::onlrStringEmission(QString log)
{
    this->ui->lrLogViewer->setText(log);
}

void TrainWindow::onlrTaskFinished()
{
    if(lrIsInitiated)
    {
        linear_regression->saveWeights();
        linear_regression->stopThread = true;
        linear_regression->quit();

        if(!linear_regression->wait(500))
        {
            linear_regression->terminate();
            linear_regression->wait();
        }

        delete linear_regression;
    }

    lrIsInitiated = false;

    this->ui->tabWidget->setEnabled(true);

    TrainWindow::changeWorkStyle1();
    this->ui->lrLogViewer->setText("");
}

void TrainWindow::onlrErrorEmission(QString err)
{
    QMessageBox::critical(this, "Linear Regression Error: ", err, QMessageBox::Ok);
}

void TrainWindow::onrfStringEmission(QString log)
{
     this->ui->rfLogViewer->setText(log);
}

void TrainWindow::onrfTaskFinished()
{
    if(rfIsInitiated)
    {
        random_forest->save_forest();
        random_forest->stopThread = true;
        random_forest->quit();

        if(!random_forest->wait(500))
        {
            random_forest->terminate();
            random_forest->wait();
        }

        delete random_forest;
    }

    rfIsInitiated = false;

    this->ui->tabWidget->setEnabled(true);

    TrainWindow::changeWorkStyle2();

    this->ui->rfLogViewer->setText("");
}

void TrainWindow::onrfErrorEmission(QString err)
{
    QMessageBox::critical(this, "Random Forest Error: ", err, QMessageBox::Ok);
}

void TrainWindow::on_lrSubmitButton_pressed()
{
    this->ui->lrSubmitButton->setFont(QFont(this->ui->lrSubmitButton->font().family(), 7, QFont::Bold));

    bool correctInput = true;

    std::string lrTrainingFilePath = this->ui->lrTrainLineEdit->text().toStdString();
    std::string lrOutFilePath = this->ui->lrOutLineEdit->text().toStdString();
    std::string lrLearnRate = this->ui->lrLearningRate->text().toStdString();
    std::string lrnumIterations = this->ui->lrNoIteration->text().toStdString();

    Common::rmvBackslash(lrTrainingFilePath);
    Common::rmvBackslash(lrOutFilePath);

    Common::bothTrim(lrTrainingFilePath);
    Common::bothTrim(lrOutFilePath);
    Common::bothTrim(lrLearnRate);
    Common::bothTrim(lrnumIterations);

    if(lrTrainingFilePath == "" || lrOutFilePath == "" || lrLearnRate == "" || lrnumIterations == "")
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "All entries should not be blank", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QFileInfo::exists(QString::fromStdString(lrTrainingFilePath)))
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "Training file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if (lrnumIterations.find(".") != std::string::npos)
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "Number of Iterations value should only be an integer.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    float learningRate;
    int numIterations;

    try
    {
        learningRate = Common::string_to<float>(lrLearnRate);
        numIterations = Common::string_to<int>(lrnumIterations);
    }
    catch (...)
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "Error while converting Learning Rate to float and Number of Iterations to int.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(learningRate <= 0 || numIterations <= 0)
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "Learing Rate and Number of Iterations should be bigger than 0.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    std::vector<std::string> elems;
    Common::split(lrOutFilePath, '/', elems);

    std::vector<std::string> elems2;

    Common::split(elems[elems.size() - 1], '.', elems2);

    if(elems2.size() > 1)
    {
        if(elems2[1] != "txt" && elems2[1] != "csv")
        {
            QMessageBox::critical(this, "Linear Regression Window Error: ", "Linear Regression Output File should only be a txt or a csv file", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        QMessageBox::critical(this, "Linear Regression Window Error: ", "Error finding Linear Regression Output File", QMessageBox::Ok);
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

        this->ui->tabWidget->setEnabled(false);

        this->ui->lrTrainButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->lrOutButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->lrTrainLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->lrOutLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->lrLearningRate->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->lrNoIteration->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->label->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_2->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_3->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_4->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->lrCheckMSE->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->lrSubmitButton->setStyleSheet("color: rgb(150, 150, 150); border: 1px solid rgb(150,150,150); background-color: rgb(64, 193, 0);");

        bool mseChecked = this->ui->lrCheckMSE->isChecked();

        linear_regression = new LinearRegression(this);
        linear_regression->setLRTrainVars(lrTrainingFilePath, lrOutFilePath, learningRate, numIterations, mseChecked);

        connect(linear_regression, SIGNAL(lrStringEmission(QString)),
                    this, SLOT(onlrStringEmission(QString)));

        connect(linear_regression, SIGNAL(lrTaskFinished()),
                    this, SLOT(onlrTaskFinished()));

        connect(linear_regression, SIGNAL(lrErrorEmission(QString)),
                    this, SLOT(onlrErrorEmission(QString)));

        lrIsInitiated = true;

        if(num_threads > 1)
        {
            linear_regression->start();
        }
        else
        {
            linear_regression->run();
        }
    }
}

void TrainWindow::on_lrSubmitButton_released()
{
    this->ui->lrSubmitButton->setFont(QFont(this->ui->lrSubmitButton->font().family(), 8, QFont::Bold));
}

void TrainWindow::on_rfTrainButton_pressed()
{
    this->ui->rfTrainButton->setFont(QFont(this->ui->rfTrainButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("RF Training File"), TrainWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TrainWindow::last_dir = QFileInfo(path).path();
        this->ui->rfTrainLineEdit->setText(path);
    }
}

void TrainWindow::on_rfTrainButton_released()
{
    this->ui->rfTrainButton->setFont(QFont(this->ui->rfTrainButton->font().family(), 8));
}

void TrainWindow::on_rfOutButton_pressed()
{
    this->ui->rfOutButton->setFont(QFont(this->ui->rfOutButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("RF Out File"), TrainWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TrainWindow::last_dir = QFileInfo(path).path();
        this->ui->rfOutLineEdit->setText(path);
    }
}

void TrainWindow::on_rfOutButton_released()
{
    this->ui->rfOutButton->setFont(QFont(this->ui->rfOutButton->font().family(), 8));
}

void TrainWindow::on_rfDefaultCheck_clicked(bool checked)
{
    if(checked)
    {
        this->ui->frame_16->setEnabled(false);
        this->ui->label_7->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_8->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_9->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_10->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_11->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->rfNoTrees->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfMaxDepth->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfNoBootStrap->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfMtry->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfIntMinSplit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
    }
    else
    {
        this->ui->frame_16->setEnabled(true);
        this->ui->label_7->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->label_8->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->label_9->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->label_10->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->label_11->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
        this->ui->rfNoTrees->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
        this->ui->rfMaxDepth->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
        this->ui->rfNoBootStrap->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
        this->ui->rfMtry->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
        this->ui->rfIntMinSplit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");

    }
}

void TrainWindow::on_rfSubmitButton_pressed()
{
    this->ui->rfSubmitButton->setFont(QFont(this->ui->rfSubmitButton->font().family(), 7, QFont::Bold));

    bool correctInput = true;

    std::string rfTrainingFilePath = this->ui->rfTrainLineEdit->text().toStdString();
    std::string rfOutFilePath = this->ui->rfOutLineEdit->text().toStdString();

    std::string rf_num_trees;
    std::string rf_max_depth;
    std::string rf_bootstrap_percent;
    std::string rf_mtry;
    std::string rf_intern_min_split;

    if (this->ui->rfDefaultCheck->checkState() == Qt::Checked)
    {
        rf_num_trees = "500";
        rf_max_depth = "9999999";
        rf_bootstrap_percent = "100";
        rf_mtry = "2";
        rf_intern_min_split = "1";
    }
    else
    {
        rf_num_trees = this->ui->rfNoTrees->text().toStdString();
        rf_max_depth = this->ui->rfMaxDepth->text().toStdString();
        rf_bootstrap_percent = this->ui->rfNoBootStrap->text().toStdString();
        rf_mtry = this->ui->rfMtry->text().toStdString();
        rf_intern_min_split = this->ui->rfIntMinSplit->text().toStdString();

        Common::bothTrim(rf_num_trees);
        Common::bothTrim(rf_max_depth);
        Common::bothTrim(rf_bootstrap_percent);
        Common::bothTrim(rf_mtry);
        Common::bothTrim(rf_intern_min_split);
    }

    Common::rmvBackslash(rfTrainingFilePath);
    Common::rmvBackslash(rfOutFilePath);

    Common::bothTrim(rfTrainingFilePath);
    Common::bothTrim(rfOutFilePath);

    if(rfTrainingFilePath == "" || rfOutFilePath == "" || rf_num_trees == "" || rf_max_depth == ""
            || rf_bootstrap_percent == "" || rf_mtry == "" || rf_intern_min_split == "")
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "All entries should not be blank", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QFileInfo::exists(QString::fromStdString(rfTrainingFilePath)))
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Training file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if (rf_num_trees.find(".") != std::string::npos)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Number of Trees value should only be an integer.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if (rf_max_depth.find(".") != std::string::npos)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Max Depth value should only be an integer.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if (rf_mtry.find(".") != std::string::npos)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Mtry value should only be an integer.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if (rf_intern_min_split.find(".") != std::string::npos)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Internal Min Split value should only be an integer.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    int numTrees;
    int maxDepth;
    float bootstrapPercent;
    int mTry;
    int internMinSplit;

    try
    {
        numTrees = Common::string_to<int>(rf_num_trees);
        maxDepth = Common::string_to<int>(rf_max_depth);
        bootstrapPercent = Common::string_to<float>(rf_bootstrap_percent);
        mTry = Common::string_to<int>(rf_mtry);
        internMinSplit = Common::string_to<int>(rf_intern_min_split);
    }
    catch (...)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Error while converting Random Forest parameters to numbers.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(numTrees <= 0 || maxDepth <= 0 || mTry <= 0 || internMinSplit <= 0)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Num of Trees, Max Depth, Mtry and Internal Min. Split should be bigger than 0.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(bootstrapPercent <= 0 || bootstrapPercent > 100)
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Percentage for Bootstrap Data should be > 0 && <= 100.", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    std::vector<std::string> elems;
    Common::split(rfOutFilePath, '/', elems);

    std::vector<std::string> elems2;

    Common::split(elems[elems.size() - 1], '.', elems2);

    if(elems2.size() > 1)
    {
        if(elems2[1] != "txt" && elems2[1] != "csv")
        {
            QMessageBox::critical(this, "Random Forest Window Error: ", "Random Forest Output File should only be a txt or a csv file", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        QMessageBox::critical(this, "Random Forest Window Error: ", "Error finding Random Forest Output File", QMessageBox::Ok);
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

        this->ui->tabWidget->setEnabled(false);

        this->ui->rfTrainButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->rfOutButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->rfTrainLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfOutLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfDefaultCheck->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->rfSubmitButton->setStyleSheet("color: rgb(150, 150, 150); border: 1px solid rgb(150,150,150); background-color: rgb(64, 193, 0);");

        this->ui->frame_16->setEnabled(false);
        this->ui->label_7->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_8->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_9->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_10->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->label_11->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->rfNoTrees->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfMaxDepth->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfNoBootStrap->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfMtry->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->rfIntMinSplit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");

        random_forest = new RandomForest(this);
        random_forest->setRFTrainVars(rfTrainingFilePath, rfOutFilePath, numTrees, maxDepth, bootstrapPercent, mTry, internMinSplit);

        connect(random_forest, SIGNAL(rfStringEmission(QString)),
                    this, SLOT(onrfStringEmission(QString)));

        connect(random_forest, SIGNAL(rfTaskFinished()),
                    this, SLOT(onrfTaskFinished()));

        connect(random_forest, SIGNAL(rfErrorEmission(QString)),
                    this, SLOT(onrfErrorEmission(QString)));

        rfIsInitiated = true;

        if(num_threads > 1)
        {
            random_forest->start();
        }
        else
        {
            random_forest->run();
        }
    }
}

void TrainWindow::on_rfSubmitButton_released()
{
    this->ui->rfSubmitButton->setFont(QFont(this->ui->rfSubmitButton->font().family(), 8, QFont::Bold));
}
