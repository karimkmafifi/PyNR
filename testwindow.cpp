#include "testwindow.h"
#include "ui_testwindow.h"

TestWindow::TestWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TestWindow)
{
    ui->setupUi(this);

    TestWindow::last_dir = QDir::currentPath();

    connect(this, SIGNAL(getTestWindowClosed()),
                this->parent(), SLOT(onGetTestWindowClosed()));

    TestWindow::changeWorkStyle();

}

TestWindow::~TestWindow()
{

    emit TestWindow::getTestWindowClosed();

    delete ui;
}

void TestWindow::changeWorkStyle()
{
    this->ui->expPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->testPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->logPushButton->setStyleSheet("color: rgb(10, 102, 204); border: 1px solid rgb(0, 0, 0);");
    this->ui->expLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->testLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->logLineEdit->setStyleSheet("border: 1px solid rgb(0, 0, 0); background-color: rgb(255, 255, 255); color: rgb(0, 0, 0);");
    this->ui->expLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->testLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->logLabel->setStyleSheet("font: bold 8.5pt; color: rgb(0, 0, 0);");
    this->ui->testSubmitButton->setStyleSheet("color: rgb(0, 0, 0); border: 1px solid rgb(0,0,0); background-color: rgb(50, 150, 0);");
}

void TestWindow::on_expPushButton_pressed()
{
    this->ui->expPushButton->setFont(QFont(this->ui->expPushButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("PDBbind Refined Data File"), TestWindow::last_dir, tr("All Files(*.*)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TestWindow::last_dir = QFileInfo(path).path();
        this->ui->expLineEdit->setText(path);
    }
}

void TestWindow::on_expPushButton_released()
{
    this->ui->expPushButton->setFont(QFont(this->ui->expPushButton->font().family(), 8));
}

void TestWindow::on_testPushButton_pressed()
{
    this->ui->testPushButton->setFont(QFont(this->ui->testPushButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("Test Data File"), TestWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TestWindow::last_dir = QFileInfo(path).path();
        this->ui->testLineEdit->setText(path);
    }
}

void TestWindow::on_testPushButton_released()
{
    this->ui->testPushButton->setFont(QFont(this->ui->testPushButton->font().family(), 8));
}

void TestWindow::on_logPushButton_pressed()
{
    this->ui->logPushButton->setFont(QFont(this->ui->logPushButton->font().family(), 7));

    QString path = "";

    path = QFileDialog::getOpenFileName(this, tr("Output File"), TestWindow::last_dir, tr("txt Files(*.txt);;csv Files(*.csv)"), 0, QFileDialog::DontUseNativeDialog);

    if(path != "")
    {
        TestWindow::last_dir = QFileInfo(path).path();
        this->ui->logLineEdit->setText(path);
    }
}

void TestWindow::on_logPushButton_released()
{
    this->ui->logPushButton->setFont(QFont(this->ui->logPushButton->font().family(), 8));
}

void TestWindow::ontestStringEmission(QString log)
{
    this->ui->testLogViewer->setText(log);
}

void TestWindow::ontestTaskFinished()
{
    if(testIsInitiated)
    {
        test_pdb_bind->stopThread = true;
        test_pdb_bind->quit();

        if(!test_pdb_bind->wait(500))
        {
            test_pdb_bind->terminate();
            test_pdb_bind->wait();
        }

        delete test_pdb_bind;
    }

    testIsInitiated = false;

    this->ui->frame->setEnabled(true);

    TestWindow::changeWorkStyle();
}

void TestWindow::ontestErrorEmission(QString err)
{
    QMessageBox::critical(this, "Test Error: ", err, QMessageBox::Ok);
}

void TestWindow::on_testSubmitButton_pressed()
{
    this->ui->testSubmitButton->setFont(QFont(this->ui->testSubmitButton->font().family(), 7, QFont::Bold));

    this->ui->testLogViewer->setText("");

    bool correctInput = true;

    std::string expDataFilePath = this->ui->expLineEdit->text().toStdString();
    std::string testDataFilePath = this->ui->testLineEdit->text().toStdString();
    std::string logFilePath = this->ui->logLineEdit->text().toStdString();

    Common::rmvBackslash(expDataFilePath);
    Common::rmvBackslash(testDataFilePath);
    Common::rmvBackslash(logFilePath);

    Common::bothTrim(expDataFilePath);
    Common::bothTrim(testDataFilePath);
    Common::bothTrim(logFilePath);

    if(expDataFilePath == "" || testDataFilePath == "" || logFilePath == "")
    {
        QMessageBox::critical(this, "Test Window Error: ", "All entries should not be blank", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QFileInfo::exists(QString::fromStdString(expDataFilePath)))
    {
        QMessageBox::critical(this, "Test Window Error: ", "PDBbind Refined data file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    if(!QFileInfo::exists(QString::fromStdString(testDataFilePath)))
    {
        QMessageBox::critical(this, "Test Window Error: ", "Test data file does not exist", QMessageBox::Ok);
        correctInput = false;
        return;
    }

    std::vector<std::string> elems;
    Common::split(logFilePath, '/', elems);

    std::vector<std::string> elems2;

    Common::split(elems[elems.size() - 1], '.', elems2);

    if(elems2.size() > 1)
    {
        if(elems2[1] != "txt" && elems2[1] != "csv")
        {
            QMessageBox::critical(this, "Test Window Error: ", "Output File should only be a txt or a csv file", QMessageBox::Ok);
            correctInput = false;
            return;
        }
    }
    else
    {
        QMessageBox::critical(this, "Test Window Error: ", "Error finding Output File", QMessageBox::Ok);
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

        this->ui->expPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->testPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->logPushButton->setStyleSheet("color: rgb(20, 129, 255); border: 1px solid rgb(150,150,150);");
        this->ui->expLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->testLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->logLineEdit->setStyleSheet("border: 1px solid rgb(150,150,150); background-color: rgb(255, 255, 255); color: rgb(150,150,150);");
        this->ui->expLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->testLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->logLabel->setStyleSheet("font: bold 8.5pt; color: rgb(150, 150, 150);");
        this->ui->testSubmitButton->setStyleSheet("color: rgb(150, 150, 150); border: 1px solid rgb(150,150,150); background-color: rgb(64, 193, 0);");

        test_pdb_bind = new TestPDBbind(this, expDataFilePath, testDataFilePath, logFilePath);

        connect(test_pdb_bind, SIGNAL(testStringEmission(QString)),
                    this, SLOT(ontestStringEmission(QString)));

        connect(test_pdb_bind, SIGNAL(testTaskFinished()),
                    this, SLOT(ontestTaskFinished()));

        connect(test_pdb_bind, SIGNAL(testErrorEmission(QString)),
                    this, SLOT(ontestErrorEmission(QString)));

        testIsInitiated = true;

        if(num_threads > 1)
        {
            test_pdb_bind->start();
        }
        else
        {
            test_pdb_bind->run();
        }
    }
}

void TestWindow::on_testSubmitButton_released()
{
    this->ui->testSubmitButton->setFont(QFont(this->ui->testSubmitButton->font().family(), 8, QFont::Bold));
}
