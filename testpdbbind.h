#ifndef TESTPDBBIND_H
#define TESTPDBBIND_H

#include "linearRegression.h"
#include "randomforest.h"

class TestPDBbind : public QThread
{
    Q_OBJECT

private:
    std::string refinedExpDataPath;
    std::string testDataPath;
    std::string logFilePath;

public:
    bool stopThread;

public:
    TestPDBbind(QObject *parent, std::string refinedExpDataPath_, std::string testDataPath_, std::string logFilePath_);
    void run();

signals:
    void testStringEmission(QString log);
    void testTaskFinished();
    void testErrorEmission(QString err);

};

#endif // TESTPDBBIND_H
