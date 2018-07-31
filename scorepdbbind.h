#ifndef SCOREPDBBIND_H
#define SCOREPDBBIND_H

#include "linearRegression.h"
#include "randomforest.h"

class ScorePDBbind : public QThread
{
    Q_OBJECT

private:
    std::string refinedExpDataPath;
    std::string complexesParentDir;
    std::string logFilePath;
    unsigned int sf_type;

public:
    bool stopThread;

public:
    ScorePDBbind(QObject *parent, std::string refinedExpDataPath_, std::string complexesParentDir_, std::string logFilePath_, unsigned int sf_type_);
    void run();

signals:
    void stringEmission(QString log);
    void taskFinished();
    void errorEmission(QString err);
};

#endif // SCOREPDBBIND_H
