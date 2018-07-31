#ifndef SCOREVS_H
#define SCOREVS_H

#include "linearRegression.h"
#include "randomforest.h"

class ScoreVS : public QThread
{
    Q_OBJECT

private:
    std::string recOrDirPath;
    std::vector<std::string> ligandsPaths;
    std::string lr_rf_FilePath;
    std::string logFilePath;
    bool isVS;
    unsigned int sfType;

public:
    bool stopThread;

public:
    ScoreVS(QObject *parent, std::string recOrDirPath_, std::vector<std::string> ligandsPaths_, std::string lr_rf_FilePath_, std::string logFilePath_, bool isVS_, unsigned int sfType_);
    void run();

signals:
    void vsStringEmission(QString log);
    void vsTaskFinished();
    void vsErrorEmission(QString err);
};

#endif // SCOREVS_H
