#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H

#include "score.h"

class LinearRegression : public QThread
{
    Q_OBJECT

private:
	std::string trainingFile;
    std::string logFile;
	unsigned int scoringFunctionType;
    std::vector<float>cost_function_terms;
    std::vector<std::pair<std::string, float>> minimum_weights;
    float learningRate;
    unsigned int noIterations;
    float min_squared_error;
	bool checkMSEIncrease;
    unsigned int num_vars;

public:
    bool stopThread;

public:
    LinearRegression(QObject *parent);
    LinearRegression();
    void setLRTrainVars(std::string trainingFile_, std::string outFile_, float learningRate_, unsigned int noIterations_, bool checkMSEIncrease_);
    void run();
    void saveWeights();
    float CalcBindingAffinity(std::vector<std::pair<std::string, float>>& terms_contributions, unsigned int scoringFunctionType_);
    void loadWeightsFile(std::string inFile_);
    ~LinearRegression();

signals:
    void lrStringEmission(QString log);
    void lrTaskFinished();
    void lrErrorEmission(QString err);

};

#endif
