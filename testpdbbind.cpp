#include "testpdbbind.h"

TestPDBbind::TestPDBbind(QObject *parent, std::string refinedExpDataPath_, std::string testDataPath_, std::string logFilePath_)
    : QThread(parent), stopThread(false), refinedExpDataPath(refinedExpDataPath_), testDataPath(testDataPath_), logFilePath(logFilePath_)
{

}

void TestPDBbind::run()
{
    try
    {
        emit TestPDBbind::testStringEmission("Reading From Refined Data File....");

        std::vector<std::pair<std::string, float>> receptor_experimental_pki;

        std::ifstream inputfile;
        inputfile.exceptions(std::ifstream::badbit);

        try
        {
            inputfile.open(TestPDBbind::refinedExpDataPath);
            std::string line;
            while (getline(inputfile, line))
            {
                if(this->stopThread)
                {
                    return;
                }

                bool found_ki_kd = false;
                std::string token;

                if (line.empty()) {}
                else if (Common::starts_with(line, "#")){}
                else if (line.find("Ki") != std::string::npos)
                {
                    found_ki_kd = true;
                    token = line.substr(0, line.find("Ki"));
                }
                else if (line.find("Kd") != std::string::npos)
                {
                    found_ki_kd = true;
                    token = line.substr(0, line.find("Kd"));
                }
                else {}

                if (found_ki_kd)
                {
                    std::string PDB_code;
                    float resolution;
                    int release_year;
                    float nlogKd_Ki;
                    std::stringstream input(token);

                    input >> PDB_code >> resolution >> release_year >> nlogKd_Ki;

                    receptor_experimental_pki.push_back(std::pair<std::string, float>(PDB_code, nlogKd_Ki));
                }
            }
        }
        catch (...)
        {
            inputfile.close();
            throw Error_report("Error reading from file " + TestPDBbind::refinedExpDataPath);
        }

        inputfile.close();

        emit TestPDBbind::testStringEmission("Reading From the Test Data File....");

        std::vector<std::pair<std::string, float>> receptor_calculated_pki;

        try
        {
            std::string line;
            unsigned int count = 0;

            inputfile.open(TestPDBbind::testDataPath);

            while (getline(inputfile, line))
            {
                if(this->stopThread)
                {
                    return;
                }

                if (count == 0)
                {
                    count = 1;

                    continue;
                }

                if (line.empty()) {}
                else
                {
                    std::vector<std::string> elems;

                    Common::split(line, '=', elems);

                    std::vector<std::string> elems2;

                    Common::split(elems[1], '|', elems2);

                    std::vector<std::string> elems3;

                    Common::split(elems2[1], ',', elems3);

                    std::string pkiValue;

                    for(int i = 0; i < elems3.size(); ++i)
                    {
                        if (elems3[i].find("Pki") != std::string::npos)
                        {
                            std::vector<std::string> elems4;
                            Common::split(elems3[i], ':', elems4);
                            pkiValue = elems4[1];
                            break;
                        }
                    }

                    elems2.clear();
                    elems3.clear();

                    Common::bothTrim(pkiValue);

                    Common::split(elems[0], ' ', elems2);

                    Common::split(elems2[0], '-', elems3);

                    receptor_calculated_pki.push_back(std::pair<std::string, float>(elems3[0], Common::string_to<float>(pkiValue)));
                }
            }

            inputfile.close();
        }
        catch (std::ifstream::failure e)
        {
            inputfile.close();
            throw Error_report(TestPDBbind::testDataPath + " " + e.what());
        }
        catch (Error_report& e) {
            inputfile.close();
            throw Error_report(TestPDBbind::testDataPath + " " + e.errorText);
        }
        catch (...)
        {
            inputfile.close();
            throw Error_report(TestPDBbind::testDataPath + " ");
        }

        emit TestPDBbind::testStringEmission("Calculating....");

        float PearsonCorrelationCof = 0.0;
        float standardDeviation = 0.0;

        float recepCalcAffinityMean = 0.0;
        float recepCalcAffinityMeanSqr = 0.0;
        float recepExpAffinityMean = 0.0;
        float recepCalcExpAffinityMean = 0.0;

        for (int i = 0; i < receptor_calculated_pki.size(); ++i)
        {
            recepCalcAffinityMean += receptor_calculated_pki[i].second;
            recepCalcAffinityMeanSqr += Common::sqr(receptor_calculated_pki[i].second);
        }

        recepCalcExpAffinityMean += recepCalcAffinityMean;

        recepCalcAffinityMean /= receptor_calculated_pki.size();
        recepCalcAffinityMeanSqr /= receptor_calculated_pki.size();

        for (int i = 0; i < receptor_experimental_pki.size(); ++i)
        {
            recepExpAffinityMean += receptor_experimental_pki[i].second;
        }

        recepCalcExpAffinityMean += recepExpAffinityMean;

        recepExpAffinityMean /= receptor_experimental_pki.size();

        recepCalcExpAffinityMean /= (receptor_calculated_pki.size() + receptor_experimental_pki.size());

        float upper_cc_term = 0.0;
        float lower_first_cc_term = 0.0;
        float lower_second_cc_term = 0.0;

        for (int i = 0; i < receptor_calculated_pki.size(); ++i)
        {
            for (int j = 0; j < receptor_experimental_pki.size(); ++j)
            {
                if (receptor_calculated_pki[i].first == receptor_experimental_pki[j].first)
                {
                    upper_cc_term += (receptor_calculated_pki[i].second - recepCalcAffinityMean) * (receptor_experimental_pki[j].second - recepExpAffinityMean);
                    lower_first_cc_term += Common::sqr(receptor_calculated_pki[i].second - recepCalcAffinityMean);
                    lower_second_cc_term += Common::sqr(receptor_experimental_pki[j].second - recepExpAffinityMean);
                }
            }
        }

        PearsonCorrelationCof = upper_cc_term / (std::sqrt(lower_first_cc_term) * std::sqrt(lower_second_cc_term));

        float b = ((recepCalcAffinityMean * recepExpAffinityMean) - recepCalcExpAffinityMean) / (Common::sqr(recepCalcAffinityMean) - recepCalcAffinityMeanSqr);
        float a = (recepExpAffinityMean - (b * recepCalcAffinityMean));

        float sd_upper_term = 0.0;

        for (int i = 0; i < receptor_calculated_pki.size(); ++i)
        {
            for (int j = 0; j < receptor_experimental_pki.size(); ++j)
            {
                if (receptor_calculated_pki[i].first == receptor_experimental_pki[j].first)
                {
                    sd_upper_term += Common::sqr((receptor_experimental_pki[j].second - (a + (b * receptor_calculated_pki[i].second))));
                }
            }
        }

        standardDeviation = std::sqrt(sd_upper_term / (receptor_calculated_pki.size() - 1));

        std::string logEmit = "Pearson Correlation: " + std::to_string(PearsonCorrelationCof) + ", Standard Deviation: " + std::to_string(standardDeviation);

        emit TestPDBbind::testStringEmission(QString::fromStdString(logEmit));

        std::ofstream logOS;
        logOS.open(TestPDBbind::logFilePath);

        logOS << logEmit;

        logOS.close();

    }
    catch (Error_report& err)
    {
        std::string tempStr = "Test Calculation Error: " + err.errorText;
        emit TestPDBbind::testErrorEmission(QString::fromStdString(tempStr));
        emit TestPDBbind::testTaskFinished();
    }
    catch (...)
    {
        emit TestPDBbind::testErrorEmission("Test Calculation Error: Internal Error, Check Input");
        emit TestPDBbind::testTaskFinished();
    }

    emit TestPDBbind::testTaskFinished();
}
