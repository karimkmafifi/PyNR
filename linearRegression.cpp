#include "linearRegression.h"

LinearRegression::LinearRegression(QObject *parent)
    : QThread(parent), stopThread(false)
{
    LinearRegression::learningRate = 0.0;
    LinearRegression::noIterations = 0;
    LinearRegression::min_squared_error = 999999.;
    LinearRegression::num_vars = 0;
}

LinearRegression::LinearRegression()
{
    LinearRegression::learningRate = 0.0;
    LinearRegression::noIterations = 0;
    LinearRegression::min_squared_error = 999999.;
    LinearRegression::num_vars = 0;
}

void LinearRegression::run()
{
    try
    {
        emit LinearRegression::lrStringEmission("Reading From Training File....");

        std::vector<std::pair<std::vector<float>, float>> all_calculated_xi;
        std::vector<std::pair<std::string, float>> weights;

        std::ifstream inputfile;
        inputfile.exceptions(std::ifstream::badbit);
        std::string line;
        unsigned int count = 0;

        try
        {
            inputfile.open(LinearRegression::trainingFile);

            while (getline(inputfile, line))
            {
                if(this->stopThread)
                {
                    return;
                }

                if (count == 0)
                {
                    if (line.find("Scoring") != std::string::npos)
                    {
                        std::vector<std::string> elems;

                        Common::split(line, '|', elems);

                        if (elems.size() > 1)
                        {
                            LinearRegression::scoringFunctionType = Common::string_to<unsigned int>(elems[1]);
                        }
                        else
                        {
                            throw Error_report("Error reading header line 1 from Xs file. ");
                        }

                        if (LinearRegression::scoringFunctionType == 0)
                        {
                            weights.push_back(std::pair<std::string, float>("Bias", 0.0));
                            weights.push_back(std::pair<std::string, float>("Van der Waals X-SCORE HC", 0.0));
                            weights.push_back(std::pair<std::string, float>("Hydrogen Bonding X-SCORE HC", 0.0));
                            weights.push_back(std::pair<std::string, float>("Hydrophobic Contact X-SCORE HC", 0.0));
                            weights.push_back(std::pair<std::string, float>("Number of Rotors X-SCORE HC", 0.0));
                        }
                        else if (LinearRegression::scoringFunctionType == 1)
                        {
                            weights.push_back(std::pair<std::string, float>("Bias", 0.0));
                            weights.push_back(std::pair<std::string, float>("Gauss 1 VINA", 0.0));
                            weights.push_back(std::pair<std::string, float>("Gauss 2 VINA", 0.0));
                            weights.push_back(std::pair<std::string, float>("Repulsion VINA", 0.0));
                            weights.push_back(std::pair<std::string, float>("Hydrophobic VINA", 0.0));
                            weights.push_back(std::pair<std::string, float>("Non Directional Hydrogen Bond VINA", 0.0));
                            weights.push_back(std::pair<std::string, float>("Num Tors VINA", 0.0));
                        }
                        else if (LinearRegression::scoringFunctionType == 2)
                        {
                            weights.push_back(std::pair<std::string, float>("Bias", 0.0));
                            weights.push_back(std::pair<std::string, float>("Gauss 1 VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Gauss 2 VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Repulsion VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Hydrophobic VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Non Directional Hydrogen Bond VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Halogen Chlorine VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Halogen Bromine VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Halogen Iodine VINAXB", 0.0));
                            weights.push_back(std::pair<std::string, float>("Num Tors VINAXB", 0.0));
                        }
                        else if (LinearRegression::scoringFunctionType == 3)
                        {
                            weights.push_back(std::pair<std::string, float>("Bias", 0.0));
                            weights.push_back(std::pair<std::string, float>("Van der Waals AutoDock", 0.0));
                            weights.push_back(std::pair<std::string, float>("Hydrogen Bonding AutoDock", 0.0));
                            weights.push_back(std::pair<std::string, float>("Electrostatic AutoDock", 0.0));
                            weights.push_back(std::pair<std::string, float>("Solvation Charge Dependent AutoDock", 0.0));
                            weights.push_back(std::pair<std::string, float>("Number of Tors AutoDock", 0.0));
                        }
                        else
                        {
                            throw Error_report("Error reading header line 1 from Xs file. Scoring function type is not valid. ");
                        }

                        for (unsigned int i = 0; i < weights.size(); ++i)
                        {
                            if(i > 0)
                            {
                                LinearRegression::minimum_weights.push_back(weights[i]);
                            }
                        }

                        count = 1;

                        continue;
                    }
                    else
                    {
                        throw Error_report("Error reading header line 1 from Xs file. ");
                    }
                }
                else if (count == 1)
                {
                    if (line.find("Receptor||folder||ligand") != std::string::npos)
                    {
                        std::vector<std::string> elems;

                        Common::split(line, '=', elems);

                        if (elems.size() <= 1)
                        {
                            throw Error_report("Error reading header line 2 from Xs file " + LinearRegression::trainingFile + " ");
                        }

                        std::vector<std::string> elems2;

                        Common::split(elems[1], ',', elems2);

                        if(elems2.size() != (weights.size() - 1))
                        {
                            throw Error_report("Linear Regression error header line 2. Number of variables is not compatible with any of the LR Scoring Functions. ");
                        }

                        LinearRegression::num_vars = elems2.size();
                    }
                    else
                    {
                        throw Error_report("Error reading header line 2 from Xs file. ");
                    }

                    count = 2;

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

                    Common::split(elems2[0], ',', elems3);

                    if (elems3.size() != (weights.size() - 1))
                    {
                        throw Error_report("Error parsing Xs file. # of Xs is not equal to: " + std::to_string(weights.size() - 1) + " ");
                    }

                    std::vector<float> elems3_dble;
                    elems3_dble.push_back(1); //bias

                    for (unsigned int i = 0; i < elems3.size(); ++i)
                    {
                        elems3_dble.push_back(Common::string_to<float>(elems3[i]));
                    }

                    all_calculated_xi.push_back(std::pair<std::vector<float>, float>(elems3_dble, Common::string_to<float>(elems2[1])));
                }
            }

            inputfile.close();
        }
        catch (std::ifstream::failure e)
        {
            inputfile.close();
            throw Error_report(LinearRegression::trainingFile + " " + e.what());
        }
        catch (Error_report& e) {
            inputfile.close();
            throw Error_report(LinearRegression::trainingFile + " " + e.errorText);
        }
        catch (...)
        {
            inputfile.close();
            throw Error_report(LinearRegression::trainingFile + " ");
        }

        emit LinearRegression::lrStringEmission("Scaling Features....");

        if (all_calculated_xi.size() > 0)
        {
            for (unsigned int i = 0; i < all_calculated_xi[0].first.size(); ++i)
            {
                if (i > 0)
                {
                    float biggest = -99999.;
                    float smallest = 99999.;
                    float sum = 0.0;

                    for (unsigned int j = 0; j < all_calculated_xi.size(); ++j)
                    {
                        if (all_calculated_xi[j].first[i] > biggest)
                        {
                            biggest = all_calculated_xi[j].first[i];
                        }

                        if (all_calculated_xi[j].first[i] < smallest)
                        {
                            smallest = all_calculated_xi[j].first[i];
                        }

                        sum += all_calculated_xi[j].first[i];
                    }

                    float mean = sum / all_calculated_xi.size();

                    for (unsigned int j = 0; j < all_calculated_xi.size(); ++j)
                    {
                        all_calculated_xi[j].first[i] = (all_calculated_xi[j].first[i] - mean) / (biggest - smallest);
                    }
                }
            }
        }
        else
        {
            throw Error_report("Error. No Xs found in file: " + LinearRegression::trainingFile + " ");
            return;
        }

        emit LinearRegression::lrStringEmission("Starting Learning....");

        for (unsigned int i = 0; i < LinearRegression::noIterations; ++i)
        {
            if(this->stopThread)
            {
                return;
            }

            std::vector<std::pair<std::vector<float>, float>> all_calculated_xi_weight;

            for (unsigned int j = 0; j < all_calculated_xi.size(); ++j)
            {
                std::vector<float> weighted_xs;

                for (unsigned int k = 0; k < all_calculated_xi[j].first.size(); ++k)
                {
                    weighted_xs.push_back(all_calculated_xi[j].first[k] * weights[k].second);
                }

                all_calculated_xi_weight.push_back(std::pair<std::vector<float>, float>(weighted_xs, all_calculated_xi[j].second));
            }

            std::vector<float> sum_of_derivative_term;

            for (unsigned int j = 0; j < weights.size(); ++j)
            {
                sum_of_derivative_term.push_back(0);
            }

            float sum_mean_sqr_error = 0.0;

            for (unsigned int j = 0; j < all_calculated_xi_weight.size(); ++j)
            {
                float hypothesis = 0.0;

                for (unsigned int k = 0; k < all_calculated_xi_weight[j].first.size(); ++k)
                {
                    hypothesis += all_calculated_xi_weight[j].first[k];
                }

                for (unsigned int k = 0; k < all_calculated_xi_weight[j].first.size(); ++k)
                {
                    if (k > 0)
                    {
                        sum_of_derivative_term[k] += (hypothesis - all_calculated_xi_weight[j].second) * all_calculated_xi[j].first[k];
                    }
                    else
                    {
                        sum_of_derivative_term[k] += hypothesis - all_calculated_xi_weight[j].second;
                    }
                }

                sum_mean_sqr_error += Common::sqr(hypothesis - all_calculated_xi_weight[j].second);
            }

            float cost_function_term = (1.0 / all_calculated_xi_weight.size()) * sum_mean_sqr_error;

            emit LinearRegression::lrStringEmission(QString::fromStdString("Interation: " + std::to_string(i + 1) + ". MSE: " + std::to_string(LinearRegression::min_squared_error)));

            if (cost_function_term < LinearRegression::min_squared_error)
            {
                LinearRegression::min_squared_error = cost_function_term;

                for (unsigned int j = 0; j < weights.size(); ++j)
                {
                    if(j > 0)
                    {
                        minimum_weights[j - 1].second = weights[j].second;
                    }
                }
            }
            else if (cost_function_term > LinearRegression::min_squared_error)
            {
                if (checkMSEIncrease)
                {
                    break;
                }
            }

            LinearRegression::cost_function_terms.push_back(cost_function_term);

            for (unsigned int j = 0; j < weights.size(); ++j)
            {
                weights[j].second = weights[j].second - (LinearRegression::learningRate * (1.0 / all_calculated_xi.size()) * sum_of_derivative_term[j]);
            }
        }

    }
    catch (Error_report& err)
    {
        std::string tempStr = "Linear Regression Calculation Error: " + err.errorText;
        emit LinearRegression::lrErrorEmission(QString::fromStdString(tempStr));
        emit LinearRegression::lrTaskFinished();
    }
    catch (...)
    {
        emit LinearRegression::lrErrorEmission("Linear Regression Calculation Error: Internal Error, Check Input");
        emit LinearRegression::lrTaskFinished();
    }

    emit LinearRegression::lrTaskFinished();
}

void LinearRegression::setLRTrainVars(std::string trainingFile_, std::string outFile_, float learningRate_, unsigned int noIterations_, bool checkMSEIncrease_)
{
    LinearRegression::trainingFile = trainingFile_;
    LinearRegression::logFile = outFile_;
    LinearRegression::learningRate = learningRate_;
    LinearRegression::noIterations = noIterations_;
    LinearRegression::checkMSEIncrease = checkMSEIncrease_;
}

void LinearRegression::loadWeightsFile(std::string inFile_)
{
    LinearRegression::minimum_weights.clear();

    std::ifstream inputfile;
    inputfile.exceptions(std::ifstream::badbit);

    unsigned int count = 0;

    try
    {
        inputfile.open(inFile_);
        std::string line;

        while (getline(inputfile, line))
        {
            if (line.empty()) {}
            else
            {
                if(count == 0)
                {
                    std::vector<std::string> elems;
                    Common::split(line, ',', elems);

                    if(elems.size() <= 0)
                    {
                         throw Error_report("Error reading header line 1 from the Linear regression weights file. ");
                    }

                    LinearRegression::scoringFunctionType = Common::string_to<unsigned int>(elems[0]);

                    count = 1;

                    continue;
                }
                else if(count == 1)
                {
                    std::vector<std::string> elems;
                    Common::split(line, ',', elems);

                    if(elems.size() <= 0)
                    {
                         throw Error_report("Error reading header line 2 from the Linear regression weights file. ");
                    }

                    for (unsigned int i = 0; i < elems.size(); ++i)
                    {
                        LinearRegression::minimum_weights.push_back(std::pair<std::string, float>("", Common::string_to<float>(elems[i])));
                    }

                    count = 2;

                    continue;
                }
            }
        }

        inputfile.close();
    }
    catch (std::ifstream::failure e)
    {
        inputfile.close();
        throw Error_report(inFile_ + " " + e.what());
    }
    catch (Error_report& e) {
        inputfile.close();
        throw Error_report(inFile_ + " " + e.errorText);
    }
    catch (...)
    {
        inputfile.close();
        throw Error_report(inFile_ + " ");
    }
}

float LinearRegression::CalcBindingAffinity(std::vector<std::pair<std::string, float>>& terms_contributions, unsigned int scoringFunctionType_)
{
    if(LinearRegression::scoringFunctionType != scoringFunctionType_)
    {
        throw Error_report("Error. Linear Regression file Scoring Function type is not the same as the entered SF. ");
    }

    if(LinearRegression::minimum_weights.size() != terms_contributions.size())
    {
        throw Error_report("Error. Inconsistent number of weights with the number of contributions. ");
    }

    float e = 0.0;

    for(unsigned int i = 0; i < terms_contributions.size(); ++i)
    {
        e += (terms_contributions[i].second * LinearRegression::minimum_weights[i].second);
    }

    return e;
}

void LinearRegression::saveWeights()
{
    std::ofstream file_log;
    file_log.open(LinearRegression::logFile);

    file_log << LinearRegression::scoringFunctionType;
    file_log << ",";
    file_log << LinearRegression::learningRate;
    file_log << ",";
    file_log << LinearRegression::noIterations;
    file_log << ",";
    file_log << LinearRegression::min_squared_error;
    file_log << ",";
    file_log << LinearRegression::checkMSEIncrease;
    file_log << ",";
    file_log << LinearRegression::num_vars;

    file_log << "\n";

    for (unsigned int i = 0; i < LinearRegression::minimum_weights.size(); ++i)
    {
        file_log << LinearRegression::minimum_weights[i].second;

        if (i < (LinearRegression::minimum_weights.size() - 1))
        {
            file_log << ",";
        }
    }

    file_log << "\n";

    file_log.close();
}

LinearRegression::~LinearRegression()
{
}
