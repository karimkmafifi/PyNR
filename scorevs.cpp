#include "scorevs.h"

ScoreVS::ScoreVS(QObject *parent, std::string recOrDirPath_, std::vector<std::string> ligandsPaths_, std::string lr_rf_FilePath_, std::string logFilePath_, bool isVS_, unsigned int sfType_)
    : QThread(parent), stopThread(false), recOrDirPath(recOrDirPath_), ligandsPaths(ligandsPaths_), lr_rf_FilePath(lr_rf_FilePath_),
      logFilePath(logFilePath_), isVS(isVS_), sfType(sfType_)
{

}

void ScoreVS::run()
{
    try
    {
        unsigned int count = 0;

        emit ScoreVS::vsStringEmission("Starting Scoring.... \n");

        Parser new_parser;
        new_parser.initialize_fconv();

        LinearRegression linear_regression;
        RandomForest random_forest;

        if(sfType < 4)
        {
            emit ScoreVS::vsStringEmission("Loading Linear Regression Weights.... \n");
            linear_regression.loadWeightsFile(lr_rf_FilePath);
        }
        else
        {
            emit ScoreVS::vsStringEmission("Loading Random Forest Trees.... \n");
            random_forest.read_trees(lr_rf_FilePath);
        }

        std::vector<std::string> logTexts;

        if(!isVS)
        {
            std::string header_ = "Receptor||folder||ligand terms contributions= .... | pki = .... \n";
            logTexts.push_back(header_);

            QDir main_directory(QString::fromStdString(ScoreVS::recOrDirPath));
            QFileInfoList folders_list = main_directory.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::NoSort);

            foreach(QFileInfo subfolder, folders_list)
            {
                if (subfolder.isDir())
                {
                    QString subfolder_path = subfolder.absoluteFilePath();
                    QDir sub_directory(subfolder_path);
                    QFileInfoList files_list = sub_directory.entryInfoList(QDir::Files | QDir::NoDotAndDotDot, QDir::NoSort);
                    bool receptor = false;
                    bool ligand = false;
                    std::string receptor_path;
                    std::string ligand_path;

                    QStringList path_string_list = subfolder_path.split("/");
                    QString folder_id = path_string_list[path_string_list.size() - 1];
                    std::string folder_name_string = folder_id.toStdString();

                    foreach(QFileInfo file, files_list)
                    {
                        if (file.isFile())
                        {
                            QString file_path = file.absoluteFilePath();

                            if (file_path.contains("protein"))
                            {
                                receptor = true;
                                receptor_path = file_path.toStdString();
                            }
                            else if (file_path.contains("ligand"))
                            {
                                ligand = true;
                                ligand_path = file_path.toStdString();
                            }

                            if (receptor == true && ligand == true)
                            {
                                try
                                {
                                    if(this->stopThread)
                                    {
                                        return;
                                    }

                                    std::string tempStr = "Scoring Complex with Folder Name= " + folder_name_string + ".... \n";

                                    emit ScoreVS::vsStringEmission(QString::fromStdString(tempStr));

                                    Score score_object(receptor_path, ligand_path, &new_parser, ScoreVS::sfType);
                                    score_object.calculate();

                                    float e = 0;

                                    if(sfType < 4)
                                    {
                                        e = linear_regression.CalcBindingAffinity(score_object.ligand->terms_contributions,  ScoreVS::sfType);
                                    }
                                    else
                                    {
                                        std::vector<float>test_xs;

                                        for(unsigned int i = 0; i < score_object.ligand->rf_inter.size(); ++i)
                                        {
                                            test_xs.push_back(score_object.ligand->rf_inter[i].second);
                                        }

                                        for(unsigned int i = 0; i < score_object.ligand->terms_contributions.size(); ++i)
                                        {
                                            test_xs.push_back(score_object.ligand->terms_contributions[i].second);
                                        }

                                        e = random_forest.predict(test_xs, ScoreVS::sfType);
                                    }

                                    std::string logText = folder_name_string + "-" + folder_name_string + " terms contributions= ";

                                    for(unsigned int i = 0; i < score_object.ligand->rf_inter.size(); ++i)
                                    {
                                        logText += ("rf_" + std::to_string(i+1) + ": " + std::to_string(score_object.ligand->rf_inter[i].second));

                                        if(i < (score_object.ligand->rf_inter.size() - 1))
                                        {
                                            logText += ", ";
                                        }
                                    }

                                    if(score_object.ligand->rf_inter.size() > 0 && score_object.ligand->terms_contributions.size() > 0)
                                    {
                                        logText += ", ";
                                    }

                                    for(unsigned int i = 0; i < score_object.ligand->terms_contributions.size(); ++i)
                                    {
                                        logText += (score_object.ligand->terms_contributions[i].first + ": " +
                                                std::to_string(score_object.ligand->terms_contributions[i].second));

                                        if(i < (score_object.ligand->terms_contributions.size() - 1))
                                        {
                                            logText += ", ";
                                        }
                                    }

                                    logText += " | Pki: " + std::to_string(e);

                                    logText += "\n";

                                    emit ScoreVS::vsStringEmission(QString::fromStdString(logText));

                                    logTexts.push_back(logText);

                                    count++;
                                }
                                catch (Error_report& err)
                                {
                                    throw Error_report(err.errorText + "| Error scoring from folder " + folder_name_string);
                                }
                                catch (...)
                                {
                                    throw Error_report("Error scoring from folder " + folder_name_string);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            std::string header_ = "Receptor||folder||ligand terms contributions= .... | pki = .... \n";
            logTexts.push_back(header_);

            std::vector<std::string> elems3;
            Common::split(recOrDirPath, '/', elems3);

            std::vector<std::string> elems4;
            Common::split(elems3[elems3.size() - 1], '.', elems4);

            for(int i = 0; i < ligandsPaths.size(); ++i)
            {
                std::vector<std::string> elems;
                Common::split(ligandsPaths[i], '/', elems);

                std::vector<std::string> elems2;

                Common::split(elems[elems.size() - 1], '.', elems2);

                try
                {
                    if(this->stopThread)
                    {
                        return;
                    }

                    std::string tempStr = "Scoring Complex with Ligand Name= " + elems[elems.size() - 1] + ".... \n";

                    emit ScoreVS::vsStringEmission(QString::fromStdString(tempStr));

                    Score score_object(recOrDirPath, ligandsPaths[i], &new_parser, ScoreVS::sfType);
                    score_object.calculate();

                    float e = 0;

                    if(sfType < 4)
                    {
                        e = linear_regression.CalcBindingAffinity(score_object.ligand->terms_contributions,  ScoreVS::sfType);
                    }
                    else
                    {
                        std::vector<float>test_xs;

                        for(unsigned int i = 0; i < score_object.ligand->rf_inter.size(); ++i)
                        {
                            test_xs.push_back(score_object.ligand->rf_inter[i].second);
                        }

                        for(unsigned int i = 0; i < score_object.ligand->terms_contributions.size(); ++i)
                        {
                            test_xs.push_back(score_object.ligand->terms_contributions[i].second);
                        }

                        e = random_forest.predict(test_xs, ScoreVS::sfType);
                    }

                    std::string logText = elems4[0] + "-" + elems2[0] + " terms contributions= ";

                    for(unsigned int i = 0; i < score_object.ligand->rf_inter.size(); ++i)
                    {
                        logText += ("rf_" + std::to_string(i+1) + ": " + std::to_string(score_object.ligand->rf_inter[i].second));

                        if(i < (score_object.ligand->rf_inter.size() - 1))
                        {
                            logText += ", ";
                        }
                    }

                    if(score_object.ligand->rf_inter.size() > 0 && score_object.ligand->terms_contributions.size() > 0)
                    {
                        logText += ", ";
                    }

                    for(unsigned int i = 0; i < score_object.ligand->terms_contributions.size(); ++i)
                    {
                        logText += (score_object.ligand->terms_contributions[i].first + ": " +
                                std::to_string(score_object.ligand->terms_contributions[i].second));

                        if(i < (score_object.ligand->terms_contributions.size() - 1))
                        {
                            logText += ", ";
                        }
                    }

                    logText += " | Pki: " + std::to_string(e);

                    logText += "\n";

                    emit ScoreVS::vsStringEmission(QString::fromStdString(logText));

                    logTexts.push_back(logText);

                    count++;
                }
                catch (Error_report& err)
                {
                    throw Error_report(err.errorText + "| Error scoring complex of ligand name: " + elems[elems.size() - 1]);
                }
                catch (...)
                {
                    throw Error_report("Error scoring complex of ligand name: " + elems[elems.size() - 1]);
                }
            }
        }

        emit ScoreVS::vsStringEmission("Writing results to file.... \n");

        std::ofstream logOS;
        logOS.open(ScoreVS::logFilePath);

        for(unsigned int i = 0; i < logTexts.size(); ++i)
        {
            logOS << logTexts[i];
        }

        logOS.close();

        emit ScoreVS::vsStringEmission("Done");
    }
    catch (Error_report& err)
    {
        std::string tempStr = "Score Calculation Error: " + err.errorText;
        emit ScoreVS::vsErrorEmission(QString::fromStdString(tempStr));
        emit ScoreVS::vsTaskFinished();
    }
    catch (...)
    {
        emit ScoreVS::vsErrorEmission("Score Calculation Error: Internal Error, Check Input");
        emit ScoreVS::vsTaskFinished();
    }

    emit ScoreVS::vsTaskFinished();
}

