#include "scorepdbbind.h"

ScorePDBbind::ScorePDBbind(QObject *parent, std::string refinedExpDataPath_, std::string complexesParentDir_, std::string logFilePath_, unsigned int sf_type_)
    : QThread(parent), stopThread(false), refinedExpDataPath(refinedExpDataPath_), complexesParentDir(complexesParentDir_), logFilePath(logFilePath_), sf_type(sf_type_)
{

}

void ScorePDBbind::run()
{
    try
    {
        emit ScorePDBbind::stringEmission("Reading From Refined Data File....");

        /*std::vector<std::pair<std::string, float>> receptor_experimental_pki;

        std::ifstream inputfile;
        inputfile.exceptions(std::ifstream::badbit);

        try
        {
            inputfile.open(ScorePDBbind::refinedExpDataPath);
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
            throw Error_report("Error reading from file " + ScorePDBbind::refinedExpDataPath);
        }

        inputfile.close();*/

        unsigned int count = 0;

        emit ScorePDBbind::stringEmission("Starting Scoring....");

        Parser new_parser;
        new_parser.initialize_fconv();

        std::vector<std::pair<std::string, std::vector<std::pair<std::string, float>>>> all_calculated_xi;
        std::vector<std::vector<std::pair<unsigned int, unsigned int>>> all_calculated_rf;

        QDir main_directory(QString::fromStdString(ScorePDBbind::complexesParentDir));
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

                        if (file_path.contains("1zkk_refine"))
                        {
                            receptor = true;
                            receptor_path = file_path.toStdString();
                        }
                        else if (file_path.contains("out"))
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

                                std::string tempStr = "Scoring Complex with Folder Name: " + folder_name_string + "....";

                                emit ScorePDBbind::stringEmission(QString::fromStdString(tempStr));

                                Score score_object(receptor_path, ligand_path, &new_parser, ScorePDBbind::sf_type);
                                score_object.calculate();

                                all_calculated_rf.push_back(score_object.ligand->rf_inter);

                                all_calculated_xi.push_back(std::pair<std::string, std::vector<std::pair<std::string, float>>>(folder_name_string, score_object.ligand->terms_contributions));
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

        if(!this->stopThread)
        {
            std::vector<unsigned int> colmnsToRemove;

            if (all_calculated_rf.size() > 0)
            {
                emit ScorePDBbind::stringEmission("Removing Random Forest Zero Columns....");

                for (unsigned int i = 0; i < all_calculated_rf[0].size(); ++i) //iterate over columns
                {
                    bool foundNonZero = false;

                    for(unsigned int j = 0; j < all_calculated_rf.size(); ++j)
                    {
                        if(all_calculated_rf[j][i].second != 0)
                        {
                            foundNonZero = true;
                            break;
                        }
                    }

                    if(!foundNonZero)
                    {
                        colmnsToRemove.push_back(i);
                    }
                }
            }

            std::vector<std::vector<std::pair<unsigned int, unsigned int>>> all_calculated_rf_edited;

            for(unsigned int i = 0; i < all_calculated_rf.size(); ++i)
            {
                std::vector<std::pair<unsigned int, unsigned int>> temp_rf;

                for(unsigned int j = 0; j < all_calculated_rf[i].size(); ++j)
                {
                    bool remove_col = false;

                    for(unsigned int k = 0; k < colmnsToRemove.size(); ++k)
                    {
                        if(colmnsToRemove[k] == j)
                        {
                            remove_col = true;
                            break;
                        }
                    }

                    if(!remove_col)
                    {
                        temp_rf.push_back(all_calculated_rf[i][j]);
                    }
                }

                all_calculated_rf_edited.push_back(temp_rf);
            }

            emit ScorePDBbind::stringEmission("Writing Results....");

            std::ofstream xi_log;
            xi_log.open(ScorePDBbind::logFilePath);

            xi_log << "Scoring function type: *X-Score HC = 0, VINA = 1, VINA-Halogen = 2, AutoDock = 3, RF-SCORE = 4, RF-SCORE & X-SCORE HC = 5, RF-SCORE & VINA = 6, RF-SCORE & VINA-Halogen = 7, RF-SCORE & AutoDock = 8* | " + std::to_string(ScorePDBbind::sf_type);
            xi_log << "\n";

            xi_log << "Receptor||folder||ligand=";

            if (all_calculated_rf_edited.size() > 0)
            {
                for (unsigned int i = 0; i < all_calculated_rf_edited[0].size(); ++i)
                {
                    xi_log << "rf_";
                    xi_log << all_calculated_rf_edited[0][i].first;

                    if (i < (all_calculated_rf_edited[0].size() - 1))
                    {
                        xi_log << ",";
                    }
                }
            }

            if (all_calculated_xi.size() > 0)
            {
                if (all_calculated_rf_edited.size() > 0)
                {
                    if (all_calculated_rf_edited[0].size() > 0 && all_calculated_xi[0].second.size() > 0)
                    {
                        xi_log << ",";
                    }
                }

                for (unsigned int i = 0; i < all_calculated_xi[0].second.size(); ++i)
                {
                    xi_log << all_calculated_xi[0].second[i].first;

                    if (i < (all_calculated_xi[0].second.size() - 1))
                    {
                        xi_log << ",";
                    }
                }
            }

            xi_log << "\n";

            for (unsigned int i = 0; i < all_calculated_xi.size(); ++i)
            {
                xi_log << all_calculated_xi[i].first;
                xi_log << "=";

                for (unsigned int j = 0; j < all_calculated_rf_edited[i].size(); ++j)
                {
                    xi_log << all_calculated_rf_edited[i][j].second;

                    if (j < (all_calculated_rf_edited[i].size() - 1))
                    {
                        xi_log << ",";
                    }
                }

                if (all_calculated_rf_edited[i].size() > 0 && all_calculated_xi[i].second.size() > 0)
                {
                    xi_log << ",";
                }

                for (unsigned int j = 0; j < all_calculated_xi[i].second.size(); ++j)
                {
                    xi_log << all_calculated_xi[i].second[j].second;

                    if (j < (all_calculated_xi[i].second.size() - 1))
                    {
                        xi_log << ",";
                    }
                }

                /*xi_log << "|";

                for (unsigned int j = 0; j < receptor_experimental_pki.size(); ++j)
                {
                    if (receptor_experimental_pki[j].first == all_calculated_xi[i].first)
                    {
                        xi_log << receptor_experimental_pki[j].second;
                        break;
                    }
                }*/

                xi_log << "\n";
            }

            xi_log.close();
        }
    }
    catch (Error_report& err)
    {
        std::string tempStr = "Get Xs Calculation Error: " + err.errorText;
        emit ScorePDBbind::errorEmission(QString::fromStdString(tempStr));
        emit ScorePDBbind::taskFinished();
    }
    catch (...)
    {
        emit ScorePDBbind::errorEmission("Get Xs Calculation Error: Internal Error, Check Input");
        emit ScorePDBbind::taskFinished();
    }

    emit ScorePDBbind::taskFinished();
}
