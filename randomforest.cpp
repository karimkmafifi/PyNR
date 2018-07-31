#include "randomforest.h"


float RandomForest::calculate_mean(std::vector<std::pair<std::vector<float>, float>>& data)
{
    float sum = 0;

    for (unsigned int i = 0; i < data.size(); ++i)
    {
        sum += data[i].second;
    }

    return sum / data.size();
}

float RandomForest::calculate_mean_mse(std::vector<float>& data)
{
    float sum = 0;

    for (unsigned int i = 0; i < data.size(); ++i)
    {
        sum += data[i];
    }

    float mean = sum / data.size();

    float mse = 0;

    for (unsigned int i = 0; i < data.size(); ++i)
    {
        mse += (data[i] - mean) * (data[i] - mean);
    }

    return mse;
}

float RandomForest::calculate_variance(float split, std::vector<float>& variable_space)
{
    std::vector<float> data1;
    std::vector<float> data2;

    for (unsigned int i = 0; i < variable_space.size(); ++i)
    {
        if (variable_space[i] <= split)
        {
            data1.push_back(variable_space[i]);
        }
        else
        {
            data2.push_back(variable_space[i]);
        }
    }

    return RandomForest::calculate_mean_mse(data1) + RandomForest::calculate_mean_mse(data2);
}

void RandomForest::delete_tree(node*& root)
{
    if (root != NULL) {
        RandomForest::delete_tree(root->left);
        RandomForest::delete_tree(root->right);
        delete root;
        root = NULL;
    }
}

void RandomForest::split_golden_ratio(std::vector<float>& variable_space, float min, float max, float previous_error)
{
    float current_error = (1 - RandomForest::r) * (max - min);

    if (current_error == previous_error)
    {
        RandomForest::error = current_error;
        return;
    }

    float min_value_phi = RandomForest::r * min;
    float max_value_phi = RandomForest::r * max;

    float x1 = max - (max_value_phi - min_value_phi);
    float x2 = min + (max_value_phi - min_value_phi);

    RandomForest::split_point = (x1 + x2) / 2;

    float fun_x1 = RandomForest::calculate_variance(x1, variable_space);
    float fun_x2 = RandomForest::calculate_variance(x2, variable_space);

    if (fun_x1 < fun_x2)
    {
        float new_max_value = x2;
        RandomForest::split_golden_ratio(variable_space, min, new_max_value, current_error);
    }
    else if(fun_x2 < fun_x1)
    {
        float new_min_value = x1;
        RandomForest::split_golden_ratio(variable_space, new_min_value, max, current_error);
    }
    else
    {
        RandomForest::error = current_error;
        return;
    }
}

node* RandomForest::create_tree(std::vector<std::pair<std::vector<float>, float>>& bag_selected_data, unsigned int depth, unsigned int caller_id, bool is_left)
{
    float mean = RandomForest::calculate_mean(bag_selected_data);

    node* current_node = new node();

    current_node->start = false;
    current_node->mean = mean;
    current_node->left = NULL;
    current_node->right = NULL;
    current_node->id = node_id;
    current_node->parent_id = caller_id;

    if (is_left)
    {
        current_node->is_left = true;
    }
    else
    {
        current_node->is_left = false;
    }

    RandomForest::node_id++;

    if(this->stopThread)
    {
        return current_node;
    }

    if (bag_selected_data.size() <= RandomForest::internalMinSplit)
    {
        return current_node;
    }

    if (depth >= RandomForest::maxDepth)
    {
        return current_node;
    }

    std::vector<unsigned int> cand_vars;

    for (unsigned int j = 0; j < RandomForest::mtry; ++j)
    {
        unsigned int random_number = std::rand() % RandomForest::num_vars;
        cand_vars.push_back(random_number);
    }

    float min_error = -1;
    float min_variance = -1;
    float min_split = -1;
    int	split_var = -1;

    for (unsigned int j = 0; j < cand_vars.size(); ++j)
    {
        std::vector<float> variable_space;

        for (unsigned int w = 0; w < bag_selected_data.size(); ++w)
        {
            variable_space.push_back(bag_selected_data[w].first[cand_vars[j]]);
        }

        float min = 99999;
        float max = -99999;

        for (unsigned int w = 0; w < variable_space.size(); ++w)
        {
            if (variable_space[w] < min)
            {
                min = variable_space[w];
            }

            if (variable_space[w] > max)
            {
                max = variable_space[w];
            }
        }

        if (min == max)
        {
            continue;
        }

        RandomForest::split_point = 0;
        RandomForest::error = 0;

        RandomForest::split_golden_ratio(variable_space, min, max, -9999);

        float variance = RandomForest::calculate_variance(RandomForest::split_point, variable_space);

        if (min_error == -1)
        {
            min_error = RandomForest::error;
            min_variance = variance;
            min_split = RandomForest::split_point;
            split_var = cand_vars[j];
        }
        else
        {
            if (RandomForest::error < min_error)
            {
                min_error = RandomForest::error;
                min_split = RandomForest::split_point;
                split_var = cand_vars[j];
            }
            else if (variance < min_variance)
            {
                min_variance = variance;
                min_split = RandomForest::split_point;
                split_var = cand_vars[j];
            }
        }

    }

    if (split_var == -1)
    {
        return current_node;
    }

    current_node->split_var = split_var;
    current_node->split_val = min_split;

    std::vector<std::pair<std::vector<float>, float>> data1;
    std::vector<std::pair<std::vector<float>, float>>	data2;

    for (unsigned int w = 0; w < bag_selected_data.size(); ++w)
    {
        if (bag_selected_data[w].first[split_var] <= min_split)
        {
            data1.push_back(bag_selected_data[w]);
        }
        else
        {
            data2.push_back(bag_selected_data[w]);
        }
    }

    current_node->left = create_tree(data1, depth + 1, current_node->id, true);
    current_node->right = create_tree(data2, depth + 1, current_node->id, false);
    return current_node;
}

void RandomForest::setRFTrainVars(std::string trainingFile_, std::string outFile_, unsigned int noTrees_, unsigned int maxDepth_, float percentOfBootstrap_, unsigned int mtry_, unsigned int internalMinSplit_)
{
    RandomForest::trainingFile = trainingFile_;
    RandomForest::logFile = outFile_;
    RandomForest::noTrees = noTrees_;
    RandomForest::maxDepth = maxDepth_;
    RandomForest::percentOfBootstrap = percentOfBootstrap_;
    RandomForest::mtry = mtry_;
    RandomForest::internalMinSplit = internalMinSplit_;
}

void RandomForest::run()
{
    try
    {
        if(RandomForest::noTrees <= 0)
        {
            throw Error_report("Random Forest number of trees should be bigger than 0. ");
        }

        if(RandomForest::maxDepth <= 0)
        {
            throw Error_report("Random Forest max depth should be bigger than 0. ");
        }

        if((RandomForest::percentOfBootstrap > 100) || (RandomForest::percentOfBootstrap <= 0))
        {
            throw Error_report("The percentage to be selected as bootstrap data should be 0 < % >= 100. ");
        }

        if(RandomForest::mtry <= 0)
        {
            throw Error_report("Random Forest mtry should be bigger than 0. ");
        }

        if(RandomForest::internalMinSplit <= 0)
        {
            throw Error_report("Random Forest internal minimum split should be bigger than 0. ");
        }

        emit RandomForest::rfStringEmission("Reading From Training File....");

        std::vector<std::pair<std::vector<float>, float>> all_data;

        std::ifstream inputfile;
        inputfile.exceptions(std::ifstream::badbit);
        std::string line;
        unsigned int count = 0;

        try
        {
            inputfile.open(RandomForest::trainingFile);

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
                            RandomForest::scoringFunctionType = Common::string_to<unsigned int>(elems[1]);
                        }
                        else
                        {
                            throw Error_report("Error reading header line 1 from Xs file. ");
                        }

                        count = 1;

                        continue;
                    }
                    else
                    {
                        throw Error_report("Error reading header line 1 from Xs file " + RandomForest::trainingFile + " ");
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
                            throw Error_report("Error reading header line 2 from Xs file. ");
                        }

                        std::vector<std::string> elems2;

                        Common::split(elems[1], ',', elems2);

                        RandomForest::num_vars = elems2.size();

                        for(unsigned int i = 0; i < elems2.size(); ++i)
                        {
                            if (elems2[i].find("rf_") != std::string::npos)
                            {
                                std::vector<std::string> elems3;
                                Common::split(elems2[i], '_', elems3);
                                RandomForest::rf_variables_used.push_back(Common::string_to<unsigned int>(elems3[1]));
                            }
                        }
                    }
                    else
                    {
                        throw Error_report("Error reading header line 2 from Xs file " + RandomForest::trainingFile + " ");
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

                    if (elems3.size() != RandomForest::num_vars)
                    {
                        throw Error_report("Error parsing Xs file. # of Xs is not equal to: " + std::to_string(RandomForest::num_vars) + " ");
                    }

                    std::vector<float> elems3_dble;

                    for (unsigned int i = 0; i < elems3.size(); ++i)
                    {
                        elems3_dble.push_back(Common::string_to<float>(elems3[i]));
                    }

                    all_data.push_back(std::pair<std::vector<float>, float>(elems3_dble, Common::string_to<float>(elems2[1])));
                }
            }

            inputfile.close();
        }
        catch (std::ifstream::failure e)
        {
            inputfile.close();
            throw Error_report(RandomForest::trainingFile + " " + e.what());
        }
        catch (Error_report& e) {
            inputfile.close();
            throw Error_report(RandomForest::trainingFile + " " + e.errorText);
        }
        catch (...)
        {
            inputfile.close();
            throw Error_report(RandomForest::trainingFile + " ");
        }

        std::srand(1);

        unsigned int num_boot = (RandomForest::percentOfBootstrap / 100) * all_data.size();

        emit RandomForest::rfStringEmission("Starting Learning....");

        for (unsigned int i = 0; i < RandomForest::noTrees; ++i)
        {
            if(this->stopThread)
            {
                return;
            }

            emit RandomForest::rfStringEmission(QString::fromStdString("Growing Tree #: " + std::to_string(i + 1) + "...."));

            std::vector<std::pair<std::vector<float>, float>> bag_selected_data;

            for (unsigned int j = 0; j < num_boot; ++j)
            {
                unsigned int random_number = std::rand() % all_data.size();
                bag_selected_data.push_back(all_data[random_number]);
            }

            RandomForest::node_id = 0;

            node* tree = RandomForest::create_tree(bag_selected_data, 0, 0, false);
            tree->start = true;

            RandomForest::trees.push_back(tree);
        }

    }
    catch (Error_report& err)
    {
        std::string tempStr = "Random Forest Calculation Error: " + err.errorText;
        emit RandomForest::rfErrorEmission(QString::fromStdString(tempStr));
        emit RandomForest::rfTaskFinished();
    }
    catch (...)
    {
        emit RandomForest::rfErrorEmission("Random Forest Calculation Error: Internal Error, Check Input");
        emit RandomForest::rfTaskFinished();
    }

    emit RandomForest::rfTaskFinished();
}

void RandomForest::find_parent(node *n, node *toadd)
{
    if (!RandomForest::found_parent_bl)
    {
        if (n == NULL)
        {
            return;
        }

        if (n->id == toadd->parent_id)
        {
            if (toadd->is_left)
            {
                n->left = toadd;
            }
            else
            {
                n->right = toadd;
            }

            RandomForest::found_parent_bl = true;
            return;
        }

        RandomForest::find_parent(n->left, toadd);
        RandomForest::find_parent(n->right, toadd);
    }
}

float RandomForest::get_tree_predict(node* root, std::vector<float>& test_xs)
{
    float predicted_mean = -9999;
    node* current_node = root;

    while (1)
    {
        if (current_node->left == NULL && current_node->right == NULL)
        {
            predicted_mean = current_node->mean;
            break;
        }
        else if (current_node->left != NULL && current_node->right != NULL)
        {
            if (test_xs[current_node->split_var] <= current_node->split_val)
            {
                current_node = current_node->left;
            }
            else
            {
                current_node = current_node->right;
            }
        }
        else if (current_node->left != NULL && current_node->right == NULL)
        {
            if (test_xs[current_node->split_var] <= current_node->split_val)
            {
                current_node = current_node->left;
            }
            else
            {
                predicted_mean = current_node->mean;
                break;
            }
        }
        else
        {
            if (test_xs[current_node->split_var] > current_node->split_val)
            {
                current_node = current_node->right;
            }
            else
            {
                predicted_mean = current_node->mean;
                break;
            }
        }
    }

    return predicted_mean;
}

float RandomForest::predict(std::vector<float>& test_xs, unsigned int scoringFunctionType_)
{
    //"Scoring function type: *X-Score HC = 0, VINA = 1, VINA_Halogen = 2, AutoDock = 3, RF-SCORE = 4, RF-SCORE & X_SCORE HC = 5, RF-SCORE & VINA = 6, RF-SCORE & VINA_Halogen = 7, RF-SCORE & AutoDock = 8* | "

    if(RandomForest::scoringFunctionType != scoringFunctionType_)
    {
        throw Error_report("Error. Random Forest file Scoring Function type is not the same as the entered SF. ");
    }

    if(RandomForest::scoringFunctionType < 4)
    {
        throw Error_report("Error. Incorrect Scoring Function type for Random Forest. ");
    }

    std::vector<float> new_test_xs;

    for(unsigned int i = 0; i < test_xs.size(); ++i)
    {
        if(i < 81)
        {
            bool rfVariableExist = false;

            for(unsigned int j = 0; j < RandomForest::rf_variables_used.size(); ++j)
            {
                if(i == RandomForest::rf_variables_used[j])
                {
                    rfVariableExist = true;
                    break;
                }
            }

            if(rfVariableExist)
            {
                new_test_xs.push_back(test_xs[i]);
            }
        }
        else
        {
            new_test_xs.push_back(test_xs[i]);
        }
    }

    std::vector<float> predicted_ys;

    for (unsigned int i = 0; i < RandomForest::trees.size(); ++i)
    {
        float tree_predicted_y;

        tree_predicted_y = RandomForest::get_tree_predict(RandomForest::trees[i], new_test_xs);

        predicted_ys.push_back(tree_predicted_y);
    }

    float sum_y = 0;

    for (unsigned int i = 0; i < predicted_ys.size(); ++i)
    {
        sum_y += predicted_ys[i];
    }

    float mean_y = sum_y / predicted_ys.size();

    return mean_y;
}

void RandomForest::read_trees(std::string treesinPath)
{
    RandomForest::delete_trees();

    std::ifstream inputfile;
    inputfile.exceptions(std::ifstream::badbit);

    unsigned int count = 0;

    try
    {
        inputfile.open(treesinPath);
        std::string line;

        while (getline(inputfile, line))
        {
            if (line.empty()) {}
            else if (Common::starts_with(line, "#"))
            {}
            else
            {
                if(count == 0)
                {
                    std::vector<std::string> elems;
                    Common::split(line, ',', elems);

                    if(elems.size() <= 0)
                    {
                         throw Error_report("Error reading header line 1 from the Random Forest file. ");
                    }

                    RandomForest::scoringFunctionType = Common::string_to<unsigned int>(elems[0]);

                    count = 1;

                    continue;
                }
                else if(count == 1)
                {
                    std::vector<std::string> elems;
                    Common::split(line, ',', elems);

                    if(elems.size() <= 0)
                    {
                         throw Error_report("Error reading header line 2 from the Random Forest file. ");
                    }

                    for (unsigned int i = 0; i < elems.size(); ++i)
                    {
                        RandomForest::rf_variables_used.push_back(Common::string_to<unsigned int>(elems[i]));
                    }

                    count = 2;

                    continue;
                }
                else
                {
                    std::vector<std::string> elems;
                    Common::split(line, ',', elems);

                    node* current_node = new node();

                    if (elems[0] == "0")
                    {
                        current_node->start = false;
                    }
                    else
                    {
                        current_node->start = true;
                    }

                    if (elems[1] == "0")
                    {
                        current_node->is_left = false;
                    }
                    else
                    {
                        current_node->is_left = true;
                    }

                    current_node->id = Common::string_to<unsigned int>(elems[2]);
                    current_node->parent_id = Common::string_to<unsigned int>(elems[3]);
                    current_node->split_var = Common::string_to<int>(elems[4]);
                    current_node->split_val = Common::string_to<float>(elems[5]);
                    current_node->mean = Common::string_to<float>(elems[6]);

                    if (!current_node->start)
                    {
                        current_node->left = NULL;
                        current_node->right = NULL;
                        RandomForest::find_parent(trees[trees.size() - 1], current_node);

                        if (RandomForest::found_parent_bl)
                        {
                            RandomForest::found_parent_bl = false;
                        }
                        else
                        {
                            throw Error_report("Could not find a node's parent for node id: " + std::to_string(current_node->id) + ". ");
                        }
                    }
                    else
                    {
                        RandomForest::trees.push_back(current_node);
                    }
                }
            }
        }
    }
    catch (std::ifstream::failure e)
    {
        inputfile.close();
        throw Error_report(treesinPath + " " + e.what());
    }
    catch (Error_report& e) {
        inputfile.close();
        throw Error_report(treesinPath + " " + e.errorText);
    }
    catch (...)
    {
        inputfile.close();
        throw Error_report(treesinPath + " ");
    }

    inputfile.close();
}

void RandomForest::iterate_pre(node *n, std::ofstream &log)
{
    if (n == NULL)
    {
        return;
    }

    log << n->start;
    log << ",";
    log << n->is_left;
    log << ",";
    log << n->id;
    log << ",";
    log << n->parent_id;
    log << ",";
    log << n->split_var;
    log << ",";
    log << n->split_val;
    log << ",";
    log << n->mean;
    log << "\n";

    RandomForest::iterate_pre(n->left, log);
    RandomForest::iterate_pre(n->right, log);
}

void RandomForest::save_forest()
{
    std::ofstream file_log;
    file_log.open(RandomForest::logFile);

    file_log << RandomForest::scoringFunctionType;
    file_log << ",";
    file_log << RandomForest::noTrees;
    file_log << ",";
    file_log << RandomForest::maxDepth;
    file_log << ",";
    file_log << RandomForest::percentOfBootstrap;
    file_log << ",";
    file_log << RandomForest::mtry;
    file_log << ",";
    file_log << RandomForest::internalMinSplit;
    file_log << ",";
    file_log << RandomForest::num_vars;
    file_log << "\n";

    for (unsigned int i = 0; i < RandomForest::rf_variables_used.size(); ++i)
    {
        file_log << RandomForest::rf_variables_used[i];

        if (i < (RandomForest::rf_variables_used.size() - 1))
        {
            file_log << ",";
        }
    }

    file_log << "\n";

    for (unsigned int i = 0; i < RandomForest::trees.size(); ++i)
    {
        RandomForest::iterate_pre(RandomForest::trees[i], file_log);

        if (i < (RandomForest::trees.size() - 1))
        {
            file_log << "#NEWTREE: \n";
        }
    }

    file_log.close();
}

RandomForest::RandomForest(QObject *parent)
    : QThread(parent), stopThread(false)
{
}

RandomForest::RandomForest()
{
}

void RandomForest::delete_trees()
{
    for (unsigned int i = 0; i < RandomForest::trees.size(); ++i)
    {
        RandomForest::delete_tree(RandomForest::trees[i]);
    }
}

RandomForest::~RandomForest()
{
    for (unsigned int i = 0; i < RandomForest::trees.size(); ++i)
    {
        RandomForest::delete_tree(RandomForest::trees[i]);
    }
}
