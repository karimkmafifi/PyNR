#ifndef RANDOMFOREST_H
#define RANDOMFOREST_H

#include "score.h"

struct node
{
public:
    bool start;
    bool is_left;
    unsigned int id;
    unsigned int parent_id;
    int split_var;
    float split_val;
    float mean;
    node* left;
    node* right;
};

class RandomForest : public QThread
{
    Q_OBJECT

private:
    std::vector<node*> trees;
    std::string trainingFile;
    std::string logFile;
    unsigned int scoringFunctionType;
    unsigned int noTrees;
    unsigned int maxDepth;
    float percentOfBootstrap;
    unsigned int mtry;
    unsigned int internalMinSplit;
    unsigned int num_vars;
    std::vector<unsigned int> rf_variables_used;
    unsigned int node_id = 0;
    float r = (std::sqrt(5) - 1) / 2;
    float split_point = 0;
    float error = 0;
    bool found_parent_bl = false;

private:
    node* create_tree(std::vector<std::pair<std::vector<float>, float>>& bag_selected_data, unsigned int depth, unsigned int caller_id, bool is_left);
    void delete_tree(node*& root);
    void split_golden_ratio(std::vector<float>& variable_space, float min, float max, float previous_error);
    float calculate_mean(std::vector<std::pair<std::vector<float>, float>>& data);
    float calculate_variance(float split, std::vector<float>& variable_space);
    float calculate_mean_mse(std::vector<float>& data);
    void iterate_pre(node *n, std::ofstream &log);
    void find_parent(node *n, node *toadd);
    float get_tree_predict(node* root, std::vector<float>& test_xs);

public:
    bool stopThread;

public:
    RandomForest(QObject *parent);
    RandomForest();
    ~RandomForest();
    void setRFTrainVars(std::string trainingFile_, std::string outFile_, unsigned int noTrees_, unsigned int maxDepth_, float percentOfBootstrap_, unsigned int mtry_, unsigned int internalMinSplit_);
    void run();
    void save_forest();
    void read_trees(std::string treesinPath);
    float predict(std::vector<float>& test_xs, unsigned int scoringFunctionType_);
    void delete_trees();

signals:
    void rfStringEmission(QString log);
    void rfTaskFinished();
    void rfErrorEmission(QString err);
};

#endif // RANDOMFOREST_H
