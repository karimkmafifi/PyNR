#ifndef SCORE_H
#define SCORE_H

#include "prepare_mol.h"

class Score
{
private:
    ScoringTerms scoring_terms;
    Parser* parser;
    unsigned int mol_type; // X-Score HC = 0, VINA = 1, VINA_Halogen = 2, AutoDock = 3, RF-SCORE = 4, RF-SCORE & X_SCORE HC = 5, RF-SCORE & VINA = 6, RF-SCORE & VINA_Halogen = 7, RF-SCORE & AutoDock = 8
    void initializeMol(const std::string& path, Molecule* mol, bool is_receptor);
    void prepare_ad_hydrogen();
    void ad_hbond_theta_eval(float& Hramp, float& racc, float& rdon, int atom_no, int recep_atom_no, std::vector<Atom*>& ligand_atoms);
    void xs_hydrogen_prepare(Molecule* mol);
    glm::vec3 xs_hbond_root(Atom* a);
    int atomnumber(const char* atomname);
public:
    Receptor* receptor;
    Ligand* ligand;
    Score(const std::string& rigid_path, const std::string& ligand_path, Parser* _parser, unsigned int type_);
    void calculate();
    ~Score();
};

#endif // SCORE_H
