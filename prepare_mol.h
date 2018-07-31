#ifndef PREPARE_MOL_H
#define PREPARE_MOL_H

#include "scoring_terms.h"
#include "ELEMENT_DATA_new_GN.h"
#include "octree_GN.hpp"
#include "mw_match_GN.h"

class Prepare_mol
{
private:
    void delete_worst_bonded(Atom* base);
    bool bonded_to_HD(Atom* a);
    bool bonded_to_heteroatom(Atom* a);
public:
    Prepare_mol();
    void get_elements(Molecule* mol, bool reset_names = true);
    void get_connections(Molecule* mol);
    void assign_types(Molecule* mol);
};

#endif // PREPARE_MOL_H
