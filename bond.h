#ifndef BOND_H
#define BOND_H

#include "common.h"

struct Bond {
public:
    int id;
    Atom* connected_atom_1;
    Atom* connected_atom_2;
    float length;
    Bond() {};
};

#endif // BOND_H
