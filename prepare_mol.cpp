#include "prepare_mol.h"

Prepare_mol::Prepare_mol()
{

}

void Prepare_mol::get_elements(Molecule* mol, bool reset_names)
{
    for (int at = 0; at < mol->atoms.size(); ++at)
    {
        mol->atoms[at]->get_element(reset_names);
    }
}

void Prepare_mol::delete_worst_bonded(Atom* base)
{
    Atom* wa;
    float wp = 1.1;
    for (int bt = 0; bt < base->atom_bonds.size(); ++bt)
    {
        Atom* connected_atom_bt = base->get_connected_atom(bt);

        float mod = 1.;
        if (connected_atom_bt->atom_bonds.size() > 4) {
            mod = 0.5;
        }
        else if (connected_atom_bt->element == "O" && connected_atom_bt->atom_bonds.size() > 2) mod = 0.5;
        float dst = std::sqrt(Common::vec_distance_sqr(base->coords, connected_atom_bt->coords));
        float bp = EDnew::EDAccess::get_max_probability(base->element, connected_atom_bt->element, dst);

        bp *= mod;

        if (bp < wp) {
            wp = bp;
            wa = connected_atom_bt;
        }
    }

    //!wa sollte jetzt das Atom mit der geringsten Bindungswahrscheinlichkeit sein
    for (int bt = 0; bt < base->atom_bonds.size(); ++bt)
    {
        if (base->get_connected_atom(bt) == wa)
        {
            base->atom_bonds.erase(base->atom_bonds.begin() + bt);
            if (wa->element != "H") base->n_heavy_bonded--;
            break;
        }
    }
    for (int bt = 0; bt < wa->atom_bonds.size(); ++bt)
    {
        if (wa->get_connected_atom(bt) == base) {
            wa->atom_bonds.erase(wa->atom_bonds.begin() + bt);
            if (base->element != "H") wa->n_heavy_bonded--;
            break;
        }
    }
}

void Prepare_mol::get_connections(Molecule* mol)
{
    EDnew::initialize();
    float dst;
    std::string key;
    int bond_id = 0;

    OCTREE<Atom>* pwp = new OCTREE<Atom>(mol->atoms, 2.5);
    for (int at = 0; at < mol->atoms.size(); ++at)
    {
        for (vector<Atom*>::iterator bt = pwp->begin(2., mol->atoms[at]->coords); bt != pwp->end(); ++bt) {
            if ((*bt)->intern_id <= mol->atoms[at]->intern_id) continue;

            dst = mol->atoms[at]->coords[0] - (*bt)->coords[0]; dst *= dst;
            if (dst > 6.75) continue;
            dst += (mol->atoms[at]->coords[1] - (*bt)->coords[1]) * (mol->atoms[at]->coords[1] - (*bt)->coords[1]);
            if (dst > 6.75) continue;
            dst += (mol->atoms[at]->coords[2] - (*bt)->coords[2]) * (mol->atoms[at]->coords[2] - (*bt)->coords[2]);
            if (dst > 6.75) continue;
            if (mol->atoms[at]->element == "H" || (*bt)->element == "H") {
                if (dst > 1.8225) continue;
            }

            if (!EDnew::EDAccess::connection_is_possible(mol->atoms[at]->element, (*bt)->element)) continue;

            dst = std::sqrt(dst);

            if (dst < EDnew::ED_min_dist || dst > EDnew::ED_max_dist) continue;

            float prob[5] = { EDnew::EDAccess::get_probability(0,dst),  // ar
                EDnew::EDAccess::get_probability(1,dst),  // 1
                EDnew::EDAccess::get_probability(2,dst),  // 2
                EDnew::EDAccess::get_probability(3,dst) }; // 3

            if (prob[1] > mol->atoms[at]->sp3) mol->atoms[at]->sp3 = prob[1];
            if (prob[1] > (*bt)->sp3) (*bt)->sp3 = prob[1];
            if (prob[2] > mol->atoms[at]->sp2) mol->atoms[at]->sp2 = prob[2];
            if (prob[2] > (*bt)->sp2) (*bt)->sp2 = prob[2];
            if (prob[0] > mol->atoms[at]->sp2) mol->atoms[at]->sp2 = prob[0];
            if (prob[0] > (*bt)->sp2) (*bt)->sp2 = prob[0];
            if (prob[3] > mol->atoms[at]->sp) mol->atoms[at]->sp = prob[3];
            if (prob[3] > (*bt)->sp) (*bt)->sp = prob[3];

            if (prob[0] > 0. || prob[1] > 0. || prob[2] > 0. || prob[3] > 0.) {

                Bond* new_bond = new Bond();
                new_bond->id = bond_id;
                new_bond->connected_atom_1 = mol->atoms[at];
                new_bond->connected_atom_2 = *bt;
                new_bond->length = dst;

                mol->bonds.push_back(new_bond);
                bond_id++;
            }
        }
    }
    delete pwp;

    for (int at = 0; at < mol->bonds.size(); ++at)
    {
        Bond* current_bond = mol->bonds[at];
        current_bond->connected_atom_1->atom_bonds.push_back(current_bond);
        current_bond->connected_atom_2->atom_bonds.push_back(current_bond);
    }

    for (int at = 0; at < mol->atoms.size(); ++at) {

        if (mol->atoms[at]->element == "H") {
            if (mol->atoms[at]->atom_bonds.size() == 0 && mol->atoms[at]->pos_metal) {
                mol->atoms[at]->type = 3;
                mol->atoms[at]->element = "Hg";
                continue;
            }
            else while (mol->atoms[at]->atom_bonds.size() > 1) { //!H mit mehr als einem gebundenen Atom
                delete_worst_bonded(mol->atoms[at]);
                continue;
            }
        }
        else if (mol->atoms[at]->element == "C" || mol->atoms[at]->element == "N" || mol->atoms[at]->element == "Si" ||
            mol->atoms[at]->element == "Se" || mol->atoms[at]->element == "As") {
            if (mol->atoms[at]->atom_bonds.size() > 4) {
                do {
                    delete_worst_bonded(mol->atoms[at]);
                } while (mol->atoms[at]->atom_bonds.size() > 4);
            }
            else if (mol->atoms[at]->atom_bonds.size() == 0) {
                if (mol->atoms[at]->pos_metal) {
                    if (mol->atoms[at]->element == "C") {
                        mol->atoms[at]->type = 3;
                        mol->atoms[at]->element = "Ca";
                    }
                    else if (mol->atoms[at]->element == "N") {
                        mol->atoms[at]->type = 3;
                        if (mol->atoms[at]->name.size() > 1) {
                            if (mol->atoms[at]->name[1] == 'I' || mol->atoms[at]->name[1] == 'i') mol->atoms[at]->element = "Ni";
                            else if (mol->atoms[at]->name[1] == 'A' || mol->atoms[at]->name[1] == 'a') mol->atoms[at]->element = "Na";
                        }
                    }
                }
            }
        }
        else if (mol->atoms[at]->element == "S") {
            if (mol->atoms[at]->atom_bonds.size() > 6) {
                do {
                    delete_worst_bonded(mol->atoms[at]);
                } while (mol->atoms[at]->atom_bonds.size() > 6);
            }
        }
        else if (mol->atoms[at]->element == "P") {
            if (mol->atoms[at]->atom_bonds.size() > 5) {
                do {
                    delete_worst_bonded(mol->atoms[at]);
                } while (mol->atoms[at]->atom_bonds.size() > 5);
            }
        }
        else if (mol->atoms[at]->element == "O") {
            while (mol->atoms[at]->atom_bonds.size() > 2) delete_worst_bonded(mol->atoms[at]);
        }
        else if (mol->atoms[at]->element == "B") {
            while (mol->atoms[at]->atom_bonds.size() > 3) delete_worst_bonded(mol->atoms[at]);
        }
        else while (mol->atoms[at]->atom_bonds.size() > 1) delete_worst_bonded(mol->atoms[at]);
    }
    for (int at = 0; at < mol->atoms.size(); ++at) {
        mol->atoms[at]->n_heavy_bonded = 0;
        for (int bt = 0; bt < mol->atoms[at]->atom_bonds.size(); ++bt) {
            if (mol->atoms[at]->get_connected_atom(bt)->element != "H")
            {
                mol->atoms[at]->n_heavy_bonded++;
            }
        }
    }
}

bool Prepare_mol::bonded_to_HD(Atom* a)
{
    for (int i = 0; i < a->atom_bonds.size(); ++i)
    {
        if (a->get_connected_atom(i)->ad == AD_TYPE_HD)
        {
            return true;
        }
    }
    return false;
}

bool Prepare_mol::bonded_to_heteroatom(Atom* a)
{
    for (int i = 0; i < a->atom_bonds.size(); ++i)
    {
        if (a->get_connected_atom(i)->is_heteroatom())
        {
            return true;
        }
    }
    return false;
}

void Prepare_mol::assign_types(Molecule* mol)
{
    for (int i = 0; i < mol->atoms.size(); ++i)
    {
        Atom* a = mol->atoms[i];
        a->assign_el();
        int& x = a->xs;

        bool acceptor = (a->ad == AD_TYPE_OA || a->ad == AD_TYPE_NA || a->ad == AD_TYPE_SA);
        bool donor_NorO = (a->el == EL_TYPE_Met || bonded_to_HD(a));

        switch (a->el) {
        case EL_TYPE_H: break;
        case EL_TYPE_C: x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H; break;
        case EL_TYPE_N: x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
        case EL_TYPE_O: x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
        case EL_TYPE_S: x = acceptor ? XS_TYPE_S_A : XS_TYPE_S_P; break;
        case EL_TYPE_P: x = XS_TYPE_P_P; break;
        case EL_TYPE_F: x = XS_TYPE_F_H; break;
        case EL_TYPE_Cl: x = XS_TYPE_Cl_H; break;
        case EL_TYPE_Br: x = XS_TYPE_Br_H; break;
        case EL_TYPE_I: x = XS_TYPE_I_H; break;
        case EL_TYPE_Met: x = XS_TYPE_Met_D; break;
        case EL_TYPE_SIZE: break;
        default: throw Error_report("Internal Error");
        }
    }
}
