#ifndef ATOM_H
#define ATOM_H

#include "atom_constants.h"

struct Atom_type {
	enum t { EL, AD, XS, SY };
    int el, ad, xs, sy;
    Atom_type() : el(EL_TYPE_SIZE), ad(AD_TYPE_SIZE), xs(XS_TYPE_SIZE), sy(SY_TYPE_SIZE) {}
    int get(t atom_typing_used) const {
		switch (atom_typing_used) {
		case EL: return el;
		case AD: return ad;
		case XS: return xs;
		case SY: return sy;
        default: Q_ASSERT(false); return max_int;
		}
	}
	bool is_hydrogen() const {
        return ad_is_hydrogen(ad);
	}
	bool is_heteroatom() const {
        return ad_is_heteroatom(ad) || xs == XS_TYPE_Met_D;
	}
	bool acceptable_type() const {
		return ad < AD_TYPE_SIZE || xs == XS_TYPE_Met_D;
	}
	void assign_el() {
        el = ad_type_to_el_type(ad);
		if (ad == AD_TYPE_SIZE && xs == XS_TYPE_Met_D)
			el = EL_TYPE_Met;
	}
    bool same_element(const Atom_type& a) const { // does not distinguish metals or unassigned types
		return el == a.el;
	}
    float covalent_radius() const {
        if (ad < AD_TYPE_SIZE) { return ad_type_property(ad).covalent_radius; }
		else if (xs == XS_TYPE_Met_D) { return metal_covalent_radius; }
		Q_ASSERT(false);          return 0; // never happens - placating the compiler
	}
    float optimal_covalent_bond_length(const Atom_type* x) const {
        return covalent_radius() + x->covalent_radius();
	}
};

inline int num_atom_types(Atom_type::t atom_typing_used) {
	switch (atom_typing_used) {
    case Atom_type::EL: return EL_TYPE_SIZE;
    case Atom_type::AD: return AD_TYPE_SIZE;
    case Atom_type::XS: return XS_TYPE_SIZE;
    case Atom_type::SY: return SY_TYPE_SIZE;
    default: Q_ASSERT(false); return max_int;
	}
}

struct Atom : Atom_type {

    unsigned int id;
    unsigned int intern_id; //interne fortlaufende Nummer (id soll spaeter ohne h_correction auskommen!)
    unsigned int branch_id;
    std::string name; //12-15 im pdb
    std::string res_name;
    std::string chain_id;
    std::string element;
    unsigned int type; //1: pdb-ATOM / 2: pdb-HETATM / 3: pdb-metal / 4: pdb-water / 5:mol2-Atom
    bool pos_metal; //!bei dem Atom koennte es sich aufgrund des namens auch um ein Metall handeln
    unsigned int n_heavy_bonded; //! 0:nicht ermittelt
    float sp; //Summe der Wahrscheinlichkeiten fuer sp-Hybridisierung
    float sp2;
    float sp3;
    glm::vec3 coords;
    float charge;
    std::vector<Bond*> atom_bonds;
    std::string intern_type;
    std::string sybyl_type; //Der mol2-Atomtype
    unsigned int res_number; // 22-25 im pdb / sub_id vom mol2
    char alt_loc_id; //16 im pdb
    unsigned int xs_donor_hydrogen_pre_count;

    Atom() : branch_id(1), n_heavy_bonded(0), charge(0), sp(0.), sp2(0.), sp3(0.), alt_loc_id(' '), pos_metal(false), intern_type("X"), type(5), sybyl_type("XX"), element("X"), xs_donor_hydrogen_pre_count(0)
    {

    }

    void get_element(bool reset_names = true)
    {
        std::string ele_help = element;
        std::ostringstream os2;
        os2 << name;
        element = os2.str();
        if (!(sybyl_type == "XX")) {
            if (sybyl_type == "H") { element = "H"; return; }
            else if (sybyl_type.size() > 1) {
                if (sybyl_type[0] == 'C' && sybyl_type[1] == '.') { element = "C"; return; }
                else if (sybyl_type[0] == 'N' && sybyl_type[1] == '.') { element = "N"; return; }
                else if (sybyl_type[0] == 'O' && sybyl_type[1] == '.') { element = "O"; return; }
                else if (sybyl_type[0] == 'S' && sybyl_type[1] == '.') { element = "S"; return; }
                else if (sybyl_type[0] == 'P' && sybyl_type[1] == '.') { element = "P"; return; }
                else if (sybyl_type[0] == 'N' && sybyl_type[1] == '.') { element = "N"; return; }
                else if (sybyl_type[1] == '.') { element = sybyl_type[0]; return; }
                else {
                    if (sybyl_type[1] < 91) sybyl_type[1] += 32;
                    if (sybyl_type == "Lp") {
                        sybyl_type = "LP";
                        element = sybyl_type;
                        return;
                    }
                    element.assign(sybyl_type, 0, 2);
                    if (element == "Cl" || element == "Br" || element == "Du") return;
                    else type = 3;
                    return;
                }
            }
            else {
                if (sybyl_type == "X") { //den hat fconv dann schonmal auf X gesetzt
                    if (name.size() > 0) {
                        if (name[0] == 'C') { element = "C"; return; }
                        if (name[0] == 'N') { element = "N"; return; }
                        if (name[0] == 'O') { element = "O"; return; }
                        if (name[0] == 'S') { element = "S"; return; }
                        if (name[0] == 'P') { element = "P"; return; }
                    }
                }
                else {
                    element = sybyl_type;
                    if (!(element == "C" || element == "N" || element == "O" || element == "S" ||
                        element == "P" || element == "F" || element == "I" || element == "B")) type = 3;
                    return;
                }
            }
        }
        else {
            if (res_name == "HOH" || res_name == "H2O") { //! neu: 10.10.2008
                element = "O"; type = 4; return;
            }

            if (ele_help != "X") {
                if (ele_help.size() == 2) {
                    if (ele_help == "C " || ele_help == "O " || ele_help == "N " ||
                        ele_help == "P " || ele_help == "S " || ele_help == "H " ||
                        ele_help == "F " || ele_help == "I " || ele_help == "B " ||
                        ele_help == "K ") {
                        element.assign(1, ele_help[0]); return;
                    }
                    else if (ele_help == " C" || ele_help == " O" || ele_help == " N" ||
                        ele_help == " P" || ele_help == " S" || ele_help == " H" ||
                        ele_help == " F" || ele_help == " I" || ele_help == " B" ||
                        ele_help == " K") {
                        element.assign(1, ele_help[1]); return;
                    }
                    else if (ele_help == "Cl" || ele_help == "Br") { element = ele_help; return; }
                    else if (ele_help == "Na" || ele_help == "Mg" || ele_help == "Ca" ||
                        ele_help == "Al" || ele_help == "Si" || ele_help == "Cr" ||
                        ele_help == "Mn" || ele_help == "Fe" || ele_help == "Co" ||
                        ele_help == "Ni" || ele_help == "Cu" || ele_help == "Zn" ||
                        ele_help == "Se" || ele_help == "Ag" || ele_help == "Au" ||
                        ele_help == "Sn" || ele_help == "Pt" || ele_help == "Hg" ||
                        ele_help == "Pb" || ele_help == "Bi") {
                        element = ele_help;
                        type = 3;
                        return;
                    }
                    else if (ele_help == "CL") { element = "Cl"; return; }
                    else if (ele_help == "BR") { element = "Br"; return; }
                    else if (ele_help == "MG") { element = "Mg"; type = 3; return; }
                    else if (ele_help == "CA") { element = "Ca"; type = 3; return; }
                    else if (ele_help == "MN") { element = "Mn"; type = 3; return; }
                    else if (ele_help == "CR") { element = "Cr"; type = 3; return; }
                    else if (ele_help == "AL") { element = "Al"; type = 3; return; }
                    else if (ele_help == "SI") { element = "Si"; type = 3; return; }
                    else if (ele_help == "NA") { element = "Na"; type = 3; return; }
                    else if (ele_help == "FE") { element = "Fe"; type = 3; return; }
                    else if (ele_help == "CO") { element = "Co"; type = 3; return; }
                    else if (ele_help == "NI") { element = "Ni"; type = 3; return; }
                    else if (ele_help == "CU") { element = "Cu"; type = 3; return; }
                    else if (ele_help == "ZN") { element = "Zn"; type = 3; return; }
                    else if (ele_help == "AG") { element = "Ag"; type = 3; return; }
                    else if (ele_help == "AU") { element = "Au"; type = 3; return; }
                    else if (ele_help == "SN") { element = "Sn"; type = 3; return; }
                    else if (ele_help == "PT") { element = "Pt"; type = 3; return; }
                    else if (ele_help == "HG") { element = "Hg"; type = 3; return; }
                    else if (ele_help == "PB") { element = "Pb"; type = 3; return; }
                    else if (ele_help == "BI") { element = "Bi"; type = 3; return; }
                    else if (ele_help == "XE") { element = "Xe"; return; }
                }
                else if (ele_help.size() == 1) { //eigentlich doch immer 2 ???
                    if (ele_help == "C" || ele_help == "O" || ele_help == "N" ||
                        ele_help == "P" || ele_help == "S" || ele_help == "H" ||
                        ele_help == "F" || ele_help == "I" || ele_help == "B" ||
                        ele_help == "K") {
                        element = ele_help; return;
                    }
                }
            }

            std::ostringstream os;

            if (name == " UNK") { element = "UNK"; return; }
            if (res_name == "NAD" || res_name == "NDP" || res_name == "GPC") {
                if (name[2] == 'H' || name[2] == 'C' || name[2] == 'N' ||
                    name[2] == 'O' || name[2] == 'S' || name[2] == 'P') {
                    if (name[2] == 'P') {
                        if (name[1] == 'O') { //!Ausnahme bei NDP
                            element.assign(1, name[1]);
                            if (reset_names) {
                                os << name[1] << name[2] << name[3];
                                name = os.str();
                            }
                            return;
                        }
                    }
                    element.assign(1, name[2]);
                    if (reset_names) {
                        os << (name[2]) << (name[3]);
                        name = os.str();
                    }
                }
                else {
                    element.assign(1, name[1]);
                    if (reset_names) {
                        os << name[1] << name[2] << name[3];
                        name = os.str();
                    }
                }
                return;
            }

            unsigned int ta = 0;

            if (name[0] == ' ' || name[0] == '1' || name[0] == '2' ||
                name[0] == '3' || name[0] == '4' || name[0] == '\'' ||
                name[0] == '"' || name[0] == '*') ta = 1;
            if (name[0 + ta] == 'H') {
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'G' || name[1 + ta] == 'g') {
                        pos_metal = true;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                element = "H";
                return;
            }
            if (name[0 + ta] == 'C') {
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'L' || name[1 + ta] == 'l') {
                        element = "Cl";
                        if (reset_names) {
                            if (name.size() > 2 + ta) os << name[0 + ta]
                                << name[1 + ta] << name[2 + ta];
                            else os << name[0 + ta] << name[1 + ta];
                            name = os.str();
                        }
                        return;
                    }
                    if (name[1 + ta] == 'A' || name[1 + ta] == 'a' ||
                        name[1 + ta] == 'U' || name[1 + ta] == 'u') {
                        pos_metal = true;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                element = "C";
                return;
            }
            if (name[0 + ta] == 'N') {
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'I' || name[1 + ta] == 'i') {
                        pos_metal = true;
                    }
                    if (name[1 + ta] == 'A' || name[1 + ta] == 'a') {
                        pos_metal = true;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                element = "N"; return;
            }
            if (name[0 + ta] == 'O') {
                element = "O";
                if (reset_names) {
                    if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                    else os << name[0 + ta] << name[1 + ta];
                    name = os.str();
                }
                return;
            }
            if (name[0 + ta] == 'P') {
                element = "P";
                if (reset_names) {
                    if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                    else os << name[0 + ta] << name[1 + ta];
                    name = os.str();
                }
                return;
            }
            if (name[0 + ta] == 'I') {
                element = "I";
                if (reset_names) {
                    if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                    else os << name[0 + ta] << name[1 + ta];
                    name = os.str();
                }
                return;
            }
            if (name[0 + ta] == 'B') {
                element = "B";
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'R' || name[1 + ta] == 'r') {
                        element = "Br";
                        if (reset_names) {
                            if (name.size() > 2 + ta) os << name[0 + ta]
                                << name[1 + ta] << name[2 + ta];
                            else os << name[0 + ta] << name[1 + ta];
                            name = os.str();
                        }
                        return;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                return;
            }
            if (name[0 + ta] == 'S') {
                element = "S";
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'I' || name[1 + ta] == 'i') {
                        pos_metal = true;
                    }
                    if (name[1 + ta] == 'E' || name[1 + ta] == 'e') {
                        pos_metal = true;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                return;
            }
            if (name[0 + ta] == 'F') {
                element = "F";
                if (name.size() > 1 + ta) {
                    if (name[1 + ta] == 'E' || name[1 + ta] == 'e') {
                        pos_metal = true;
                    }
                    if (reset_names) {
                        if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                        else os << name[0 + ta] << name[1 + ta];
                        name = os.str();
                    }
                }
                return;
            }

            element.assign(1, name[0 + ta]);
            if (reset_names) {
                if (name.size() > 2 + ta) os << name[0 + ta] << name[1 + ta] << name[2 + ta];
                else os << name[0 + ta] << name[1 + ta];
                name = os.str();
            }

            pos_metal = true; //!wenn sonst nicht zuzuordnen ist ein metall recht wahrscheinlich
        }
    }

    inline Atom* get_connected_atom(unsigned int bond_number)
    {
        if (this->atom_bonds[bond_number]->connected_atom_1 == this)
        {
            return this->atom_bonds[bond_number]->connected_atom_2;
        }
        else
        {
            return this->atom_bonds[bond_number]->connected_atom_1;
        }
    }
};


#endif
