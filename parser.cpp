#include "parser.h"

void Parser::initialize_fconv() {
    if (Parser::vdW_map.size() > 0) return;
    for (int i = 0; i<number_of_elements; ++i)
    {
        Parser::vdW_map[vdW_elements[i]] = vdW_radii[i];
        Parser::square_vdW_map[vdW_elements[i]] = vdW_radii[i] * vdW_radii[i];
        Parser::clash_map[vdW_elements[i]] = clash_radii[i];
        Parser::square_clash_map[vdW_elements[i]] = clash_radii[i] * clash_radii[i];
        Parser::mol_weight[vdW_elements[i]] = mweights[i];
    }
    for (int i = 0; i<31; ++i) Parser::known_aacids.insert(normal_aacids[i]);
    for (int i = 0; i<11; ++i) Parser::known_mod_aacids.insert(modified_aacids[i]);
    for (int i = 0; i<11; ++i) Parser::known_nucleic_acids.insert(nucleic_acids[i]);
}

//eine HETATM-Line einer pdb-Datei auswerten
bool Parser::parse_pdb_hetatm_string(Atom* atom_, const std::string& str, int atom_count, bool parse_hydrogens, bool parse_ligands, bool parse_water, bool& is_pdbqt)
{
    std::istringstream iss;

    if (!parse_water) { //kein Wasser mit einlesen
        if (str[13] == 'O') {
            if ((str[17] == 'H') && (str[18] == 'O') && (str[19] == 'H')) return false;
            if ((str[17] == 'H') && (str[18] == '2') && (str[19] == 'O')) return false;
        }
    }
    if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
        if ((str[13] == 'H') || (str[12] == 'H')) {
            return false;
        }
    }

    //!18.08.06:
    std::string s2(str, 12, 4);
    if (s2 == " DUM" || (s2[0] == ' ' && s2[1] == 'Q')) return false;

    atom_->type = 2;
    atom_->intern_id = atom_count; //fortlaufende id fr interne Verwaltung

    std::string s1(str, 6, 5), s4(str, 17, 3), s6(str, 22, 4), //id / name / res_name / res_number
        s7(str, 30, 8), s8(str, 38, 8), s9(str, 46, 8); //Koordinaten

    //! 07.10.09: weil modifizierte AS als HETATM abgelegt werden: =======
    if (!parse_ligands) {
        if (Parser::known_mod_aacids.find(s4) == Parser::known_mod_aacids.end()) {

            return false;
        }
        else {

            bool valid_atom = parse_pdb_atom_string(atom_, str, atom_count, parse_hydrogens, parse_ligands, parse_water, is_pdbqt);

            if (valid_atom == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else {
        if (Parser::known_mod_aacids.find(s4) != Parser::known_mod_aacids.end()) {

            bool valid_atom = parse_pdb_atom_string(atom_, str, atom_count, parse_hydrogens, parse_ligands, parse_water, is_pdbqt);

            if (valid_atom == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if (Parser::known_nucleic_acids.find(s4) != Parser::known_nucleic_acids.end()) {

            bool valid_atom = parse_pdb_atom_string(atom_, str, atom_count, parse_hydrogens, parse_ligands, parse_water, is_pdbqt);

            if (valid_atom == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    //!===================================================================

    std::string s10, s11, s12, s13;
    try { //siehe ATOM
        if (is_pdbqt)
        {
            s10.assign(str, 54, 6); s11.assign(str, 60, 6); s12.assign(str, 77, 2); s13.assign(str, 70, 6); //occup. / b_factor / element / charge
        }
        else
        {
            s10.assign(str, 54, 6); s11.assign(str, 60, 6); s12.assign(str, 76, 2); s13.assign(str, 78, 6); //occup. / b_factor / element / charge
        }
    }
    catch (...) {
        throw Error_report("Error: Parsing pdb line: " + str);
        return false;
    };

    iss.clear(); //streamobjekt zuruecksetzen
    iss.str(s1); //string dem streamobjekt zuweisen
    iss >> atom_->id;     //Die "Seriennummer" des Atoms

    atom_->name = s2;  //Atom-name  --  z.B.  CA  fr C-alpha
    Common::rmvWhiteSpaces(atom_->name);
    atom_->alt_loc_id = str[16]; //alternative location fr das Atom?
    atom_->res_name = s4; //residue-name
    atom_->chain_id = str[21]; //chain_id
    iss.clear(); iss.str(s6); iss >> atom_->res_number; //integer
    iss.clear(); iss.str(s7); iss >> atom_->coords[0]; //x
    iss.clear(); iss.str(s8); iss >> atom_->coords[1]; //y
    iss.clear(); iss.str(s9); iss >> atom_->coords[2]; //z

    if (is_pdbqt)
    {
        std::string ad_string = s12;
        Common::rmvWhiteSpaces(ad_string);

        int ad = string_to_ad_type(ad_string);

        atom_->ad = ad;

        if (is_non_ad_metal_name(ad_string))
        {
            atom_->xs = XS_TYPE_Met_D;
        }

        if (!(atom_->acceptable_type()))
        {
            throw Error_report(" ( " + ad_string + " ) is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
        }
    }

    try {
        iss.clear();

        if (is_pdbqt)
        {
            if (s12 == "Cl" || s12 == "Br") atom_->element = s12;
            else if (s12[0] == 'A' || s12[0] == 'C') atom_->element = "C";
            else if (s12[0] == 'O') atom_->element = "O";
            else if (s12[0] == 'N') atom_->element = "N";
            else if (s12[0] == 'S') atom_->element = "S";
            else if (s12[0] == 'P') atom_->element = "P";
            else if (s12[0] == 'H') atom_->element = "H";
            else if (s12[0] == 'F') atom_->element = "F";
            else if (s12[0] == 'I') atom_->element = "I";
            else
            {
                atom_->element.assign(1, s12[0]);
                if (atom_->element.size() > 1)
                {
                    atom_->element = s12;
                    if (s12[1] < 97) atom_->element[1] += 32;
                }
            }
        }
        else
        {
            atom_->element = s12; //element
        }

        if (s13.size() > 0 && s13 != "  ") atom_->charge = Common::string_to<float>(s13);
    }
    catch (...) {
        throw Error_report("Error: Parsing pdb line: " + str);
        return false;
    };

    for (int i = 0; i < n_metal_list; ++i)
    {
        if (atom_->res_name == metal_list_ra[i]) { //handelt es sich um ein Metall, dass in elements.hpp bercksichtigt ist?
            atom_->type = 3;
        }
    }

    // jetzt kommt eigentlich nur noch kleines Ion oder Wasser in frage
    if (str[13] == 'O' || str[13] == 'H')
    {
        if (((str[17] == 'H') && (str[18] == 'O') && (str[19] == 'H')) ||
            ((str[17] == 'H') && (str[18] == '2') && (str[19] == 'O'))) {
            atom_->type = 4;
        }
    }

    return true;
}


bool Parser::parse_pdb_atom_string(Atom* atom_, const std::string& str, int atom_count, bool parse_hydrogens, bool parse_ligands, bool parse_water, bool& is_pdbqt)
{

    std::istringstream is;

    //! wenn kein H gelesen wird muessen noch alle anderen id's angepasst werden
    //! => noch Variable anlegen, die die H's zaehlt und den Wert von folgenden id's abzieht
    //!    (auch im CONECT Teil nicht vergessen!)
    //! (erst beim Schreiben bercksichtigen?)
    //! gleiches muss eventuell auch noch fuer HOH implementiert werden!

    //!11.04.07:  Ketten nach einem HETATM trotzdem weiterparsen, obwohl dann auch Peptide mit geparsed werden
    //!if (hetatoms && !(parse_ligands)) return;

    if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
        //!GN 03.05.06:  auch row[12]
        if ((str[13] == 'H') || (str[12] == 'H')) { //kann in diesem Fall nur Wasserstoff sein
            return false;
        }
    }

    //!18.08.06:
    std::string s2(str, 12, 4);
    if ((s2 == " DUM") || (s2[0] == ' ' && s2[1] == 'Q')) return false;

    atom_->type = 1;
    atom_->intern_id = atom_count; //fortlaufende id fuer interne Verwaltung

    //line in einzelne strings zerlegen
    std::string s1(str, 6, 5), s4(str, 17, 3), s6(str, 22, 4), //id / name / res_name / res_number
        s7(str, 30, 8), s8(str, 38, 8), s9(str, 46, 8); //Koordinaten

    //! 07.10.09: Fuer den Fall das die Kristallographen wieder vergessen haben HETATMs auch
    //! als solche zu kennzeichnen
    if (Parser::known_aacids.find(s4) == Parser::known_aacids.end())
    {
        if (Parser::known_mod_aacids.find(s4) == Parser::known_mod_aacids.end())
        {
            if (Parser::known_nucleic_acids.find(s4) == Parser::known_nucleic_acids.end())
            {
                bool valid_atom = parse_pdb_hetatm_string(atom_, str, atom_count, parse_hydrogens, parse_ligands, parse_water, is_pdbqt);

                if (valid_atom == true)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
    }

    std::string s10, s11, s12, s13;

    try {

        if (is_pdbqt)
        {
            s10.assign(str, 54, 6); s11.assign(str, 60, 6); s12.assign(str, 77, 2); s13.assign(str, 70, 6); //occup. / b_factor / element / charge
        }
        else
        {
            s10.assign(str, 54, 6); s11.assign(str, 60, 6); s12.assign(str, 76, 2); s13.assign(str, 78, 6); //occup. / b_factor / element / charge
        }
    }
    catch (...) {
        throw Error_report("Error: Parsing pdb line: " + str);
        return false;
    };

    is.clear();
    is.str(s1);

    is >> atom_->id;    //Die "Seriennummer" des Atoms
    atom_->name = s2;  //Atom-name  --  z.B.  CA  fuer C-alpha
    Common::rmvWhiteSpaces(atom_->name);
    atom_->alt_loc_id = str[16]; //alternative location fuer das Atom?
    atom_->res_name = s4; //residue-name
    atom_->chain_id = str[21]; //chain_id
    is.clear(); is.str(s6); is >> atom_->res_number; //integer

    //!Kommentar: Es gibt Probleme, wenn sich im pdb 2 chains (z.B. A und B) staendig abwechseln und keine TER-Eintraege
    //!           vorhanden sind und in der Kette AS fehlen

    is.clear(); is.str(s7); is >> atom_->coords[0]; //x
    is.clear(); is.str(s8); is >> atom_->coords[1]; //y
    is.clear(); is.str(s9); is >> atom_->coords[2]; //z

    if (is_pdbqt)
    {
        std::string ad_string = s12;
        Common::rmvWhiteSpaces(ad_string);

        int ad = string_to_ad_type(ad_string);

        atom_->ad = ad;

        if (is_non_ad_metal_name(ad_string))
        {
            atom_->xs = XS_TYPE_Met_D;
        }

        if (!(atom_->acceptable_type()))
        {
            throw Error_report(" ( " + ad_string + " ) is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
        }
    }

    try { //siehe oben
        is.clear();

        if (is_pdbqt)
        {
            if (s12 == "Cl" || s12 == "Br") atom_->element = s12;
            else if (s12[0] == 'A' || s12[0] == 'C') atom_->element = "C";
            else if (s12[0] == 'O') atom_->element = "O";
            else if (s12[0] == 'N') atom_->element = "N";
            else if (s12[0] == 'S') atom_->element = "S";
            else if (s12[0] == 'P') atom_->element = "P";
            else if (s12[0] == 'H') atom_->element = "H";
            else if (s12[0] == 'F') atom_->element = "F";
            else if (s12[0] == 'I') atom_->element = "I";
            else
            {
                atom_->element.assign(1, s12[0]);
                if (atom_->element.size() > 1)
                {
                    atom_->element = s12;
                    if (s12[1] < 97) atom_->element[1] += 32;
                }
            }
        }
        else
        {
            atom_->element = s12; //element
        }

        if (s13.size() > 0 && s13 != "  ") atom_->charge = Common::string_to<float>(s13); //!so gut wie nie vorhanden
    }
    catch (...)
    {
        throw Error_report("Error: Parsing pdb line: " + str);
        return false;
    };


    return true;
}

void Parser::parse_pdb(const std::string& file_path, Molecule* mol, bool is_pdbqt) {

    std::ifstream inputfile;
    inputfile.exceptions(std::ifstream::badbit);
    unsigned int count = 0;
    unsigned int atom_count = 0;
    unsigned int num_of_chains_ter = 1;
    unsigned int branch_id = 1;

    try
    {
        inputfile.open(file_path);
        std::string line;
        while (getline(inputfile, line))
        {
            count++;

            if (line.empty()) {}
            else if (Common::starts_with(line, "TER")) { num_of_chains_ter++; }
            else if (Common::starts_with(line, "ATOM"))
            {
                try
                {
                    Atom* new_atom = new Atom();
                    bool correct_atom = parse_pdb_atom_string(new_atom, line, atom_count, true, true, true, is_pdbqt);
                    new_atom->branch_id = branch_id;

                    if (correct_atom == true)
                    {
                        if(new_atom->type == 4)
                        {
                            mol->waters.push_back(new_atom);
                        }
                        else
                        {
                            mol->atoms.push_back(new_atom);
                        }

                        atom_count++;
                    }
                }
                catch (Error_report& err)
                {
                    throw Error_report("Error: Line " + std::to_string(count) + ": pdb ATOM syntax incorrect: " + err.errorText);
                }
                catch (...) {

                    throw Error_report("Error: Line " + std::to_string(count) + ": pdb ATOM syntax incorrect");
                }
            }
            else if (Common::starts_with(line, "HETATM"))
            {
                try
                {
                    Atom* new_atom = new Atom();
                    bool correct_atom = parse_pdb_hetatm_string(new_atom, line, atom_count, true, true, true, is_pdbqt);
                    new_atom->branch_id = branch_id;

                    if (correct_atom == true)
                    {
                        if(new_atom->type == 4)
                        {
                            mol->waters.push_back(new_atom);
                        }
                        else
                        {
                            mol->atoms.push_back(new_atom);
                        }

                        atom_count++;
                    }
                }
                catch (Error_report& err)
                {
                    throw Error_report("Error: Line " + std::to_string(count) + ": pdb ATOM syntax incorrect: " + err.errorText);
                }
                catch (...) {

                    throw Error_report("Error: Line " + std::to_string(count) + ": pdb ATOM syntax incorrect");
                }
            }
            else if (Common::starts_with(line, "BRANCH"))
            {
                branch_id++;
            }
            else if (Common::starts_with(line, "MODEL")) { throw Error_report("Error: Line " + std::to_string(count) + ": Unexpected multi-MODEL input. Split models first?"); }
            else {}
        }

        inputfile.close();
    }
    catch (std::ifstream::failure e)
    {
        inputfile.close();
        throw Error_report(file_path + " " + e.what());
    }
    catch (Error_report& e) {
        inputfile.close();
        throw Error_report(file_path + " " + e.errorText);
    }
    catch (...)
    {
        inputfile.close();
        throw Error_report(file_path + " ");
    }
}


void Parser::parse_sdf(const std::string& file_path, Molecule* mol)
{
    std::istringstream is;

    std::ifstream inputfile;
    inputfile.exceptions(std::ifstream::badbit);

    try {

        inputfile.open(file_path);

        int n_parsed = 0;
        std::string row, key;
        int n_atoms = 0;

        getline(inputfile, row);
        if (inputfile.eof()) return;

        getline(inputfile, row);
        if (inputfile.eof()) return;
        getline(inputfile, row);
        if (inputfile.eof()) return;

        int a_intern_id = 0;
        getline(inputfile, row);


        std::string s1(row, 0, 3);
        is.clear(); is.str(s1);
        n_atoms = 0;
        is >> n_atoms;

        //jetzt die Atome einlesen:
        for (int i = 0; i < n_atoms; ++i)
        {
            getline(inputfile, row);
            if (inputfile.eof()) {
                return;
            }

            Atom* atm = new Atom();
            std::string sx(row, 0, 10), sy(row, 10, 10), sz(row, 20, 10), se(row, 31, 3);
            is.clear(); is.str(sx);
            is >> atm->coords[0];
            is.clear(); is.str(sy);
            is >> atm->coords[1];
            is.clear(); is.str(sz);
            is >> atm->coords[2];
            is.clear(); is.str(se);
            is >> atm->element;

            if (atm->element.size() > 1) {
                if (atm->element[1] < 97) atm->element[1] += 32;
            }

            atm->name = atm->element;
            atm->type = 5;
            atm->intern_id = a_intern_id; a_intern_id++;
            atm->id = atm->intern_id;
            atm->res_number = 0;
            mol->atoms.push_back(atm);
        }

        inputfile.close();
    }
    catch (std::ifstream::failure e)
    {
        inputfile.close();
        throw Error_report(file_path + " " + e.what());
    }
    catch (Error_report& e) {
        inputfile.close();
        throw Error_report(file_path + " " + e.errorText);
    }
    catch (...)
    {
        inputfile.close();
        throw Error_report(file_path + " ");
    }
}

void Parser::reset_mol2_flags(bool& atom_reading, bool& bond_reading)
{
    atom_reading = false;
    bond_reading = false;
}

void Parser::parse_mol2(const std::string& file_path, Molecule* mol, bool atoms_only, std::vector<std::pair<std::pair<int, int>, std::string>>& mol2_bonds)
{
    std::ifstream inputfile;
    inputfile.exceptions(std::ifstream::badbit);
    int count = 0;
    std::istringstream is;
    std::istringstream is_tmp;
    std::string dummy;
    std::map<int, int> id_index_atom_map;

    try
    {
        inputfile.open(file_path);
        std::string line;
        bool atom_reading = false;
        bool bond_reading = false;
        bool molecule_found = false;
        bool atom_found = false;
        bool bond_found = false;
        int atom_count = 0;
        int bond_count = 0;

        while (getline(inputfile, line))
        {
            count++;

            if (line.empty()) {}
            else if (Common::starts_with(line, "@<TRIPOS>MOLECULE"))
            {
                if (molecule_found == true)
                {
                    throw Error_report("Error: Line " + std::to_string(count) + ": Unexpected mol2 multi-@<TRIPOS>MOLECULE input. Split molecules?");
                }

                molecule_found = true;
                continue;
            }
            else if (Common::starts_with(line, "@<TRIPOS>ATOM"))
            {
                if (atom_found == true)
                {
                    throw Error_report("line " + std::to_string(count) + ": Unexpected mol2 multi-@<TRIPOS>ATOM input. Split molecules?");
                }

                atom_found = true;
                atom_reading = true;
                bond_reading = false;
                continue;
            }
            else if (Common::starts_with(line, "@<TRIPOS>BOND"))
            {
                if (!atoms_only)
                {
                    if (bond_found)
                    {
                        throw Error_report("line " + std::to_string(count) + ": Unexpected mol2 multi-@<TRIPOS>BOND input. Split molecules?");
                    }

                    bond_found = true;
                    bond_reading = true;
                }

                atom_reading = false;
                continue;
            }
            else if (Common::starts_with(line, "@<TRIPOS>ALT_TYPE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>ANCHOR_ATOM")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>ASSOCIATED_ANNOTATION")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>CENTER_OF_MASS")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>CENTROID")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>COMMENT")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>CRYSIN")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>DATA_FILE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>DICT")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>EXTENSION_POINT")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FF_PBC")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FFCON_ANGLE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FFCON_DIST")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FFCON_MULTI")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FFCON_RANGE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>FFCON_TORSION")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>LINE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>LSPLANE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>NORMAL")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>QSAR_ALIGN_RULE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>RING_CLOSURE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>ROTATABLE_BOND")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>SEARCH_DIST")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>SEARCH_OPTS")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>SET")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>SUBSTRUCTURE")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>U_FEAT")) {}
            else if (Common::starts_with(line, "@<TRIPOS>UNITY_ATOM_ATTR")) { reset_mol2_flags(atom_reading, bond_reading); }
            else if (Common::starts_with(line, "@<TRIPOS>UNITY_BOND_ATTR")) { reset_mol2_flags(atom_reading, bond_reading); }

            if (atom_reading)
            {
                try
                {
                    if (line == "" || line[0] == '#') continue;
                    is.clear();
                    is.str(line);

                    Atom* new_atom = new Atom();

                    int adder = 0;

                    new_atom->type = 5;
                    is >> new_atom->id >> new_atom->name;

                    new_atom->intern_id = atom_count; //fortlaufende id fuer interne Verwaltung
                    atom_count += 1; //fortlaufende Nummer hochzaehlen

                    is >> new_atom->coords[0] >> new_atom->coords[1] >> new_atom->coords[2];
                    is >> new_atom->sybyl_type;

                    //!NEU: 11.03.2010:  Lonepairs NICHT mit einlesen
                    if (new_atom->sybyl_type == "LP" || new_atom->sybyl_type == "Lp") {
                        if(atom_count > 0)
                        {
                            atom_count -= 1;
                        }
                        continue;
                    }

                    if (!is.eof()) {
                        //    is >> atm->res_number; //als int gelesen
                        //! 14.8.09: Leider steht da nicht immer ein int => kleiner Umstand
                        is >> dummy;
                        is_tmp.clear(); is_tmp.str(dummy);
                        is_tmp >> new_atom->res_number;

                        if (new_atom->res_number == 0) adder = 1;
                        new_atom->res_number += adder;
                    }

                    if (!is.eof()) {
                        is >> new_atom->res_name;

                        //! neu: 08.04.09: Laenge des resname auf 3 beschraenken (fuer Sequenzalignment):
                        //! --> NEIN! : lieber nur dort (im Alignment) beschneiden !!!
                    }
                    if (!is.eof()) is >> new_atom->charge; //als string gelesen

                    id_index_atom_map.insert(std::pair<int, int>(new_atom->id, mol->atoms.size()));
                    mol->atoms.push_back(new_atom);
                }
                catch (Error_report& err)
                {
                    throw Error_report("line " + std::to_string(count) + ": Mol2 ATOM syntax incorrect: " + err.errorText);
                }
                catch (...) {

                    throw Error_report("line " + std::to_string(count) + ": Mol2 ATOM syntax incorrect");
                }

            }


            if (bond_reading)
            {
                try
                {
                    int from, to;
                    if (line == "" || line[0] == '#') continue;
                    is.clear(); //streamobjekt zurcksetzen
                    is.str(line);

                    int bond_id;
                    std::string bond_type;

                    is >> bond_id >> from >> to >> bond_type;

                    if (from == to)
                    {
                        continue;
                    }

                    std::map<int, int>::iterator it;

                    it = id_index_atom_map.find(from);
                    if (it == id_index_atom_map.end())
                    {
                        continue;
                    }

                    it = id_index_atom_map.find(to);
                    if (it == id_index_atom_map.end())
                    {
                        continue;
                    }

                    mol2_bonds.push_back(std::pair<std::pair<int, int>, std::string>(std::pair<int, int>(id_index_atom_map.find(from)->second, id_index_atom_map.find(to)->second), bond_type));

                    bond_count++;
                }
                catch (Error_report& err)
                {
                    throw Error_report("line " + std::to_string(count) + ": Mol2 BOND syntax incorrect: " + err.errorText);
                }
                catch (...) {

                    throw Error_report("line " + std::to_string(count) + ": Mol2 BOND syntax incorrect");
                }

            }
        }

        inputfile.close();
    }
    catch (std::ifstream::failure e)
    {
        inputfile.close();
        throw Error_report(file_path + " " + e.what());
    }
    catch (Error_report& e) {
        inputfile.close();
        throw Error_report(file_path + " " + e.errorText);
    }
    catch (...)
    {
        inputfile.close();
        throw Error_report(file_path + " ");
    }
}

void Parser::add_mol2_bonds(Molecule* mol, std::vector<std::pair<std::pair<int, int>, std::string>>& mol2_bonds)
{
    int bond_count = 0;

    for(int i = 0; i < mol2_bonds.size(); ++i)
    {
        Bond* new_bond = new Bond();
        new_bond->connected_atom_1 = mol->atoms[mol2_bonds[i].first.first];
        new_bond->connected_atom_2 = mol->atoms[mol2_bonds[i].first.second];
        float r2 = Common::vec_distance_sqr(new_bond->connected_atom_1->coords, new_bond->connected_atom_2->coords);
        new_bond->length = std::sqrt(r2);
        new_bond->id = bond_count;
        mol->bonds.push_back(new_bond);
        bond_count++;
    }
}
