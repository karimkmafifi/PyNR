#ifndef PARSER_H
#define PARSER_H
#include "filesystem.h"
#include "molecule.h"


class Parser
{
private:
    std::unordered_map<std::string, float> vdW_map; //!vdW-Radien
    std::unordered_map<std::string, float> square_vdW_map;
    std::unordered_map<std::string, float> clash_map; //!Atom(/Ionen)-Radien
    std::unordered_map<std::string, float> square_clash_map;
    std::unordered_map<std::string, float> mol_weight;

    std::unordered_set<std::string> known_aacids;
    std::unordered_set<std::string> known_mod_aacids;
    std::unordered_set<std::string> known_nucleic_acids;
private:
    bool parse_pdb_hetatm_string(Atom* atom_, const std::string& str, int atom_count, bool parse_hydrogens, bool parse_ligands, bool parse_water, bool& is_pdbqt);
    bool parse_pdb_atom_string(Atom* atom_, const std::string& str, int atom_count, bool parse_hydrogens, bool parse_ligands, bool parse_water, bool& is_pdbqt);
    void reset_mol2_flags(bool& atom_reading, bool& bond_reading);
public:
    void initialize_fconv();
    void parse_pdb(const std::string& file_path, Molecule* mol, bool is_pdbqt);
    void parse_sdf(const std::string& file_path, Molecule* mol);
    void parse_mol2(const std::string& file_path, Molecule* mol, bool atoms_only, std::vector<std::pair<std::pair<int, int>, std::string>>& mol2_bonds);
    void add_mol2_bonds(Molecule* mol, std::vector<std::pair<std::pair<int, int>, std::string>>& mol2_bonds);
};

#endif
