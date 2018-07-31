#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
#include "filesystem.h"

struct Molecule
{
public:
    std::string path;
    std::string format;
    std::vector<Atom*> atoms;
    std::vector<Atom*> waters;
    std::vector<Bond*> bonds;

    ~Molecule()
    {
        for(unsigned int i = 0; i < atoms.size(); ++i)
        {
            delete atoms[i];
        }

        for(unsigned int i = 0; i < waters.size(); ++i)
        {
            delete waters[i];
        }

        for(unsigned int i = 0; i < bonds.size(); ++i)
        {
            delete bonds[i];
        }
    }
};

struct Receptor : Molecule
{
public:
  std::vector<int> rexp;
  std::vector<glm::vec3> rvector;
  std::vector<glm::vec3> rvector2;
};

struct Ligand : Molecule
{
public:
  std::vector<std::pair<std::string, float>> terms_contributions;
  std::vector<std::pair<unsigned int, unsigned int>> rf_inter;
};

#endif // MOLECULE_H
