
//============================================================================
// mw_match_GN.h -*- C++ -*-; maximum weighted match for Graphs where no
//                            nodes with degree > 3 exist
//
// Copyright (C) 2010, 2011 Gerd Neudert
//
// This file is part of fconv.
// fconv is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------------
//
// author:      Gerd Neudert
//              gneudert(place_at_here)web.de
// supervisor:  Prof. Dr. Gerhard Klebe
//              klebe(place_at_here)staff.uni-marburg.de
// Department of Pharmaceutical Chemistry, Philipps-University of Marburg
//
// This library is part of the program fconv. Please see the
// documentation in the file fconv.cpp first!
// The matching only works for graphs where no node has a degree higher
// than 3. This is sufficient in case of the optimization of a network
// of alternating double bonds.
//============================================================================


#ifndef MW_MATCHGN_H
#define MW_MATCHGN_H

#include<vector>


using namespace std;


class MW_EDGE {
public:
    int node_i;
    int node_j;
    float weight;
    bool stacked;
    bool labeled;
    vector<MW_EDGE*> adjacent; // max. 4 benachbarte Kanten
    vector<bool> label;
    MW_EDGE(int const& i,int const& j,float const& w);
    ~MW_EDGE();
    void set_max_label(unsigned int const& ml);
    void sort_adjacent();
};


class MW_MATCH {
private:
    vector<MW_EDGE*> max_match;
    vector<MW_EDGE*> non_max;
    char* is_best;
public:
    MW_MATCH();
    ~MW_MATCH();
    void solve(int const& n_vertices,vector<MW_EDGE*>& edges);
    vector<MW_EDGE*> const& get_max_match();
    vector<MW_EDGE*> const& get_non_max();
    bool is_in_max(int const& id);
};

#endif
