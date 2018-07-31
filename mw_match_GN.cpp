
//============================================================================
// mw_match_GN.cpp -*- C++ -*-; maximum weighted match for Graphs where no
//                              nodes with degree > 3 exist
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


#include"mw_match_GN.h"
#include<iostream>
#include<stack>
#include<string.h>
#include<list>
#include<algorithm>
#include<unordered_set>


bool edge_sort(MW_EDGE* const l,MW_EDGE* const r) {return (l->weight < r->weight);}


MW_EDGE::MW_EDGE(int const& i,int const& j,float const& w):node_i(i),node_j(j),weight(w),
                                                           stacked(false),labeled(false) {}


MW_EDGE::~MW_EDGE() {}


void MW_EDGE::set_max_label(unsigned int const& ml) {
    if (ml > label.size()) {
        int diff = ml-label.size();
        for (int i=0; i<diff; ++i) label.push_back(false);
    }
}


void MW_EDGE::sort_adjacent() {
    sort(adjacent.begin(),adjacent.end(),edge_sort);
}


MW_MATCH::MW_MATCH():is_best(0) {}


MW_MATCH::~MW_MATCH() {
    if (is_best) delete[] is_best;
}


void MW_MATCH::solve(int const& n_vertices,vector<MW_EDGE*>& edges) {
    if (is_best) delete[] is_best;
    is_best = new char[n_vertices];
    memset(is_best,0,n_vertices);
    
    for (vector<MW_EDGE*>::iterator it=edges.begin(); it!=edges.end(); ++it) {
        for (vector<MW_EDGE*>::iterator jt=it+1; jt!=edges.end(); ++jt) {
            if ((*it)->node_i == (*jt)->node_i || (*it)->node_i == (*jt)->node_j ||
                (*it)->node_j == (*jt)->node_j || (*it)->node_j == (*jt)->node_i) {
                (*it)->adjacent.push_back(*jt);
                (*jt)->adjacent.push_back(*it);
            }
        }
    }

    for (vector<MW_EDGE*>::iterator it=edges.begin(); it!=edges.end(); ++it) {
        (*it)->set_max_label(2);
        (*it)->sort_adjacent();
    }

    for (vector<MW_EDGE*>::iterator vit=edges.begin(); vit!=edges.end(); ++vit) {
        if ((*vit)->labeled) continue;

        vector<MW_EDGE*> processed_edges;


//        unsigned int max_list_size = 0; // DEBUG
//        unsigned int max_end_set = 0; // DEBUG


        int max_label = 0;
        list<int> label_list;
        label_list.push_back(0);

        list<MW_EDGE*> pq;
        pq.push_back(*vit);
        (*vit)->stacked = true;

        vector<float> label_weights;
        label_weights.push_back(0.);

        while (!pq.empty()) 
		{
            MW_EDGE* curr = pq.back(); pq.pop_back();

//            cerr << "curr = " << curr->node_i << "," << curr->node_j << "   max_label = "
//                 << max_label << "   list_size = " << label_list.size() << endl;

//            if (label_list.size() > max_list_size) max_list_size = label_list.size(); // DEBUG

            for (unsigned int i=0; i<curr->adjacent.size(); ++i) curr->adjacent[i]->set_max_label(max_label+1);
            curr->set_max_label(max_label+1);

            int loc_max = max_label;
            for (list<int>::iterator l=label_list.begin(); l!=label_list.end(); ++l) {
                if (curr->label[*l]) continue;

                bool conti = false;

                for (unsigned int i=0; i<curr->adjacent.size(); ++i) {
                    if (!curr->adjacent[i]->label[*l]) continue;
                    conti = true;
                    if (curr->weight >= curr->adjacent[i]->weight) {
                        ++loc_max;
                        curr->set_max_label(loc_max+1);
                        label_weights.push_back(curr->weight);
                        curr->label[loc_max] = true;
                        curr->labeled = true;

//                        cerr << "  -> curr has better " << *l << " than " << curr->adjacent[i]->node_i << "," << curr->adjacent[i]->node_j
//                             << " : setting " << loc_max << endl;

                        for (vector<MW_EDGE*>::iterator pe=processed_edges.begin(); pe!=processed_edges.end(); ++pe) {
							bool pe_is_equal = false;
							for (int adj_index = 0; adj_index < curr->adjacent.size(); ++adj_index)
							{
								if (adj_index > 3)
								{
									break;
								}
								else
								{
									if (*pe == curr->adjacent[adj_index])
									{
										pe_is_equal = true;
										break;
									}
								}
							}
                            if(pe_is_equal) continue;
                            (*pe)->set_max_label(loc_max+1);
                            if ((*pe)->label[*l]) {
                                (*pe)->label[loc_max] = true;
                                label_weights[loc_max] += (*pe)->weight;
                            }
                        }
                        break;
                    }
                }
                for (unsigned int i=0; i<curr->adjacent.size(); ++i) {
                    if (!curr->adjacent[i]->stacked) {
                        pq.push_back(curr->adjacent[i]);
                        curr->adjacent[i]->stacked = true;
                    }
                }
                if (conti) continue;

//                cerr << " -> setting label " << *l << endl;

                curr->label[*l] = true;
                curr->labeled = true;
                label_weights[*l] += curr->weight;
            }
            if (!curr->labeled) {
                ++loc_max;

//                cerr << " -> setting new label " << loc_max << endl;
                
                curr->set_max_label(loc_max+1);
                label_weights.push_back(curr->weight);
                curr->label[loc_max] = true;
                curr->labeled = true;
            }

            for (int nl=max_label+1; nl<=loc_max; ++nl) label_list.push_back(nl);
            max_label = loc_max;
            processed_edges.push_back(curr);

            if (pq.size() < 22) {
                unordered_set<MW_EDGE*> end_set;
                list<MW_EDGE*>::iterator next_pq;
                int best_n_labeled = -1;
                for (list<MW_EDGE*>::iterator pt=pq.begin(); pt!=pq.end(); ++pt) {
                    int n_labeled = 0;
                    for (vector<MW_EDGE*>::iterator adt=(*pt)->adjacent.begin(); adt!=(*pt)->adjacent.end(); ++adt) {
                        if ((*adt)->labeled) {
                            end_set.insert(*adt);
                            ++n_labeled;
                        }
                    }
                    if (n_labeled > best_n_labeled) {
                        best_n_labeled = n_labeled;
                        next_pq = pt;
                    }
                }
                if (best_n_labeled > -1) {
                    // Sicherstellen, dass als naechstes die Kante prozessiert wird, an die die
                    // meisten gelabelten Kanten angrenzen => es werden dann weniger Terminale
                    // Kanten im naechsten Durchlauf!
                    pq.push_back(*next_pq);
                    pq.erase(next_pq);
                }


//                if (end_set.size() > max_end_set) max_end_set = end_set.size(); // DEBUG


                if (end_set.size() < 22) {
                    vector<MW_EDGE*> end_edges;
                    for (unordered_set<MW_EDGE*>::iterator et=end_set.begin(); et!=end_set.end(); ++et) {
                        (*et)->set_max_label(max_label+1);
                        end_edges.push_back(*et);
                    }

                    char* best_labels = new char[max_label+1];
                    memset(best_labels,0,max_label+1);
                    //! Fuer alle moeglichen Kombinationen aus Ende als DB/nicht-DB das beste Label bestimmen
                    int n_comb = 1;
                    n_comb <<= end_edges.size(); // Anzhal der Kombinationen
                    vector<int> comb_best(n_comb,-1);
                    vector<float> comb_weight(n_comb,-999999.);
                    
                    //! Label der entsprechenden Kombination zuordnen:
                    for (list<int>::iterator l=label_list.begin(); l!=label_list.end(); ++l) {
                        int index = 0;
                        int mul = 1;
                        for (unsigned int pt=0; pt<end_edges.size(); ++pt) {
                            if (end_edges[pt]->label[*l]) index += mul;
                            mul <<= 1;
                        }
                        if (label_weights[*l] > comb_weight[index]) {
                            comb_best[index] = *l;
                            comb_weight[index] = label_weights[*l];
                        }
                    }

                    //! In jeder Kombi das beste Label bestimmen:
                    for (int pt=0; pt<n_comb; ++pt) {
                        if (comb_best[pt] < 0) continue;
                        else best_labels[comb_best[pt]] = 1;
                    }

                    if ((max_label - int(label_list.size())) > 100) {
                        vector<int> old_best;
                        vector<float> old_weight;
                        for (list<int>::iterator l=label_list.begin(); l!=label_list.end(); ++l) {
                            if (best_labels[*l]) {
                                old_best.push_back(*l);
                                old_weight.push_back(label_weights[*l]);
                            }
                        }
                        for (vector<MW_EDGE*>::iterator pe=processed_edges.begin(); pe!=processed_edges.end(); ++pe) {
                            (*pe)->set_max_label(max_label+1);
                            for (int nl=0; nl<int(old_best.size()); ++nl) {
                                if ((*pe)->label[old_best[nl]]) (*pe)->label[nl] = true;
                                else (*pe)->label[nl] = false;
                            }
                            for (int ll=int(old_best.size()); ll<=max_label; ++ll) (*pe)->label[ll] = false;
                        }
                        label_list.clear();
                        for (int i=0; i<int(old_best.size()); ++i) label_list.push_back(i);
                        max_label = int(old_best.size()) - 1;
                        label_weights.clear();
                        for (int i=0; i<int(old_best.size()); ++i) label_weights.push_back(old_weight[i]);
                    } 
					else
					{
						int count_label_list = -1;
						for (list<int>::iterator l = label_list.begin(); l != label_list.end(); ++l) 
						{
							count_label_list++;
							if (best_labels[*l]) continue;

							bool break_from_loop = false;

							if (count_label_list == (label_list.size() - 1))
							{
								break_from_loop = true;
							}

							label_list.erase(l++);

							if (break_from_loop)
							{
								break;
							}
						}
					}

                    delete[] best_labels;
                }
            }
        }

        // Label mit hoechstem Gewicht bestimmen und entsprechende Kanten dem max_match zufuegen:
        int best_label = 0;
        float best_weight = 0.;
        for (int i=0; i<=max_label; ++i) {
            if (label_weights[i] > best_weight) {
                best_weight = label_weights[i];
                best_label = i;
            }
        }

//        cerr << "best_label = " << best_label << "  weight = " << best_weight
//             << "   max_terminal = " << max_end_set << "   max_labels = " << max_list_size << endl;

        for (vector<MW_EDGE*>::iterator pe=processed_edges.begin(); pe!=processed_edges.end(); ++pe) {
            if (int((*pe)->label.size()) < (best_label+1)) {
                non_max.push_back(*pe);
            } else if ((*pe)->label[best_label]) {
                max_match.push_back(*pe);
                is_best[(*pe)->node_i] = 1;
                is_best[(*pe)->node_j] = 1;
            } else non_max.push_back(*pe);
        }
    }
}


vector<MW_EDGE*> const& MW_MATCH::get_max_match() {
    return max_match;
}


vector<MW_EDGE*> const& MW_MATCH::get_non_max() {
    return non_max;
}


bool MW_MATCH::is_in_max(int const& id) {
    return is_best[id];
}
