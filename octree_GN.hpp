
//============================================================================
// octree_GN.h -*- C++ -*-; Octree for objects with 'coords' attribute
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
// This library provides a generic octree for all objects with an
// attribute 'glm::vec3 coords'
//----------------------------------------------------------------------------


#ifndef __OCTREE_GN
#define __OCTREE_GN

#include<stdlib.h>
#include<iostream>
#include<vector>
#include"stl_ptr_GN.hpp"


using namespace std;


//! Der Typ T MUSS ein Attribut 'glm::vec3 coords' haben !!!
template<class T>
class OCTREE {
    public:
        vector<T*> objects;
        OCTREE<T>* rel_tree;
        typename vector<T*>::iterator begin_iterator;
        typename vector<T*>::iterator end_iterator;
        vector<OCTREE<T>*> leafs;
        float x_min; float x_max; float x_half;
        float y_min; float y_max; float y_half;
        float z_min; float z_max; float z_half;
        float max_edge;
        float range;
        OCTREE<T>* o1; // x_min -> x_half; y_min -> y_half; z_min -> z_half
        OCTREE<T>* o2; // x_half -> x_max; y_min -> y_half; z_min -> z_half
        OCTREE<T>* o3; // x_min -> x_half; y_half -> y_max; z_min -> z_half
        OCTREE<T>* o4; // x_half -> x_max; y_half -> y_max; z_min -> z_half
        OCTREE<T>* o5; // x_min -> x_half; y_min -> y_half; z_half -> z_max
        OCTREE<T>* o6; // x_half -> x_max; y_min -> y_half; z_half -> z_max
        OCTREE<T>* o7; // x_min -> x_half; y_half -> y_max; z_half -> z_max
        OCTREE<T>* o8; // x_half -> x_max; y_half -> y_max; z_half -> z_max
        OCTREE<T>* parent;

        OCTREE(vector<stl_ptr<T> > const& objs,float const& r,float const& extra_range=-1.,OCTREE<T>* pp=0) :
               x_min(999999.),x_max(-999999.),y_min(999999.),y_max(-999999.),
               z_min(999999.),z_max(-999999.),range(r),o1(0),o2(0),
               o3(0),o4(0),o5(0),o6(0),o7(0),o8(0),parent(pp) {
            for (typename vector<stl_ptr<T> >::const_iterator at=objs.begin(); at!=objs.end(); ++at) {
                if ((*at)->coords[0] < x_min) x_min = (*at)->coords[0];
                if ((*at)->coords[0] > x_max) x_max = (*at)->coords[0];
                if ((*at)->coords[1] < y_min) y_min = (*at)->coords[1];
                if ((*at)->coords[1] > y_max) y_max = (*at)->coords[1];
                if ((*at)->coords[2] < z_min) z_min = (*at)->coords[2];
                if ((*at)->coords[2] > z_max) z_max = (*at)->coords[2];
                objects.push_back(at->get_pointer());
            }
            x_min -= range; x_max += range;
            y_min -= range; y_max += range;
            z_min -= range; z_max += range;
            if (extra_range > 0.) {
                x_min -= extra_range; x_max += extra_range;
                y_min -= extra_range; y_max += extra_range;
                z_min -= extra_range; z_max += extra_range;
            }
            float xedge = x_max - x_min; max_edge = xedge;
            float yedge = y_max - y_min; if (yedge > max_edge) max_edge = yedge;
            float zedge = z_max - z_min; if (zedge > max_edge) max_edge = zedge;
            x_half = x_min + (xedge / 2.);
            y_half = y_min + (yedge / 2.);
            z_half = z_min + (zedge / 2.);
        }


        OCTREE(vector<T*> const& objs,float const& r,float const& extra_range=-1.,OCTREE<T>* pp=0) :
               x_min(999999.),x_max(-999999.),y_min(999999.),y_max(-999999.),
               z_min(999999.),z_max(-999999.),range(r),o1(0),o2(0),
               o3(0),o4(0),o5(0),o6(0),o7(0),o8(0),parent(pp) {
            for (typename vector<T*>::const_iterator at=objs.begin(); at!=objs.end(); ++at) {
                if ((*at)->coords[0] < x_min) x_min = (*at)->coords[0];
                if ((*at)->coords[0] > x_max) x_max = (*at)->coords[0];
                if ((*at)->coords[1] < y_min) y_min = (*at)->coords[1];
                if ((*at)->coords[1] > y_max) y_max = (*at)->coords[1];
                if ((*at)->coords[2] < z_min) z_min = (*at)->coords[2];
                if ((*at)->coords[2] > z_max) z_max = (*at)->coords[2];
                objects.push_back(*at);
            }
            x_min -= range; x_max += range;
            y_min -= range; y_max += range;
            z_min -= range; z_max += range;
            if (extra_range > 0.) {
                x_min -= extra_range; x_max += extra_range;
                y_min -= extra_range; y_max += extra_range;
                z_min -= extra_range; z_max += extra_range;
            }
            float xedge = x_max - x_min; max_edge = xedge;
            float yedge = y_max - y_min; if (yedge > max_edge) max_edge = yedge;
            float zedge = z_max - z_min; if (zedge > max_edge) max_edge = zedge;
            x_half = x_min + (xedge / 2.);
            y_half = y_min + (yedge / 2.);
            z_half = z_min + (zedge / 2.);
        }


        OCTREE(vector<T*> const& objs,float const& r,float const& xmi,float const& xma,
               float const& ymi,float const& yma,float const& zmi,float const& zma,OCTREE<T>* pp=0) :
               x_min(xmi),x_max(xma),y_min(ymi),y_max(yma),z_min(zmi),z_max(zma),range(r),
               o1(0),o2(0),o3(0),o4(0),o5(0),o6(0),o7(0),o8(0),parent(pp) {
            float xc_min = x_min - range; float xc_max = x_max + range;
            float yc_min = y_min - range; float yc_max = y_max + range;
            float zc_min = z_min - range; float zc_max = z_max + range;
            for (typename vector<T*>::const_iterator at=objs.begin(); at!=objs.end(); ++at) {
                if ((*at)->coords[0] < xc_min) continue;
                else if ((*at)->coords[0] > xc_max) continue;
                else if ((*at)->coords[1] < yc_min) continue;
                else if ((*at)->coords[1] > yc_max) continue;
                else if ((*at)->coords[2] < zc_min) continue;
                else if ((*at)->coords[2] > zc_max) continue;
                else objects.push_back(*at);
            }
            float xedge = x_max - x_min; max_edge = xedge;
            float yedge = y_max - y_min; if (yedge > max_edge) max_edge = yedge;
            float zedge = z_max - z_min; if (zedge > max_edge) max_edge = zedge;
            x_half = x_min + (xedge / 2.);
            y_half = y_min + (yedge / 2.);
            z_half = z_min + (zedge / 2.);
        }


        ~OCTREE() {
            if (o1) delete o1;
            if (o2) delete o2;
            if (o3) delete o3;
            if (o4) delete o4;
            if (o5) delete o5;
            if (o6) delete o6;
            if (o7) delete o7;
            if (o8) delete o8;
        }


        bool is_inside(glm::vec3 const& coords) {
            //! Diese Funktion NUR fuer den Wurzelknoten aufrufen!
            if (coords[0] < x_min) return false;
            else if (coords[0] > x_max) return false;
            else if (coords[1] < y_min) return false;
            else if (coords[1] > y_max) return false;
            else if (coords[2] < z_min) return false;
            else if (coords[2] > z_max) return false;
            else return true;
        }


        vector<T*>& get_rel_atoms(float const& max_spacing,glm::vec3 const& coords,
                                  unsigned int const& min_objs = 1,unsigned int const& max_objs = 9999999) {
            if ((max_spacing < max_edge && objects.size() > min_objs) || objects.size() > max_objs) {
                if (coords[0] < x_half) { // o1, o3, o5, o7
                    if (coords[1] < y_half) { // o1, o5
                        if (coords[2] < z_half) { // o1
                            if (o1 == 0) {
                                o1 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_min,z_half,this);
                            }
                            return o1->get_rel_atoms(max_spacing,coords);
                        } else { // o5
                            if (o5 == 0) {
                                o5 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_half,z_max,this);
                            }
                            return o5->get_rel_atoms(max_spacing,coords);
                        }
                    } else { // o3, o7
                        if (coords[2] < z_half) { // o3
                            if (o3 == 0) {
                                o3 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_min,z_half,this);
                            }
                            return o3->get_rel_atoms(max_spacing,coords);
                        } else { // o7
                            if (o7 == 0) {
                                o7 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_half,z_max,this);
                            }
                            return o7->get_rel_atoms(max_spacing,coords);
                        }
                    }
                } else { // o2, o4, o6, o8
                    if (coords[1] < y_half) { // o2, o6
                        if (coords[2] < z_half) { // o2
                            if (o2 == 0) {
                                o2 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_min,z_half,this);
                            }
                            return o2->get_rel_atoms(max_spacing,coords);
                        } else { // o6
                            if (o6 == 0) {
                                o6 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_half,z_max,this);
                            }
                            return o6->get_rel_atoms(max_spacing,coords);
                        }
                    } else { // o4, o8
                        if (coords[2] < z_half) { // o4
                            if (o4 == 0) {
                                o4 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_min,z_half,this);
                            }
                            return o4->get_rel_atoms(max_spacing,coords);
                        } else { // o8
                            if (o8 == 0) {
                                o8 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_half,z_max,this);
                            }
                            return o8->get_rel_atoms(max_spacing,coords);
                        }
                    }
                }
            } else return objects;
        }


        OCTREE<T>* p_get_rel_atoms(float const& max_spacing,glm::vec3 const& coords,
                                   unsigned int const& min_objs = 1,unsigned int const& max_objs = 9999999) {
            if ((max_spacing < max_edge && objects.size() > min_objs) || objects.size() > max_objs) {
                if (coords[0] < x_half) { // o1, o3, o5, o7
                    if (coords[1] < y_half) { // o1, o5
                        if (coords[2] < z_half) { // o1
                            if (o1 == 0) {
                                o1 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_min,z_half,this);
                            }
                            return o1->p_get_rel_atoms(max_spacing,coords);
                        } else { // o5
                            if (o5 == 0) {
                                o5 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_half,z_max,this);
                            }
                            return o5->p_get_rel_atoms(max_spacing,coords);
                        }
                    } else { // o3, o7
                        if (coords[2] < z_half) { // o3
                            if (o3 == 0) {
                                o3 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_min,z_half,this);
                            }
                            return o3->p_get_rel_atoms(max_spacing,coords);
                        } else { // o7
                            if (o7 == 0) {
                                o7 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_half,z_max,this);
                            }
                            return o7->p_get_rel_atoms(max_spacing,coords);
                        }
                    }
                } else { // o2, o4, o6, o8
                    if (coords[1] < y_half) { // o2, o6
                        if (coords[2] < z_half) { // o2
                            if (o2 == 0) {
                                o2 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_min,z_half,this);
                            }
                            return o2->p_get_rel_atoms(max_spacing,coords);
                        } else { // o6
                            if (o6 == 0) {
                                o6 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_half,z_max,this);
                            }
                            return o6->p_get_rel_atoms(max_spacing,coords);
                        }
                    } else { // o4, o8
                        if (coords[2] < z_half) { // o4
                            if (o4 == 0) {
                                o4 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_min,z_half,this);
                            }
                            return o4->p_get_rel_atoms(max_spacing,coords);
                        } else { // o8
                            if (o8 == 0) {
                                o8 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_half,z_max,this);
                            }
                            return o8->p_get_rel_atoms(max_spacing,coords);
                        }
                    }
                }
            } else return this;
        }


        typename vector<T*>::iterator const& begin(float const& max_spacing,glm::vec3 const& coords,
                                                   unsigned int const& min_objs = 1,
                                                   unsigned int const& max_objs = 9999999) {
            if (is_inside(coords)) {
                rel_tree = p_get_rel_atoms(max_spacing,coords,min_objs,max_objs);
                begin_iterator = rel_tree->objects.begin();
                end_iterator = rel_tree->objects.end();
                return begin_iterator;
            } else {
                end_iterator = objects.end();
                return end_iterator;
            }
        }


        typename vector<T*>::iterator const& end() {
            return end_iterator;
        }


        void full_subdivide() {
            if (o1 == 0) o1 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_min,z_half,this);
            if (o2 == 0) o2 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_min,z_half,this);
            if (o3 == 0) o3 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_min,z_half,this);
            if (o4 == 0) o4 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_min,z_half,this);
            if (o5 == 0) o5 = new OCTREE<T>(objects,range,x_min,x_half,y_min,y_half,z_half,z_max,this);
            if (o6 == 0) o6 = new OCTREE<T>(objects,range,x_half,x_max,y_min,y_half,z_half,z_max,this);
            if (o7 == 0) o7 = new OCTREE<T>(objects,range,x_min,x_half,y_half,y_max,z_half,z_max,this);
            if (o8 == 0) o8 = new OCTREE<T>(objects,range,x_half,x_max,y_half,y_max,z_half,z_max,this);
        }


        void generate_full_tree(float const& max_spacing,unsigned int const& min_objs = 1,
                                unsigned int const& max_objs = 9999999) {
            vector<OCTREE<T>*> curr_leafs;
            vector<OCTREE<T>*> new_leafs;
            curr_leafs.push_back(this);
            while (curr_leafs.size() > 0) {
                for (typename vector<OCTREE<T>*>::iterator it=curr_leafs.begin(); it!=curr_leafs.end(); ++it) {
                    if (max_spacing >= (*it)->max_edge || objects.size() < min_objs) {
                        if ((*it)->objects.size() > 0) leafs.push_back(*it);
                    } else {
                        (*it)->full_subdivide();
                        if ((max_spacing < (*it)->o1->max_edge && (*it)->o1->objects.size() > min_objs) ||
                             (*it)->o1->objects.size() > max_objs) new_leafs.push_back((*it)->o1);
                        else if ((*it)->o1->objects.size() > 0) leafs.push_back((*it)->o1);
                        if ((max_spacing < (*it)->o2->max_edge && (*it)->o2->objects.size() > min_objs) ||
                             (*it)->o2->objects.size() > max_objs) new_leafs.push_back((*it)->o2);
                        else if ((*it)->o2->objects.size() > 0) leafs.push_back((*it)->o2);
                        if ((max_spacing < (*it)->o3->max_edge && (*it)->o3->objects.size() > min_objs) ||
                             (*it)->o3->objects.size() > max_objs) new_leafs.push_back((*it)->o3);
                        else if ((*it)->o3->objects.size() > 0) leafs.push_back((*it)->o3);
                        if ((max_spacing < (*it)->o4->max_edge && (*it)->o4->objects.size() > min_objs) ||
                             (*it)->o4->objects.size() > max_objs) new_leafs.push_back((*it)->o4);
                        else if ((*it)->o4->objects.size() > 0) leafs.push_back((*it)->o4);
                        if ((max_spacing < (*it)->o5->max_edge && (*it)->o5->objects.size() > min_objs) ||
                             (*it)->o5->objects.size() > max_objs) new_leafs.push_back((*it)->o5);
                        else if ((*it)->o5->objects.size() > 0) leafs.push_back((*it)->o5);
                        if ((max_spacing < (*it)->o6->max_edge && (*it)->o6->objects.size() > min_objs) ||
                             (*it)->o6->objects.size() > max_objs) new_leafs.push_back((*it)->o6);
                        else if ((*it)->o6->objects.size() > 0) leafs.push_back((*it)->o6);
                        if ((max_spacing < (*it)->o7->max_edge && (*it)->o7->objects.size() > min_objs) ||
                             (*it)->o7->objects.size() > max_objs) new_leafs.push_back((*it)->o7);
                        else if ((*it)->o7->objects.size() > 0) leafs.push_back((*it)->o7);
                        if ((max_spacing < (*it)->o8->max_edge && (*it)->o8->objects.size() > min_objs) ||
                             (*it)->o8->objects.size() > max_objs) new_leafs.push_back((*it)->o8);
                        else if ((*it)->o8->objects.size() > 0) leafs.push_back((*it)->o8);
                    }
                }
                curr_leafs.clear();
                for (typename vector<OCTREE<T>*>::iterator it=new_leafs.begin(); it!=new_leafs.end(); ++it) {
                    curr_leafs.push_back(*it);
                }
                new_leafs.clear();
            }
        }
};


#endif
