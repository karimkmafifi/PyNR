// This file was automatically created by 'derive_geometries.cpp'
// It supports distance data from the CSD

#ifndef __ELEMENT_DATA_new_GN
#define __ELEMENT_DATA_new_GN

#include<iostream>
#include<unordered_map>
#include<string>

using namespace std;

namespace EDnew {
template<typename T,int N1,int N2>
class ED2DTuple {
public:
    T t[N1][N2];
    ED2DTuple() {for (int i=0; i<N1; ++i) for (int j=0; j<N2; ++j) t[i][j] = T(0);}
    ED2DTuple(T const& v) {for (int i=0; i<N1; ++i) for (int j=0; j<N2; ++j) t[i][j] = v;}
    void set_element(int const& i1,int const& i2,T const& v) {t[i1][i2] = v;}
    T const& get_element(int const& i1,int const& i2) const {return t[i1][i2];}
};

extern float bond_data[7020];
extern unordered_map<string,ED2DTuple<int,5,3> > position_map; // 'ar','1','2','3','am' -> start_pos,start_bin,end_bin
extern float const ED_dummy;
extern float const ED_max_dist;
extern float const ED_min_dist;
extern float const ED_bin_size;
extern float angle_data[441];
extern unordered_map<string,ED2DTuple<int,3,3> > angle_position_map; // sp,sp2,sp3 -> start_pos,start_bin,end_bin
extern float const ED_max_angle;
extern float const ED_min_angle;
extern float const ED_angle_bin_size;
extern bool is_initialized;
void initialize();

class EDAccess {
private:
    static unordered_map<string,ED2DTuple<int,5,3> >::const_iterator last_questioned;
    static unordered_map<string,ED2DTuple<int,3,3> >::const_iterator angle_last_questioned;
    static string last_e1;
    static string last_e2;
    static string angle_last_e;
public:
    static bool connection_is_possible(string const& ele1,string const& ele2);
    static bool angle_exists(string const& ele);
    static float const& get_probability(int const& bond_order,float const& dist); // bond_order muss 0,1,2,3,4 sein <- Keine Pruefung!
    static float const& get_angle_probability(int const& hyb,float const& angle);
    static float const& get_probability(string const& ele1,string const& ele2,int const& bond_order,float const& dist); // bond_order muss 0,1,2,3,4 sein <- Keine Pruefung!
    static float const& get_angle_probability(string const& ele,int const& hyb,float const& angle);
    static float get_max_probability(float const& dist);
    static float get_max_angle_probability(float const& angle);
    static float get_max_probability(string const& ele1,string const& ele2,float const& dist);
    static float get_max_angle_probability(string const& ele,float const& angle);
};
}
#endif
