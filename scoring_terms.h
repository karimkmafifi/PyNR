#ifndef SCORING_TERMS_H
#define SCORING_TERMS_H
#include "parser.h"

class ScoringTerms
{
private:

    inline float optimal_distance(int xs_t1, int xs_t2)
    {
        return xs_radius(xs_t1) + xs_radius(xs_t2);
    }

   float gaussian(float x, float width);
   float slope_step(float x_bad, float x_good, float x);
   unsigned int atom_rotors(Atom* atm);

   template<unsigned n, unsigned m>
   inline void find_vdw_hb_coefficients(float position, float depth, float& c_n, float& c_m)
   {
       Q_ASSERT(n != m);

       c_n = pow(position, n) * depth * m / (float(n) - float(m));
       c_m = pow(position, m) * depth * n / (float(m) - float(n));
   }

   inline float calc_ddd_Mehler_Solmajer(float distance);
   inline float solvation_parameter(Atom_type* a);
   inline float volume(Atom_type* a);
   inline float xs_hbond_distance_piece_wise(float optimal_distance, float distance);
   inline float xs_hbond_theta_piece_wise(float theta);
   inline float xs_hydrophobic_contact_piece_wise(float optimal_distance, float distance);

public:
    ScoringTerms();
    float gaussVINA_eval(Atom* a, Atom* b, float r, float offset, float width);
    float repulsionVINA_eval(Atom* a, Atom* b, float r, float offset);
    float hydrophobicVINA_eval(Atom* a, Atom* b, float r, float good, float bad);
    float non_dir_h_bondVINA_eval(Atom* a, Atom* b, float r, float good, float bad);
    float halogen_cl_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float bad);
    float halogen_br_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float good_3, float bad);
    float halogen_i_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float good_3, float bad);
    void calculate_conf_independent(Molecule* mol, float& num_tors, float& num_rotors);

    template<unsigned i, unsigned j>
    float vdw_eval(Atom* a, Atom* b, float r, float smoothing, float cap, bool is_autodock)
    {
        bool calculate_flag = false;

        if (is_autodock == true)
        {
            if (!(ad_h_bond_possible(a->ad, b->ad)))
            {
                calculate_flag = true;
            }
        }
        else
        {
            calculate_flag = true;
        }

        if (calculate_flag == true)
        {
            float d0 = 0;
            float depth = 1;

            if (is_autodock == false)
            {
                d0 = ScoringTerms::optimal_distance(a->xs, b->xs);
            }
            else if (is_autodock == true)
            {
                d0 = (ad_type_property(a->ad).radius + ad_type_property(b->ad).radius);
                depth = sqrt(ad_type_property(a->ad).depth * ad_type_property(b->ad).depth);
            }

            float c_i = 0;
            float c_j = 0;
            find_vdw_hb_coefficients<i, j>(d0, depth, c_i, c_j);

            if (r > d0 + smoothing)
            {
                r -= smoothing;
            }
            else if (r < d0 - smoothing)
            {
                r += smoothing;
            }
            else
            {
                r = d0;
            }


            float r_i = pow(r, i);
            float r_j = pow(r, j);


            if (r_i > epsilon_float && r_j > epsilon_float)
            {
                return (std::min)(cap, c_i / r_i + c_j / r_j);
            }
            else
            {
                return cap;
            }
        }
        else
        {
            return 0.;
        }
    }

    float hydrogen_bondAD_eval(Atom* a, Atom* b, float r, float smoothing, float cap);
    float electrostaticAD_eval(Atom* a, Atom* b, float r);
    float solvationAD_eval(Atom* a, Atom* b, float r, float desolvation_sigma, float solvation_q, bool charge_dependent);

    float eval_dir_h_bond_XS(Atom* a, Atom* b, float r, float theta_1, float theta_2);
    float evalHydrophobicContactXS(Atom* a, Atom* b, float r);

    ~ScoringTerms();
};

#endif
