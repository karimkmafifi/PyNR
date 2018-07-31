#include "scoring_terms.h"

ScoringTerms::ScoringTerms()
{

}

float ScoringTerms::gaussian(float x, float width)
{
    return std::exp(-Common::sqr(x / width));
}

float ScoringTerms::slope_step(float x_bad, float x_good, float x)
{
    if (x_bad < x_good) {
        if (x <= x_bad) return 0;
        if (x >= x_good) return 1;
    }
    else {
        if (x >= x_bad) return 0;
        if (x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}


float ScoringTerms::gaussVINA_eval(Atom* a, Atom* b, float r, float offset, float width)
{
    float op_dst = ScoringTerms::optimal_distance(a->xs, b->xs);
    float x = r - (op_dst + offset);
    return gaussian(x, width);
}


float ScoringTerms::repulsionVINA_eval(Atom* a, Atom* b, float r, float offset)
{
    float d = r - (ScoringTerms::optimal_distance(a->xs, b->xs) + offset);
    if (d > 0)
        return 0;
    return d*d;
}


float ScoringTerms::hydrophobicVINA_eval(Atom* a, Atom* b, float r, float good, float bad)
{
    if (xs_is_hydrophobic(a->xs) && xs_is_hydrophobic(b->xs))
    {
        return slope_step(bad, good, r - ScoringTerms::optimal_distance(a->xs, b->xs));
    }
    else
    {
        return 0;
    }
}

float ScoringTerms::non_dir_h_bondVINA_eval(Atom* a, Atom* b, float r, float good, float bad)
{
    if (xs_h_bond_possible(a->xs, b->xs))
    {
        return slope_step(bad, good, r - ScoringTerms::optimal_distance(a->xs, b->xs));
    }

    return 0;
}


float ScoringTerms::halogen_cl_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float bad)
{
    float good;
    if (theta > 165)
    {
        good = good_1;
    }
    else
    {
        good = good_2;
    }
    float ClA = 0.14;
    float ClB = -0.016;
    float Clv = 2.46;
    float ClAngle = (ClA*cos(((180 - theta)*Clv)* pi / 180.0) + ClB) / (ClA + ClB); //max is 1
    ClAngle = ClAngle < 0 ? 0 : ClAngle;
    ClAngle = theta < 140 ? 0 : ClAngle;
    if (xs_hal_cl_bond_possible(a->xs, b->xs))
        return slope_step(bad, good, r - ScoringTerms::optimal_distance(a->xs, b->xs)) * ClAngle;
    return 0;
}

float ScoringTerms::halogen_br_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float good_3, float bad) {
    float good;
    if (theta > 165) {
        good = good_1;
    }
    else if (theta > 150) {
        good = good_2;
    }
    else {
        good = good_3;
    }
    float BrA = 0.23;
    float BrB = 0.15;
    float Brv = 2.42;
    float BrAngle = (BrA*cos(((180 - theta)*Brv)* pi / 180.0) + BrB) / (BrA + BrB); //max is 1
    BrAngle = BrAngle < 0 ? 0 : BrAngle;
    BrAngle = theta < 120 ? 0 : BrAngle;
    if (xs_hal_br_bond_possible(a->xs, b->xs))
        return slope_step(bad, good, r - ScoringTerms::optimal_distance(a->xs, b->xs)) * BrAngle;
    return 0;
}

float ScoringTerms::halogen_i_bondVINAXB_eval(Atom* a, Atom* b, float r, float theta, float good_1, float good_2, float good_3, float bad) {

    float good;
    if (theta > 165) {
        good = good_1;
    }
    else if (theta > 150) {
        good = good_2;
    }
    else {
        good = good_3;
    }
    float IA = 0.46;
    float IB = 0.29;
    float Iv = 2.23;
    float IAngle = (IA*cos(((180 - theta)*Iv)*pi / 180.0) + IB) / (IA + IB); //max is 1
    IAngle = IAngle < 0 ? 0 : IAngle;
    IAngle = theta < 120 ? 0 : IAngle;
    if (xs_hal_i_bond_possible(a->xs, b->xs)) {
        return slope_step(bad, good, r - ScoringTerms::optimal_distance(a->xs, b->xs)) * IAngle;
    }
    return 0;
}

float ScoringTerms::hydrogen_bondAD_eval(Atom* a, Atom* b, float r, float smoothing, float cap)
{
    if (ad_h_bond_possible(a->ad, b->ad))
    {
        float d0 = 0;
        float depth = 1;

        if (ad_is_acceptor(a->ad))
        {
            d0 = ad_type_property(a->ad).hbond_radius;
            depth = ad_type_property(a->ad).hbond_depth;
        }
        else
        {
            d0 = ad_type_property(b->ad).hbond_radius;
            depth = ad_type_property(b->ad).hbond_depth;
        }

        float c_i = 0;
        float c_j = 0;

        find_vdw_hb_coefficients<10, 12>(d0, depth, c_i, c_j);

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

        float r_i = pow(r, 10);
        float r_j = pow(r, 12);

        if (r_i > epsilon_float && r_j > epsilon_float)
            return (std::min)(cap, c_i / r_i + c_j / r_j);
        else
            return cap;
    }
    return 0;
}

inline float ScoringTerms::calc_ddd_Mehler_Solmajer(float distance)
{
    //____________________________________________________________________________
    //Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
    //___________________________________________________________________________
    float epsilon = 1.0L;
    float lambda = 0.003627L;
    float epsilon0 = 78.4L;
    float A = -8.5525L;
    float B;
    B = epsilon0 - A; //8.69525e1
    float rk = 7.7839L;
    float lambda_B;
    lambda_B = -lambda * B; //-0.3153767175

    epsilon = A + B / (1.0L + rk*exp(lambda_B * distance));

    if (epsilon < APPROX_ZERO) {
        epsilon = 1.0L;
    }
    return epsilon;
}

float ScoringTerms::electrostaticAD_eval(Atom* a, Atom* b, float r)
{
    float q1q2 = a->charge * b->charge;

    if (r == 0.0)
    {
        r = 0.0001;
    }

    float epsilon = calc_ddd_Mehler_Solmajer(r);

    float r_epsilon = r * epsilon;
    r_epsilon = 1.0 / r_epsilon;

    float elec = q1q2 * r_epsilon;

    return elec;
}

inline float ScoringTerms::solvation_parameter(Atom_type* a)
{
    if (a->ad < AD_TYPE_SIZE)
    {
        return ad_type_property(a->ad).solvation;
    }
    else if (a->xs == XS_TYPE_Met_D)
    {
        return metal_solvation_parameter;
    }

    throw Error_report("Internal error. Can't find Solvation Parameter");

    return 0; // placating the compiler
}

inline float ScoringTerms::volume(Atom_type* a)
{
    if (a->ad < AD_TYPE_SIZE)
    {
        return ad_type_property(a->ad).volume;
    }
    else if (a->xs < XS_TYPE_SIZE)
    {
        return 4 * pi / 3 * pow(xs_radius(a->xs), 3);
    }

    throw Error_report("Internal error. Can't find Volume");

    return 0; // placating the compiler
}

float ScoringTerms::solvationAD_eval(Atom* a, Atom* b, float r, float desolvation_sigma, float solvation_q, bool charge_dependent)
{
    float q1 = a->charge;
    float q2 = b->charge;

    if (!(q1 < 0.1 * max_float)) { throw Error_report("Solvation AD Error"); };
    if (!(q2 < 0.1 * max_float)) { throw Error_report("Solvation AD Error"); };

    float solv1 = solvation_parameter(a);
    float solv2 = solvation_parameter(b);

    float volume1 = volume(a);
    float volume2 = volume(b);

    float my_solv = charge_dependent ? solvation_q : 0;

    float tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 +
        (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-Common::sqr(r / (2 * desolvation_sigma)));

    if (!(tmp < 0.1 * max_float)) { throw Error_report("Solvation AD Error"); };

    return tmp;
}

inline float ScoringTerms::xs_hbond_distance_piece_wise(float optimal_distance, float distance)
{
    if (optimal_distance <= (optimal_distance - 0.7))
    {
        return 1.0;
    }
    else if ((optimal_distance - 0.7 < distance) && (distance <= optimal_distance))
    {
        return (1 / 0.7) * (optimal_distance - distance);
    }
    else if (distance > optimal_distance)
    {
        return 0.0;
    }
    else
    {
        return 1.0;
    }

    throw Error_report("X-Score Hydrogen Bonding Distance Piece Wise Error");

    return 1.0;
}

inline float ScoringTerms::xs_hbond_theta_piece_wise(float theta)
{
    if (theta >= 120)
    {
        return 1.0;
    }
    else if ((120 > theta) && (theta >= 60))
    {
        return ((1. / 60) * (theta - 60));
    }
    else if (theta < 60)
    {
        return 0.0;
    }

    throw Error_report("X-Score Hydrogen Bonding Theta Piece Wise Error");

    return 1.0;
}

float ScoringTerms::eval_dir_h_bond_XS(Atom* a, Atom* b, float r, float theta_1, float theta_2)
{
    float distance_term = ScoringTerms::xs_hbond_distance_piece_wise(ScoringTerms::optimal_distance(a->xs, b->xs), r);
    float theta_1_term = ScoringTerms::xs_hbond_theta_piece_wise(theta_1);
    float theta_2_term = ScoringTerms::xs_hbond_theta_piece_wise(theta_2);
    return  distance_term * theta_1_term * theta_2_term;
}

inline float ScoringTerms::xs_hydrophobic_contact_piece_wise(float optimal_distance, float distance)
{
    if (distance <= (optimal_distance + 0.5))
    {
        return 1.0;
    }
    else if ((optimal_distance + 0.5 < distance) && (distance <= optimal_distance + 2.0))
    {
        return (1 / 1.5) * ((optimal_distance + 2.0) - distance);
    }
    else if (distance > optimal_distance + 2.0)
    {
        return 0.0;
    }

    throw Error_report("X-Score Hydrophobic Contact Piece Wise Error");

    return 1.0;
}

float ScoringTerms::evalHydrophobicContactXS(Atom* a, Atom* b, float r)
{
    if (xs_is_hydrophobic(a->xs) && xs_is_hydrophobic(b->xs))
    {
        return ScoringTerms::xs_hydrophobic_contact_piece_wise(ScoringTerms::optimal_distance(a->xs, b->xs), r);
    }
    else
    {
        return 0.0;
    }
}

unsigned int ScoringTerms::atom_rotors(Atom* atm)
{
    unsigned int acc = 0;

    for (int j = 0; j < atm->atom_bonds.size(); ++j)
    {
        Atom* a = atm->get_connected_atom(j);

        if ((atm->branch_id != a->branch_id) && !a->is_hydrogen() && a->n_heavy_bonded > 1) // not counting CH_3, etc
        {
            ++acc;
        }
    }

    return acc;
}

void ScoringTerms::calculate_conf_independent(Molecule* mol, float& num_tors, float& num_rotors)
{
    for (int i = 0; i < mol->atoms.size(); ++i)
    {
        Atom* a = mol->atoms[i];
        if (a->el != EL_TYPE_H)
        {
            unsigned int ar = atom_rotors(a);

            num_tors += 0.5 * ar;

            if (ar > 2)
            {
                num_rotors += 0.5;
            }
            else
            {
                num_rotors += 0.5 * ar;
            }
        }
    }
}

ScoringTerms::~ScoringTerms()
{

}
