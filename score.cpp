#include "score.h"

void Score::prepare_ad_hydrogen()
{
    Score::receptor->rexp.reserve(Score::receptor->atoms.size());
    Score::receptor->rvector.reserve(Score::receptor->atoms.size());
    Score::receptor->rvector2.reserve(Score::receptor->atoms.size());

    for (unsigned int i = 0; i < Score::receptor->atoms.size(); ++i)
    {
        Score::receptor->rexp.push_back(0);
        Score::receptor->rvector.push_back(glm::vec3(0.0, 0.0, 0.0));
        Score::receptor->rvector2.push_back(glm::vec3(0.0, 0.0, 0.0));
    }

    for (int ia = 0; ia < Score::receptor->atoms.size(); ++ia)
    {
        int from = ((ia - 20) > (0)) ? (ia - 20) : (0);
        int to = ((ia + 20) < (Score::receptor->atoms.size() - 1)) ? (ia + 20) : (Score::receptor->atoms.size() - 1);

        Atom* ia_atom = Score::receptor->atoms[ia];

        if (ad_type_property(ia_atom->ad).hbond == 2)
        {
            for (int ib = from; ib < to; ++ib)
            {
                if (ib != ia)
                {
                    Atom* ib_atom = Score::receptor->atoms[ib];
                    glm::vec3 d = ia_atom->coords - ib_atom->coords;
                    float rd2 = Common::sqr(d[0]) + Common::sqr(d[1]) + Common::sqr(d[2]);

                    if (rd2 < 1.90)
                    {
                        if (rd2 < APPROX_ZERO)
                        {
                            rd2 = APPROX_ZERO;
                        }

                        float inv_rd = 1. / sqrt(rd2);

                        if ((ib_atom->ad != AD_TYPE_OA) && (ib_atom->ad != AD_TYPE_SA))
                        {
                            Score::receptor->rexp[ia] = 2;
                        }

                        if ((ib_atom->ad == AD_TYPE_OA) || (ib_atom->ad == AD_TYPE_SA))
                        {
                            Score::receptor->rexp[ia] = 4;
                        }

                        Score::receptor->rvector[ia] = d * inv_rd;

                        break;
                    }
                }
            }
        }
        else if (ad_type_property(ia_atom->ad).hbond == 5)
        {
            int nbond = 0;
            int i1 = 0, i2 = 0;

            for (int ib = from; ib < to; ++ib)
            {
                if (ib != ia)
                {
                    Atom* ib_atom = Score::receptor->atoms[ib];
                    float rd2 = 0.;

                    glm::vec3 dc = ia_atom->coords - ib_atom->coords;
                    rd2 = Common::sqr(dc[0]) + Common::sqr(dc[1]) + Common::sqr(dc[2]);

                    if (((rd2 < 3.61) && ((ib_atom->ad != AD_TYPE_HD) && (ib_atom->ad != AD_TYPE_H))) ||
                        ((rd2 < 1.69) && ((ib_atom->ad == AD_TYPE_HD) || (ib_atom->ad == AD_TYPE_H)))) {

                        if (nbond == 2) {
                            continue;
                            //throw Error_report("Error: Found an H-bonding atom with three bonded atoms, atom serial number " + (ia + 1));
                        }
                        if (nbond == 1) {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0) {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }
            }

            if (nbond == 0) {
                //throw Error_report("Error: Oxygen atom found with no bonded atoms, atom serial number " + (ia + 1) + ad_type_property(ia_atom->ad).name);
                continue;
            }

            if (nbond == 1)
            {
                Atom* i1_atom = Score::receptor->atoms[i1];

                float rd2 = 0.;

                Score::receptor->rvector[ia] = ia_atom->coords - i1_atom->coords;
                rd2 = Common::sqr(Score::receptor->rvector[ia][0]) + Common::sqr(Score::receptor->rvector[ia][1]) + Common::sqr(Score::receptor->rvector[ia][2]);

                if (rd2 < APPROX_ZERO) {
                    rd2 = APPROX_ZERO;
                }

                float inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector[ia] *= inv_rd;

                for (int i2 = from; i2 < to; ++i2)
                {
                    if ((i2 != i1) && (i2 != ia)) {

                        Atom* i2_atom = Score::receptor->atoms.at(i2);

                        rd2 = 0.;

                        glm::vec3 dc = i1_atom->coords - i2_atom->coords;
                        rd2 = Common::sqr(dc[0]) + Common::sqr(dc[1]) + Common::sqr(dc[2]);

                        if (((rd2 < 2.89) && (i2_atom->ad != AD_TYPE_HD)) ||
                            ((rd2 < 1.69) && (i2_atom->ad == AD_TYPE_HD))) {

                            rd2 = 0.;

                            glm::vec3 d = i2_atom->coords - i1_atom->coords;
                            rd2 = Common::sqr(d[0]) + Common::sqr(d[1]) + Common::sqr(d[2]);

                            if (rd2 < APPROX_ZERO) {
                                rd2 = APPROX_ZERO;
                            }

                            inv_rd = 1. / sqrt(rd2);

                            d *= inv_rd;

                            Score::receptor->rvector2[ia][0] = Score::receptor->rvector[ia][1] * d[2] - Score::receptor->rvector[ia][2] * d[1];
                            Score::receptor->rvector2[ia][1] = Score::receptor->rvector[ia][2] * d[0] - Score::receptor->rvector[ia][0] * d[2];
                            Score::receptor->rvector2[ia][2] = Score::receptor->rvector[ia][0] * d[1] - Score::receptor->rvector[ia][1] * d[0];

                            rd2 = 0.;

                            rd2 = Common::sqr(Score::receptor->rvector2[ia][0]) + Common::sqr(Score::receptor->rvector2[ia][1]) + Common::sqr(Score::receptor->rvector2[ia][2]);

                            if (rd2 < APPROX_ZERO) {
                                rd2 = APPROX_ZERO;
                            }

                            inv_rd = 1. / sqrt(rd2);

                            Score::receptor->rvector2[ia] *= inv_rd;

                        }
                    }
                }
            }

            if (nbond == 2)
            {
                Atom* i1_atom = Score::receptor->atoms[i1];

                Atom* i2_atom = Score::receptor->atoms[i2];

                float rd2 = 0.;

                Score::receptor->rvector2[ia] = i2_atom->coords - i1_atom->coords;
                rd2 = Common::sqr(Score::receptor->rvector2[ia][0]) + Common::sqr(Score::receptor->rvector2[ia][1]) + Common::sqr(Score::receptor->rvector2[ia][2]);

                if (rd2 < APPROX_ZERO) {
                    rd2 = APPROX_ZERO;
                }

                float inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector2[ia] *= inv_rd;

                float rdot = 0.;

                for (int i = 0; i < 3; ++i) {
                    rdot += (ia_atom->coords[i] - i1_atom->coords[i]) * Score::receptor->rvector2[ia][i];
                }

                rd2 = 0.;

                Score::receptor->rvector[ia] = ia_atom->coords - ((rdot*Score::receptor->rvector2[ia]) + i1_atom->coords);

                rd2 = Common::sqr(Score::receptor->rvector[ia][0]) + Common::sqr(Score::receptor->rvector[ia][1]) + Common::sqr(Score::receptor->rvector[ia][2]);

                if (rd2 < APPROX_ZERO) {
                    rd2 = APPROX_ZERO;
                }

                inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector[ia] *= inv_rd;

            }
        }
        else if (ad_type_property(ia_atom->ad).hbond == 4)
        {
            int nbond = 0;
            int i1 = 0, i2 = 0, i3 = 0;

            for (int ib = from; ib < to; ++ib)
            {
                if (ib != ia)
                {
                    Atom* ib_atom = Score::receptor->atoms[ib];
                    float rd2 = 0.;

                    glm::vec3 dc = ia_atom->coords - ib_atom->coords;
                    rd2 = Common::sqr(dc[0]) + Common::sqr(dc[1]) + Common::sqr(dc[2]);

                    if (((rd2 < 2.89) && ((ib_atom->ad != AD_TYPE_HD) && (ib_atom->ad != AD_TYPE_H))) ||
                        ((rd2 < 1.69) && ((ib_atom->ad == AD_TYPE_HD) || (ib_atom->ad == AD_TYPE_H)))) {

                        if (nbond == 2) {
                            nbond = 3;
                            i3 = ib;
                        }
                        if (nbond == 1) {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0) {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }
            }

            if (nbond == 0) {
                continue;
                //throw Error_report("Error: Nitrogen atom found with no bonded atoms, atom serial number " + (ia + 1) + ad_type_property(ia_atom->ad).name);
            }


            if (nbond == 1) {

                Atom* i1_atom = Score::receptor->atoms[i1];

                float rd2 = 0.;

                Score::receptor->rvector[ia] = ia_atom->coords - i1_atom->coords;
                rd2 = Common::sqr(Score::receptor->rvector[ia][0]) + Common::sqr(Score::receptor->rvector[ia][1]) + Common::sqr(Score::receptor->rvector[ia][2]);

                if (rd2 < APPROX_ZERO) {
                    rd2 = APPROX_ZERO;
                }

                float inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector[ia] *= inv_rd;

            }


            if (nbond == 2) {

                Atom* i1_atom = Score::receptor->atoms[i1];
                Atom* i2_atom = Score::receptor->atoms[i2];

                float rd2 = 0.;

                for (int i = 0; i < 3; ++i) {
                    Score::receptor->rvector[ia][i] = ia_atom->coords[i] - (i2_atom->coords[i] + i1_atom->coords[i]) / 2.;
                    rd2 += Common::sqr(Score::receptor->rvector[ia][i]);
                }

                if (rd2 < APPROX_ZERO) {
                    rd2 = APPROX_ZERO;
                }

                float inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector[ia] *= inv_rd;
            }

            if (nbond == 3) {

                Atom* i1_atom = Score::receptor->atoms[i1];
                Atom* i2_atom = Score::receptor->atoms[i2];
                Atom* i3_atom = Score::receptor->atoms[i3];

                float rd2 = 0.;

                for (int i = 0; i < 3; ++i) {
                    Score::receptor->rvector[ia][i] = ia_atom->coords[i] - (i1_atom->coords[i] + i2_atom->coords[i] + i3_atom->coords[i]) / 3.;
                    rd2 += Common::sqr(Score::receptor->rvector[ia][i]);
                }

                if (rd2 < APPROX_ZERO)
                {
                    rd2 = APPROX_ZERO;
                }

                float inv_rd = 1. / sqrt(rd2);

                Score::receptor->rvector[ia] *= inv_rd;

            }
        }
    }
}

void Score::ad_hbond_theta_eval(float& Hramp, float& racc, float& rdon, int atom_no, int recep_atom_no, std::vector<Atom*>& ligand_atoms)
{
    int closestH = 0;
    float rmin = 999999.;

    for (int ia = 0; ia < Score::receptor->atoms.size(); ++ia)
    {
        Atom* ia_atom = Score::receptor->atoms[ia];

        if ((ad_type_property(ia_atom->ad).hbond == 1) || (ad_type_property(ia_atom->ad).hbond == 2)) {

            glm::vec3 d = ia_atom->coords - ligand_atoms[atom_no]->coords;

            float r = sqrt(Common::sqr(d[0]) + Common::sqr(d[1]) + Common::sqr(d[2]));
            if (r < rmin) {
                rmin = r;
                closestH = ia;
            }
        }
    }

    Atom* ia_atom = Score::receptor->atoms[recep_atom_no];

    glm::vec3 d = ia_atom->coords - ligand_atoms[atom_no]->coords;

    float r = sqrt(Common::sqr(d[0]) + Common::sqr(d[1]) + Common::sqr(d[2]));

    if (r < APPROX_ZERO) {
        r = APPROX_ZERO;
    }

    float inv_r = 1. / r;
    float inv_rmax = 1. / (r > 0.5) ? r : 0.5;

    d *= inv_r;

    racc = 1.;
    rdon = 1.;
    // NEW2 Hramp ramps in Hbond acceptor probes
    Hramp = 1.;
    // END NEW2 Hramp ramps in Hbond acceptor probes

    if (ad_type_property(ia_atom->ad).hbond == 2)
    {
        float cos_theta = 0.;

        for (int i = 0; i < 3; ++i)
        {
            cos_theta -= d[i] * Score::receptor->rvector[recep_atom_no][i];
        }

        if (cos_theta <= 0.)
        {
            racc = 0.;
        }
        else {
            switch (Score::receptor->rexp[recep_atom_no]) {
            case 1:
            default:
                racc = cos_theta;
                break;
            case 2:
                racc = cos_theta*cos_theta;
                break;
            case 4:
                float tmp = cos_theta*cos_theta;
                racc = tmp*tmp;
                break;
            }

            // racc = pow( cos_theta, (float)rexp[ia]);
            // NEW2 calculate dot product of bond vector with bond vector of best hbond
            if (recep_atom_no == closestH) {
                Hramp = 1.;
            }
            else
            {
                cos_theta = 0.;
                for (int i = 0; i < 3; ++i) {
                    cos_theta += Score::receptor->rvector[closestH][i] * Score::receptor->rvector[recep_atom_no][i];
                }
                cos_theta = (cos_theta < 1.0) ? cos_theta : 1.0;
                cos_theta = (cos_theta > -1.0) ? cos_theta : -1.0;
                float theta = acos(cos_theta);
                Hramp = 0.5 - 0.5*cos(theta * 120. / 90.);
            } // ia test for closestH
        }
    }
    else if (ad_type_property(ia_atom->ad).hbond == 4)
    {
        // NEW Directional N acceptor
        //
        //  ia-th macromolecule atom = Nitrogen ( 4 = H )
        //  calculate rdon for H-bond Donor PROBES at this grid pt.
        //            ====     ======================
        //
        float cos_theta = 0.;
        //
        //  d[] = Unit vector from current grid pt to ia_th m/m atom.
        //  cos_theta = d dot rvector == cos(angle) subtended.
        //
        for (int i = 0; i < 3; ++i) {
            cos_theta -= d[i] * Score::receptor->rvector[recep_atom_no][i];
        }

        if (cos_theta <= 0.) {
            //
            //  H->current-grid-pt vector >= 90 degrees from
            //  X->N vector,
            //
            rdon = 0.;
        }
        else {
            //
            //  racc = [cos(theta)]^2.0 for H->N
            //
            rdon = cos_theta*cos_theta;
        }
        // endif (atom_type[ia] == nitrogen)
        // end NEW Directional N acceptor
    }
    else if (ad_type_property(ia_atom->ad).hbond == 5) {//A2// ia-th receptor atom = Oxygen  => receptor H-bond acceptor, oxygen.

        rdon = 0.;

        // check to see that probe is in front of oxygen, not behind
        float cos_theta = 0.;
        for (int i = 0; i < 3; ++i) {
            cos_theta -= d[i] * Score::receptor->rvector[recep_atom_no][i];
        }
        //
        // t0 is the angle out of the lone pair plane, calculated
        // as 90 deg - acos (vector to grid point DOT lone pair
        // plane normal)

        float t0 = 0.;
        for (int i = 0; i < 3; ++i) {
            t0 += d[i] * Score::receptor->rvector2[recep_atom_no][i];
        }
        if (t0 > 1.) {
            t0 = 1.;
            //(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", t0);
            //print_error( logFile, WARNING, message );Feb2012
        }
        else if (t0 < -1.) {
            t0 = -1.;
            //(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", t0);
            //print_error( logFile, WARNING, message );Feb2012
        }
        t0 = PI_halved - acos(t0);
        //
        // ti is the angle in the lone pair plane, away from the
        // vector between the lone pairs,
        // calculated as (grid vector CROSS lone pair plane normal)
        // DOT C=O vector - 90 deg

        glm::vec3 cross;

        cross[0] = d[1] * Score::receptor->rvector2[recep_atom_no][2] - d[2] * Score::receptor->rvector2[recep_atom_no][1];
        cross[1] = d[2] * Score::receptor->rvector2[recep_atom_no][0] - d[0] * Score::receptor->rvector2[recep_atom_no][2];
        cross[2] = d[0] * Score::receptor->rvector2[recep_atom_no][1] - d[1] * Score::receptor->rvector2[recep_atom_no][0];
        float rd2 = Common::sqr(cross[0]) + Common::sqr(cross[1]) + Common::sqr(cross[2]);
        if (rd2 < APPROX_ZERO) {
            rd2 = APPROX_ZERO;
        }
        float inv_rd = 1. / sqrt(rd2);
        float ti = 0.;
        for (int i = 0; i < 3; ++i) {
            ti += cross[i] * inv_rd * Score::receptor->rvector[recep_atom_no][i];
        }

        // rdon expressions from Goodford
        rdon = 0.;
        if (cos_theta >= 0.) {
            if (ti > 1.) {
                ti = 1.;
                //(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", ti);
                //print_error( logFile, WARNING, message );Feb2012
            }
            else if (ti < -1.) {
                ti = -1.;
                //(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", ti);
                //print_error( logFile, WARNING, message );Feb2012
            }
            ti = acos(ti) - PI_halved;
            if (ti < 0.) {
                ti = -ti;
            }
            // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
            rdon = (0.9 + 0.1*sin(ti + ti))*cos(t0);
        }
        else if (cos_theta >= -0.34202) {
            // 0.34202 = cos (100 deg)
            rdon = 562.25*pow(0.116978 - Common::sqr(cos_theta), 3.)*cos(t0);
        }
        // endif atom_type == oxygen, not disordered
    }
}

void Score::xs_hydrogen_prepare(Molecule* mol)
{
    for (int i = 0; i < mol->atoms.size(); ++i)
    {
        Atom* a = mol->atoms[i];

        if (xs_is_donor(a->xs))
        {
            for (int j = 0; j < a->atom_bonds.size(); ++j)
            {
                Atom* connected_atom_1 = a->get_connected_atom(j);

                if (connected_atom_1->ad == AD_TYPE_HD)
                {
                    a->xs_donor_hydrogen_pre_count++;
                }
            }
        }
    }
}

void Score::initializeMol(const std::string& path, Molecule* mol, bool is_receptor)
{
    Prepare_mol prepare;
    std::string file_name = Common::get_file_name(path);

    std::vector<std::string> elems2;
    Common::split(file_name, '.', elems2);

    std::string file_type =  elems2[1];

    if(Score::mol_type == 0 || Score::mol_type == 1 || Score::mol_type == 2 || Score::mol_type == 3 || Score::mol_type == 5 || Score::mol_type == 6 || Score::mol_type == 7 || Score::mol_type == 8) // X-Score HC = 0, VINA = 1, VINA_Halogen = 2, AutoDock = 3, RF-SCORE = 4, RF-SCORE & X_SCORE HC = 5, RF-SCORE & VINA = 6, RF-SCORE & VINA_Halogen = 7, RF-SCORE & AutoDock = 8
    {
        if(file_type != "pdbqt")
        {
            throw Error_report("Error: PDBQT is the only file type allowed for X-Score HC, VINA, VINA_Halogen, AutoDock, RF-SCORE & X_SCORE HC, RF-SCORE & VINA, RF-SCORE & VINA_Halogen and RF-SCORE & AutoDock. ");
        }
    }
    else if (Score::mol_type == 4) //RF-SCORE = 4
    {
        if (file_type != "pdb" && file_type != "sdf")
        {
            throw Error_report("Error: PDB Receptor & SDF Ligand are the only file types allowed for RF-SCORE. ");
        }

        if (is_receptor)
        {
            if (file_type != "pdb")
            {
                throw Error_report("Error: PDB is the only file type allowed for receptors in RF-SCORE. ");
            }
        }
        else
        {
            if (file_type != "sdf")
            {
                throw Error_report("Error: SDF is the only file type allowed for ligands in RF-SCORE. ");
            }
        }
    }
    else
    {
        throw Error_report("Error: Unknown Scoring Function Type. ");
    }

    mol->path = path;
    mol->format = file_type;

    if(file_type == "pdbqt")
    {
        parser->parse_pdb(path, mol, true);
        prepare.get_elements(mol);
        prepare.get_connections(mol);
        prepare.assign_types(mol);
    }
    else if (file_type == "pdb") //only for RF-SCORE so no need for connections or assigning types
    {
        parser->parse_pdb(path, mol, false);
        prepare.get_elements(mol);
    }
    else if (file_type == "sdf") //only for RF-SCORE so no need for connections or assigning types
    {
        parser->parse_sdf(path, mol);
    }
    else
    {
        throw Error_report("Error: Unknown file format");
    }
}

glm::vec3 Score::xs_hbond_root(Atom* a)
{
    glm::vec3 root = { 0, 0, 0 };

    unsigned int count = 0;

    if (xs_is_acceptor(a->xs) || xs_is_donor(a->xs))
    {
        for (int j = 0; j < a->atom_bonds.size(); ++j)
        {
            Atom* connected_atom = a->get_connected_atom(j);

            if (connected_atom->ad != AD_TYPE_HD)
            {
                root += connected_atom->coords;
                count++;
            }
        }
    }

    if (count == 0)
    {
        return root;
    }
    else
    {
        return root /= count;
    }
}

int Score::atomnumber(const char* atomname)
{
    int atnum = -9999;
    //All commented atom types did not occur once in the pdbbind refined and core set. Uncomment and
    //re-size the rf_inter 2D dimensions

    if (!strncmp(atomname, "C", 4)) atnum = 0;
    else if (!strncmp(atomname, "CA", 4)) atnum = 0;
    else if (!strncmp(atomname, "CB", 4)) atnum = 0;
    else if (!strncmp(atomname, "CD", 4)) atnum = 0;
    else if (!strncmp(atomname, "CD1", 4)) atnum = 0;
    else if (!strncmp(atomname, "CD2", 4)) atnum = 0;
    else if (!strncmp(atomname, "CE", 4)) atnum = 0;
    else if (!strncmp(atomname, "CE1", 4)) atnum = 0;
    else if (!strncmp(atomname, "CE2", 4)) atnum = 0;
    else if (!strncmp(atomname, "CE3", 4)) atnum = 0;
    else if (!strncmp(atomname, "CG", 4)) atnum = 0;
    else if (!strncmp(atomname, "CG1", 4)) atnum = 0;
    else if (!strncmp(atomname, "CG2", 4)) atnum = 0;
    else if (!strncmp(atomname, "CH2", 4)) atnum = 0;
    else if (!strncmp(atomname, "CZ", 4)) atnum = 0;
    else if (!strncmp(atomname, "CZ2", 4)) atnum = 0;
    else if (!strncmp(atomname, "CZ3", 4)) atnum = 0;
    else if (!strncmp(atomname, "O", 4)) atnum = 1;
    else if (!strncmp(atomname, "OD1", 4)) atnum = 1;
    else if (!strncmp(atomname, "OD2", 4)) atnum = 1;
    else if (!strncmp(atomname, "OE1", 4)) atnum = 1;
    else if (!strncmp(atomname, "OE1A", 4)) atnum = 1;
    else if (!strncmp(atomname, "OE1B", 4)) atnum = 1;
    else if (!strncmp(atomname, "OE2", 4)) atnum = 1;
    else if (!strncmp(atomname, "OG", 4)) atnum = 1;
    else if (!strncmp(atomname, "OG1", 4)) atnum = 1;
    else if (!strncmp(atomname, "OH", 4)) atnum = 1;
    else if (!strncmp(atomname, "OXT", 4)) atnum = 1;
    else if (!strncmp(atomname, "N", 4)) atnum = 2;
    else if (!strncmp(atomname, "NE", 4)) atnum = 2;
    else if (!strncmp(atomname, "NE1", 4)) atnum = 2;
    else if (!strncmp(atomname, "NE2", 4)) atnum = 2;
    else if (!strncmp(atomname, "NE2A", 4)) atnum = 2;
    else if (!strncmp(atomname, "NE2B", 4)) atnum = 2;
    else if (!strncmp(atomname, "ND1", 4)) atnum = 2;
    else if (!strncmp(atomname, "ND2", 4)) atnum = 2;
    else if (!strncmp(atomname, "NH1", 4)) atnum = 2;
    else if (!strncmp(atomname, "NH2", 4)) atnum = 2;
    else if (!strncmp(atomname, "NZ", 4)) atnum = 2;
    else if (!strncmp(atomname, "P", 4)) atnum = 3;
    else if (!strncmp(atomname, "S", 4)) atnum = 4;
    else if (!strncmp(atomname, "SD", 4)) atnum = 4;
    else if (!strncmp(atomname, "SG", 4)) atnum = 4;
    else if (!strncmp(atomname, "Cl", 4)) atnum = 5;
    else if (!strncmp(atomname, "H", 4)) atnum = 6;
    else if (!strncmp(atomname, "11HH", 4)) atnum = 6;
    else if (!strncmp(atomname, "12HH", 4)) atnum = 6;
    else if (!strncmp(atomname, "21HH", 4)) atnum = 6;
    else if (!strncmp(atomname, "22HH", 4)) atnum = 6;
    else if (!strncmp(atomname, "21HD", 4)) atnum = 6;
    else if (!strncmp(atomname, "22HD", 4)) atnum = 6;
    else if (!strncmp(atomname, "21HE", 4)) atnum = 6;
    else if (!strncmp(atomname, "22HE", 4)) atnum = 6;
    else if (!strncmp(atomname, "HO", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN3", 4)) atnum = 6;
    else if (!strncmp(atomname, "HA", 4)) atnum = 6;
    else if (!strncmp(atomname, "HA1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HA2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HB", 4)) atnum = 6;
    else if (!strncmp(atomname, "HB1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HB2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HB3", 4)) atnum = 6;
    else if (!strncmp(atomname, "HD1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HD2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HE", 4)) atnum = 6;
    else if (!strncmp(atomname, "HE1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HE2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HE3", 4)) atnum = 6;
    else if (!strncmp(atomname, "HG", 4)) atnum = 6;
    else if (!strncmp(atomname, "HG1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HG2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HH2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HN3", 4)) atnum = 6;
    else if (!strncmp(atomname, "HH", 4)) atnum = 6;
    else if (!strncmp(atomname, "HZ", 4)) atnum = 6;
    else if (!strncmp(atomname, "HZ1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HZ2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HZ3", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HD1", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HD1", 4)) atnum = 6;
    else if (!strncmp(atomname, "3HD1", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HD2", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HD2", 4)) atnum = 6;
    else if (!strncmp(atomname, "3HD2", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HE2", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HE2", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HG1", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HG1", 4)) atnum = 6;
    else if (!strncmp(atomname, "3HG1", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HG2", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HG2", 4)) atnum = 6;
    else if (!strncmp(atomname, "3HG2", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HH1", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HH1", 4)) atnum = 6;
    else if (!strncmp(atomname, "1HH2", 4)) atnum = 6;
    else if (!strncmp(atomname, "2HH2", 4)) atnum = 6;
    else if (!strncmp(atomname, "HH1", 4)) atnum = 6;
    else if (!strncmp(atomname, "HXT", 4)) atnum = 6;

    if (atnum == -9999)
    {
        if (!strncmp(atomname, "F", 4)) atnum = 7;
        else if (!strncmp(atomname, "Br", 4)) atnum = 8;
        else if (!strncmp(atomname, "I", 4)) atnum = 9;
        /*else if (!strncmp(atomname, "B", 4)) atnum = 10;
        else if (!strncmp(atomname, "Na", 4)) atnum = 11;
        else if (!strncmp(atomname, "Mg", 4)) atnum = 12;
        else if (!strncmp(atomname, "Al", 4)) atnum = 13;
        else if (!strncmp(atomname, "K", 4)) atnum = 14;
        else if (!strncmp(atomname, "Ca", 4)) atnum = 15;
        else if (!strncmp(atomname, "V", 4)) atnum = 16;
        else if (!strncmp(atomname, "M ", 4)) atnum = 17;
        else if (!strncmp(atomname, "Fe", 4)) atnum = 18;
        else if (!strncmp(atomname, "Co", 4)) atnum = 19;
        else if (!strncmp(atomname, "Ni", 4)) atnum = 20;
        else if (!strncmp(atomname, "Cu", 4)) atnum = 21;
        else if (!strncmp(atomname, "Zn", 4)) atnum = 22;
        else if (!strncmp(atomname, "As", 4)) atnum = 23;
        else if (!strncmp(atomname, "Se", 4)) atnum = 24;
        else if (!strncmp(atomname, "Cd", 4)) atnum = 25;
        else if (!strncmp(atomname, "Yb", 4)) atnum = 26;
        else if (!strncmp(atomname, "Au", 4)) atnum = 27;
        else if (!strncmp(atomname, "Hg", 4)) atnum = 28;
        else if (!strncmp(atomname, "Pb", 4)) atnum = 29;
        else if (!strncmp(atomname, "U", 4)) atnum = 30;
        else if (!strncmp(atomname, "Mn", 4)) atnum = 31;*/
        else
        {
            //This atom type will be ignored later. No need to throw an error.
            //throw Error_report("RF-SCORE Error: Atomname " + std::string(atomname) + " not recognised");
        }
    }

    return atnum;
}

void Score::calculate()
{
    if (Score::mol_type > 3 && Score::mol_type < 9)
    {
        std::vector<int> receptor_atomNumbers;

        for (unsigned int i = 0; i < Score::receptor->atoms.size(); ++i)
        {
            receptor_atomNumbers.push_back(Score::atomnumber(Score::receptor->atoms[i]->name.c_str()));
        }

        unsigned int** rf_inter_array;

        rf_inter_array = new unsigned int*[9];

        for (unsigned int j = 0; j < 9; ++j)
        {
            rf_inter_array[j] = new unsigned int[9];
        }

        for (unsigned int j = 0; j < 9; ++j)
        {
            for (unsigned int k = 0; k < 9; ++k)
            {
                rf_inter_array[j][k] = 0;
            }
        }

        for (unsigned int j = 0; j < Score::ligand->atoms.size(); ++j)
        {
            Atom* a = Score::ligand->atoms[j];
            int ligandAtomNumber = Score::atomnumber(a->name.c_str());

            for (unsigned int k = 0; k < Score::receptor->atoms.size(); ++k)
            {
                Atom* b = Score::receptor->atoms[k];

                if ((receptor_atomNumbers[k] >= 0) && (ligandAtomNumber >= 0))
                {
                    float ddum = std::sqrt(Common::vec_distance_sqr(a->coords, b->coords));

                    if (ddum < 12)
                    {
                        int mmax = rf_inter_array[receptor_atomNumbers[k]][ligandAtomNumber];
                        mmax++;
                        rf_inter_array[receptor_atomNumbers[k]][ligandAtomNumber] = mmax;
                    }
                }
            }
        }

        unsigned int rfXID = 0;

        for (unsigned int j = 0; j < 9; ++j)
        {
            for (unsigned int k = 0; k < 9; ++k)
            {
                Score::ligand->rf_inter.push_back(std::pair<unsigned int, unsigned int>(rfXID, rf_inter_array[j][k]));
                rfXID++;
            }
        }

        for (unsigned int j = 0; j < 9; ++j)
        {
            delete[] rf_inter_array[j];
        }

        delete[] rf_inter_array;

    }

    if(Score::mol_type == 0 || Score::mol_type == 5) // X-SCORE HC || RF-SCORE & X_SCORE HC
    {
        Score::xs_hydrogen_prepare(Score::receptor);

        int nat = num_atom_types(Atom_type::XS);

        std::vector<unsigned int> xs_recep_don_h_num_bonds;
        std::vector<unsigned int> xs_recep_acc_h_num_bonds;
        std::vector<unsigned int> xs_ligand_don_h_num_bonds;
        std::vector<unsigned int> xs_ligand_acc_h_num_bonds;

        for(int j = 0; j < Score::receptor->atoms.size(); ++j)
        {
            xs_recep_don_h_num_bonds.push_back(0);
            xs_recep_acc_h_num_bonds.push_back(0);
        }

        for(int j = 0; j < Score::ligand->atoms.size(); ++j)
        {
            xs_ligand_don_h_num_bonds.push_back(0);
            xs_ligand_acc_h_num_bonds.push_back(0);
        }

        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Van der Waals X-SCORE HC", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Hydrogen Bonding X-SCORE HC", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Hydrophobic Contact X-SCORE HC", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Number of Rotors X-SCORE HC", 0.0));

        Score::xs_hydrogen_prepare(Score::ligand);

        for(int j = 0; j < Score::ligand->atoms.size(); ++j)
        {
            Atom* a = Score::ligand->atoms[j];

            int t1 = a->get(Atom_type::XS);
            if (t1 >= nat) continue;

            for (int k = 0; k < Score::receptor->atoms.size(); ++k)
            {
                Atom* b = Score::receptor->atoms[k];

                int t2 = b->get(Atom_type::XS);
                if (t2 >= nat) continue;

                float r2 = Common::vec_distance_sqr(a->coords, b->coords);
                float r = std::sqrt(r2);

                if (r2 < 64)
                {
                    Score::ligand->terms_contributions[0].second += scoring_terms.vdw_eval<4, 8>(a, b, r, 0, 100, false);

                    bool xs_hydrogen_bond = false;

                    if (xs_h_bond_possible(a->xs, b->xs))
                    {

                        if (xs_is_donor(a->xs) && xs_is_acceptor(b->xs))
                        {
                            if ((xs_ligand_don_h_num_bonds[j] < a->xs_donor_hydrogen_pre_count)
                                && xs_recep_acc_h_num_bonds[k] < ad_type_property(b->ad).lone_pairs)
                            {
                                xs_ligand_don_h_num_bonds[j]++;
                                xs_recep_acc_h_num_bonds[k]++;
                                xs_hydrogen_bond = true;
                            }
                        }
                        else if (xs_is_acceptor(a->xs) && xs_is_donor(b->xs))
                        {
                            if ((xs_recep_don_h_num_bonds[k] < b->xs_donor_hydrogen_pre_count)
                                && xs_ligand_acc_h_num_bonds[j] < ad_type_property(a->ad).lone_pairs)
                            {
                                xs_recep_don_h_num_bonds[k]++;
                                xs_ligand_acc_h_num_bonds[j]++;
                                xs_hydrogen_bond = true;
                            }
                        }

                        if (xs_hydrogen_bond)
                        {
                            glm::vec3 root_1 = Score::xs_hbond_root(a);
                            glm::vec3 root_2 = Score::xs_hbond_root(b);

                            glm::vec3 vector_11 = a->coords - root_1;
                            glm::vec3 vector_12 = a->coords - b->coords;

                            glm::vec3 vector_21 = b->coords - root_2;
                            glm::vec3 vector_22 = b->coords - a->coords;

                            float theta_1 = Common::angle_radian(vector_11, vector_12) * 180.0 / pi;
                            float theta_2 = Common::angle_radian(vector_21, vector_22) * 180.0 / pi;

                            Score::ligand->terms_contributions[1].second += scoring_terms.eval_dir_h_bond_XS(a, b, r, theta_1, theta_2);
                        }
                    }

                    Score::ligand->terms_contributions[2].second += scoring_terms.evalHydrophobicContactXS(a, b, r);
                }
            }
        }

        float num_tors = 0;
        float num_rotors = 0;
        scoring_terms.calculate_conf_independent(Score::ligand, num_tors, num_rotors);

        Score::ligand->terms_contributions[Score::ligand->terms_contributions.size() - 1].second = num_rotors;
    }
    else if(Score::mol_type == 1 || Score::mol_type == 2 || Score::mol_type == 6 || Score::mol_type == 7) // VINA || VINA_Halogen || RF-SCORE & VINA || RF-SCORE & VINA_Halogen
    {
        int nat = num_atom_types(Atom_type::XS);

        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Gauss 1 VINA", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Gauss 2 VINA", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Repulsion VINA", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Hydrophobic VINA", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Non Directional Hydrogen Bond VINA", 0.0));

        if(Score::mol_type == 2 || Score::mol_type == 7)
        {
            Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Halogen Chlorine VINAXB", 0.0));
            Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Halogen Bromine VINAXB", 0.0));
            Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Halogen Iodine VINAXB", 0.0));
        }

        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Number of Tors VINA", 0.0));

        for(int j = 0; j < Score::ligand->atoms.size(); ++j)
        {
            Atom* a = Score::ligand->atoms[j];

            int t1 = a->get(Atom_type::XS);
            if (t1 >= nat) continue;

            for(int k = 0; k < Score::receptor->atoms.size(); ++k)
            {
                Atom* b = Score::receptor->atoms[k];

                int t2 = b->get(Atom_type::XS);
                if (t2 >= nat) continue;

                float r2 = Common::vec_distance_sqr(a->coords, b->coords);
                float r = std::sqrt(r2);

                if (r2 < 64)
                {
                    Score::ligand->terms_contributions[0].second += scoring_terms.gaussVINA_eval(a, b, r, 0, 0.5);
                    Score::ligand->terms_contributions[1].second += scoring_terms.gaussVINA_eval(a, b, r, 3, 2.0);
                    Score::ligand->terms_contributions[2].second += scoring_terms.repulsionVINA_eval(a, b, r, 0.0);
                    Score::ligand->terms_contributions[3].second += scoring_terms.hydrophobicVINA_eval(a, b, r, 0.5, 1.5);
                    Score::ligand->terms_contributions[4].second += scoring_terms.non_dir_h_bondVINA_eval(a, b, r, -0.7, 0);

                    if(Score::mol_type == 2 || Score::mol_type == 7)
                    {
                        float halogen_theta = 180;

                        if (xs_hal_any_bond_possible(a->xs, b->xs)) {     //if there is a halogen bond possible, we need the angle
                            if (xs_is_halogen(a->xs)) {                  //a is the halogen, b is the acceptor
                                for (int k = 0; k < a->atom_bonds.size(); ++k)
                                {
                                    Atom* c = a->get_connected_atom(k);
                                    if (c->el == EL_TYPE_C) {              //the halogen needs to be connected to a carbon
                                        glm::vec3 v1 = a->coords - c->coords;  //vector between halogen and carbon
                                        glm::vec3 v2 = a->coords - b->coords;  //vector between halogen and acceptor
                                        halogen_theta = Common::angle_radian(v1, v2) * 180.0 / pi; //convert to theta
                                    }
                                }
                            }
                            else {                                  //b is the halogen, a is the acceptor (might never happen?)
                                for (int k = 0; k < b->atom_bonds.size(); ++k)
                                {
                                    Atom* c = b->get_connected_atom(k);
                                    if (c->el == EL_TYPE_C) {              //the halogen needs to be connected to a carbon
                                        glm::vec3 v1 = b->coords - c->coords;  //vector between halogen and carbon
                                        glm::vec3 v2 = b->coords - a->coords;  //vector between halogen and acceptor
                                        halogen_theta = Common::angle_radian(v1, v2) * 180.0 / pi; //convert to theta
                                    }
                                }
                            }
                        }

                        Score::ligand->terms_contributions[5].second += scoring_terms.halogen_cl_bondVINAXB_eval(a, b, r, halogen_theta, -0.25, -0.15, 0);
                        Score::ligand->terms_contributions[6].second += scoring_terms.halogen_br_bondVINAXB_eval(a, b, r, halogen_theta, -0.45, -0.35, -0.25, 0);
                        Score::ligand->terms_contributions[7].second += scoring_terms.halogen_i_bondVINAXB_eval(a, b, r, halogen_theta, -0.55, -0.45, -0.35, 0);
                    }
                }
            }
        }

        float num_tors = 0;
        float num_rotors = 0;
        scoring_terms.calculate_conf_independent(Score::ligand, num_tors, num_rotors);

        Score::ligand->terms_contributions[Score::ligand->terms_contributions.size() - 1].second = 1 + num_tors / 5.0;

    }
    else if(Score::mol_type == 3 || Score::mol_type == 8) // AutoDock || RF-SCORE & AutoDock
    {
        Score::prepare_ad_hydrogen();

        int nat = num_atom_types(Atom_type::AD);

        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Van der Waals AutoDock", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Hydrogen Bonding AutoDock", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Electrostatic AutoDock", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Solvation Charge Dependent AutoDock", 0.0));
        Score::ligand->terms_contributions.push_back(std::pair<std::string, float>("Number of Tors AutoDock", 0.0));

        for(int j = 0; j < Score::ligand->atoms.size(); ++j)
        {
            Atom* a = Score::ligand->atoms[j];

            int t1 = a->get(Atom_type::AD);
            if (t1 >= nat) continue;

            for(int k = 0; k < Score::receptor->atoms.size(); ++k)
            {
                Atom* b = Score::receptor->atoms[k];

                int t2 = b->get(Atom_type::AD);
                if (t2 >= nat) continue;

                float r2 = Common::vec_distance_sqr(a->coords, b->coords);
                float r = std::sqrt(r2);

                if (r2 < 64)
                {

                    Score::ligand->terms_contributions[0].second += scoring_terms.vdw_eval<6, 12>(a, b, r, 0.5, 100, true);
                    float e_hbond = scoring_terms.hydrogen_bondAD_eval(a, b, r, 0.5, 100);

                    float racc = 1.;
                    float rdon = 1.;
                    float Hramp = 1.;

                    Score::ad_hbond_theta_eval(Hramp, racc, rdon, j, k, Score::ligand->atoms);

                    float rsph = e_hbond / 100.;
                    rsph = ((rsph > 0.) ? rsph : 0.);
                    rsph = ((rsph < 1.) ? rsph : 1.);

                    if ((ad_type_property(a->ad).hbond == 3 || ad_type_property(a->ad).hbond == 5)
                        && (ad_type_property(b->ad).hbond == 1 || ad_type_property(b->ad).hbond == 2))
                    {
                        e_hbond *= Hramp * (racc + (1. - racc)*rsph);
                    }
                    else if ((ad_type_property(a->ad).hbond == 4)
                        && (ad_type_property(b->ad).hbond == 1 || ad_type_property(b->ad).hbond == 2))
                    {
                        e_hbond *= (racc + (1. - racc)*rsph);
                    }
                    else if ((ad_type_property(b->ad).hbond > 2)
                        && (ad_type_property(a->ad).hbond == 1 || ad_type_property(a->ad).hbond == 2))
                    {
                        e_hbond *= (rdon + (1. - rdon)*rsph);
                    }

                    Score::ligand->terms_contributions[1].second += e_hbond;
                }

                ////no cutoff
                Score::ligand->terms_contributions[2].second += scoring_terms.electrostaticAD_eval(a, b, r);
                Score::ligand->terms_contributions[3].second += scoring_terms.solvationAD_eval(a, b, r, 3.6, 0.01097, true);
            }
        }

        float num_tors = 0;
        float num_rotors = 0;
        scoring_terms.calculate_conf_independent(Score::ligand, num_tors, num_rotors);

        Score::ligand->terms_contributions[Score::ligand->terms_contributions.size() - 1].second = num_tors;

        Score::receptor->rexp.clear();
        Score::receptor->rvector.clear();
        Score::receptor->rvector2.clear();
    }
    else if (Score::mol_type == 4) // RF-SCORE "To not throw an error for RF-SCORE"
    {

    }
    else
    {
        throw Error_report("Error: Unknown Scoring Function Type. ");
    }
}

Score::Score(const std::string& rigid_path, const std::string& ligand_path, Parser* _parser, unsigned int type_)
{
    Score::mol_type = type_;
    Score::parser = _parser;

    Score::receptor = new Receptor();
    Score::initializeMol(rigid_path, Score::receptor, true);

    Score::ligand = new Ligand();
    Score::initializeMol(ligand_path, Score::ligand, false);
}

Score::~Score()
{
    delete Score::receptor;
    delete Score::ligand;
}
