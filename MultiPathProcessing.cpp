#include "MultiPathProcessing.h"
#include "functions.h"

// Global Data
PathInfoQueue_t  path_info_queue;
PathDataVector_t multi_paths_data;

///////////////////////////////////////////////////////////////////////////////
/// VARIABLES
///////////////////////////////////////////////////////////////////////////////
extern const gsl_rng_type * TT;
extern gsl_rng * rr;

extern int N_bath;
extern int N_slice;
extern int Njump;
extern int N_slice;
extern int SS0;
extern int SS1;
extern int SS2;
extern int SS3;
extern int counter;
extern double TSLICE;

double abs_d;

extern double *RR;
extern double *PP;
extern double *R1;
extern double *v;
extern double *m;

extern complex<double> z[4];
extern int S[4];

double (*www[2][4][4])(double cosa, double sina, double de, double Pdotdhat);
extern double (*dens_init[4])(double*, double*);
extern double (*obs[4])(double*, double*);
extern double (*obs1[4])(double*, double*);

extern double *abszsum1;
extern double *argzsum1;
extern double *habszsum1;
extern double *hargzsum1;


// =============================================================================
// Path Processing
// =============================================================================

void process_path(PathInfo& path_info) {

    ///////////////////////////////////////////////////////////////////////////////
    /// VARIABLES
    ///////////////////////////////////////////////////////////////////////////////
    int signPdotdhat;
    double phase0 = 0.0;
    double p0;
    double p1;
    double p2;
    double p3;
    double ap0;
    double ap1;
    double ap2;
    double ap3;
    double dn2;
    double prob0;
    double prob1;
    int change;
    double (*phi)(double*, double*);
    complex<double> I(0,1);
    complex<double> initd(0,0);
    initd = dens_init[SS3](R1,v);
    double *dhat;
    double *Pperp;
    dhat = new double[N_bath];
    Pperp = new double[N_bath];

    double de;
    double sina;
    double cosa;
    double Pdotdhat;

    /*!< Giving the surfaces starting phase factor*/
    for (int i = 0; i < N_PATHS; ++i){
        if (i == path_info.surface){ /*!< Current surface is given stored phase factor */
            z[i]=path_info.probability;
        }
        else{ /*!< Otherwise set to 1(dummy) */
            z[i] = {1.0,1.0};
        }
    }

    change = 0; /*!< Condition to check if jump completed or not, if so, exiting for loop*/
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Process Path " << path_info.id << " (level " << path_info.level << " )" << endl;

    counter = path_info.clock; /*!< Counter is given the current clock time */

    /*!< Ancestor of path segment is read */
    long ancestor_id = path_info.parent_id;
    while (ancestor_id >= 0) { // root path id is 0, parent of root path id is -1
        cout << "ancestor: " << ancestor_id << endl;
        ancestor_id = multi_paths_data[ancestor_id].parent_id; // next ancestor
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// PHASE SPACE CALCULATION
    ///////////////////////////////////////////////////////////////////////////////


    while (change == 0 && counter < N_slice){ /*!< While loop iterates until jump, or time runs out */
        //cout << "Stored Inital Surface Value " <<path_info.surface << endl;
        //SS0 = path_info.surface; /*!< Initial surface value is pulled from PathInfo */

        //cout << "Counter: " << counter << endl;

        ///////////////////////////////////////////////////////////////////////////////
        /// ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        phase0 = U(RR, PP, SS0, TSLICE*0.5); // exp(iLd/2) (before jump)
        z[SS0] *= exp(I * phase0);



        ///////////////////////////////////////////////////////////////////////////////
        /// NON-ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        dd(dhat, RR); /*!< calculating non-adiabatic coupling matrix */
        de = dE(RR); /*!< calculating energy */
        double alpha = 0.0;
        Pdotdhat = 0;
        for (int i = 0; i < N_bath; ++i) {
            Pdotdhat += PP[i] * dhat[i]; //parallel component of dhat to momentum
            alpha += PP[i] * dhat[i] / m[i];
        }
        alpha *= abs_d * TSLICE;

        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if neg, 1 if pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = PP[i] - signPdotdhat * Pdotdhat * dhat[i]; // perp component of dhat
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        /*!< Importance Sampling - non-adiabatic coupling matrix gives probabilities */
        ap0 = fabs(p0 = ((www[1][SS0][0](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][0](cosa, sina, de, Pdotdhat)));
        ap1 = fabs(p1 = ((www[1][SS0][1](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][1](cosa, sina, de, Pdotdhat)));
        ap2 = fabs(p2 = ((www[1][SS0][2](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][2](cosa, sina, de, Pdotdhat)));
        ap3 = fabs(p3 = ((www[1][SS0][3](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][3](cosa, sina, de, Pdotdhat)));
        dn2 = ap0 + ap1 + ap2 + ap3;
        double xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        //alpha goes to 0, pdotdhat very small, matrix becomes identiy and prob of jumping goes to 0
        //cout << "Prob:" << "ap0: " << ap0 <<" ap1: "<< ap1 <<" ap2: " << ap2 <<" ap3: "<< ap3 << endl;
        SS2 = SS0;

       /*!< Probability Weighting - probability of various surfaces calculated */
        if (SS0 == 0){
            prob0 = ap0/dn2;
            prob1 = (ap1 + ap2 + ap3)/dn2;
        }
        else if (SS0 == 1){
            prob0 = ap1/dn2;
            prob1 = (ap0 + ap2 + ap3)/dn2;
        }
        else if (SS0 == 2){
            prob0 = ap2/dn2;
            prob1 = (ap0 + ap1 + ap3)/dn2;
        }
        else {
            prob0 = ap3/dn2;
            prob1 = (ap0 + ap1 + ap2)/dn2;
        }

        ///////////////////////////////////////////////////////////////////////////////
        /// DECISION - NOT OR JUMP
        ///////////////////////////////////////////////////////////////////////////////

        // ==========================================================================
        // ==========================================================================
        // ==========================================================================

        ///////////////////////////////////////////////////////////////////////////////
        /// NO JUMP - PROPOGATES ADIABATICALLY
        ///////////////////////////////////////////////////////////////////////////////

        if (xx < prob0){ /*!< if random number less than probability of propagating adiabatically, stay on surface */
            SS1 = SS0;
           // cout << "Propogating Adiabatically" << endl;
            S[SS1] = SS1;
            /*!< Assign value of phase factor to 'new' surface */
            if (SS1 == 0){
                z[SS1] *= p0;
            }
            else if (SS1 == 1){
                z[SS1] *= p1;
            }
            else if (SS1 == 2){
                z[SS1] *= p2;
            }
            else {
                z[SS1] *= p3;
            }

            path_info.surface = SS1; /*!< Set 'new' surface value into the surface spot of PathInfo */

            if (www[1][SS0][SS1](cosa, sina, de, Pdotdhat) != 9999.0) /*! updating momentum values */
                for (int i = 0; i < N_bath; ++i)
                    PP[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1](cosa, sina, de, Pdotdhat) * dhat[i];


            ///////////////////////////////////////////////////////////////////////////////
            /// ADIABATIC PROPAGATOR
            ///////////////////////////////////////////////////////////////////////////////
            phase0 = U(RR,PP,SS1,TSLICE*0.5); // exp(iLd/2) (after jump)
            z[SS1] *= exp(I*phase0);


            ///////////////////////////////////////////////////////////////////////////////
            /// CALCULATING SUM VALUES FOR EVERY TIME INTERVAL (Solving Integral for Observable and Initial Density)
            ///////////////////////////////////////////////////////////////////////////////
            phi = obs[SS1];  /*!< Observable 1 function */
            abszsum1[counter]  += real(z[SS1]*phi(RR,PP)*initd);
            argzsum1[counter]  += imag(z[SS1]*phi(RR,PP)*initd);

            phi = obs1[SS1]; /*!< Observable 2 (Hamiltonian) function */
            habszsum1[counter]  += real(z[SS1]*phi(RR,PP)*initd);
            hargzsum1[counter]  += imag(z[SS1]*phi(RR,PP)*initd);
            //cout << endl;
            }

            ///////////////////////////////////////////////////////////////////////////////
            /// POSSIBLE JUMP
            ///////////////////////////////////////////////////////////////////////////////
        else{
            change = 1; /*!< Change exit condition of for loop - path segment ended */
            Njump++; /*!< Increase jump counter */
           // cout << "Jump Possible" << endl;
          //  cout << "Probability Adiabatic: " << prob0 << ", Probability Jump: " << prob1 << endl;

            /*!< Depending on which initial surface had been chosen, phase factor of all surfaces updated accordingly*/
            if (SS0 == 0){
                S[0] = 0;
                S[1] = 1;
                S[2] = 2;
                S[3] = 3;
                z[0] *= p0;
                z[1] *= p1/prob1;
                z[2] *= p2/prob1;
                z[3] *= p3/prob1;
            }
            else if (SS0 == 1){
                S[0] = 0;
                S[1] = 1;
                S[2] = 2;
                S[3] = 3;
                z[0] *= p0/prob1;
                z[1] *= p1;
                z[2] *= p2/prob1;
                z[3] *= p3/prob1;
            }
            else if (SS0 == 2){
                S[0] = 0;
                S[1] = 1;
                S[2] = 2;
                S[3] = 3;
                z[0] *= p0/prob1;
                z[1] *= p1/prob1;
                z[2] *= p2;
                z[3] *= p3/prob1;
            }
            else {
                S[0] = 0;
                S[1] = 1;
                S[2] = 2;
                S[3] = 3;
                z[0] *= p0/prob1;
                z[1] *= p1/prob1;
                z[2] *= p2/prob1;
                z[3] *= p3;
            }


            ///////////////////////////////////////////////////////////////////////////////
            /// IF JUMP: PROPAGATE ALONG NEW SURFACE AND CALCULATE SUM
            ///////////////////////////////////////////////////////////////////////////////
            for (int i = 0; i < N_PATHS; ++i) {
                if (www[1][SS0][S[i]](cosa, sina, de, Pdotdhat) != 9999.0){ /*! updating momentum values */
                    for (int j = 0; j < N_bath; ++j){
                        PP[j] = Pperp[j] + signPdotdhat * www[1][SS0][S[i]](cosa, sina, de, Pdotdhat) * dhat[j];
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////
                /// ADIABATIC PROPAGATOR
                ///////////////////////////////////////////////////////////////////////////////
                phase0 = U(RR, PP, S[i], TSLICE * 0.5); // exp(iLd/2) (after jump)
                z[S[i]] *= exp(I * phase0);


                ///////////////////////////////////////////////////////////////////////////////
                /// CALCULATING SUM VALUES FOR EVERY TIME INTERVAL (Solving Integral for Observable and Initial Density)
                ///////////////////////////////////////////////////////////////////////////////
                phi = obs[SS1]; /*!< Observable 1 function */
                abszsum1[counter]  += real(z[SS1]*phi(RR,PP)*initd);
                argzsum1[counter]  += imag(z[SS1]*phi(RR,PP)*initd);

                phi = obs1[SS1]; /*!< Observable 2 (Hamiltonian) function */
                habszsum1[counter]  += real(z[SS1]*phi(RR,PP)*initd);
                hargzsum1[counter]  += imag(z[SS1]*phi(RR,PP)*initd);
                }
            counter++; /*!< increase time interval counter */
            break;
        }

        counter++;
    } ;

    ///////////////////////////////////////////////////////////////////////////////
    /// WRITE PATH SEGMENT DATA INTO PathData
    ///////////////////////////////////////////////////////////////////////////////

    cout << "writing" << endl;

    path_info.clock = counter; /*!< Place current time interval into clock counter*/

    multi_paths_data[path_info.id].valid = true; /*!< Set path segment as having been completed */
    multi_paths_data[path_info.id].parent_id = path_info.parent_id; /*!< Assign parent's ID */
    for (int i = 0; i < N_PATHS; ++i){
        multi_paths_data[path_info.id].probability[i] = z[i]; /*!< Place calculated phase factors into vector */
        multi_paths_data[path_info.id].surface[i] = S[i]; /*!< Place surface values into vector */
    }

    for (int i = 0; i < multi_paths_data[path_info.id].n_data1D; ++i) { /*!< Set new version of sum1 into Path's matrix*/
        multi_paths_data[path_info.id].abszsum1[i] = abszsum1[i];
        multi_paths_data[path_info.id].argzsum1[i] = argzsum1[i];
        multi_paths_data[path_info.id].habszsum1[i] = habszsum1[i];
        multi_paths_data[path_info.id].hargzsum1[i] = hargzsum1[i];
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// CREATING CHILDREN
    ///////////////////////////////////////////////////////////////////////////////

    if ((path_info.Njump < N_JUMPS) &&  (counter < N_slice)) { /*!< If jump counter or time counter not at cut-off, spawn and create children paths(ID's) */

        /*!< Calculate children's path ID's */
        long id_min_level = 0;
        for (int l=0; l<path_info.level; ++l)
            id_min_level += pow(N_PATHS, l);
        long id_min_next_level = id_min_level + pow(N_PATHS, path_info.level); // lowest path ID in next level
        long path_id = id_min_next_level + (path_info.id - id_min_level)*N_PATHS; // first path ID of following paths in next level

        cout << "Children paths created:"<< endl;
        for (int p = 0; p < N_PATHS; ++p){ /*!< Place children's PathInfo into queue to be processed */
            cout << path_id + p << endl;
            if (S[p] == SS0){ /*!< If adiabatically propagated surface after jump, stored jump counter does not increase */
                path_info_queue.emplace(PathInfo(path_info.id, path_id + p, S[p], z[p], path_info.level + 1, Njump, counter));
            }
            else{ /*!< If a different surface than initial, stored jump counter increased */
                path_info_queue.emplace(PathInfo(path_info.id, path_id + p, S[p], z[p], path_info.level + 1, Njump + 1, counter));
            }
        }
    }

    delete [] dhat; delete [] Pperp;
}