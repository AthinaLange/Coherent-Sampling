#include "MultiPathProcessing.h"
#include "functions.h"

// Global Data
PathInfoQueue_t  path_info_queue;
PathDataVector_t multi_paths_data;

// =============================================================================
// VARIABLES
// =============================================================================
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
extern double de;
extern double alpha;
extern double sina;
extern double cosa;
extern double Pdotdhat;
extern double abs_d;


extern double *RR;
extern double *PP;
extern double *dhat;
extern double *Pperp;
extern double *R1;
extern double *v;

extern complex<double> z[4];
extern int S[4];

extern double (*www[2][4][4])();
extern double (*phi)(double*, double*);
extern double (*dens_init[4])(double*, double*);
extern double (*obs[4])(double*, double*);
extern double (*obs1[4])(double*, double*);

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

complex<double> I(0,1);

double **abszsum1;
double **argzsum1;
double **habszsum1;
double **hargzsum1;


// =============================================================================
// Path Processing
// =============================================================================

void process_path(PathInfo& path_info) {

    abszsum1 = new double*[N_PATHS];
    argzsum1 = new double*[N_PATHS];
    habszsum1 = new double*[N_PATHS];
    hargzsum1 = new double*[N_PATHS];
    for (int i = 0; i < N_PATHS; ++i){
        abszsum1[i] = new double[N_slice];
        argzsum1[i] = new double[N_slice];
        habszsum1[i] = new double[N_slice];
        hargzsum1[i] = new double[N_slice];
        for (int j = 0; j < N_slice; ++j){
            abszsum1[i][j] = 0.0;
            argzsum1[i][j] = 0.0;
            habszsum1[i][j] = 0.0;
            hargzsum1[i][j] = 0.0;
        }

        if (i == path_info.surface){
            z[i]=path_info.probability;
        }
        else{
            z[i] = {1.0,1.0};
        }
    }

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Process Path " << path_info.id << " (level " << path_info.level << " )" << endl;
    change = 0; // start loop calculation
    counter = path_info.clock;

    // Read ancestors PathData
    long ancestor_id = path_info.parent_id;
    while (ancestor_id >= 0) { // root path id is 0, parent of root path id is -1
        cout << "ancestor: " << ancestor_id << endl;
        // How to access PathData
        // multi_paths_data[ancestor_id].valid ... always valid
        // multi_paths_data[ancestor_id].parent_id
        // multi_paths_data[ancestor_id].data1D[i]
        // multi_paths_data[ancestor_id].data2D[i][j]
        ancestor_id = multi_paths_data[ancestor_id].parent_id; // next ancestor
    }


    // =========================================================================
    // Phase Space Calculation
    // =========================================================================

    while (change == 0 && counter < N_slice){
        cout << "Stored Inital Surface Value " <<path_info.surface << endl;
        SS0 = path_info.surface; // Put new surface as 'original' surface

        cout << "Counter: " << counter << endl;

        cout << "Starting z" << endl;
        for (int i = 0; i < N_PATHS;++i){
            cout << z[i] << endl;
        }
        // Trotter-Suziki Approx. from exp(iLd/2) =========================================
        phase0 = U(RR, PP, SS0, TSLICE*0.5);
        z[SS0] *= exp(I * phase0);

        cout << "1st propogator: surface: " << SS0 << endl;
        for (int i = 0; i < N_PATHS;++i){
            cout << z[i] << endl;
        }

        // Trotter-Suziki Approx. from exp(iJd) ===========================================
        dd(RR); // non-adiabatic coupling matrix
        de = dE(RR); // energy
        alpha = 0.0;
        Pdotdhat = 0;
        for (int i = 0; i < N_bath; ++i) {
            Pdotdhat += PP[i] * dhat[i]; //parallel component of dhat to momentum
        }
        alpha = Pdotdhat * abs_d * TSLICE;

        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if neg, 1 if pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = PP[i] - signPdotdhat * Pdotdhat * dhat[i]; // perp component of dhat
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        // Probability calculation ============================================
        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;
        double xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        //alpha goes to 0, pdotdhat very small, matrix becomes identiy and prob of jumping goes to 0

        cout << "Prob:" << "ap0: " << ap0 <<" ap1: "<< ap1 <<" ap2: " << ap2 <<" ap3: "<< ap3 << endl;
        SS2 = SS0;

        // Probability Weighting ==============================================
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

        // DECISION ============================================================================

        // ========================================================================
        // NO JUMP
        // ========================================================================

        if (xx < prob0){
            SS1 = SS0;
            cout << "Propogating Adiabatically" << endl;
            S[SS1] = SS1;
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
            cout << "after jump propogator: SS1: " << SS1 << endl;
            for (int i = 0; i < N_PATHS;++i){
                cout << z[i] << endl;
            }
            path_info.surface = SS1;

            // Changing values for transition matrices ============================

            if (www[1][SS0][SS1]() != 9999.0)
                for (int i = 0; i < N_bath; ++i)
                    PP[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1]() * dhat[i];

            // Trotter-Suziki Approx. from exp(iLd/2) =============================
            phase0 = U(RR,PP,SS1,TSLICE*0.5); // exp(iLd/2) (after jump)
            z[SS1] *= exp(I*phase0);
            cout << "after 3rd propogator: " << endl;
            for (int i = 0; i < N_PATHS;++i){
                cout << z[i] << endl;
            }
         /*   for (int i = 0; i < N_PATHS; ++i){
                cout << "Probabilities for each surface "<< i << ": " << z[i] << endl;
            }*/
            // Calculating new phase space points =================================
            phi = obs[SS1];
            abszsum1[SS1][counter]  = real(z[SS1]*phi(RR,PP)*dens_init[SS3](R1,v));
            argzsum1[SS1][counter]  = imag(z[SS1]*phi(RR,PP)*dens_init[SS3](R1,v));

            phi = obs1[SS1];
            habszsum1[SS1][counter]  = real(z[SS1]*phi(RR,PP)*dens_init[SS3](R1,v));
            hargzsum1[SS1][counter]  = imag(z[SS1]*phi(RR,PP)*dens_init[SS3](R1,v));
            cout << endl;
        }


            // ========================================================================
            // POSSIBLE JUMP
            // ========================================================================
        else{
            change = 1; // Has changed surface
            Njump++;
            cout << "Jump Possible" << endl;
            cout << "Propability Adiabatic: " << prob0 << ", Probability Jump: " << prob1 << endl;

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
            cout << "after jump propogator" << endl;
            for (int i = 0; i < N_PATHS;++i){
                cout << z[i] << endl;
            }

        /*    for (int k = 0; k < N_PATHS; ++k){
                cout << "Probabilities for each surface "<< k << ": " << z[k] << endl;
            }*/
            // need to fix probabilities index

            cout << "after 3rd propogator" << endl;
            for (int i = 0; i < N_PATHS; ++i) {
                if (www[1][SS0][S[i]]() != 9999.0){
                    for (int j = 0; j < N_bath; ++j){
                        PP[j] = Pperp[j] + signPdotdhat * www[1][SS0][S[i]]() * dhat[j];
                    }
                }

                // Trotter-Suziki Approx. from exp(iLd/2) ========================================
                phase0 = U(RR, PP, S[i], TSLICE * 0.5); // exp(iLd/2) (after jump)
                z[S[i]] *= exp(I * phase0);


                cout << z[i] << endl;

                // Calculating new phase space points =================================
                phi = obs[S[i]];
                abszsum1[S[i]][counter] = real(z[S[i]] * phi(RR, PP) * dens_init[SS3](R1, v));
                argzsum1[S[i]][counter] = imag(z[S[i]] * phi(RR, PP) * dens_init[SS3](R1, v));


                phi = obs1[S[i]];
                habszsum1[S[i]][counter] = real(z[S[i]] * phi(RR, PP) * dens_init[SS3](R1, v));
                hargzsum1[S[i]][counter] = imag(z[S[i]] * phi(RR, PP) * dens_init[SS3](R1, v));
            }
            counter++;
            break;
        }

        counter++;
    } ;

    // =======================================================================================================

    // Write PathData
    cout << "writing" << endl;

    path_info.clock = counter;

    multi_paths_data[path_info.id].valid = true;
    multi_paths_data[path_info.id].parent_id = path_info.parent_id;
    for (int i = 0; i < N_PATHS; ++i){
        multi_paths_data[path_info.id].probability[i] = z[i];
        multi_paths_data[path_info.id].surface[i] = S[i];
    }

    /* for (long i=0; i< multi_paths_data[path_info.id].n_data1D; ++i) {
         multi_paths_data[path_info.id].data1D[i] = path_info.id;
     }*/
    for (int i = 0; i < multi_paths_data[path_info.id].n_data2D_2; ++i) {
        multi_paths_data[path_info.id].abszsum1[0][i] = abszsum1[0][i];
        multi_paths_data[path_info.id].argzsum1[0][i] = argzsum1[0][i];
        multi_paths_data[path_info.id].habszsum1[0][i] = habszsum1[0][i];
        multi_paths_data[path_info.id].hargzsum1[0][i] = hargzsum1[0][i];

        multi_paths_data[path_info.id].abszsum1[1][i] = abszsum1[1][i];
        multi_paths_data[path_info.id].argzsum1[1][i] = argzsum1[1][i];
        multi_paths_data[path_info.id].habszsum1[1][i] = habszsum1[1][i];
        multi_paths_data[path_info.id].hargzsum1[1][i] = hargzsum1[1][i];

        multi_paths_data[path_info.id].abszsum1[2][i] = abszsum1[2][i];
        multi_paths_data[path_info.id].argzsum1[2][i] = argzsum1[2][i];
        multi_paths_data[path_info.id].habszsum1[2][i] = habszsum1[2][i];
        multi_paths_data[path_info.id].hargzsum1[2][i] = hargzsum1[2][i];

        multi_paths_data[path_info.id].abszsum1[3][i] = abszsum1[3][i];
        multi_paths_data[path_info.id].argzsum1[3][i] = argzsum1[3][i];
        multi_paths_data[path_info.id].habszsum1[3][i] = habszsum1[3][i];
        multi_paths_data[path_info.id].hargzsum1[3][i] = hargzsum1[3][i];

    }

    if ((path_info.Njump < N_JUMPS) &&  (counter < N_slice)) {

        // Calculate the following paths ids using a formula (vs. using a shared counter)
        // to avoid synchronization between parallel executions of process_path().
        // lowest path id in current level
        long id_min_level = 0;
        for (int l=0; l<path_info.level; ++l) id_min_level += pow(N_PATHS, l);
        // lowest path id in next level
        long id_min_next_level = id_min_level + pow(N_PATHS, path_info.level);
        // first path id of following paths in next level
        long path_id = id_min_next_level + (path_info.id - id_min_level)*N_PATHS;


        // Enqueue path following paths information
        //for (long p=0; p<(N_PATHS); ++p) {=
        // Random Number Generator
        // new random states based on different seeds create with random state from current path
        //unsigned long seed = path_info.random_state.uniform_int(0, path_info.random_state.MAX_INT);
        //cout << "Child Seed " << p << ": " << seed << endl;
        //RandomState random_state = RandomState(seed);
        //path_info_queue.emplace(PathInfo(path_info.id, path_id + p, S[p], z[p] path_info.level + 1, Njump+1, counter/*, random_state*/));
        // (parent_id, id, level, clock, random_state)
        //}
        // Pass on random state from current path to one of the following paths
        // after the random state has been used to generate new seeds for the following paths
        // to make sure that different seeds are generated in the following path.
        //path_info_queue.emplace(PathInfo(path_info.id, path_id + (N_PATHS-1), path_info.level+1, Njump, counter, path_info.random_state));
        // (parent_id, id, level, clock, random_state)
        cout << "Children paths created:"<< endl;
        for (int p = 0; p < 4; ++p){
            cout << path_id + p << endl;
            if (S[p] == SS0){
                path_info_queue.emplace(PathInfo(path_info.id, path_id + p, S[p], z[p], path_info.level + 1, Njump, counter));
            }
            else{
                path_info_queue.emplace(PathInfo(path_info.id, path_id + p, S[p], z[p], path_info.level + 1, Njump + 1, counter));
            }

        }
    }

    for (int j = 0; j < N_PATHS; ++j){
        delete [] abszsum1[j];
        delete [] argzsum1[j];
        delete [] habszsum1[j];
        delete [] hargzsum1[j];
    }
    delete [] abszsum1;
    delete [] argzsum1;
    delete [] habszsum1;
    delete [] hargzsum1;


}