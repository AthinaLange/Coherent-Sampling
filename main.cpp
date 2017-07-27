/*!
 *  \mainpage
 *  \brief Using a Trotter Approximation to calculate a quantum-classical non-adiabatic approximation to
 * separate the quantum and the quasi-classical degrees of freedom to allow a surface-hopping scheme to
 * be implemented that creates a tree-like structure.
 *
 * This code follows one of the possible paths through the tree structure. The non-adiabatic propagator
 * determines the likelihood of a jump based on importance sampling
 *
 * \author Donal MacKernan, Athina Lange and Philip McGrath.
 * \date 24.7.17
 */
#include <iostream>
#include <complex>
#include <fstream>

#include "Global.h"
#include "MultiPathProcessing.h"

#include "random.h"
#include "functions.h"

// Global Data
extern PathInfoQueue_t  path_info_queue;
extern PathDataVector_t multi_paths_data;

using namespace std;
#include <gsl/gsl_rng.h>

///////////////////////////////////////////////////////////////////////////////
/// VARIABLES
///////////////////////////////////////////////////////////////////////////////

char  datafilename[80];
//FILE   *stream;

const gsl_rng_type * TT; /*!< Random Number Generator seed based on Gaussian Distribution */
gsl_rng * rr;


int N_bath; /*!< Size of bath */
int N_slice; /*!< Number of time intervals */
int Ncut; /*!< Truncation parameter */
double timestep; /*!< Size of time interval */
int Nsample; /*!< Sample Size (No. of trees calculated) */
double delta; /*!< MD Integrating Timestep \f$ (\delta) \f$*/
double ppower; /*!< */

double ddd; /*!<  \f$ \delta^2 \f$ */
double ddd4; /*!< \f$ \delta^2/4 \f$ */
double *m; /*!< Mass of particles */
double *c; /*!< */
double *w; /*!< */
double *f; /*!< Force on particles */

double *dgam;
double *mww;
double *sig; /*!< Sigma/Variance */
double *RR;
double *PP;
double *R1;
double *v;
double *mu;
double TSLICE;

double Njump;
int SS0;
int SS1;
int SS2;
int SS3;
int counter;
complex<double> z[4] = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
int S[4] = {0,0,0,0};

double (*dens_init[4])(double*, double*); /*!< Initial Density Matrix*/
double (*obs[4])(double*, double*); /*!< Observable Matrix*/
double (*obs1[4])(double*, double*); /*!< Another Observable Matrix*/
void (*force[4])(double *); /*!< Hellman-Feynman Forces*/

extern double ranVector[10001];

double *abszsum1;
double *argzsum1;
double *habszsum1;
double *hargzsum1;

///////////////////////////////////////////////////////////////////////////////
/// MULTI-PATH PROCESSING PROGRAM - BREATH FIRST
///////////////////////////////////////////////////////////////////////////////


int main() {

    ///////////////////////////////////////////////////////////////////////////////
    /// VARIABLES
    //////////////////////////////////////////////////////////////////////////////////
    double w_max;
    double eta;
    double T;
    double beta; /*!< Inverse Temperature */
    int init_seed;


    ///////////////////////////////////////////////////////////////////////////////
    /// SYSTEM INPUT
    ///////////////////////////////////////////////////////////////////////////////
    /* Initialize Stream - Scope 0
    cout << "Print information about new stream: "<< endl;
    cout << "Input datafilename, N_bath, N_slice, Ncut" << endl;
    cout << "timestep, T, init_seed, Nsample" << endl;
    cout << "w_max, eta, beta, delta, power" << endl;
    cin >> datafilename >> N_bath >> N_slice >> Ncut >> timestep >> T >> init_seed >> Nsample >> w_max >> eta >> beta >> delta >> ppower;
    */
    cout << " Print information about new stream:" << endl;
    cout << "Input datafilename" << endl;
    cin >> datafilename;
    N_bath = 200;
    N_slice = 20;
    Ncut = 10;
    timestep = 0.05;
    T = 15;
    init_seed = 0;
    Nsample = 100;
    w_max = 3;
    eta = 0.13;
    beta = 25;
    delta = 0.8;
    ppower = 100000;

    ///////////////////////////////////////////////////////////////////////////////
    /// ALLOCATING MEMORY
    ///////////////////////////////////////////////////////////////////////////////
    mww = new double[N_bath];
    mu = new double[N_bath];
    sig =  new double[2*N_bath];
    dgam = new double[N_bath];
    R1 = new double[N_bath];
    v = new double[N_bath];
    f = new double[N_bath];
    c = new double[N_bath];
    m = new double[N_bath];
    w = new double[N_bath];
    RR = new double[N_bath];
    PP = new double[N_bath];
    abszsum1 = new double[N_slice];
    argzsum1 = new double[N_slice];
    habszsum1 = new double[N_slice];
    hargzsum1 = new double[N_slice];
    /*!< Initializing sum1 counters*/
    for (int i = 0; i < N_slice; ++i) {
        abszsum1[i] = 0.0;
        argzsum1[i] = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;

    }

    /*!< Sets up the use of a Gaussian Random Number Generator from GSL */
    gsl_rng_env_setup();
    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);

    /*!< Setting dimensions of PathData */
    long n_data1D= N_slice;
    unsigned long n_paths = pow(N_PATHS, (N_LEVELS+1.0)) - 1; /*!< n_paths: maximum possible no. of path segements */
    cout << "Max no. of paths: " << n_paths << endl;
    multi_paths_data.resize(n_paths, PathData(n_data1D));

    ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZATION OF SYSTEM
    ///////////////////////////////////////////////////////////////////////////////

    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;

    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;

    bath_para(eta,w_max);       /*!< Defining Bath Parameters */

    for (int i = 0; i < N_bath; ++i)
        mu[i] = beta*w[i]*0.5;
    for (int i = 0; i < N_bath; ++i){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
    }

    for (int i = 0; i < N_bath; ++i)
        sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));
    /*!< Defining force field */
    force[0] = F1;
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;

    /*! Defining non-adiabatic coupling matrix */
    setwww();




    ///////////////////////////////////////////////////////////////////////////////
    /// 1st OUTPUT
    ///////////////////////////////////////////////////////////////////////////////
    ofstream outputFile;
    outputFile.open("Filename.txt");
    outputFile << "w_max: " << w_max << ", eta: " << eta << ", beta: " << beta << ", delta: " << delta << ", ppower: " << ppower << ", N_bath: " << N_bath << ", N_slice: " << N_slice << endl;
    outputFile.close();


    ///////////////////////////////////////////////////////////////////////////////
    /// PROCESSING TREE - MPI DIRECTIVE - Breath-First
    ///////////////////////////////////////////////////////////////////////////////
    // for loop needed iteratating 10,000 times (entire tree calculated for every particle)

    // might be able to keep multi_paths_data internal to each particle, and final end result stored in seperate folder
    // known to everyone using MPI to parallize for loop

    for (int k = 0; k < Nsample; ++k) { /*!< Iterates over the Nsample size */
        cout << endl;
        cout << "===========================================================" << endl;
        cout << "===========================================================" << endl;
        cout << "tree no: " << k << endl;
        ///////////////////////////////////////////////////////////////////////////////
        /// INITIALIZATION OF INITIAL SURFACE
        ///////////////////////////////////////////////////////////////////////////////

        gauss_init_W(R1, v); /*!< Creates initial random sampling */
        double yy = 4.0 * (gsl_rng_uniform(rr));
        /*!< Setting initial surface value */
        if (yy < 1.0)
            SS3 = (SS0 = 0);
        else if (yy < 2.0) {
            SS0 = 1;
            SS3 = 2;
        } else if (yy < 3.0) {
            SS0 = 2;
            SS3 = 1;
        } else {
            SS3 = (SS0 = 3);
        }

        /*!< Allocating values for position and momentum from Gaussian */
        for (int l = 0; l < N_bath; ++l) {
            RR[l] = R1[l];
            PP[l] = v[l];
        }
        SS1 = SS0;
        z[SS0] = 4.0; /*!< Setting phase factor of initial surface to 4 */
        S[SS0] = SS1;
        counter = 0; /*!< Setting interval counter to 0 */
        Njump = 0;


        /*!< Place root path segment into processing queue */
        path_info_queue.emplace(PathInfo(-1, 0, S[SS0], z[SS0], 0, 0, 0));
        // (parent_id -1 for no parent, id, surface, phase, level, jump counter, clock)

        /*!< all path segments that are in the queue can be processed in parallel! */

        ///////////////////////////////////////////////////////////////////////////////
        /// ITERATION FOR TREE: CALCULATE EACH PATH SEGMENT - OPENMP DIRECTIVE
        ///////////////////////////////////////////////////////////////////////////////

        while (!path_info_queue.empty()) {



            /*!< Give thread PathInfo for path segment to be processed */
            PathInfo path_info = path_info_queue.front();
            path_info_queue.pop();

            /*!< Process path segment */
            process_path(path_info);
        }

        // SAMPLE OUTPUT

        //long path = n_paths-1; // choose any one between 0 and n_paths-1

/*    cout << endl;
    cout << "path " << path << endl;

    cout << "valid: " << multi_paths_data[path].valid << endl;
    cout << "parent_id: " << multi_paths_data[path].parent_id << endl;
/*    for (long i=0; i< multi_paths_data[path].n_data1D; ++i) {
        cout << "multi_paths_data[" << path << "].data1D[" << i << "]: " << multi_paths_data[path].data1D[i] << endl;
    }
    for (long i=0; i< multi_paths_data[path].n_data2D_1; ++i) {
        for (long j=0; j < multi_paths_data[path].n_data2D_2; ++j) {
            cout << "multi_paths_data[" << path << "].data2D[" << i << "][" << j << "]: " << multi_paths_data[path].data2D[i][j] << endl;
        }
    }*/
        /*       for (int i = 0; i < n_paths; ++i) {
               multi_paths_data[path_info.id].abszsum1[i] = abszsum1[i];
               multi_paths_data[path_info.id].argzsum1[i] = argzsum1[i];
               multi_paths_data[path_info.id].habszsum1[i] = habszsum1[i];
               multi_paths_data[path_info.id].hargzsum1[i] = hargzsum1[i];
           }*/

    }

    ///////////////////////////////////////////////////////////////////////////////
    /// 3rd OUTPUT
    ///////////////////////////////////////////////////////////////////////////////
    outputFile.open("Filename.txt");
    outputFile << "timestep: " << timestep << ", T: " << T << ", Nsample: " << Nsample << endl;
    for (int i = 0; i < N_slice; ++i){
        outputFile << "Njump: " << i <<  ", " << TSLICE*(i+1) <<"     " << abszsum1[i] << "   " << argzsum1[i] << "   " << habszsum1[i] <<"   " << hargzsum1[i] << endl;
        cout << "Njump: " << i <<  ", " << TSLICE*(i+1) << "     " << abszsum1[i] << "   " << argzsum1[i] << "   " << habszsum1[i] <<"   " << hargzsum1[i] << endl;

    }
    outputFile.close();


    ///////////////////////////////////////////////////////////////////////////////
    /// DEALLOCATING MEMEORY
    ///////////////////////////////////////////////////////////////////////////////
    delete [] mww;
    delete [] mu;
    delete [] sig;
    delete [] dgam;
    delete [] R1;
    delete [] v;
    delete [] f;
    delete [] c;
    delete [] m;
    delete [] w;
    delete [] RR;
    delete [] PP;

    delete [] abszsum1;
    delete [] argzsum1;
    delete [] habszsum1;
    delete [] hargzsum1;

    return 0;
}