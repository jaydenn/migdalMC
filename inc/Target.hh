
#ifndef TARGET_HH
#define TARGET_HH

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <iostream>
#include <sys/time.h>
#include <assert.h>

class Target {
    public:
        
        int migType;
        int migOptimize;
        double maxE;

        //target properties
        double Mt;
        double Leff;

        int elementZ;
        int isoN;
        std::vector<double> isoTable;
        std::vector<double> isoFracTable;
        std::vector<std::vector<double>> nl_energy;
        std::vector<int> Nlmax;
        
        //interpolation objects
        gsl_spline *Znl_spline[6][3];
        gsl_interp_accel *Znl_accel[6][3];
        gsl_spline *Znl_int_spline[6][3];
        gsl_interp_accel *Znl_int_accel[6][3];
        gsl_spline *invZnl_int_spline[6][3];
        gsl_interp_accel *invZnl_int_accel[6][3];

        double Znl_y_max[6][3];
        double Znl_x_max[6][3];
        double Znl_x_min[6][3];
        double maxCumulativeProb=1;
        
        //random number generator
        const gsl_rng_type * T;
        gsl_rng * r;
        struct timeval tv;

        //methods
        double Z_nl(int N, int l, double Qe, double Ee);
        double Z_nl_max(int N, int l, double Qe);
        double invZ_nl_integrated(int N, int l, double Qe, double Znl);
        double Z_nl_integrated(int N, int l, double Qe, double max_Ee);
        double totalMigProb(double maxE, double ENR, double MT);
        std::vector<double> rand_migdalE(double ERnr);
        double rand_isotope(double ENR);
        double ERmax(double Ein, double EM);
        double ERmax(double Ein, double EM, double MT);
        double EMmax(double Ein);
        double EMmax(double Ein, double MT);
        double EMmaxER(double Ein, double ER);
        double EMmaxER(double Ein, double ER, double MT);
        
        int init();
        
        Target(int Z, int type, double E, int op){
            elementZ = Z;
            migType = type;
            maxE = E;
            migOptimize = op;
            
            //initialize random number generator
            gsl_rng_env_setup();
            gettimeofday(&tv,0);     //use time as seed
            unsigned long mySeed = tv.tv_sec + tv.tv_usec;
            T = gsl_rng_default; // Generator setup
            r = gsl_rng_alloc (T);
            gsl_rng_set(r, mySeed);
            
            //initialize migdal calc
            if(init()<0)
                assert(0);
        };
};

#endif /* TARGET_HH */

