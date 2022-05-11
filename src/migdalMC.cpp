using namespace std;

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
#include <iostream>
#include <array>

#include "Target.hh"

//const gsl_rng_type * T;
//gsl_rng * r;
//struct timeval tv;


int main(int argc, char** argv) 
{

    
    //initial set up variables
    int targetZ = -1;
    cout << "Enter target element Z = ";
    cin >> targetZ; 

    double maxE = -1;       //max ingoing particle's kinetic energy in keV
    cout << "Enter incoming neutron energy (keV) = ";
    cin >> maxE; 

    int migType = 1;        // 1 = mono-energetic heavy particle
    
    int migdalOptimize = 1; //set to 1 to scale Migdal probabilities
    cout << "Optimize Migdal probability? (1=yes,0=no)\n";
    cin >> migdalOptimize;
    
    //Create target object and initialize the Migdal calculation 
    Target myTarget(targetZ, migType, maxE, migdalOptimize);
    
    //output theoretical spectrum
    //TargetcalcMigdalSpectrum(&(spec.spline_spectrum_prep));

    //set up vector to store output and choose recoil energy
    vector<double> migdalE = {0,0,0,0};
    double ENR=-1; //in keV

    cout << "Enter nuclear recoil energy < Emax = " <<  myTarget.ERmax(maxE, 0) << " keV:\n";
    cin >> ENR; 
    
    //do calculation
    migdalE = myTarget.rand_migdalE(ENR);  // Returns tuple of [electron energy, binding energy of shell, shell index, mass number of recoiling atom]
    
    cout << "Electron energy (keV)\tBinding energy (keV)\tShell (n)\n";
    cout << migdalE[0] << "\t\t\t" << migdalE[1] << "\t\t\t" << migdalE[2] << endl; 
    
    return 1;

}
