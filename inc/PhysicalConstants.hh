
#ifndef PHYSICALCONSTANTS_HH
#define PHYSICALCONSTANTS_HH

/////////////////////////////////////////////////
//     Define some physical constants          //
/////////////////////////////////////////////////
const double Mn = 939565.42;       //neutron mass in keV
const double AMU = Mn/1.00727;     //1 atomic mass unit in keV
const double Me = 510.999;         //electron mass

//Atomic numbers for targets
const int XENON = 54;
const int ARGON = 18;

//Migdal type
const int NEUTRON = 1;
const int NEUTRINO = 2;

//mass numbers for targets
double MtXe = 131.29*AMU;    //RAM of xenon in keV
double MtAr = 39.948*AMU;    //RAM of argon in keV

//isotope data
const int isoNXe=7;
const vector<double> isoTableXe = {128,129,130,131,132,134,136}; //mass numbers of isotopes
const vector<double> isoFracTableXe = {0.0191421,0.283724,0.324514,0.536981,0.806574,0.911205,1}; //Cumulative abundance of naturally occuring isotopes of xenon

const int isoNAr=1;
const vector<double> isoTableAr = {40}; //atomic numbers of isotopes
const vector<double> isoFracTableAr = {1}; //naturally occuring isotopic abundances of argon

//atomic binding energies
const vector<vector<double>> nl_energyXe = {{    -1,    -1,    -1,-1}, //the atomic binding energies of shells indexed by [n][l]
                          {    35,    -1,    -1,-1},    
                          {   5.4,   4.9,    -1,-1},
                          {   1.1,  0.93,  0.66,-1},
                          {   0.2,  0.14,6.1e-2,-1},
                          {2.1e-2,9.8e-3,    -1,-1},
}; 
const vector<int> NlmaxXe = {-1,1,2,3,3,2};

const vector<vector<double>> nl_energyAr = {{  -1,    -1,    -1,-1}, //the atomic binding energies of shells indexed by [n][l]
                          {   3.2,    -1,    -1,-1},    
                          {   .32,   .24,    -1,-1},
                          {  .027,  .013,    -1,-1},
                          {    -1,    -1,    -1,-1},
                          {    -1,    -1,    -1,-1},
}; 
const vector<int> NlmaxAr = {-1,1,2,2,-1,-1}; 


double LeffXe = 0.15;          //approximate quenching factors
double LeffAr = 0.25;


/////////////////////////////////////////////////

#endif /* PHYSICALCONSTANTS_HH */

