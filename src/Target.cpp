
#include "Target.hh"
#include "PhysicalConstants.hh"

using namespace std;
using std::vector;


//reduced mass definition
double mu(double x, double y)
{
    return x*y/(x+y);
}

//electron recoil momentum
double qe(double ERnr, double MT)
{
    return Me*sqrt(2 * ERnr / MT);
}

double Target::ERmax(double Ein, double EM, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return pow(mu(MT,Mn),2)/(MT*Mn)*Ein*(pow(1-sqrt(1-EM/(mu(MT,Mn)*Ein/Mn)),2)+4*sqrt(1-EM/(mu(MT,Mn)*Ein/Mn)));
        case NEUTRINO:
            return pow(2.*Ein - EM,2)/(2.*(2.*Ein+MT));;
        default:
            return -1;
    }
    
}
double Target::ERmax(double Ein, double EM)
{
    return ERmax( Ein, EM, Mt);
}

double Target::EMmax(double Ein, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return mu(Mn,Mt)*Ein/Mn;
        case NEUTRINO:
            return Ein;
        default:
            return -1;
    }
    
}
double Target::EMmax(double Ein)
{
    return EMmax( Ein, Mt);
}
    
double Target::EMmaxER(double Ein, double ENR, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return 2.*sqrt(Ein*ENR*MT/Mn) - ENR*(MT+Mn)/Mn;
        case NEUTRINO:
            return Ein-ENR;
        default:
            return -1;
    }
    
}
double Target::EMmaxER(double Ein, double ENR)
{
    return EMmaxER( Ein, ENR, Mt);
}
    
//returns total prob below max_Ee
double Target::Z_nl_integrated(int N, int l, double Qe, double max_Ee)
{ 
    if ( max_Ee > 69.9 )
    {
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_int_spline[N][l], 69.9, Znl_int_accel[N][l]);
    }
    else if ( max_Ee < .001 )
        return 0;
    else
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_int_spline[N][l], max_Ee, Znl_int_accel[N][l]);
}

//returns total prob below max_Ee
double Target::invZ_nl_integrated(int N, int l, double Qe, double Znl)
{ 
    if(Znl>Z_nl_integrated( N, l, Qe, 70))
    {    cout << "Prob. too high result may be skewed" << endl; return 70; }
    return gsl_spline_eval(invZnl_int_spline[N][l], Znl/pow((Qe / 0.001),2), invZnl_int_accel[N][l]);
}

//returns maximum diff probability of ionization for a given level
double Target::Z_nl_max(int N, int l, double Qe)
{
    return pow((Qe / 0.001),2) * Znl_y_max[N][l];
}


int Target::init()
{

    std::fstream RFF;
    //Set target element
    if (elementZ == XENON)
    {
        Mt = MtXe;
        isoN = isoNXe;
        isoTable = isoTableXe;
        isoFracTable = isoFracTableXe;
        nl_energy = nl_energyXe;
        Nlmax = NlmaxXe;
        Leff = LeffXe;
        RFF.open("src/Xe_new.dat");
    }
    else if (elementZ == ARGON)
    {
        Mt = MtAr;
        isoN = isoNAr;
        isoTable = isoTableAr;
        isoFracTable = isoFracTableAr;
        nl_energy = nl_energyAr;
        Nlmax = NlmaxAr;
        Leff = LeffAr;
        RFF.open("src/Ar.dat");
    }
    else
    {
        cout << "target not found\n";
        return -1;
    }

    if(!RFF.is_open())
    {
        cout << "error opening migdal data file" << endl;
        return -1;
    }
    
    double*** EMnl = new double**[6];
    double*** Znl = new double**[6];
    
    for (int ni=1; ni < 6; ni++)
    {
        EMnl[ni] = new double*[3]; 
        Znl[ni] = new double*[3];
    }
    
    int n,l;
    string line;
    while (getline(RFF,line))
    {
        stringstream ss(line);
        ss >> n >> l;
        EMnl[n][l] = new double[251];
        Znl[n][l] = new double[251];
        Znl_y_max[n][l] = 0;
        int i = 0;
        while (i<251)
        {
            getline(RFF,line);
            stringstream ss(line);
            ss >> EMnl[n][l][i] >> Znl[n][l][i];
            EMnl[n][l][i]/=1000;   //convert units
            Znl[n][l][i]*=1000/2/M_PI;
            if ( Znl[n][l][i] > Znl_y_max[n][l] )
                Znl_y_max[n][l] = Znl[n][l][i];
            i++;         
        }
        Znl_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
        Znl_accel[n][l] = gsl_interp_accel_alloc();
        gsl_spline_init(Znl_spline[n][l], EMnl[n][l], Znl[n][l], 251);
        Znl_x_min[n][l] = EMnl[n][l][0];
        Znl_x_max[n][l] = EMnl[n][l][250];
        
    }
    
    //create integrated probabilities to a maximum Ee
    double Znl_int[251];
    for(n=1; n<6; n++)
    {
        for(l=0; l<Nlmax[n]; l++)
        {
            for(int j=0; j<251; j++)
            {
                Znl_int[j] = gsl_spline_eval_integ(Znl_spline[n][l], Znl_x_min[n][l], EMnl[n][l][j], Znl_accel[n][l]);
            }
            Znl_int_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
            Znl_int_accel[n][l] = gsl_interp_accel_alloc();
            gsl_spline_init(Znl_int_spline[n][l], EMnl[n][l], Znl_int, 251);
            
            invZnl_int_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
            invZnl_int_accel[n][l] = gsl_interp_accel_alloc();
            gsl_spline_init(invZnl_int_spline[n][l], Znl_int, EMnl[n][l], 251);
        }
    }
    
    if(RFF.bad())
        perror("error while reading file ");
    RFF.close();
    
    double maxENR = 0;
    maxENR = ERmax(maxE,0,isoTable[0]*AMU);
   
    //find maximum probability and normalize to it
    if(migOptimize==1)
    {
        cout << "optimizing migdal probabilities..";
        double maxProb = 0;
        double nr = .001;
        while (totalMigProb(maxE, nr, AMU*isoTable[0]) > maxProb && nr < maxENR)
        {
            maxProb = totalMigProb(maxE, nr, AMU*isoTable[0]);
            nr*=1.001;
        }
        maxCumulativeProb = maxProb;
        cout << " Migdal probabilities *~= " << 1/maxCumulativeProb << endl;
    }
    
    return 0;
}


//calculate the total migdal probability for a given recoil
double Target::totalMigProb(double maxE, double ENR, double MT)
{
    double maxEM=0, maxEe=0;
    maxEM  = EMmaxER(maxE, ENR);

    double prob=0;
    for(int N=1; N<6; N++)
    {
        for (int L=0; L<Nlmax[N]; L++)
        {
            maxEe = maxEM - nl_energy[N][L];
            if (maxEe < 0)
                continue;
            prob += Z_nl_integrated(N,L,qe(ENR,MT),maxEe);
        }
    }
    return prob;
}

//function for the ionization probabilities
double Target::Z_nl(int N, int l, double Qe, double Ee)
{
    if ( Ee > 70 )
    {
        return 0;
    }
    else if ( Ee < .001 )
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_spline[N][l], 0.001, Znl_accel[N][l]);
    else
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_spline[N][l], Ee, Znl_accel[N][l]);
}


//eject a random electron? (proportional to prob) - main Migdal MC function
vector<double> Target::rand_migdalE(double ENR)
{
    //select random isotope
    double MT = -1;
    if (isoN>1)
        MT = AMU*rand_isotope(ENR);
    else
        MT = Mt;
        
    double maxEM=0;
    maxEM = EMmaxER( maxE, ENR, MT);

    double Ee,EeMax;

    double t,prob=0;
    double maxProb=totalMigProb(maxE, ENR, MT)/maxCumulativeProb;

    t = gsl_rng_uniform(r);
    if ( t < maxProb )
    {
        for(int N=1; N<6; N++)
        {
            for (int L=0; L<Nlmax[N]; L++)
            {
                EeMax = maxEM - nl_energy[N][L];
                if (EeMax < 0)
                    continue;
                prob+=Z_nl_integrated(N,L,qe(ENR,MT),EeMax)/maxCumulativeProb;
                if( t < prob )
                {
                    t = gsl_rng_uniform(r)*Z_nl_integrated(N,L,qe(ENR,MT),EeMax);
                    Ee = invZ_nl_integrated(N, L, qe(ENR,MT), t);
                    return {Ee,nl_energy[N][L],(double)N,MT/AMU};
                }
            }
        }
    }
    return {0,0,0,0};

}


//returns a randomly selected isotope consistent with ENR
double Target::rand_isotope(double ENR)
{
    isoLoop:
        double u = gsl_rng_uniform(r);
        double maxER=0;
        for (int i=0;i<isoN;i++)
        {   
            if (u < isoFracTable[i])
            {
                maxER = ERmax(maxE, 0, AMU*isoTable[i]);
                if( ENR < maxER ) //ensures kinematic consistency
                    return isoTable[i];
                else
                    goto isoLoop;
            }
        }
        return Mt; //in case something goes wrong
}

