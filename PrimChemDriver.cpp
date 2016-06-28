
#include <stdio.h>
#include <iostream>
#include <iomanip>
using namespace std;

#include "PrimChemDriver.H"
using namespace PrimChemNS;
//#include <stdlib.h>
#include <math.h>

#include "mol_data.h"
#include "spline.cpp"


void PChem_Burner(REAL dt, REAL* eosIn, REAL* xIn, REAL* xOut,REAL* aion,REAL* zion,REAL* bion,REAL* gamma,int chem_cc_case, REAL* RadParams, REAL dl, REAL divv, REAL* temptab, REAL *rate1);

FILE *outfile;

//void o2ChemDrive(REAL dt, REAL den, REAL Temp)
int main(int argc, char **argv)
{
    outfile=fopen("output.dat","w");
    int i, count;
    
    double dl = 1.e14; double divv = 0;
    
    // Initialize arrays (xIN will be read in, actually)
    REAL xIn[NSPECIES] = {0.0};
    REAL xOut[NSPECIES] = {0.0};
    REAL eosIn[6];  // Density, Temperature, EngDen, Pressure, Metal Fraction, gamma
    
    //Initialization parameters
    //and other parameters that might become read in
    int chem_cc_case = 1; // =0 Calls the optically thin cooling table (Osterbrock case A), where ionizing photons emitted during recombination are lost to the system. =1 calls the optically thick cooling table (Osterbrock case B), where the photons are locally absorbed, reducing the recombination rates.
    REAL fssh2 = 1.0; // Allow for H_2 dissociation from radiation (also requires j21>0)
    REAL fsshd = 1.0; // Allow for HD dissociation from radiation (also requires j21>0)
    REAL J21 = 0;//1.0;  // J21 radiation
    REAL RadParams[3] = {J21,fssh2,fsshd};
    
    // Debugging Routine
    SetUp(xIn,eosIn);
    //REAL dt = 5.25e7;
    REAL dt = 2.e9;
    
    
    
    // Initialize some code values (This was the contents of Chem_init.F90)
    // It be nice if these could become part of the namespace, but nothing compiles
    // correctly when I try... :/
    //int idx[NSPECIES];
    REAL aion[NSPECIES];
    REAL zion[NSPECIES];
    REAL bion[NSPECIES];
    REAL gamma[NSPECIES];
    initChem(aion,zion,bion,gamma);
    
    
    double *temptab = (double*)malloc(sizeof(double) * (NMD+1));
    double nmd, i_doub;
    double tmin=1.e0, tmax=2.e8, dtlog;
    
    nmd = NMD;
    dtlog = log10(tmax) / (nmd - 1.);
    for (i=0; i<NMD; i++) {i_doub = (double) i; temptab[i] = pow(10.0, i_doub*dtlog); }
    
    double h2_temp[NH2DATA];
    for (i=0; i<NH2DATA; i++) {i_doub = (double) i; h2_temp[i] = pow(10.0, 2.e0 + 5.e-2*(i_doub)); }
    

    int tcount=0;
    double epst = 1.e-2, t_curr = 0.0, tstep = 1.e5, tff, tcool, dt1, dt2, gam, e, x;
    double GN = 6.67e-8, PROTONMASS = 1.6726e-24, PI = 3.14159;
    double den, xn=1.0, xnold, temp, temp_old, Tdot, ndot=1.e-20, nchange;
    double mu, mu_inv, h2frac, hefrac, hfrac, muh2, muh2in, dTdt_ad = 1.e-20;
    
    //while(t_curr < dt)
    while(xn < 1.e12)
    {
        // Charge neutrality
        //xIn[iELEC] = (xIn[iHP]/aion[iHP] + xIn[iHEP]/aion[iHEP] + xIn[iDP]/aion[iDP]) * aion[iELEC];
        xIn[iELEC] = (xIn[iHP]/aion[iHP] + xIn[iHEP]/aion[iHEP] + xIn[iDP]/aion[iDP]
                      + xIn[iHEPP]/aion[iHEPP] + xIn[iHDP]/aion[iHDP]
                      + xIn[iH2P]/aion[iH2P] + xIn[iD2P]/aion[iD2P]
                      - xIn[iHM]/aion[iHM] - xIn[iDM]/aion[iDM]) * aion[iELEC] / 1.3158;
        
        // Load equation of state stuff (some read in, some calculated in EOS)
        den = eosIn[0];
        eosIn[2] = 0.0; //ei;
        eosIn[3] = 0.0; //P;
        Chem_EOS(1,eosIn,xIn,aion,gamma,NSPECIES);
 
    
        mu=1.2195;
        h2frac=xIn[iH2];
        hefrac = xIn[iHE];
        hfrac = 1. - h2frac - hefrac;
        
        mu_inv = hfrac + h2frac/2 + hefrac/4;
        mu = 1./mu_inv;

        xn = xnold = eosIn[0] / (mu * PROTONMASS);
        temp = temp_old = eosIn[1];
        
        e=2.7182818;
        x=6100.0/temp;
        //gam=1.0 + (1.0/mu)*pow((1.5*0.24/4.0)+(1.5*0.76*(1.0-h2frac))
        //      + ((2.5+pow(x,2.0)*pow(e,x)/pow((pow(e,x)-1.0),2.0))*(.76/2.0)*h2frac), -1.e0);
        
        //gam = eosIn[5];
        
        // Call BURNER
        PChem_Burner(tstep,eosIn,xIn,xOut,aion,zion,bion,gamma,chem_cc_case,RadParams,dl,divv,temptab,rate1);
        Chem_EOS(1,eosIn,xOut,aion,gamma,NSPECIES);
    
        // We're all done, let's return the results and go home.
    
        // eosIN and xOut must then be returned
        //for(i=0;i<NSPECIES;i++) cout << i+1 << ": " << xOut[i] << endl;
        cout << setprecision(15);
    
    
        //for(i=0;i<NSPECIES;i++) cout << "i = " << i << " , " << xOut[i] << endl;
        cout << endl;
        for(i=0;i<6;i++) cout << eosIn[i] << endl;
  
        
        tff = sqrt(3.e0*PI/(32.e0*GN*den));
        nchange = xn/tff;
        xn = xn + (nchange*tstep);
        ndot=(xn-xnold)/tstep;
        
        tcool = temp/Tdot;
        dTdt_ad = 0.01*(eosIn[5] - 1)*temp*ndot/xnold;
        temp = temp + dTdt_ad*tstep;
        Tdot = (temp - temp_old)/tstep;
        
        dt1=epst*temp/fabs(Tdot);
        dt2=epst*xn/fabs(ndot);
        
         
        //tstep = 1.e7;
        tstep = min(dt1, dt2);
        
        for(i=0;i<NSPECIES;i++) xIn[i] = xOut[i];
        eosIn[0] = xn * mu * PROTONMASS;
        eosIn[1] = temp;

        double dummy=0.0;
        
        printf("tcount = %d, xn = %lg, mu = %lg, tff = %lg, tcool = %lg, tstep = %lg, t_curr = %lg \n", tcount, xn, mu, tff, tcool, tstep, t_curr);
        t_curr = t_curr + tstep;
        tcount++;
        fprintf(outfile, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", xn, xOut[iH2], xOut[iHD], xOut[iDP], xOut[iHEP], xOut[iHP], dummy, dummy, dummy, dummy, temp, dummy, dummy, dummy, dummy, dummy);
    }
    
    free(temptab);
    //free(rate1);
    
    fclose(outfile);
    return 0;
    
}
