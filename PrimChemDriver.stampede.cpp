
#include <stdio.h>
#include <iostream>
#include <iomanip>
using namespace std;

#include "PrimChemDriver.H"
using namespace PrimChemNS;
#include <stdlib.h>
#include <math.h>

#include "mol_data.h"
#include "spline.cpp"
#define  NH2DATA 41
#define  NMD 10000

void PChem_Burner(REAL dt, REAL* eosIn, REAL* xIn, REAL* xOut,REAL* aion,REAL* zion,REAL* bion,REAL* gamma,int chem_cc_case, REAL* RadParams, REAL dl, REAL divv, REAL* temptab, REAL *rate1);

//void PrimChemDriver(double dt_in, double &den_in, double &ei_in, double* x_fromOrion);

//void o2ChemDrive(REAL dt, REAL den, REAL Temp)

void PrimChemDriver(double dt_in, double den_in, double &ei_in, double gamma_here, double* x_fromOrion, double dl, double divv)
{
    int i, use_orion=1;
    double random, rfac = 1.e-6;
    random = rand();

    //printf("Hello we are in PrimChemDriver\n");
        
    // Initialize arrays (xIN will be read in, actually)
    REAL xIn[NSPECIES] = {0.0};
    REAL xOut[NSPECIES] = {0.0};
    REAL eosIn[6];  // Density, Temperature, EngDen, Pressure, Metal Fraction, gamma

    //Initialization parameters
    //and other parameters that might become read in
    int chem_cc_case = 1;
    REAL fssh2 = 1.0;
    REAL fsshd = 1.0;
    REAL J21 = 0.0;  // J21 radiation
    REAL RadParams[3] = {J21,fssh2,fsshd};
    
    // Debugging Routine
    SetUp(xIn,eosIn);
    //REAL dt = 1e15;
    
    // Load from Orion2
    if(use_orion == 1 || use_orion == 2) for(i=0;i<NSPECIES;i++) xIn[i] = x_fromOrion[i];

    if(random < rfac*RAND_MAX)
      for(i=0;i<NSPECIES;i++) printf("i = %d, xIn = %lg, x_fromOrion = %lg\n", i, xIn[i], x_fromOrion[i]);

    REAL dt;
    dt = 1.e7;
    if(use_orion == 1 || use_orion == 2) dt = dt_in;

    if(random < rfac*RAND_MAX) printf("dt = %lg, dt_in = %lg\n", dt, dt_in);
    if(dt < 1.e1) dt = 1.e1;

    if(use_orion == 1 || use_orion == 2) 
       {
       if(use_orion == 1) eosIn[0] = den_in;
       eosIn[1] = 0.0; //Temp
       eosIn[2] = ei_in;
       eosIn[3] = 0.0; //P
       eosIn[4] = 0.0; // mass frac
       eosIn[5] = gamma_here;
       } 

    // Initialize some code values (This was the contents of Chem_init.F90)
    // It be nice if these could become part of the namespace, but nothing compiles
    // correctly when I try... :/
    //int idx[NSPECIES];
    REAL aion[NSPECIES];
    REAL zion[NSPECIES];
    REAL bion[NSPECIES];
    REAL gamma[NSPECIES];
    initChem(aion,zion,bion,gamma);
  
    if(random < rfac*RAND_MAX) 
      for(i=0;i<NSPECIES;i++) printf("i = %d, aion = %lg, bion = %lg, gamma = %lg\n", i, aion[i], bion[i], gamma[i]);
 
    // Charge neutrality
    //xIn[iELEC] = (xIn[iHP]/aion[iHP] + xIn[iHEP]/aion[iHEP] + xIn[iDP]/aion[iDP]) * aion[iELEC];
    xIn[iELEC] = (xIn[iHP]/aion[iHP] + xIn[iHEP]/aion[iHEP] + xIn[iDP]/aion[iDP] 
                 + xIn[iHEPP]/aion[iHEPP] + xIn[iHDP]/aion[iHDP] 
                 + xIn[iH2P]/aion[iH2P] + xIn[iD2P]/aion[iD2P] 
                 - xIn[iHM]/aion[iHM] - xIn[iDM]/aion[iDM]) * aion[iELEC];   

    if(random < rfac*RAND_MAX) for(i=0;i<6;i++) 
      printf("1 i = %d, eos = %lg\n", i, eosIn[i]);
 
    // Load equation of state stuff (some read in, some calculated in EOS)
    if(use_orion == 1) 
      {
      Chem_EOS(2,eosIn,xIn,aion,gamma,NSPECIES);
      double tempK_min = 50.0;
      if(eosIn[1] < tempK_min)
        {  
        eosIn[1] = tempK_min;
        Chem_EOS(1,eosIn,xIn,aion,gamma,NSPECIES);
        printf("Driver, We have hit the temperature floor! nh = %lg, dt = %lg\n", eosIn[0]/1.67e-24, dt);
        for(i=0;i<NSPECIES;i++) printf("floor i = %d, xIn = %lg, xOut = %lg\n", i, xIn[i], xOut[i]);
        }
      }
    else
      {
      eosIn[2] = 1.0; //ei;
      eosIn[3] = 1.0; //P;
      Chem_EOS(1,eosIn,xIn,aion,gamma,NSPECIES);
      }

       if(random < rfac*RAND_MAX)
         {
         for(i=0;i<6;i++) printf("2 i = %d, eos = %lg\n", i, eosIn[i]);
         for(i=0;i<NSPECIES;i++) printf("2 i = %d, xIn = %lg, xOut = %lg\n", i, xIn[i], xOut[i]);
         }    

    double temptab[NMD], nmd, i_doub;  
    double rate1[NMD], tmin=1.e0, tmax=2.e8, dtlog; 
    nmd = NMD;   
    dtlog = log10(tmax) / (nmd - 1.);
    for (i=0; i<NMD; i++) {i_doub = (double) i; temptab[i] = pow(10.0, i_doub*dtlog); }
    
    double h2_temp[NH2DATA];
    for (i=0; i<NH2DATA; i++) {i_doub = (double) i; h2_temp[i] = pow(10.0, 2.e0 + 5.e-2*(i_doub)); }
    
    spline_eval(NH2DATA, h2_temp, h2_lte, NMD, temptab, rate1);

    // Call BURNER
    PChem_Burner(dt,eosIn,xIn,xOut,aion,zion,bion,gamma,chem_cc_case,RadParams,dl,divv,temptab,rate1);
    Chem_EOS(1,eosIn,xOut,aion,gamma,NSPECIES);
   
       if(random < rfac*RAND_MAX)
         {
         for(i=0;i<6;i++) printf("3 i = %d, eos = %20.15g\n", i, eosIn[i]);
         for(i=0;i<NSPECIES;i++) printf("3 i = %d, xIn = %lg, xOut = %20.15g\n", i, xIn[i], xOut[i]);
         for (i=0; i<NH2DATA; i++) printf("i = %d, h2_temp = %lg\n", i, h2_temp[i]);
         for (i=0; i<NH2DATA; i++) printf("i = %d, h2_lte = %lg\n", i, h2_lte[i]);
         for (i=0; i<NMD; i++) printf("i = %d, temptab = %lg\n", i, temptab[i]);
         for (i=0; i<NMD; i++) printf("i = %d, rate1 = %lg\n", i, rate1[i]);
         } 

    // We're all done, let's return the results and go home.
    
    // eosIN and xOut must then be returned
    //for(i=0;i<NSPECIES;i++) cout << i+1 << ": " << xOut[i] << endl;
    //cout << setprecision(15);
    //for(i=0;i<NSPECIES;i++) cout << xOut[i] << endl;
    //cout << endl;
    //for(i=0;i<5;i++) cout << eosIn[i] << endl;
    //for(i=0;i<5;i++) cout << "eosIn[" << i << "]: " << eosIn[i] << endl;
    
    for(i=0;i<NSPECIES;i++) x_fromOrion[i] = xOut[i];
    den_in = eosIn[0];
    ei_in = eosIn[2];    

}
