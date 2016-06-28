
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#define  NMD 10000
using namespace std;

#include "PrimChemBurner.H"

//#include "PrimChemGlobals.H"
//using namespace PrimChemNS;


// From other files
//void PrimChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* ymass,REAL tol, void (*f)(REAL,REAL,REAL*,REAL*,REAL*,bool,REAL*));

void PrimChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass,int* oOutputs, void (*f)(REAL,REAL*,REAL* ,REAL* ),REAL* RadParams,REAL* aion,REAL* RatRaw,int &jcounts);


void PChem_Burner(REAL dt, REAL* eosIn, REAL* xIn, REAL* xOut,REAL* aion,REAL* zion,REAL* bion,REAL* gamma, int chem_cc_case, REAL* RadParams, REAL dl, REAL divv, REAL *temptab, REAL* rate1)
{
    // Integration/Debug stuff
    int i, j;
 
    double avo = 6.0221367e23, ev2erg = 1.602e-12;
    double conv = ev2erg*1.0e0*avo;
    //double tol = 1.0e-4, odescal = 1.0e-6;
    
    // Cooling?
    int Do_cool = 1;
    int Do_cool_metal = 0;
    
    REAL tmin  = dt*1.e-20, ttot = dt;
    
    // Writing to separate file
    //bool WriteOut = 1;
    //char fOut[] = "Orion2ChemOut.txt";
    //ofstream myFile;
    //myFile.open(fOut);
    
    // Define xmass and ymass arrays ready
    REAL xmass[NSPECIES],ymass[NSPECIES], ysave[NSPECIES];

    for(i=0;i<NSPECIES;i++)
    {
        xmass[i] = xIn[i];
        ymass[i] = xmass[i]/aion[i] * 1.3158;     //(1.0 + 4.0 * ABHE) = 1.3158
        if(ymass[i]< smallx ) ymass[i] = smallx;
        ysave[i] = ymass[i];
    }
   
    //for(i=0;i<NSPECIES;i++) printf("1 i = %d, ymass = %lg\n", i, ymass[i]);
 
    //If writing, save to file
    /*
    if(WriteOut)
    {
        myFile << 1.0 << " ";
        for(i=0;i<NSPECIES;i++) myFile << ymass[i]*aion[i] << " " ;
        myFile << eosIn[1] << endl;
    }
    */
    
    // Prepare the network rates (inputs: temp, density, array)
    //REAL RatRaw[NREACTION];
    double *RatRaw = (double*)malloc(sizeof(double) * NREACTION);

    //ARS TEST -- use constant temperature at high density to test *only* chemistry network
//    if(eosIn[0] < 1.e7/1.e-24)
      GetNetworkRates(eosIn[1],eosIn[0],ymass,aion,RatRaw,chem_cc_case,RadParams);
//    else
//      GetNetworkRates(1000.0,eosIn[0],ymass,aion,RatRaw,chem_cc_case,RadParams);

    //for(i=0;i<NSPECIES;i++) printf("2 i = %d, ymass = %lg\n", i, ymass[i]);
    //printf("1a dt = %lg\n", dt);
    
    // Get the values of derivatives for the ODE solver
    REAL ddyddx[NSPECIES];
    GetNetwork(dt,ddyddx,RatRaw,ymass);


    //for(i=0;i<NSPECIES;i++) printf("1b i = %d, y = %lg dydx = %lg RatRaw = %lg\n", i, ymass[i], ddyddx[i], RatRaw[2*i]);
     //printf("1c dt = %lg\n", dt);

    // Perhaps some other outputs we want to hang onto
    int oOutputs[5];
    oOutputs[0] = chem_cc_case;
    
    
    // Given these, we can estimate the first timestep to try
    double firstDT = dt;
    double tempDT = 0.0;
    double eps = 0.1;
    double ei_old, T_old, Tdot, tempK_min = 50.0;
    int isSmall = -1;
    for(i=0;i<NSPECIES-1;i++)
    {
        tempDT = fabs(eps*(ymass[i]+0.01*ymass[iHP])/ddyddx[i]);
        //cout << "tempDT[" << i << "] = " << tempDT << endl;
        if(tempDT < firstDT) if(ymass[i] > 1.e2*smallx)
        {
            firstDT = tempDT;
            isSmall = i;
        }
    }
  
    //////ARS adding estimate of cooling rate and timescale
    T_old = eosIn[1]; ei_old = eosIn[2];
    Cool_function(eosIn, ymass, aion, gamma, RatRaw, dl, divv, temptab, rate1, chem_cc_case, Do_cool_metal,firstDT);
    Tdot = fabs(T_old - eosIn[1]) / firstDT;
    tempDT = T_old / Tdot;
    if(tempDT < firstDT) firstDT = tempDT; 
    eosIn[1] = T_old; eosIn[2] = ei_old;
    ////
     
    int counts = 0, jcounts = 0;
    double advT = 0.0;
    double endT;
    while(dt>0) // Main integration loop
    {
        counts = counts + 1;
        
        // We have now set the time to integrate before updating photochem
        
        tmin = firstDT*1.0e-20;
        endT = firstDT;
        
        //for(i=0;i<NSPECIES;i++) printf("4a i = %d, ymass = %lg\n", i, ymass[i]);
        //printf("2 dt = %lg, firstDT = %lg\n", dt, firstDT);

        //if(eosIn[0] > 4.4e-16)
        if(jcounts < 1 && firstDT < 1.e0) 
          for(i=0;i<6;i++)
            printf("0 PrimChem NEG i = %d, eos = %lg, firstDT = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, eosIn[i], firstDT, dt, endT, advT, jcounts);

        //if(eosIn[0] > 4.4e-16)
        if(jcounts < 1 && firstDT < 1.e0)
          for(i=0;i<NSPECIES;i++)
            printf("0 PrimChem NEG i = %d, ymass = %lg, ddyddx = %lg firstDT = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, ymass[i], ddyddx[i], firstDT, dt, endT, advT, jcounts);

        ymass[iELEC] = ymass[iHP] + ymass[iDP] + ymass[iHEP] + 2.0*ymass[iHEPP] + ymass[iHDP] + ymass[iH2P] + ymass[iD2P] - ymass[iHM] - ymass[iDM];

        //ymass[iHM] = 1.e-7 * ymass[iHP];
        //ymass[iH2P] =  1.e-9 * ymass[iHP];
        ymass[iH] = 1.e0 - 2.e0 * ymass[iH2]  - ymass[iHP]  - ymass[iHD] - 2.*ymass[iH2P] - ymass[iHM];
        ymass[iD]  = abundD - ymass[iDP] - ymass[iHD];
        ymass[iHE] = abhe - ymass[iHEP] + ymass[iHEPP];
        ymass[iDM] = ymass[iHDP] = ymass[iD2P] = ymass[iD2] = smallx; 

        PrimChemIntegrate(0,firstDT,tmin,endT,eosIn,ymass,oOutputs, GetNetwork,RadParams,aion,RatRaw, jcounts);
        
        // Charge neutrality
        ymass[iELEC] = ymass[iHP] + ymass[iDP] + ymass[iHEP] + 2.0*ymass[iHEPP] + ymass[iHDP] + ymass[iH2P] + ymass[iD2P] - ymass[iHM] - ymass[iDM];
    
        //more neutrality -- ARS test only TO BE REMOVED LATER
        //ymass[iHM] = 1.e-7 * ymass[iHP];
        //ymass[iH2P] =  1.e-9 * ymass[iHP];
        ymass[iH] = 1.e0 - 2.e0 * ymass[iH2]  - ymass[iHP]  - ymass[iHD] - 2.*ymass[iH2P] - ymass[iHM]; 
        ymass[iD]  = abundD - ymass[iDP] - ymass[iHD];
        ymass[iHE] = abhe - ymass[iHEP] + ymass[iHEPP];
        ymass[iDM] = ymass[iHDP] = ymass[iD2P] = ymass[iD2] = smallx; 
 
        // Prevent negative values
        for(i=0;i<NSPECIES;i++) if(ymass[i] < smallx ) ymass[i] = smallx;
        
        //for(i=0;i<NSPECIES;i++) printf("4b i = %d, ymass = %lg\n", i, ymass[i]);        

        // Determines how much ionization occured from collisions vs. photoionizations
        double sdotRate=0.0;
        double photoc[NSPECIES] = {0.0};
        Chemistry_photo(ymass,RatRaw,photoc);
        for(i=0;i<NSPECIES;i++) sdotRate = sdotRate + (ymass[i]-ysave[i])*bion[i]*photoc[i];
       
        //for(i=0;i<NSPECIES;i++) printf("4c i = %d, ymass = %lg\n", i, ymass[i]);
 
        // Update energy density from this photochem / Update chemEOS from there
        eosIn[2] = eosIn[2] + sdotRate*conv;
        if(eosIn[2]<0) cout << "ERROR: energy density is negative!!!" << endl;
        chemEOS(ymass,eosIn,aion,gamma,Do_cool);
       
        //for(i=0;i<NSPECIES;i++) printf("4d i = %d, ymass = %lg\n", i, ymass[i]);
 
        //if(eosIn[0] > 4.4e-16)
        if(jcounts < 1 && firstDT < 1.e0)
          for(i=0;i<6;i++)
            printf("1 PrimChem NEG i = %d, eos = %lg, firstDT = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, eosIn[i], firstDT, dt, endT, advT, jcounts);

        // Cooling
        if(Do_cool)
        {
            
            if(eosIn[1] < tempK_min)
              {
              eosIn[1] = tempK_min;
              Chem_EOSB(1,eosIn,xIn,aion,gamma,NSPECIES);
              //printf("We have hit the temperature floor! nh = %lg\n", eosIn[0]/1.67e-24);
              }

            ei_old = eosIn[2]; 

             
            //cout << "Entering cooling " << eosIn[2] << endl;
            Cool_function(eosIn, ymass,aion,gamma, RatRaw, dl, divv, temptab, rate1, chem_cc_case, Do_cool_metal,endT); //updates e, T
            chemEOS(ymass,eosIn,aion,gamma,Do_cool);

            if(eosIn[2] < 0.5*ei_old)
              {
              eosIn[2] = 0.5*ei_old;
              Chem_EOSB(2,eosIn,xIn,aion,gamma,NSPECIES);
              eps = eps*0.9;
              printf("We are cooling too fast! eps =%lg, nh = %lg\n", eps, eosIn[0]/1.67e-24);
              }

            //ARS applyting temperature and energy floor
            if(eosIn[1] < tempK_min)
              {
              eosIn[1] = tempK_min;
              Chem_EOSB(1,eosIn,xIn,aion,gamma,NSPECIES);
              //printf("We have hit the temperature floor! nh = %lg\n", eosIn[0]/1.67e-24);
              }
 
        }

    
        // Re-evaluate rates and network using new state values and species abundances
   //     if(eosIn[0] < 1.e7/1.e-24)
          GetNetworkRates(eosIn[1],eosIn[0],ymass,aion,RatRaw,chem_cc_case,RadParams);
   //     else
   //       GetNetworkRates(1000.0,eosIn[0],ymass,aion,RatRaw,chem_cc_case,RadParams);

        GetNetwork(0.0,ddyddx,RatRaw,ymass); // first entry time, not actually used
       
        //for(i=0;i<NSPECIES;i++) printf("5 i = %d, ymass = %lg\n", i, ymass[i]);
 
        // We have now successfully advanced 'endT' worth of time
        // Update time left
        advT = advT + endT;
        dt = dt - endT;
        firstDT = dt;
        
        // Re-evaluate firstDT
        for(i=0;i<NSPECIES-1;i++)
        {
            tempDT = fabs(eps*(ymass[i]+0.01*ymass[iHP])/ddyddx[i]);
            //cout << "tempDT[" << i << "] = " << tempDT << endl;
            if(tempDT < firstDT) if(ymass[i] > 1.e2*smallx)
            {
                firstDT = tempDT;
                isSmall = i;
            }
        }
        
        
        ////ARS adding estimate of cooling rate and timescale
        T_old = eosIn[1]; ei_old = eosIn[2];
        Cool_function(eosIn, ymass,aion,gamma, RatRaw, dl, divv, temptab, rate1, chem_cc_case, Do_cool_metal, firstDT);
        Tdot = fabs(T_old - eosIn[1]) / firstDT;
        tempDT = T_old / Tdot;
        if(tempDT < firstDT) firstDT = tempDT;
        eosIn[1] = T_old; eosIn[2] = ei_old;

        
        if(firstDT > dt) firstDT = dt; // should never happen, but why not check.
        
        // Save for photoc derivative
        for(i=0;i<NSPECIES;i++) ysave[i] = ymass[i];
      
        //if(eosIn[0] > 4.4e-16) 
        if(jcounts < 1 && firstDT < 1.e0)
          for(i=0;i<6;i++)
            printf("2 PrimChem NEG i = %d, eos = %lg, firstDT = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, eosIn[i], firstDT, dt, endT, advT, jcounts);

         if(jcounts%10000 == 0)
          printf("3 PrimChem i = %d, Temp = %lg, firstDT = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, eosIn[1], firstDT, dt, endT, advT, jcounts);
        
        //If writing, save to file
        /*
        if(WriteOut)
        {
        myFile << advT << " ";
        for(i=0;i<NSPECIES;i++) myFile << ymass[i]*aion[i] << " " ;
        myFile << eosIn[1] << endl;
        }
        */
        //free(photoc);
    }

    free(RatRaw);
    
    // All done!
    for(i=0;i<NSPECIES;i++)
    {
        xOut[i] = ymass[i]*aion[i] / 1.3158;
        if(xOut[i] < smallx) xOut[i] = smallx;
    }

    //if(jcounts > 1000)
      //{
      printf("Took a little while!  jcounts = %d, nh = %lg, dt = %lg, advT = %lg, ttot = %lg\n", jcounts, eosIn[0]/1.67e-24, dt, advT, ttot);
      for(i=0;i<6;i++) printf("4 i = %d, eos = %20.15g\n", i, eosIn[i]);
      for(i=0;i<NSPECIES;i++) printf("4 i = %d, xIn = %lg, xOut = %20.15g\n", i, xIn[i], xOut[i]);
      //}

    for(i=0;i<NSPECIES;i++)
      if(xOut[i] != xOut[i]) 
        printf("3 PrimChem NEG i = %d, xout = %lg, ymass = %lg, ddyddx = %lg, aion = %lg\n", i, xOut[i], ymass[i], ddyddx[i], aion[i]);
    //for(i=0;i<6;i++)
      //if(eosIn[i] != eosIn[i])
        //printf("3 PrimChem NEG i = %d, eos = %lg, dt = %lg, endT = %lg, advT = %lg, jcounts = %d\n", i, eosIn[i], dt, endT, advT, jcounts);
 
    //cout << "End of Burner: counts = " << counts << " , jcounts = " << jcounts << endl;
    //myFile.close();
};






