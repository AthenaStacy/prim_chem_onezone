/* Primordial Chemistry Routine
 
 Based on Fortran code from Will Gray, re-written for C++ by Aaron Lee, 2014
 
 */



#include "PrimChemGlobals.H"

using namespace PrimChemNS;


void SetUp(REAL* xIn, REAL* eosIn)
{
    
    xIn[iHP] = 4.e-8; //0.0667; //1.1e-10; //                 ## H+
    xIn[iH] = 0.76;// 0.0667; //1.99e-10; //                 ## H
    xIn[iHM] = 5.e-16;// 0.0667; //1.1e-10; //                 ## H-
    xIn[iH2P] = 1.e-18;// 0.0667; //1.1e-10; //                 ## H_2+
    xIn[iH2] = 7.e-4;// 0.0667; //0.75178; //                 ## H_2
    xIn[iDP] = 1.e-13;// 0.0667; //1.1e-10; //                 ## D+
    xIn[iD] = 1.e-5;// 0.0667; //3.90962e-5; //              ## D
    xIn[iDM] = 2.e-20;// 0.0667; //1.1e-10; //                 ## D-
    xIn[iHDP] = 2.e-20;// 0.0667; //1.1e-10; //                 ## HD+
    xIn[iHD] = 3.e-8;// 0.0667; //1.1e-10; //                 ## HD
    xIn[iD2P] = 4.e-20;// 0.0667; //1.1e-10; //                 ## D_2+
    xIn[iD2] = 1.e-11;// 0.0667; //1.1e-10; //                 ## D_2
    xIn[iHEP] = 4.e-20;// 0.0667; //1.1e-10; //                 ## HE+
    xIn[iHE] = 0.24;// 0.0667; //0.24818; //                 ## HE
    xIn[iHEPP] = 4.e-20;// 0.0662; //1.1e-10; //                 ## HE++
    xIn[iELEC] = 1.-12;// 0.01; //1.1e-10; //                 ## elec
    
    eosIn[0] = 1.e-17; //Rho
    eosIn[1] = 1000.0; // Temp
    eosIn[4] = 0.0; // Metal Fraction

}


void initChem(double* aion,double* zion,double* bion,double* gamma)
//void initChem(int* idx,double* zion,double* bion,double* gamma)
{
    
    /* Index array (might not be needed)
    idx[iHP] = iHP;
    idx[iH] = iH;
    idx[iHM] = iHM;
    idx[iH2P] = iH2P;
    idx[iH2] = iH2;
    idx[iDP] = iDP;
    idx[iD] = iD;
    idx[iDM] = iDM;
    idx[iHDP] = iHDP;
    idx[iHD] = iHD;
    idx[iD2P] = iD2P;
    idx[iD2] = iD2;
    idx[iHEP] = iHEP;
    idx[iHE] = iHE;
    idx[iHEPP] = iHEPP;
    idx[iELEC] = iELEC;
    */

    //Total atomic mass
    //REAL aion[NSPECIES+1];
    aion[iELEC] = 0.000549;
    aion[iHP]   = 0.999451;
    aion[iH]    = 1.0;
    aion[iHM]  = 1.000549;
    aion[iH2P]  = 1.999451;
    aion[iH2]   = 2.0;
    aion[iDP]   = 1.999451;
    aion[iD]    = 2.0;
    aion[iDM ]    = 2.000549;
    aion[iHDP ]  = 2.999451;
    aion[iHD ]   = 3.0;
    aion[iD2P ]  = 3.999451;
    aion[iD2 ]   = 4.0;
    aion[iHEP ]  = 3.999451;
    aion[iHE ]   = 4.0;
    aion[iHEPP ] = 3.998902;
    
    //Protons
    zion[iELEC] = 0.0;
    zion[iHP]   = 1.0;
    zion[iH]    = 1.0;
    zion[iHM]   = 1.0;
    zion[iH2P]  = 2.0;
    zion[iH2]   = 2.0;
    zion[iDP]   = 1.0;
    zion[iD]    = 1.0;
    zion[iDM]   = 1.0;
    zion[iHDP]  = 2.0;
    zion[iHD]   = 2.0;
    zion[iD2P]  = 2.0;
    zion[iD2]   = 2.0;
    zion[iHEP]  = 2.0;
    zion[iHE]   = 2.0;
    zion[iHEPP] = 2.0;
    
    //Binding Energy (eV)
    bion[iELEC] = 0.00;
    bion[iHP]   = -13.59844;
    bion[iH]    = 0.0000000;
    bion[iHM]   = 0.7500000;
    bion[iH2P]  = -10.94710;
    bion[iH2]   = 4.4778800;
    bion[iDP]   = -13.59844;
    bion[iD]    = 0.0;
    bion[iDM]   = 0.75;
    bion[iHDP]  = -10.9471;
    bion[iHD]   = 4.47788;
    bion[iD2P]  = -10.9471;
    bion[iD2]   = 0.0;
    bion[iHEP]  = -24.58741;
    bion[iHE]   = 0.0;
    bion[iHEPP] = -79.00510;
    
    //Gamma
    gamma[iELEC] = 1.6667;
    gamma[iHP]   = 1.6667;
    gamma[iH]    = 1.6667;
    gamma[iHM]   = 1.6667;
    gamma[iH2P]  = 1.4;
    gamma[iH2]   = 1.4;
    gamma[iDP]   = 1.6667;
    gamma[iD]    = 1.6667;
    gamma[iDM]   = 1.6667;
    gamma[iHDP]  = 1.4;
    gamma[iHD]   = 1.4;
    gamma[iD2P]  = 1.4;
    gamma[iD2]   = 1.4;
    gamma[iHEP]  = 1.6667;
    gamma[iHE]   = 1.6667;
    gamma[iHEPP] = 1.6667;

}

const void eos_getabar(double& a,double* w,double* aion,int n)
{
    double val=0.0;
    
    for(int i=0;i<n;i++)
    {
        val = val + w[i]/aion[i];
    }
    
    a = val;
}

const void Chem_EOS(int mode, double* eosIn, double* x, double* aion,double* gam, int nSpec)
{
    double gasconst = 8.314472e7; // [erg/Kelvin]
    
    double den,temp,ei,pres;
    den=eosIn[0];
    temp=eosIn[1];
    ei=eosIn[2];
    pres=eosIn[3];
    
    double gamma1j[nSpec];
    for(int i=0;i<nSpec;i++)
    {
        gamma1j[i] = 1.0 / (gam[i]-1.0);
        //cout << "gam: " << gam[i] << " invgam: " << gamma1j[i] << endl;
    }
    
    
    double weight[nSpec];
    for(int i=0;i<nSpec;i++)
    {
        weight[i] = x[i]*gamma1j[i];
    }
    
    double abar=0.0;
    double abarinv=0.0;
    eos_getabar(abarinv,x,aion,nSpec);
    abar = 1.0/abarinv;
    
    double rt;
    eos_getabar(rt,weight,aion,nSpec);
    
    double gc = 1.0+1.0/(rt*abar);
    cout << "Effective gamma = " << gc << endl;
    eosIn[5] = gc;
    
    if(mode==1) // determines ei and presure from rho and T
    {
        eosIn[0]=den;
        eosIn[1]=temp;
        eosIn[2]=gasconst*temp/abar/(gc-1.0);
        eosIn[3]=gasconst*den*temp/abar;
    }
    else if(mode==2) //determines T and P from rho and ei
    {
        eosIn[0]=den;
        eosIn[1]=ei*(gc-1.0)*abar/gasconst;
        eosIn[2]=ei;
        eosIn[3]=den*ei*(gc-1.0);
    }
    else if(mode==3) // determines T and ei from rho and P
    {
        eosIn[0]=den;
        eosIn[2]=pres/den/(gc-1.0);
        eosIn[1]=eosIn[2]*(gc-1.0)*abar/gasconst;
        eosIn[3]=pres;
    }
    else // doesn't change anything... not good.
    {
        
        
    }
}


