#include <iostream>
#include <stdlib.h>
#include <iomanip>
using namespace std;
#include <math.h>


#include "PrimChemGlobals.H"
using namespace PrimChemNS;


//#include "Jacobian.H"
void Jacobian(REAL x,REAL* ys,REAL (*dfdy)[NSPECIES], REAL* ratraw);






void leqs(REAL (*a)[NCALC],REAL* b)
{
    int n1 = NCALC-1;
    int i,j,k,l;
    
    REAL c,r;
    
 
    
    
    for(i=0;i<NCALC;i++)
    {
        r = fabs(a[i][0]);
        for(j=1;j<NCALC;j++)
        {
            c=fabs(a[i][j]);
            if(r<c) r=c;
        }
        
        for(j=0;j<NCALC;j++) a[i][j] = a[i][j]/r;
        b[i] = b[i]/r;
    }
    
    
    
    for(j=0;j<n1;j++)
    {
        l = j+1;
        for(i=l;i<NCALC;i++)
        {
            r = -a[i][j];
            if(r==0.0) continue; //goto onward;
            r=r/(a[j][j]);
            for(k=l;k<NCALC;k++) a[i][k] = a[i][k] + r*a[j][k];
            b[i] = b[i] + r*b[j];
        //onward:
            continue;
        }
    }
    
    
    b[n1] = b[n1] / a[n1][n1]; //last entry
    
    for(l=0;l<n1;l++)
    {
        i= NCALC-2-l;
        r=0;
        int imax = i+1;
        for(j=imax;j<NCALC;j++)
        {
            int jj=i+NCALC-j;
            r = r+(a[i][jj])*b[jj];
        }
        b[i]= (b[i]-r)/(a[i][i]);
    }
  
    //cout << "LEQS start" << endl;
    //for(i=0;i<NSPECIES;i++) cout << setprecision(15) << b[i] << endl;
    //cout << "LEQS end" << endl;
};


void StepMe(REAL* y,REAL* dydx,REAL &x,REAL htry,REAL &hdid,REAL &hnext,REAL* yscal, void (*derivs)(REAL,REAL*,REAL*,REAL*), int &jcounts, REAL* ratraw)
{

    
    int i,j,jtry;
    
    int TryMax = 30000000;
    
    
    REAL a21 = 2.0e0, a31 = 48.0e0/25.0e0, a32 = 6.0e0/25.0e0, c21 = -8.0e0, c31 = 372.0e0/25.0e0, c32 = 12.e0/5.0e0, c41 = -112.0e0/125.0e0,   c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,   c1x = 1.0e0/2.0e0, c2x = -3.0e0/2.0e0,  c3x = 121.0e0/50.0e0, c4x = 29.0e0/250.0e0, a2x = 1.0e0, a3x = 3.0e0/5.0e0;
    
    REAL eps = 1.0e-4;
    
    
    // Original values (needed if step fails)
    REAL xsav,ysav[NSPECIES], dysav[NSPECIES];
    for(i=0;i<NSPECIES;i++) ysav[i] = y[i];
    for(i=0;i<NSPECIES;i++) dysav[i] = dydx[i];
    xsav = x;
    
    
    // Get the dense Jacobian
    REAL dfdy[NSPECIES][NSPECIES];
    Jacobian(xsav,ysav,dfdy,ratraw);
    REAL h = htry;
    
    
    // Needed for matrix solve
    REAL dmat[NCALC][NCALC];
    REAL av[NCALC][NCALC];
    REAL bb[NCALC];
    
    
    REAL gamf = 0.5;
    // Solve equations loop
    for(jtry=1;jtry<TryMax;jtry++)
    {
        
        
        // Form the matrix
        REAL xx = 1.0/h/gamf;
        
        for(j=0;j<NCALC;j++)
            for(i=0;i<NCALC;i++)
                dmat[i][j] = -dfdy[i][j];
        for(i=0;i<NCALC;i++) dmat[i][i] = xx + dmat[i][i];
        
        
        // Set up and solve RHS for g1
        REAL g1[NCALC];
        for(i=0;i<NCALC;i++) g1[i] = dysav[i];
        for(j=0;j<NCALC;j++)
            for(i=0;i<NCALC;i++)
                av[i][j] = dmat[i][j];
        for(i=0;i<NCALC;i++) bb[i] = g1[i];
        
        
        leqs(av,bb);
        for(j=0;j<NCALC;j++) g1[j] = bb[j];
        if(iELEC < NCALC) g1[iELEC] = g1[iHP] + g1[iDP] + 2.0*g1[iHEPP] + g1[iHEP] + g1[iHDP] + g1[iH2P] + g1[iD2P] - g1[iHM] - g1[iDM];
        
        for(i=0;i<NCALC;i++) y[i] = ysav[i] + a21 * g1[i];
        if(iELEC < NCALC)                   
          {
          //g1[iH] = 1.e0 - 2.e0 * g1[iH2]  - g1[iHP]  - g1[iHD] - 2.*g1[iH2P] - g1[iHM];
          //g1[iD]  = abundD - g1[iDP] - g1[iHD];
          //g1[iHE] = abhe - g1[iHEP] + g1[iHEPP];
          //g1[iDM] = g1[iHDP] = g1[iD2P] = g1[iD2] = smallx;
          }

        
        // Set up and solve RHS for g2
        x = xsav + a2x * h ;
        REAL g2[NCALC];
        derivs(x,dydx,ratraw,y);
        for(i=0;i<NCALC;i++) g2[i] = dydx[i] + c21*g1[i]/h;
        for(j=0;j<NCALC;j++)
            for(i=0;i<NCALC;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NCALC;i++) bb[i] = g2[i];
        leqs(av,bb);
        for(j=0;j<NCALC;j++) g2[j] = bb[j];
        if(iELEC < NCALC) g2[iELEC] = g2[iHP] + g2[iDP] + 2.0*g2[iHEPP] + g2[iHEP] + g2[iHDP] + g2[iH2P] + g2[iD2P] - g2[iHM] - g2[iDM];
        if(iELEC < NCALC)
          {
          //g2[iH] = 1.e0 - 2.e0 * g2[iH2]  - g2[iHP]  - g2[iHD] - 2.*g2[iH2P] - g2[iHM];
          //g2[iD]  = abundD - g2[iDP] - g2[iHD];
          //g2[iHE] = abhe - g2[iHEP] + g2[iHEPP];
          //g2[iDM] = g2[iHDP] = g2[iD2P] = g2[iD2] = smallx;
          }                

 
        for(i=0;i<NCALC;i++) y[i] = ysav[i] + a31*g1[i] + a32*g2[i];
        
        
        // Set up and solve RHS for g3
        x = xsav + a3x*h;
        REAL g3[NCALC];
        derivs(x,dydx,ratraw,y);
        for(i=0;i<NCALC;i++) g3[i] = dydx[i] + (c31*g1[i] + c32*g2[i])/h;
        for(j=0;j<NCALC;j++)
            for(i=0;i<NCALC;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NCALC;i++) bb[i] = g3[i];
        leqs(av,bb);
        for(j=0;j<NCALC;j++) g3[j] = bb[j];
        if(iELEC < NCALC) g3[iELEC] = g3[iHP] + g3[iDP] + 2.0*g3[iHEPP] + g3[iHEP] + g3[iHDP] + g3[iH2P] + g3[iD2P] - g3[iHM] - g3[iDM];
        if(iELEC < NCALC)                   
          {
          //g3[iH] = 1.e0 - 2.e0 * g3[iH2]  - g3[iHP]  - g3[iHD] - 2.*g3[iH2P] - g3[iHM];
          //g3[iD]  = abundD - g3[iDP] - g3[iHD];
          //g3[iHE] = abhe - g3[iHEP] + g3[iHEPP];
          //g3[iDM] = g3[iHDP] = g3[iD2P] = g3[iD2] = smallx;
          }
       
 
        // And finally, g4
        REAL g4[NCALC];
        for(i=0;i<NCALC;i++) g4[i]= dydx[i]+(c41*g1[i]+c42*g2[i]+c43*g3[i])/h;
        for(j=0;j<NCALC;j++)
            for(i=0;i<NSPECIES;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NCALC;i++) bb[i] = g4[i];
        leqs(av,bb);
        for(j=0;j<NCALC;j++) g4[j] = bb[j];
        if(iELEC < NCALC) g4[iELEC] = g4[iHP] + g4[iDP] + 2.0*g4[iHEPP] + g4[iHEP] + g4[iHDP] + g4[iH2P] + g4[iD2P] - g4[iHM] - g4[iDM];
        if(iELEC < NCALC)                   
          {
          //g4[iH] = 1.e0 - 2.e0 * g4[iH2]  - g4[iHP]  - g4[iHD] - 2.*g4[iH2P] - g4[iHM];
          //g4[iD]  = abundD - g4[iDP] - g4[iHD];
          //g4[iHE] = abhe - g4[iHEP] + g4[iHEPP];
          //g4[iDM] = g4[iHDP] = g4[iD2P] = g4[iD2] = smallx;
          }
       
 
        
        // Compute 3rd and 4th order estimates for y
        REAL err[NCALC];
        double b1 = 19.0e0/9.0e0, b2 = 1.0e0/2.0e0,
        b3 = 25.0e0/108.0e0, b4 = 125.0e0/108.0e0;
        double e1 = 17.0e0/54.0e0, e2 = 7.0e0/36.0e0,
        e3 = 0.0, e4 = 125.0e0/108.0e0;
        
        for(i=0;i<NCALC;i++)
        {
            y[i] = ysav[i] + b1*g1[i] + b2*g2[i] + b3*g3[i] + b4*g4[i];
            err[i] = e1*g1[i] + e2*g2[i] + e3*g3[i] + e4*g4[i];
        }
        
        
        
        // Time advances
        x = xsav + h;
        
        if(x==xsav)
        {
            // Something has gone wrong
            cout << "ERROR: Uh oh, h=0" << endl;
        }
        
        
        // Determine the scaled accuracy
        REAL errmax = 0.0;
        for(i=0;i<NCALC;i++) errmax = max(errmax,fabs(err[i]/yscal[i]));
        errmax = errmax/eps;
        
        
        REAL safety = 0.9, grow= 1.5, pgrow= -0.25, shrnk= 0.5, pshrnk= -1.0/3.0;
        REAL errcon = 0.1296;
        if(errmax <= 1.0) // All went well
        {
            hdid = h;
            if(errmax > errcon) hnext = safety*h*pow(errmax,pgrow);
            else hnext = grow*h;
            jcounts = jcounts+jtry;
            
            if(jcounts>400000 && jcounts%100000 == 0)
            {
                // Why is this taking so long?!
                cout << "Jcounts is too high (jcounts="<<jcounts<<")" << endl;
            }

            return;
            
        }
        else // Error too large, gots to try again
        {
            hnext = safety*h*pow(errmax,pshrnk);
            double hsign = 1.0;
            if(h<0) hsign = -1.0;
            if(h<0) cout << "ERROR: hsign is < 0?!" << endl; // should never happen
            h = hsign*max(fabs(hnext),shrnk*fabs(h));
            //cout << "Failed step " << h << " " << hnext << " " << hdid << endl;
            
            
        }
        
       
        
        
    } // End of solve equation loop

    // should never get here
    cout << "ERROR: Too many step tries (jtry ="<<jtry<<")" <<endl;
};

