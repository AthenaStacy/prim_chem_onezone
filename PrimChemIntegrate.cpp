#include "PrimChemIntegrate.H"

//#include "PrimChemGlobals.H"
//using namespace PrimChemNS;


// GetNetworkRates(REAL Temp,REAL Den,REAL* ys, REAL* aion, REAL* RatRaw, bool chem_cc_case,REAL* RadParams)

void PrimChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass,int* oOutputs, void (*derivs)(REAL,REAL*,REAL*,REAL* ),REAL* RadParams,REAL* aion, REAL* RatRaw, int &jcounts)
{
    // REAL start = start time (0)
    // REAL firstDT = first attempted dt
    // REAL tmin = min allowed dt
    // REAL dt = end time
    // REAL* ymass = initial conditions of all species
    // REAL* oOutputs = other information we want want to store and outside of integrate
    //
    int i;
    
    //REAL dxsav = 1;
    //int kmax = 1;
    REAL odescal = 1.0e-6;
    int stepMax = 1000000, stepCur = 0;
    

    int chem_cc_case = oOutputs[0];
    
    int nok = 0;
    int nbad = 0;
    
    // Store boundary conditions
    REAL y[NSPECIES];
    REAL yscal[NSPECIES];
    for(i=0;i<NSPECIES;i++) y[i] = ymass[i];
    
    REAL dydx[NSPECIES];
    
    REAL x = start;
    REAL h = firstDT, hdid = 0, hnext=0;
    
    for(stepCur=1;stepCur < stepMax ; stepCur++)
    {
     

        y[iELEC] = y[iHP] + y[iDP] + y[iHEP] + 2.0*y[iHEPP] + y[iHDP] + y[iH2P] + y[iD2P] - y[iHM] - y[iDM];
        y[iH] = 1.e0 - 2.e0 * y[iH2]  - y[iHP]  - y[iHD] - 2.*y[iH2P] - y[iHM];
        y[iD]  = abundD - y[iDP] - y[iHD];
        y[iHE] = abhe - y[iHEP] + y[iHEPP];
        y[iDM] = y[iHDP] = y[iD2P] = y[iD2] = smallx;
        
        // Yet another check that everything is positive
        for(i=0;i<NSPECIES;i++) y[i] = max(y[i],1.0e-20);
        
        // Update the derivatives
        //cout << "About to enter loop derivs" << endl;
        derivs(h,dydx,RatRaw,y);
        
        // Used in error checking
        for(i=0;i<NSPECIES;i++) yscal[i] = max(odescal,fabs(y[i]));
        
        // Check to make sure we're not overshooting the timestep
        if( (x+h-dt)*(x+h-start) > 0.0 ) h = dt-x;
        
        hdid = 0; hnext = 0;
        
        // Step away! (y,x,hdid,hnext,jcounts updated)
        StepMe(y,dydx,x,h,hdid,hnext,yscal,derivs,jcounts,RatRaw);
        
            
        // Tally up
        if( hdid == h ) nok = nok+1;
        else nbad = nbad + 1;
        
        
        // Checks
        if(fabs(h) < tmin)
        {
            // print an error message
            cout << "ERROR: with fabs(h) < tmin" << endl;
        }
        
        if(stepCur==stepMax-1)
        {
            // Exceed max number of steps...
            cout << "ERROR: Exceeded max number of steps" << endl;
        }
        
        // Are we done?
        if( stepCur==stepMax-1 || (x-dt)*(dt-start) >= 0)
        {
            for(i=0;i<NSPECIES;i++) ymass[i] = y[i];
            //cout << "Exiting Integrate (stepCur = " << stepCur <<")" << endl;
            break;
        }
        
        
        // If not, get ready to start over
        h = hnext;
        
        
    } // end of stepping loop
    
    //cout << "jcount = " << jcounts << " nok = " << nok << " and nbad = " << nbad << endl;
    return;
    
};
