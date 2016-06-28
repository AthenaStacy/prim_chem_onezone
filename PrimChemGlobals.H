#ifndef PRIMCHEMGLOBALS_H
#define PRIMCHEMGLOBALS_H

#define NSPECIES 16
#define NCALC    16
#define NREACTION 124
#define real double
#define REAL double
#define abhe   0.078947 
#define abundD 2.6e-5
#define smallx 1.e-40

namespace PrimChemNS {
    
    enum Species
    {
        iHP   = 0,
        iHM   = 1,
        iH2P  = 2,
        iH2   = 3,
        iDP   = 4,
        iHD   = 5,
        iHEP  = 6,
        iHEPP = 7,
        
        
        iH    = 8,
        iD    = 9,
        iDM   = 10,
        iHDP  = 11,
        iD2P  = 12,
        iD2   = 13,
        iHE   = 14,
        iELEC = 15
    };
    
    enum Reactions{
        
        iR001 = 0, iR002, iR003, iR004, iR005, iR006, iR007, iR008, iR009, iR010,
        iR011, iR012, iR013, iR014, iR015, iR016, iR017, iR018, iR019, iR020,
        iR021, iR022, iR023, iR024, iR025, iR026, iR027, iR028, iR029, iR030,
        iR031, iR032, iR033, iR034, iR035, iR036, iR037, iR038, iR039, iR040,
        iR041, iR042, iR043, iR044, iR045, iR046, iR047, iR048, iR049, iR050,
        iR051, iR052, iR053, iR054, iR055, iR056, iR057, iR058, iR059, iR060,
        iR061, iR062, iR063, iR064, iR065, iR066, iR067, iR068, iR069, iR070,
        iR071, iR072, iR073, iR074, iR075, iR076, iR077, iR078, iR079, iR080,
        iR081, iR082, iR083, iR084, iR085, iR086, iR087, iR088, iR089, iR090,
        iR091, iR092, iR093, iR094, iR095, iR096, iR097, iR098, iR099, iR100,
        iR101, iR102, iR103, iR104, iR105, iR106, iR107, iR108, iR109, iR110,
        iR111, iR112, iR113, iR114, iR115, iR116, iR117, iR118, iR119, iR120,
        iR121, iR122, iR123, iR124
    };
    
}

#endif



