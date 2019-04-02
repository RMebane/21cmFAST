#ifndef _SOURCES_
#define _SOURCES_

// Function definitions for source populations
// Rick Mebane 2018

#include "INIT_PARAMS.H"
#include "ANAL_PARAMS.H"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// more general, allows for halos with multiple populations
typedef struct{
    double (*fesc)(double, double); //escape fraction as a function of halo mass and redshift
    double (*Nion)(double, double); //Nion, generally just set to 4e4 for Pop III, 5e3 for Pop II, and weighted average if a halo has both
    double (*fstar)(double, double); //star formation efficiency, defined as fraction of a halo's gas in stars. combine both Pop II and Pop III here, differences set in Nion and fPopIII
                                     //NOTE: fstar here is CUMULATIVE, which is very important for Pop III... make sure to add up all the stellar mass that has been in the halo, not that is there currently
    double (*fx)(double, double); //boost to xray emissivity, fiducially 1 for Pop II and 100 for Pop III in Yuxiang's 21cmFAST w/ minihalos
    double (*minMass)(double); //minimum halo mass for star formation (both Pop II and Pop III)
    double (*fPopIII)(double); //mass fraction of Pop III halos at redshift z. If no mixed populations, just set to 0 or 1. Needed for spectral emissivity calculation
} sources;

// atomic cooling threshold
// used for min mass comparisons
double minMassATC(double z)
{
    double mu;
    mu = 1.22;
    return (double) TtoM(z, 1.0e4, mu);
}

// initialize source struct
sources setSources(double (*fesc)(double, double),
                   double (*Nion)(double, double),
                   double (*fstar)(double, double),
                   double (*fx)(double, double),
                   double (*minMass)(double),
                   double (*fPopIII)(double))
{
    sources s;
    s.fesc = fesc;
    s.Nion = Nion;
    s.fstar = fstar;
    s.minMass = minMass;
    s.fx = fx;
    s.fPopIII = fPopIII;
    return s;
}

double zeta_at_M(sources s, double z, double M)
{
    return 1.22 * s.fstar(z, M) * s.fesc(z, M) * s.Nion(z, M);
}

// BEGIN SPECIFIC SOURCE FUNCTIONS

double fstar_jordan(double z, double M)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    double table[26*61];

    const char *s = " ";
    char *token = NULL;

    int i = 0;
    int j;

    fp = fopen("../Parameter_files/smhm_uvlf.txt", "r");
    if(fp == NULL)
    {
        printf("Error opening");
        exit(EXIT_FAILURE);
    }

    while ((read=getline(&line, &len, fp)) != -1)
    {
        token = strtok(line, s);
        while(token != NULL)
        {
            table[i] = atof(token);
            token=strtok(NULL, s);
            i++;
        }
    }

    fclose(fp);

    // get fstar from table
    // Not general, but just wanted to get this working quick
    int zind, Mind;
    zind = (int) round(z-5);
    Mind = (int) round(((log10(M) - 8) / 5) * 50 + 10);
    return table[Mind * 26 + zind] * .28 / .046 * 2.0;

}

// Eq. 1 from Feng & Holder 2018
// Not really a source pop function, but putting it here for now...
// Tcmb given at z=0 
double T_background(double Tcmb, double excess_frac, double z)
{
    double beta, Tr;
    beta = -2.62;
    Tr = 1.19;
    return (Tcmb + excess_frac * Tr * pow(1.428 / (1.0 + z), beta)) * (1.0 + z);
}

double Mmin_jordan(double z)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    double table[26*61];

    const char *s = " ";
    char *token = NULL;

    int i = 0;
    int j;

    fp = fopen("../Parameter_files/smhm_uvlf.txt", "r");
    if(fp == NULL)
    {
        printf("Error opening");
        exit(EXIT_FAILURE);
    }

    while ((read=getline(&line, &len, fp)) != -1)
    {
        token = strtok(line, s);
        while(token != NULL)
        {
            table[i] = atof(token);
            token=strtok(NULL, s);
            i++;
        }
    }

    fclose(fp);

    double logM = 7;
    int zind, Mind;
    while(logM < 10)
    {
        Mind = (int) round(((logM - 8) / 5) * 50 + 10);
        zind = (int) round(z-5);

        if(table[Mind * 26 + zind] > 0) 
        {
            return pow(10.0, logM);
        }
        logM += 0.1;
    }
    return 1.0e8;
    
}

// just returns 1 if res is > 1
// useful for fstar and fesc functions
double maxOne(double res)
{
    if(res > 1.0) return 1.0;
    return res;
}

// Weighted by the exponential cutoff
// Check to make sure fstar isn't above one BEFORE applying cutoff
double fstar_21cmfid(double z, double M)
{
    double res;
    res = maxOne(0.05 * pow(M / 1.0e10, 0.5));
    return res * exp(-5e8 / M);
}

double fesc_21cmfid(double z, double M)
{
    double res;
    res = 0.1 * pow(M / 1.0e10, -0.5);
    return maxOne(res);
}

double fxDef(double z, double M)
{
    if(M > minMassATC(z)) return 1.0;
    return 100.0;
}

double minMassDef(double z)
{
    //return 1.0e4;
    double res;
    res = pow(10.0, -3.0/44.0 * z + 8.4); // quick linear fit to fiducial min mass in log space for testing
                                           // will do this right once we make sure everything works
    if(res > minMassATC(z)) return minMassATC(z);
    return res;
}

double minMassDefMom(double z)
{
    double res;
    res = pow(10.0, -1.0/14.0 * z + 60.0/7.0);
    if(res > minMassATC(z)) return minMassATC(z);
    return res;
}

double NionDef(double z, double M)
{
    if(M > minMassATC(z)) return 5000.0;
    return 40000.0;
}

double fstarDef(double z, double M)
{
    if(M > minMassATC(z)) return maxOne(0.03 * pow(M / 1.0e10, 0.4));
    return 0.0005;
}

double fstarTest(double z, double M)
{
    if(M > 1.0e8) return 0.05;
    return 0.005;
}

double fescDef(double z, double M)
{
    if(M > minMassATC(z)) return 0.1;
    return 1.0;
}

double etaSteveEnergy(double z, double M)
{
    return pow(pow(10, 11.5) / M, 2.0/3.0) * 9.0 / (1.0 + z);
}

double fstarSteveEnergy(double z, double M)
{
    return maxOne(1.0 / (1.0 + etaSteveEnergy(z, M)));
}

double etaSteveMom(double z, double M)
{
    return 0.2 * pow(pow(10, 11.5) / M, 1.0/3.0) * pow(9.0 / (1.0 + z), 1.0/2.0);
}

double fstarSteveMom(double z, double M)
{
    return maxOne(1.0 / (1.0 + etaSteveMom(z, M)));
}

double fstarFid(double z, double M)
{
    if (M > minMassATC(z)) return fstarSteveEnergy(z, M);
    return 0.005;
}

double fstarFidMom(double z, double M)
{
    if (M > minMassATC(z)) return fstarSteveMom(z, M);
    return 0.005;
}

double noPopIII(double z)
{
    return 0;
}

double allPopIII(double z)
{
    return 1;
}

double PopIII_frac_test(double z)
{
    if(z < 20) return 0;
    if(z > 25) return 1.0;
    return 0.5;
}

double minMass_21cmfid(double z) {return 5e8 / 50.0;}

// source functions for comparing with sharp cutoff model
double minMassCutoff(double z) {return 5e8;}
double NionCutoff(double z, double M) {return 5000;}
double fescCutoff(double z, double M) {return 0.1;}
double fstarCutoff(double z, double M) {return 0.05;}
double fxCutoff(double z, double M) {return 1.0;}

sources MomPopIII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarFidMom, fxDef, minMassDefMom, noPopIII);
    return s;
}

sources EnergyPopIII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarFid, fxDef, minMassDef, noPopIII);
    return s;
}

sources MomPopII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarSteveMom, fxDef, minMassATC, noPopIII);
    return s;
}

sources EnergyPopII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarSteveEnergy, fxDef, minMassATC, noPopIII);
    return s;
}

sources fid_21cmFAST()
{
    sources s;
    s = setSources(fesc_21cmfid, NionCutoff, fstar_21cmfid, fxCutoff, minMass_21cmfid, noPopIII);
    return s;
}

// set the sources here
sources defaultSources()
{
    sources s;
    s = MomPopIII();
    return s;
}

#endif
