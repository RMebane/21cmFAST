#ifndef _SOURCES_
#define _SOURCES_

// Function definitions for source populations
// This might go somewhere else in the code, but doing it all here for now
// Rick Mebane 2017

#include "INIT_PARAMS.H"
#include "ANAL_PARAMS.H"

// struct containing function pointers for source populations
typedef struct{
    double (*fesc)(double, double);
    double (*fescIII)(double, double);
    double (*Nion)(double, double);
    double (*NionIII)(double, double);
    double (*fstar)(double, double);
    double (*fstarIII)(double, double);
    double (*minMass)(double);
    double (*minMassIII)(double);
} sources;

// params struct for integration
struct parameters_gsl_ionEff_int_{
    sources s;
    double z;
};

// atomic cooling threshold
// used for min mass comparisons


double minMassATC(double z)
{
    double mu;
    mu = 1.22;
    return (double) TtoM(z, 1.0e4, mu);
}

// initialize source struct
// if not using Pop III stars, make minMassIII(z) = -1
sources setSources(double (*fesc)(double, double),
                   double (*fescIII)(double, double),
                   double (*Nion)(double, double),
                   double (*NionIII)(double, double),
                   double (*fstar)(double, double),
                   double (*fstarIII)(double, double),
                   double (*minMass)(double),
                   double (*minMassIII)(double))
{
    sources s;
    s.fesc = fesc;
    s.fescIII = fescIII;
    s.Nion = Nion;
    s.NionIII = NionIII;
    s.fstar = fstar;
    s.fstarIII = fstarIII;
    s.minMass = minMass;
    s.minMassIII = minMassIII;
    return s;
}

//integrand for ionEff
// returns dn/dlnM * fesc * fstar * Nion * A_He
double ionEffIntegrand(double lnM, void *params)
{
    sources s;
    double A_He, z, M;
    struct parameters_gsl_ionEff_int_ vals = 
        *(struct parameters_gsl_ionEff_int_ *)params;
    s = vals.s;
    z = vals.z;
    M = exp(lnM);
    A_He = 1.22;
    // this is a Pop II star forming halo
    if(M >= s.minMass(z))
        return A_He * s.fstar(z, M) * s.fesc(z, M) * s.Nion(z, M) 
               * dNdM_st(z, M) * M;
    // this is a Pop III star forming halo
    else if (M >= s.minMassIII(z) && s.minMassIII(z) > 0)
        return A_He * s.fstarIII(z, M) * s.fescIII(z, M) * s.NionIII(z, M) 
               * dNdM_st(z, M) * M;
    // this halo is below the minimum mass
    return 0.0;
}
    
// integrand for integral over mass function
double massFuncInt(double lnM, void *params)
{
    double z, M;
    z = *(double *)params;
    M = exp(lnM);
    return dNdM_st(z, M) * M;
}


// integrate over the mass function to calculate average ionizing efficiency
// for the given source popultion
double ionEff(double z, sources s)
{
    double mMinII, mMinIII, A_He, rel_tol, lower_limit, upper_limit;
    double resultEff, resultMf, error, test;
    gsl_function F;

    rel_tol = 0.001;
    mMinII = s.minMass(z);
    mMinIII = s.minMassIII(z);

    if(mMinIII < 0 || mMinIII >= mMinII) lower_limit = log(mMinII);
    else lower_limit = log(mMinIII);
    upper_limit = log(FMAX(1.0e16, exp(lower_limit)));

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    struct parameters_gsl_ionEff_int_ parameters_gsl_ionEff = {
        .s = s,
        .z = z,
    };

    F.function = &ionEffIntegrand;
    F.params = &parameters_gsl_ionEff;

    gsl_integration_qag(&F, lower_limit, upper_limit, 0, rel_tol, 1000, 
                        GSL_INTEG_GAUSS61, w, &resultEff, &error);

    F.function = &massFuncInt;
    F.params = &z;


    gsl_integration_qag(&F, lower_limit, upper_limit, 0, rel_tol, 1000,
                        GSL_INTEG_GAUSS61, w, &resultMf, &error);

    gsl_integration_workspace_free(w);

    return resultEff / resultMf;
}

// BEGIN SPECIFIC SOURCE FUNCTIONS

double fescDef(double z, double M) {return 0.1 * pow(M / 1.0e10, -0.4);}
double fescIIIDef(double z, double M) {return 0.1;}
double NionDef(double z, double M) {return 5000.0;}
double NionIIIDef(double z, double M) {return 40000.0;}
double fstarDef(double z, double M) {return 0.03 * pow(M / 1.0e10, 0.4);}
double fstarIIIDef(double z, double M) {return 1.0e-4;}
double minMassIIIDef(double z) {return minMassATC(z) / 10.0;}

sources defaultSources()
{
    sources s;
    s = setSources(fescDef, fescIIIDef, NionDef, NionIIIDef, 
                   fstarDef, fstarIIIDef, minMassATC, minMassIIIDef);
    return s;
}


   

#endif
