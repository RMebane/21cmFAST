/*
NOTE:  We encourage the user to perform their own convergance tests if modifying the code significantly from the fiducial model.  Preliminary convergance studies were done around the fiducial model however.
 */

#ifndef _HEAT_PARAMS_H_
#define _HEAT_PARAMS_H_


/*
  Flag specifying whether one wants to calculate the spin temperature or just assume (Ts>>Tk)
*/
#define USE_TS_IN_21CM (int) (0)


/*
  Soft band X-ray luminosity, L_X. Used for determining the number of X-ray photons produced per stellar baryon.
  In v1.4, the old ZETA_X is replaced with a soft-band luminosity L_X. The soft band luminosity now determines the normalisation, 
 and is available as a free parameter. This is defined as the <2keV luminosity, hence NU_X_BAND_MAX = 2000.00.
  Defined as log10(L_X). E.g. 40 = log10(10^40)
*/
#define L_X (double) (40.5)
#define NU_X_BAND_MAX (double) (2000*NU_over_EV)


/*
  Smallest frequency X-ray photon contributing to IGM heating and ionization.
  This is directly related to the mean HI column density of the first galaxies.
  Note default 0.5 keV corresponds to ~ <N_HI> = 10^41.5 cm^-2 (c.f Das et al. 2017)
  There might be wierd behaviour if you set this to less than the HeII ionization
  threshold, <~54.4 eV
*/
#define NU_X_THRESH (double) (500*NU_over_EV)


/*
  X-Ray spectral index at frequencies higher than NU_X_THRESH
 */
#define X_RAY_SPEC_INDEX (double) (1.0)


/*
Filter choice used for smoothing the linear density field to obtain the collapsed fraction in Ts.c
0 = use real space top hat filter
1 = use sharp k-space filter
2 = use gaussian filter
*/
#define HEAT_FILTER (int) 0


#define CLUMPING_FACTOR (double) (2) /* sub grid scale.  note that if you want to run-down from a very high redshift (>50), you should set this to one.. */


/*
  Maximum redshift used in the Tk and x_e evolution equations.
  Temperature and x_e are assumed to be homogenous at higher redshifts.
*/
#define Z_HEAT_MAX (float) (30)

/*
  Maximum radius of influence for computing Xray and Lya pumping in cMpc. This
  should be larger than the mean free path of the relevant photons.  If increasing,
  you might want to adjust the z'' loop in Ts.c to skip integrating from early times,
  to increase speed (i.e. the code integrates as far back as ZPP_MAX, where R_XLy_MAX
  is the comoving separation between ZPP_MAX and Z_HEAT_MAX).
*/
#define R_XLy_MAX (float) (500)

/*
  Number of spherical annuli used to compute df_coll/dz' in the simulation box.
  The spherical annuli are evenly spaced in log R, ranging from the cell size to the box size.
  Ts.c will create this many boxes of size HII_DIM, so be wary of memory usage if values are high.
*/
#define NUM_FILTER_STEPS_FOR_Ts (int) (40)

/*
  Redshift step-size used in the z' integral.  Logarithmic dz.
*/
#define ZPRIME_STEP_FACTOR (float) (1.02)


/*
  If the following are >0, then the user chooses to overwrite the default boundary conditions
  at Z_HEAT_MAX obtained from RECFAST, and use his/her own.
 */
#define TK_at_Z_HEAT_MAX (double) (-1)
#define XION_at_Z_HEAT_MAX (double) (-1)



/*
  Set this to 1, if you want the spin temperature field printed at eary step of the evolution.
  This would be useful if investigating redshift evolution, so Ts.c needs to be called only once
  at the low redshift, and then one has all of the higher redshift Ts boxes upon completion.
 */
#define Ts_verbose (int) (1)


/*
 Filenames of the appropriate output from RECFAST to be used as boundary conditions in Ts.c
 as well as other tables used to compute the spin temperature
*/
#define RECFAST_FILENAME (const char *) "../External_tables/recfast_LCDM.dat"
#define STELLAR_SPECTRA_FILENAME (const char *) "../External_tables/stellar_spectra.dat"
#define KAPPA_EH_FILENAME (const char *) "../External_tables/kappa_eH_table.dat"
#define KAPPA_PH_FILENAME (const char *) "../External_tables/kappa_pH_table.dat"


/*
  Stellar Population responsible for early heating
  Pop == 2 Pop2 stars heat the early universe
  Pop == 3 Pop3 stars heat the early universe
*/
#define Pop (int) (2)


/*
  Number of ionizing photons per baryon of the two stellar species
*/
#define Pop2_ion (float) (4361)
#define Pop3_ion (float) (44021)

/*
  Flag to turn off or on verbose status messages in Ts.  The GSL libraries are very finicky,
  and this is useful when to help issolate why Ts crashed, if it did...
 */
#define DEBUG_ON (int) (0)

#endif
