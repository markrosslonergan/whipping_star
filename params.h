#ifndef params_H_
#define params_H_


/**************************
 * Bin and Spectra Parameters
 **************************/
#define N_m_bins 19
#define N_e_bins 11

#define N_e_spectra 7
#define N_m_spectra 2

#define N_dets 3
#define N_anti 2


/**************************
 * Physical Masses
 **************************/

#define MPROTON  0.938
#define MPION   0.13957
#define MKAON   0.493
#define MSIGMA  1.189


/**************************
 *   Fit specific params
 ***************************/

#define psmear  0.05
#define pismear  0.05
#define EMsmear  0.15
#define MUsmear  0.06

#define p_thresh  0.021
#define pip_thresh  0.00
#define pim_thresh  0.00
#define vertex_thresh  0.05
#define EM_thresh  0.200

/**************************
 *    Mode selection tools
 ***************************/

#define APP_ONLY 0
#define DIS_ONLY 1
#define BOTH_ONLY 2
#define WIERD_ONLY 3 	// This is combined, but with e_disapearance off.

#define NU_MODE 0
#define NU_NUBAR_MODE 1

/**************************
 *  	Detector Parameters
 **************************/
#define DET_SBND 0
#define DET_UBOONE 1
#define DET_ICARUS 2

//Metric to imperial ton, tonne, tonnes
#define MET2IMP 1 //0.9071847



#endif
