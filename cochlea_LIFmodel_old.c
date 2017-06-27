// [vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel_old(iSyn, Fe)
//
// Reduced and simplified version of FilterX.c from Jan Simon : 
// http://www.mathworks.com/matlabcentral/fileexchange/32261-filterm
//
// iSyn             : Input synaptic current (output of rectifier)      - double vector
// Fe               : Sampling rate of iSyn                             - scalar
// vNeuron          : Neuron potential                                  - double vector
// spikeTimeInd     : List of spikes temporal indices
// spikeNeuronInd   : List of spikes neuron indices (neuron which fired)
// Yrect: output vector, filtered and recitfier signal      - double or single vector

// Compilation:
//   mex -O cochlea_LIFmodel_old.c
// Consider C99 comments on Linux with GCC:
//   mex -O CFLAGS="\$CFLAGS -std=c99" cochlea_LIFmodel_old.c

// Headers: --------------------------------------------------------------------
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Define neurone caracteristics
#define V_REST      -70.
#define V_THRESH    -50.
#define V_SPIKE     35.
#define TAU         0.001
#define T_REFRACT   0.003
#define RM          300.

// Some macros to reduce the confusion:
#define iSyn_in         prhs[0]
#define fe_in           prhs[1]
#define vNeuron_out     plhs[0]
#define spkTimeInd_out  plhs[1]
#define spkChanInd_out  plhs[2]

// For message errors
#define ERR_ID   "cochlea_LIFmodel_old:"
#define ERR_HEAD "cochlea_LIFmodel_old[mex]: "


// Main function ----------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *iSyn, *vNeuron, *spkTimeInd, *spkChanInd, tauFe, denomInv, fe;
    mwSize nPnts, nChan;
    mwSize iTime, iChan, iTimeChan, tRefractInd, *lastSpikeInd;
    const mwSize *iSynDims;

    // Check number of inputs
    if (nrhs!=2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "Wrong number of inputs : [vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel(iSyn,Fe)");
    }
    
    // Get dimensions of inputs:
    nPnts       = mxGetM(iSyn_in);    // nPnts
    nChan       = mxGetN(iSyn_in);    // nChan
    iSynDims    = mxGetDimensions(iSyn_in);
    if (mxGetNumberOfDimensions(iSyn_in)>2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "iSyn input must be a matrix or vector (2D max)");
    }
    // make sure the Fe input argument is scalar
    if( !mxIsDouble(fe_in) || mxIsComplex(fe_in) || mxGetNumberOfElements(fe_in)!=1 ) {
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "fe argument must be a scalar");
    }
    iSyn            = mxGetPr(iSyn_in);
    fe              = mxGetScalar(fe_in);
    
    // Output variables
	vNeuron_out     = mxCreateNumericArray(2, iSynDims, mxDOUBLE_CLASS, mxREAL);
    spkTimeInd_out  = mxCreateNumericArray(2, iSynDims, mxINT32_CLASS,  mxREAL);
    spkChanInd_out  = mxCreateNumericArray(2, iSynDims, mxINT32_CLASS,  mxREAL);
    
    vNeuron         = mxGetPr(vNeuron_out);
    spkTimeInd      = mxGetPr(spkTimeInd_out);
    spkChanInd      = mxGetPr(spkChanInd_out);
    
    vNeuron[0]      = V_REST;
    tauFe           = TAU*fe;
    denomInv        = 1/(1+tauFe);
    tRefractInd     = (mwSize) (T_REFRACT*fe);
    
    // Allocate memory for lastSpikeInd and 
    if ((lastSpikeInd = mxCalloc(nChan, sizeof(mwSize))) == NULL) {
        mexErrMsgIdAndTxt(ERR_ID   "NoMemory", ERR_HEAD "No memory for lastSpikeInd");
    }
    // initalize lastSpikeInd and vNeuron at time 0
    for (iChan=0;iChan<nChan;iChan++){
        lastSpikeInd[iChan]     = -2*tRefractInd;        
        vNeuron[iChan*nPnts]    = V_REST; 
    }

//     FOR ONE CHANNEL
//     lastSpikeInd    = -2*tRefractInd;
//     for (iTime=1;iTime<nPnts;iTime++){
//         if ( (iTime-lastSpikeInd)<tRefractInd ){  // Neuron in refractory period, V=V_REST
//            vNeuron[iTime]   = V_REST;
//         }
//         else{
//             vNeuron[iTime]  = denomInv*(V_REST+tauFe*vNeuron[iTime-1]+RM*iSyn[iTime]);
//             if ( vNeuron[iTime]>V_THRESH ){ // Spike
//                 vNeuron[iTime]  = V_SPIKE;
//                 lastSpikeInd    = iTime;
//             }
//         }
//     }
    
    for (iTime=1;iTime<nPnts;iTime++){
        for (iChan=0;iChan<nChan;iChan++){
            iTimeChan = iChan*nPnts+iTime; 
            if ( (iTime-lastSpikeInd[iChan])<tRefractInd ){  // Neuron in refractory period, V=V_REST
               vNeuron[iTimeChan]   = V_REST;
            }
            else{
                vNeuron[iTimeChan]  = denomInv*(V_REST+tauFe*vNeuron[iTimeChan-1]+RM*iSyn[iTimeChan]);
                if ( vNeuron[iTimeChan]>V_THRESH ){ // Spike
                    vNeuron[iTimeChan]  = V_SPIKE;
                    lastSpikeInd[iChan] = iTime;
                }
            }
        }
    }
    
    // Cleanup
    mxFree(lastSpikeInd);
    
    return;
}


