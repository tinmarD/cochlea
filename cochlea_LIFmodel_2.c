// [vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel_2(iSyn, Fe)
//
// NEW : no more refractory period but instead voltage is reset to a negative value V_RESET
// (this way the time to go back to rest potential depends on the input amplitude) (a very little bit slower)
//
// iSyn             : Input synaptic current (output of rectifier)      - double vector
// Fe               : Sampling rate of iSyn                             - scalar
// vNeuron          : Neuron potential                                  - double vector
// spikeTimeInd     : List of spikes temporal indices
// spikeNeuronInd   : List of spikes neuron indices (neuron which fired)
// Yrect: output vector, filtered and recitfier signal      - double or single vector

// Compilation:
//   mex -O cochlea_LIFmodel_2.c
// Consider C99 comments on Linux with GCC:
//   mex -O CFLAGS="\$CFLAGS -std=c99" cochlea_LIFmodel_2.c

// Headers: --------------------------------------------------------------------
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Define neurone caracteristics
#define V_THRESH    20.
#define V_SPIKE     35.
#define TAU         0.01
#define V_RESET  	-30

// Some macros to reduce the confusion:
#define iSyn_in         prhs[0]
#define fe_in           prhs[1]
#define vNeuron_out     plhs[0]
#define spkTimeInd_out  plhs[1]
#define spkChanInd_out  plhs[2]

// For message errors
#define ERR_ID   "cochlea_LIFmodel:"
#define ERR_HEAD "cochlea_LIFmodel[mex]: "


// Main function ----------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *iSyn, *vNeuron, tauFe, vNeuronMult, iMult, fe;
    mwSize nPnts, nChan, nSpk=0;
    mwSize i, spkCount=0, iTime, iChan, iTimeChan, iChanTime;
    double *spkTimeInd, *spkChanInd, *tempSpkTimeInd, *tempSpkChanInd;
    const mwSize *iSynDims;
    mwSize *outSpkListDims;

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
    
    tempSpkTimeInd  = mxCalloc(nPnts*nChan, sizeof(double));
    tempSpkChanInd  = mxCalloc(nPnts*nChan, sizeof(double));
    
    // Output variables
	vNeuron_out     = mxCreateNumericArray(2, iSynDims, mxDOUBLE_CLASS, mxREAL);
    vNeuron         = mxGetPr(vNeuron_out);
  
    vNeuron[0]      = 0;
    vNeuronMult     = (1-1/(fe*TAU));
    iMult           = 1/(fe*TAU);
    
    // initalize vNeuron at time 0
    for (iChan=0;iChan<nChan;iChan++){       
        vNeuron[iChan*nPnts]    = 0; 
    }
    
    nSpk = 0;
    for (iTime=1;iTime<nPnts;iTime++){
        for (iChan=0;iChan<nChan;iChan++){
            iTimeChan = iChan*nPnts+iTime;
            iChanTime = iTime*nChan+iChan;
            if ( vNeuron[iTimeChan-1]==V_SPIKE ){  // Neuron in refractory period, set potential to V_RESET
               vNeuron[iTimeChan]   = V_RESET;
            }
            else{
                vNeuron[iTimeChan]  = vNeuron[iTimeChan-1]*vNeuronMult+iSyn[iTimeChan-1]*iMult;
                if ( vNeuron[iTimeChan]>V_THRESH ){ // Spike
                    vNeuron[iTimeChan]          = V_SPIKE;
                    tempSpkTimeInd[iChanTime]   = (double)iTime;
                    tempSpkChanInd[iChanTime]   = (double)iChan;
                    nSpk++;
                }
            }
        }
    }
    
    outSpkListDims      = mxCalloc(2, sizeof(mwSize));
    outSpkListDims[0]   = nSpk;
    outSpkListDims[1]   = 1;
    spkTimeInd_out      = mxCreateNumericArray(2, outSpkListDims, mxDOUBLE_CLASS,  mxREAL);
    spkChanInd_out      = mxCreateNumericArray(2, outSpkListDims, mxDOUBLE_CLASS,  mxREAL);
    spkTimeInd          = mxGetPr(spkTimeInd_out);
    spkChanInd          = mxGetPr(spkChanInd_out);
    spkCount            = 0;
    for (i=0; i<(mwSize)nPnts*nChan; i++){
        if (tempSpkTimeInd[i]!=0){
            spkTimeInd[spkCount] = tempSpkTimeInd[i];
            spkChanInd[spkCount] = tempSpkChanInd[i];
            spkCount++;
        }
    }
    
    // Cleanup
    mxFree(tempSpkTimeInd);
    mxFree(tempSpkChanInd);
    mxFree(outSpkListDims);
    
    return;
}


