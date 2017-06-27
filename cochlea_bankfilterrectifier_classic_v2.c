// [Y,Yrect] = cochlea_bankfilterrectifier_classic_v2(b,a,bLPfilt,aLPfilt,X)
//
// New in v2 : The cutoff freq of LP filtering after rectification is not constant
// as in previous version but follow the characteristic frequency of each bandpass filter
// 
// Reduced and simplified version of FilterX.c from Jan Simon : 
// http://www.mathworks.com/matlabcentral/fileexchange/32261-filterm
//
// b,a  : filter coefficients                               - double vector
// X    : input vector                                      - double or single vector
// Y    : output vector, filtered signal                    - double or single vector
// Yrect: output vector, filtered and recitfier signal      - double or single vector
// 
// Compilation:
//   mex -O cochlea_bankfilterrectifier_classic_v2.c
// Consider C99 comments on Linux with GCC:
//   mex -O CFLAGS="\$CFLAGS -std=c99" cochlea_bankfilterrectifier_classic_v2.c

// Headers: --------------------------------------------------------------------
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Some macros to reduce the confusion:
#define b_in        prhs[0]
#define a_in        prhs[1]
#define bLP_in      prhs[2]
#define aLP_in      prhs[3]
#define X_in        prhs[4]
#define Y_out       plhs[0]
#define Yrect_out   plhs[1]

// For message errors
#define ERR_ID   "cochlea_bankfilterrectifier_classic_v2:"
#define ERR_HEAD "cochlea_bankfilterrectifier_classic_v2[mex]: "

// Prototypes: -----------------------------------------------------------------
void CoreDoubleN(double *X, mwSize MX, double *a, double *b, mwSize order, double *Y);
void CoreDouble3(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y);
void CoreDouble5(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y);
void CoreDouble7(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y);
// LP filter
void CoreDouble2LPFilt(double *X, mwSize MX, mwSize NX, double *a, double *b, double *Y);
void CoreDouble3LPFilt(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y);
// Rectifier
void rectifierDouble(double *Y, double *Yrect, mwSize MY, mwSize nFilters, double *periodChar);
void rectifierDouble_v2(double *Y, double *Yrect, mwSize MY, mwSize nFilters, double *periodChar);


// Main function ----------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *a, *b, *aLP, *bLP, *X, *Y, *Yrect, *periodChar; 
//     mwSize *periodChar;
    mwSize order, orderLP, nParam, MX, NX, Xndims;
    mwSize Mb,MbLP,Nb,NbLP,Ma,MaLP,Na,NaLP,i,bnDims,anDims,nFilters;
    bool   allocate_ba = false;
    char   Rev[2];
    const mwSize *Xdims;
    mwSize Ydims[2];
        
    // Check number of inputs
    if (nrhs!=5){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "Wrong number of inputs : [Y,Yrect]= cochlea_bankfilterrectifier_classic_v2(b,a,X,periodChar)");
    }
    
    // Check input data type
    if (!mxIsDouble(b_in) || !mxIsDouble(a_in)){
        mexErrMsgIdAndTxt(ERR_ID   "TypeInput",ERR_HEAD "Filter parameters must be DOUBLES.");
    }
    // Check input data type
    if (!mxIsDouble(bLP_in) || !mxIsDouble(aLP_in)){
        mexErrMsgIdAndTxt(ERR_ID   "TypeInput",ERR_HEAD "Filter parameters must be DOUBLES.");
    }
    if (!mxIsDouble(X_in) && !mxIsSingle(X_in)){
        mexErrMsgIdAndTxt(ERR_ID   "TypeInput",ERR_HEAD "Input signal X must be DOUBLE or SINGLE.");
    }

    // Get dimensions of inputs:
    MX = mxGetM(X_in);    // First dimension
    NX = mxGetN(X_in);    // Product of trailing dimensions (?)
    
    Xndims = mxGetNumberOfDimensions(X_in);
    Xdims  = mxGetDimensions(X_in);
    if (Xndims!=2 || NX!=1){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "input vector X must be a 1-D column vector");
    }
    
    bnDims = mxGetNumberOfDimensions(b_in);
    anDims = mxGetNumberOfDimensions(a_in);
    Mb      = mxGetM(b_in);     // length of filters
    Nb      = mxGetN(b_in);     // number of filters
    Ma      = mxGetM(a_in);     // length of filters
    Na      = mxGetN(a_in);     // number of filters
    MbLP    = mxGetM(bLP_in);
    MaLP    = mxGetM(aLP_in);
    NbLP    = mxGetN(bLP_in);
    NaLP    = mxGetN(aLP_in);

    nFilters = Nb;   
    if (bnDims!=2 || anDims!=2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "b and a must be 2D matrices or 1D vectors");
    }
    if (Mb!=Ma || Nb!=Na){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "b and a must be of the same size");
    }
    if (Mb<2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "filter order must be a least 1 - bad filter coefficients");
    }
    if (MbLP!=MaLP || NbLP!=NaLP){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "bLP and aLP must be of the same size");
    }
    if (MbLP<2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "filter order must be a least 1 - bad filter coefficients (for LP filtering)");
    }
    //mexPrintf("%d filters of length %d (order %d)\n",Nb,Mb,Mb-1);
    b   = mxGetPr(b_in);
    a   = mxGetPr(a_in);
    bLP = mxGetPr(bLP_in);
    aLP = mxGetPr(aLP_in);

    // Check normalization of filters : a[0,:] == 0
    for (i=0;i<Nb;i++){
        if (a[i*Mb]!=1){
            mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "a[0] must be equal to 1");
        }
    }
    order   = Mb-1;
    orderLP = MbLP-1;
      
    // Y_out must be of size (length(x),numberOfFilters)
    Ydims[0] = MX;
    Ydims[1] = nFilters;
    //mexPrintf("Ydims: %d,%d\n",Ydims[0],Ydims[1]);
    // Filtering
    if (mxIsDouble(X_in)){
        Y_out = mxCreateNumericArray(2, (const) Ydims, mxDOUBLE_CLASS, mxREAL);
        Y     = mxGetPr(Y_out);
        X     = mxGetPr(X_in);
        // Filtering (band-pass filter bank)
        switch (order) {
            case 2:   CoreDouble3(X, MX, nFilters, a, b, Y);  break;
            case 4:   CoreDouble5(X, MX, nFilters, a, b, Y);  break;
            case 6:   CoreDouble7(X, MX, nFilters, a, b, Y);  break;
            default:  mexPrintf("Wrong filter order");        break;
        }
        // Rectifier
        Yrect_out   = mxCreateNumericArray(2, (const) Ydims, mxDOUBLE_CLASS, mxREAL);
        Yrect       = mxGetPr(Yrect_out);
        // Half-wave rectification
        rectifierDouble_v2(Y,Yrect,MX,nFilters,periodChar);
        // Filtering (band-pass filter bank)
        switch (orderLP) {
            case 2:   CoreDouble3LPFilt(Yrect, MX, nFilters, aLP, bLP, Yrect);  break;
            default:  mexPrintf("Wrong LP filter order"); break;
        }
    }
    
    return;
}


// Rectifier: ---------------------------------------------------------------
void rectifierDouble_v2(double *Y, double *Yrect, mwSize MY, mwSize nFilters, double *periodChar)
{
    mwSize i, nElements;
    nElements = MY*nFilters;
    
    // Half-wave rectification
    for (i=0;i<nElements;i++){
        if (Y[i]<0){
            Yrect[i] = 0;
        }
        else{
            Yrect[i] = Y[i];
        }        
    }
}


// Filtering functions: ------------------------------------------------------

void CoreDouble2LPFilt(double *X, mwSize MX, mwSize NX, double *a, double *b, double *Y)
{
    // Filter with loop unrolled for 2 parameters (filter order 1).
    // Same input as the CoreDoubleN, but ommited [order], because it is 1.

    double Xi, Yi, z0, a1 = a[1];
    mwSize i = 0, C;

    while (NX--) {
        z0 = 0;
        C  = i + MX;
        while (i < C) {
            Xi = X[i];
            Yi = b[0] * Xi + z0;
            z0 = b[1] * Xi - a1 * Yi;
            Y[i++] = Yi;//pow(Yi,1/3.);
        }
    }

    return;
}

void CoreDouble3LPFilt(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y)
{
	double Xi, Yi, z0, z1, a1 = a[1], a2 = a[2];
	mwSize iX=0, iY=0, iFilt=0, iSig=0;
    
    while (nFilters--){
        iX = 0;
        a1 = a[iFilt+1];
        a2 = a[iFilt+2];
        z0 = 0;
        z1 = 0;
        while (iX < MX){
            Xi = X[iSig+iX++];
            Yi = b[iFilt]   * Xi + z0;
            z0 = b[iFilt+1] * Xi + z1 - a1 * Yi;
            z1 = b[iFilt+2] * Xi      - a2 * Yi;
            Y[iY++] = Yi;
        }
        iFilt = iFilt+3;
        iSig  = iSig+MX;
    }

    return;
}



// Filtering of Input X
void CoreDouble3(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y)
{
	double Xi, Yi, z0, z1, a1 = a[1], a2 = a[2];
	mwSize iX=0, iY=0, iFilt=0;
    
    while (nFilters--){
        iX = 0;
        a1 = a[iFilt+1];
        a2 = a[iFilt+2];
        z0 = 0;
        z1 = 0;
        while (iX < MX){
            Xi = X[iX++];
            Yi = b[iFilt]   * Xi + z0;
            z0 = b[iFilt+1] * Xi + z1 - a1 * Yi;
            z1 = b[iFilt+2] * Xi      - a2 * Yi;
            Y[iY++] = Yi;
        }
        iFilt = iFilt+3;
    }

    return;
}

void CoreDouble5(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y)
{
    // Still 33% faster than the loop method.
    double Xi, Yi, z0, z1, z2, z3,
    a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6];
    mwSize iX=0,iY=0,iFilt=0;

    while (nFilters--) {
        iX = 0;
        a1 = a[iFilt+1];
        a2 = a[iFilt+2];
        a3 = a[iFilt+3];
        a4 = a[iFilt+4];
        z0 = 0;
        z1 = 0;
        z2 = 0;
        z3 = 0;
        while (iX < MX) {
            Xi = X[iX++];
            Yi = b[iFilt]   * Xi + z0;
            z0 = b[iFilt+1] * Xi + z1 - a1 * Yi;
            z1 = b[iFilt+2] * Xi + z2 - a2 * Yi;
            z2 = b[iFilt+3] * Xi + z3 - a3 * Yi;
            z3 = b[iFilt+4] * Xi      - a4 * Yi;
            Y[iY++] = Yi;
        }
        iFilt = iFilt+5;
    }
    return;
}

void CoreDouble7(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Y)
{
    // Still 33% faster than the loop method.
    double Xi, Yi, z0, z1, z2, z3, z4, z5,
    a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6];
    mwSize iX=0,iY=0,iFilt=0;

    while (nFilters--) {
        iX = 0;
        a1 = a[iFilt+1];
        a2 = a[iFilt+2];
        a3 = a[iFilt+3];
        a4 = a[iFilt+4];
        a5 = a[iFilt+5];
        a6 = a[iFilt+6];
        z0 = 0;
        z1 = 0;
        z2 = 0;
        z3 = 0;
        z4 = 0;
        z5 = 0;
        while (iX < MX) {
            Xi = X[iX++];
            Yi = b[iFilt] * Xi + z0;
            z0 = b[iFilt+1] * Xi + z1 - a1 * Yi;
            z1 = b[iFilt+2] * Xi + z2 - a2 * Yi;
            z2 = b[iFilt+3] * Xi + z3 - a3 * Yi;
            z3 = b[iFilt+4] * Xi + z4 - a4 * Yi;
            z4 = b[iFilt+5] * Xi + z5 - a5 * Yi;
            z5 = b[iFilt+6] * Xi      - a6 * Yi;
            Y[iY++] = Yi;
        }
        iFilt = iFilt+7;
    }
    return;
}

