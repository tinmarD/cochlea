// [Y,Yrect] = cochlea_bankfilterrectifier_hybrid(b,a,X,periodChar)
//
// Reduced and simplified version of FilterX.c from Jan Simon : 
// http://www.mathworks.com/matlabcentral/fileexchange/32261-filterm
//
// b,a  : filter coefficients                               - double vector
// X    : input vector                                      - double or single vector
// Y    : output vector, filtered signal                    - double or single vector
// Yrect: output vector, filtered and recitfier signal      - double or single vector

// Compilation:
//   mex -O cochlea_bankfilterrectifier_hybrid.c
// Consider C99 comments on Linux with GCC:
//   mex -O CFLAGS="\$CFLAGS -std=c99" cochlea_bankfilterrectifier_hybrid.c

// Headers: --------------------------------------------------------------------
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Some macros to reduce the confusion:
#define b_in        prhs[0]
#define a_in        prhs[1]
#define X_in        prhs[2]
#define pc_in       prhs[3]
#define Y_out       plhs[0]
#define Yrect_out   plhs[1]

// For message errors
#define ERR_ID   "cochlea_bankfilterrectifier_hybrid:"
#define ERR_HEAD "cochlea_bankfilterrectifier_hybrid[mex]: "

// Prototypes: -----------------------------------------------------------------
void CoreDoubleN(double *X, mwSize MX, double *a, double *b, mwSize order, double *Z, double *Y);
void CoreDouble3(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y);
void CoreDouble5(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y);
void CoreDouble7(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y);
// Rectifier
void rectifierDouble(double *Y, double *Yrect, mwSize MY, mwSize nFilters, double *periodChar);


// Main function ----------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *a, *b, *X, *Y, *Yrect, *Z, *periodChar; 
//     mwSize *periodChar;
    mwSize order, nParam, MX, NX, Xndims;
    mwSize Mb,Nb,Ma,Na,Mpc,Npc,i,bnDims,anDims,nFilters;
    bool   allocate_ba = false;
    char   Rev[2];
    const mwSize *Xdims;
    mwSize Ydims[2];
        
    // Check number of inputs
    if (nrhs!=4){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "Wrong number of inputs : [Y,Yrect]= cochlea_bankfilterrectifier(b,a,X,periodChar)");
    }
    
    // Check input data type
    if (!mxIsDouble(b_in) || !mxIsDouble(a_in)){
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
    Mb = mxGetM(b_in); // length of filters
    Nb = mxGetN(b_in); // number of filters
    Ma = mxGetM(a_in); // length of filters
    Na = mxGetN(a_in); // number of filters
    Mpc = mxGetM(pc_in); // number of filters
    Npc = mxGetN(pc_in); // 1
    nFilters = Nb;
    
    if (Npc!=1){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "periodChar must be a column vector");
    }
    if (Mpc!=nFilters){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "periodChar length must equal the number of filters");
    }    
    if (bnDims!=2 || anDims!=2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "b and a must be 2D matrices or 1D vectors");
    }
    if (Mb!=Ma || Nb!=Na){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "b and a must be of the same size");
    }
    if (Mb<2){
        mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "filter order must be a least 1 - bad filter coefficients");
    }
    mexPrintf("%d filters of length %d (order %d)\n",Nb,Mb,Mb-1);
    
    b  = mxGetPr(b_in);
    a  = mxGetPr(a_in);
    periodChar = mxGetPr(pc_in);

    // Check normalization of filters : a[0,:] == 0
    for (i=0;i<Nb;i++){
        if (a[i*Mb]!=1){
            mexErrMsgIdAndTxt(ERR_ID "input",ERR_HEAD "a[0] must be equal to 1");
        }
    }
    order = Mb-1;

    // Allocate memory for Z (initial and final? condition)
    if ((Z = mxCalloc(order, sizeof(double))) == NULL){
        mexErrMsgIdAndTxt(ERR_ID   "NoMemory", ERR_HEAD "No memory for initial conditions.");
    }
    
    
    // Y_out must be of size (length(x),numberOfFilters)
    Ydims[0] = MX;
    Ydims[1] = nFilters;
    mexPrintf("Ydims: %d,%d\n",Ydims[0],Ydims[1]);
    // Filtering
    if (mxIsDouble(X_in)){
        Y_out = mxCreateNumericArray(2, (const) Ydims, mxDOUBLE_CLASS, mxREAL);
        Y     = mxGetPr(Y_out);
        X     = mxGetPr(X_in);
        // Filtering (band-pass filter bank)
        switch (order) {
            case 2:   CoreDouble3(X, MX, nFilters, a, b, Z, Y);  break;
            case 4:   CoreDouble5(X, MX, nFilters, a, b, Z, Y);  break;
            case 6:   CoreDouble7(X, MX, nFilters, a, b, Z, Y);  break;
        }
        // Rectifier
        Yrect_out   = mxCreateNumericArray(2, (const) Ydims, mxDOUBLE_CLASS, mxREAL);
        Yrect       = mxGetPr(Yrect_out);
        rectifierDouble(Y,Yrect,MX,nFilters,periodChar);
        
    }
    
    // Cleanup:
    mxFree(Z);
    
    return;
}


// Rectifier: ---------------------------------------------------------------
void rectifierDouble(double *Y, double *Yrect, mwSize MY, mwSize nFilters, double *periodChar)
{
    mwSize i, iY=0, iFilt, maxX=0, period;
    double maxY=0;
    
    Yrect[0] = Y[0];
    
    for (iFilt=0;iFilt<nFilters;iFilt++)
    {
        period          = (mwSize) periodChar[iFilt];
        maxX            = 0;
        maxY            = 0;
        
        for (i=0;i<MY;i++) 
        {
            if ( Y[iY]>maxY || (i-maxX)>period ) // Update when increase of Y or last update was long ago
            {
                if (Y[iY]<0){
                    Yrect[iY] = 0;
                }
                else{
                    maxX = i;
                    maxY = Y[iY];
                    Yrect[iY] = Y[iY]; // Initialize Yrect with Y at the beginning to reduce the time
                }
            } 
            else // Maintien when decrease of Y
            { 
                Yrect[iY] = maxY;
            }
            iY++;
        }
    }
}


// Filtering functions: ------------------------------------------------------

void CoreDouble3(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y)
{
	double Xi, Yi, z0, z1, a1 = a[1], a2 = a[2];
	mwSize iX=0, iY=0, iFilt=0;
    
    while (nFilters--){
        iX = 0;
        a1 = a[iFilt+1];
        a2 = a[iFilt+2];
        z0 = Z[0];
        z1 = Z[1];
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

void CoreDouble5(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y)
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
        z0 = Z[0];
        z1 = Z[1];
        z2 = Z[2];
        z3 = Z[3];
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

void CoreDouble7(double *X, mwSize MX, mwSize nFilters, double *a, double *b, double *Z, double *Y)
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
        z0 = Z[0];
        z1 = Z[1];
        z2 = Z[2];
        z3 = Z[3];
        z4 = Z[4];
        z5 = Z[5];
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

