
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void computeQ_aQuasi(const mxArray *img, const mxArray *guidance, mwSize nPixels, mwSize imWidth, mwSize imHeight,
					int aQuasiWindow, double aQuasiSigma, double aQuasiQuantile, 
					mxUint32 *pointerRows, mxUint32 *pointerColse);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
 {
	// check for proper number of input and output arguments
	if(nrhs!=5){
		mexErrMsgTxt("5 input arguments are required!");
	}
		
	if(nlhs<2 || nlhs >3){
		mexErrMsgTxt("Only two or three output arguments are allowed!");
	}
		
	const mxArray *img = prhs[0];
	const mxArray *guidance = prhs[1];
	int aQuasiWindow = *(mxGetPr(prhs[2]));
	double aQuasiSigma = *(mxGetPr(prhs[3]));
	double aQuasiQuantile = *(mxGetPr(prhs[4]));

	mwSize nDimImg = mxGetNumberOfDimensions(img);
	const mwSize* imDimImg = mxGetDimensions(img);	
	mwSize nDimGuide = mxGetNumberOfDimensions(guidance);
	const mwSize* imDimGuide = mxGetDimensions(guidance); 

	if (nDimImg != nDimGuide){
		mexErrMsgTxt("The number of dimensions of the image and the guidance does not fit!");
	}

	if (nDimImg != 2){
		mexErrMsgTxt("Method requires a 2-D image!");
	}

	for (int i = 0; i < nDimImg; i++){
		if(imDimImg[i] != imDimGuide[i]){
			mexErrMsgTxt("The size of the image and the guidance does not fit!");
		}
	}

	int imWidth,imHeight,nPixels;
	imWidth = imDimImg[1];
	imHeight = imDimImg[0];
	nPixels = imWidth*imHeight;	
	mxArray *lookUpRows = mxCreateNumericMatrix(nPixels,1,mxUINT32_CLASS,mxREAL);
	mxArray *lookUpCols = mxCreateNumericMatrix(nPixels,1,mxUINT32_CLASS,mxREAL);
	mxArray *result = mxCreateNumericMatrix(imHeight,imWidth,mxUINT32_CLASS,mxREAL);
	mxArray *Q = mxCreateNumericMatrix(imHeight,imWidth,mxUINT32_CLASS,mxREAL);

	mxUint32 *pointerRows = mxGetUint32s(lookUpRows);
	mxUint32 *pointerCols = mxGetUint32s(lookUpCols);

	computeQ_aQuasi(img, guidance, nPixels,imWidth,imHeight, aQuasiWindow, aQuasiSigma, aQuasiQuantile, pointerRows,pointerCols);



	plhs[0] = lookUpRows;
	plhs[1] = lookUpCols;

}

void computeQ_aQuasi(const mxArray *img, const mxArray *guidance, mwSize nPixels, mwSize imWidth,mwSize imHeight, 
			int aQuasiWindow, double aQuasiSigma, double aQuasiQuantile, 
			mxUint32 *pointerRows, mxUint32 *pointerCols)
{
	
	int *position = new int[aQuasiWindow*aQuasiWindow];;
	double *weights = new double[aQuasiWindow*aQuasiWindow];
	double *intensity = new double[aQuasiWindow*aQuasiWindow];
	int *permutation = new int[aQuasiWindow*aQuasiWindow];

	int halfWindow = aQuasiWindow/2;

	double *pI = mxGetPr(img);
	double *pG = mxGetPr(guidance);

	mxUint32 iter = 0;

	for(mwSize x = 0; x < imWidth; x++) {	
		for(mwSize y = 0;  y < imHeight; y++) {

			pointerRows[iter] = iter+1;				

			int offset = (x)*imHeight + (y);	

			int count = 0;
			double weights_total = 0;	

			for (int wx = -halfWindow; wx <= halfWindow; wx++) {					
				if (x+wx < 0) continue;
				if (x+wx >= imWidth) continue;	

				for (int wy = -halfWindow; wy <= halfWindow; wy++) {				
					if (y+wy < 0) continue;
					if (y+wy >= imHeight) continue;

					permutation[count] = count;
					
					position[count] = (x+wx)*imHeight + (y+wy);					
					intensity[count] = pI[position[count]];	
					
					weights[count] = exp(-
						((pG[offset] - pG[position[count]])*(pG[offset] - pG[position[count]]))
						/ (2*aQuasiSigma*aQuasiSigma));
					weights_total += weights[count];

					count++;
				}
				
			}

			
			std::sort(permutation, permutation+count, [&](int l, int r)
			{
				return intensity[l] < intensity[r];
			});

			//qsort(permutation, count, sizeof(int), aQuasiCmp);

			double weightQuantile = weights_total*aQuasiQuantile;
			double weights_cumulative = 0;
			// for (auto w_position = 0; weights_cumulative < weightQuantile; ++w_position) {
			int w_position = -1;

			while (weights_cumulative <= weightQuantile) {
				w_position++;
				weights_cumulative += weights[permutation[w_position]];
			}
			
			int position_of_quantile = position[permutation[w_position]];

			pointerCols[iter] = position_of_quantile+1;			
			iter++;		
			
			
		}
	
	}



	delete[] position;
	delete[] weights;
	delete[] intensity;
	delete[] permutation;



}
 
