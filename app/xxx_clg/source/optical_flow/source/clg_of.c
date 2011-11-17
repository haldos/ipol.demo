#include "iio.c"
#include "nrutil.c"
#include <math.h> // only for "exp"
#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"

#define EPS			1E-12	// Precision threshold for linear system solving.
#define NPRE 			2	// Number of relaxation sweeps before : : :
#define NPOST 			2	// : : : and after the coarse-grid correction is computed.
#define JROWS 			4	// Size of array for storing J coeffs.
#define JCOLS 			4

double lin2by2det(double a,double b,double c,double d);
void nrMatrixSmooth(double **nrMatrix, int rows, int cols, int filterSide, double filterSigma);
void relax(double **u, double **v, double **J[JROWS][JCOLS], int rows, int cols, double alpha);
int writeflow(double* uflow, double* vflow, int* pm, int* pn);
int calc_clg(double* prevImage, double* currImage, double* uflow, double* vflow, int pCols, int pRows, int numIterations, double alpha, double rho, double sigma);

inline double lin2by2det(double a, double b, double c, double d) {
	return a*d - b*c;
}

void computeNRImageDerivatives(double **prevFrame, double **currFrame, double **J[JROWS][JCOLS], int numRows, int numCols) {

	int i=1, j=1;

	// --- Image derivatives -----------------------------------------------------------
	double **dfdx, **dfdy, **dfdt;
	//dfdx = new double*[numRows];
	//dfdy = new double*[numRows];
	//dfdt = new double*[numRows];
	dfdx = (double **) malloc((size_t) (numRows * sizeof(double)));
	dfdy = (double **) malloc((size_t) (numRows * sizeof(double)));
	dfdt = (double **) malloc((size_t) (numRows * sizeof(double)));

	double h = 1.0 / (numCols-1);  // TO DO: check cols vs. numCols...
	double imDerStencil[] = {1.0/12/h, -8.0/12/h, 0.0/12/h, 8.0/12/h, -1.0/12/h};
	int imDerStencilLength = 5;

	for (i=1; i<numRows; i++) {
		dfdx[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
		dfdy[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
		dfdt[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
		for (j=1; j<numCols; j++) {
			dfdx[i][j] = 0.0;
			dfdy[i][j] = 0.0;

			if(i>2 && i<numRows-2) {
				int count = 0;
				for(count=0; count<imDerStencilLength; count++) {
					dfdx[i][j] += imDerStencil[count]
						//*((uchar*)(currentFrame1C->imageData + currentFrame1C->widthStep*(i+count-3)))[j-1];
						*currFrame[i+count-2][j];
				}
			}
			if (j>2 && j<numCols-2) {
				int count = 0;
				for(count=0; count<imDerStencilLength; count++) {
					dfdy[i][j] += imDerStencil[count]
						//*((uchar*)(currentFrame1C->imageData + currentFrame1C->widthStep*(i-1)))[j+count-3];
						*currFrame[i][j+count-2];
				}
			}
			//dfdt[i][j]=((uchar*)(currentFrame1C->imageData + currentFrame1C->widthStep*(i-1)))[j-1]
			//-((uchar*)(previousFrame1C->imageData + previousFrame1C->widthStep*(i-1)))[j-1];
			dfdt[i][j] = currFrame[i][j] - prevFrame[i][j];
		}
	}
	//for(i=1;i<3;i++)
	i = 2;
	for(j=1; j<numCols; j++) {
		dfdx[i][j]         = (currFrame[3][j] - currFrame[1][j]) / 2.0 / h;
		dfdx[numRows-i][j] = (currFrame[numRows-1][j] - currFrame[numRows-3][j]) / 2.0 / h;//dfdx[numRows-4][j];
	}
	i = 1;
	for (j=1; j<numCols; j++) {
		dfdx[i][j] = dfdx[2][j];
		dfdx[numRows-i][j] = dfdx[numRows-2][j];
	}

	for(i=1;i<numRows;i++)
		//for(j=1;j<3;j++)
	{
		j = 2;
		dfdy[i][j] = (currFrame[i][3]-currFrame[i][1]) / 2.0 / h;//dfdy[i][3];
		dfdy[i][numCols-j] = (currFrame[i][numCols-1] - currFrame[i][numCols-3]) / 2.0 / h;//dfdy[i][numCols-4];

		j = 1;
		dfdy[i][j] = dfdy[i][2];
		dfdy[i][numCols-j] = dfdy[i][numCols-2];
	}
	//}

	//if(!loadComputedDerivatives) {
	int k=0, l=0;
	for (i=1; i<numRows; i++)
		for (j=1; j<numCols; j++) {
			k=1; l=1;
			J[k][l][i][j] = dfdx[i][j] * dfdx[i][j];
			k=1; l=2;
			J[k][l][i][j] = dfdx[i][j] * dfdy[i][j];
			k=1; l=3;
			J[k][l][i][j] = dfdx[i][j] * dfdt[i][j];

			//k=2; l=1;
			//J[k][l][i][j]=J[1][2][i][j];
			k=2; l=2;
			J[k][l][i][j] = dfdy[i][j] * dfdy[i][j];
			k=2; l=3;
			J[k][l][i][j] = dfdy[i][j] * dfdt[i][j];
			k=3; l=3;
			J[k][l][i][j] = dfdt[i][j] * dfdt[i][j];
		}

	// memory is freed
	for (i=1; i<numRows; i++) {
		free(dfdx[i]);
		free(dfdy[i]);
		free(dfdt[i]);
	}
	free(dfdx);
	free(dfdy);
	free(dfdt);

}

void computeNRImageDerivativesHS(double **prevFrame, double **currFrame, double **J[JROWS][JCOLS], int numRows, int numCols) {

	int i=1, j=1;

	// --- Image derivatives
	double **dfdx, **dfdy, **dfdt;
	//dfdx = new double*[numRows];
	//dfdy = new double*[numRows];
	//dfdt = new double*[numRows];
	dfdx = (double **) malloc((size_t) (numRows * sizeof(double)));
	dfdy = (double **) malloc((size_t) (numRows * sizeof(double)));
	dfdt = (double **) malloc((size_t) (numRows * sizeof(double)));

	double h = 1.0 / (numCols-1);

	for (i=1; i<numRows; i++) {
		//dfdx[i] = new double[numCols];
		//dfdy[i] = new double[numCols];
		//dfdt[i] = new double[numCols];
		dfdx[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
		dfdy[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
		dfdt[i] = (double *) malloc((size_t) (numCols * sizeof(double)));
	}

	for (i=1; i<numRows-1; i++) {
		for (j=1; j<numCols-1; j++) {

			dfdx[i][j] = 0.25 / h * (prevFrame[i+1][j] - prevFrame[i][j]
								+ prevFrame[i+1][j+1]  - prevFrame[i][j+1]
								+ currFrame[i+1][j]    - currFrame[i][j]
								+ currFrame[i+1][j+1]  - currFrame[i][j+1]);
			dfdy[i][j] = 0.25 / h * (prevFrame[i][j+1] - prevFrame[i][j]
								+ prevFrame[i+1][j+1]  - prevFrame[i+1][j]
								+ currFrame[i][j+1]    - currFrame[i][j]
								+ currFrame[i+1][j+1]  - currFrame[i+1][j]);
			dfdt[i][j] = 0.25 * ( currFrame[i][j]	- prevFrame[i][j]
								+ currFrame[i+1][j]    - prevFrame[i+1][j]
								+ currFrame[i][j+1]    - prevFrame[i][j+1]
								+ currFrame[i+1][j+1]  - prevFrame[i+1][j+1]);
		}
	}

	//for(j=1; j<numCols; j++) {
	//	dfdx[numRows-1][j]=dfdx[numRows-2][j];
	//	dfdy[numRows-1][j]=dfdy[numRows-2][j];
	//	dfdt[numRows-1][j]=dfdt[numRows-2][j];
	//}

	//for(i=1; i<numRows-1; i++) {
	//	dfdx[i][numCols-1]=dfdx[i][numCols-2];
	//	dfdy[i][numCols-1]=dfdy[i][numCols-2];
	//	dfdt[i][numCols-1]=dfdt[i][numCols-2];
	//}

	for (j=1; j<numCols; j++) {
		dfdx[numRows-1][j] = 0.0;
		dfdy[numRows-1][j] = 0.0;
		dfdt[numRows-1][j] = 0.0;
	}

	for (i=1; i<numRows-1; i++) {
		dfdx[i][numCols-1] = 0.0;
		dfdy[i][numCols-1] = 0.0;
		dfdt[i][numCols-1] = 0.0;
	}

	int k=0, l=0;
	for (i=1; i<numRows; i++)
		for (j=1; j<numCols; j++) {
			k=1; l=1;
			J[k][l][i][j] = dfdx[i][j] * dfdx[i][j];
			k = 1; l = 2;
			J[k][l][i][j] = dfdx[i][j] * dfdy[i][j];
			k=1; l=3;
			J[k][l][i][j] = dfdx[i][j] * dfdt[i][j];

			//k=2; l=1;
			//J[k][l][i][j]=J[1][2][i][j];
			k=2; l=2;
			J[k][l][i][j] = dfdy[i][j]*dfdy[i][j];
			k=2; l=3;
			J[k][l][i][j] = dfdy[i][j]*dfdt[i][j];

			k=3; l=3;
			J[k][l][i][j] = dfdt[i][j]*dfdt[i][j];
		}

	// memory is free
	for(i=1; i<numRows; i++) {
		free(dfdx[i]);
		free(dfdy[i]);
		free(dfdt[i]);
	}
	free(dfdx);
	free(dfdy);
	free(dfdt);
}


void ofclg2d_coefs(double **prevFrame, double **currFrame, double **J[JROWS][JCOLS], int numRows, int numCols, double gaussSigma, double gaussRho) {

	const int nu1 = NPRE, nu2 = NPOST;

	int i = 1, j = 1;

	// --- Gaussian pre smoothing
	int filterSide = (int) (2.0 * (int)(2.5 * gaussSigma) + 1.0);

	if (filterSide >= 2) {
		fprintf(stderr, "gaussian smoothing for each frame\n");
		nrMatrixSmooth(prevFrame, numRows, numCols, filterSide, gaussSigma);
		nrMatrixSmooth(currFrame, numRows, numCols, filterSide, gaussSigma);
	}

	// --- J coeffs obtention
	int k = 0, l = 0;

	computeNRImageDerivatives(prevFrame, currFrame, J, numRows, numCols);

	//filterSide=(int)(2.0*(int)(2.5*gaussRho)+1.0);
	filterSide = 2 * ((int) (2.5 * gaussRho)) + 1;  // for rho>=0, filterSigma=2*floor(2.5*rho)+1

	if (filterSide < 3)
		filterSide = 3;

	fprintf(stderr, "local spatio temporal smoothing\n");
	if (gaussRho > 1e-7 && filterSide >= 2) {
		for (k=1; k<JROWS; k++)
			for (l=k; l<JCOLS; l++)
				nrMatrixSmooth(J[k][l], numRows, numCols, filterSide, gaussRho);
	}

}


int calc_clg(double* prevImage, double* currImage, double* uflow, double* vflow, int pCols, int pRows, int numIterations, double alpha, double rho, double sigma) {

	fprintf(stderr, "calc_clg setup variables\n");
	int i=0, j=0, numCols = pCols+1, numRows = pRows+1;

	// --- IDL variable to NR variable --------
	double **prevFrame, **currFrame;
	prevFrame = dmatrix(1, numRows-1, 1, numCols-1);
	currFrame = dmatrix(1, numRows-1, 1, numCols-1);

	for (j=1; j<numRows; j++) {
		for(i=1; i<numCols; i++) {
			prevFrame[j][i] = prevImage[i-1+(j-1) * (numCols-1)];
			currFrame[j][i] = currImage[i-1+(j-1) * (numCols-1)];
		}
	}
	fprintf(stderr, "allocating memory for u, v\n");
	double **u, **v;
	u = dmatrix(1, numRows-1, 1, numCols-1);
	v = dmatrix(1, numRows-1, 1, numCols-1);

	for (j=1; j<numRows; j++) {
		for (i=1; i<numCols; i++) {
			u[j][i] = uflow[i-1+(j-1) * (numCols-1)];
			v[j][i] = vflow[i-1+(j-1) * (numCols-1)];
		}
	}

	fprintf(stderr, "allocating memory for J matrices (derivatives)\n");
	double **J[JROWS][JCOLS];
	int k=1, l=1;
	for (k=1; k<JROWS; k++) {
		for (l=1; l<JCOLS; l++) {
			J[k][l] = dmatrix(1, numRows-1, 1, numCols-1);
		}
	}

	fprintf(stderr, "going to ofclg2d_coefs\n");
	ofclg2d_coefs(prevFrame, currFrame, J, numRows, numCols, sigma, rho);

	if (numIterations == 0)
		numIterations = (int)((numRows-1.0) * (numCols-1.0) / 8.0);

	int count = 0;
	fprintf(stderr, "goinf to relax \n");
	for (count=0; count<numIterations; count++)
		relax(u, v, J, numRows-1, numCols-1, alpha);


	// Optical flow scaling (image supposed to be rx1.0)
	fprintf(stderr, "scaling u,v \n");
	for (i=1; i<numRows; i++) {
		for (j=1; j<numCols; j++) {
			u[i][j] = (numCols - 1.0) * u[i][j];
			v[i][j] = (numCols - 1.0) * v[i][j];
		}
	}

	// --- NR results to IDL variable -------------------------------------------------
	fprintf(stderr, "post formatting\n");
	for (j=1; j<numRows; j++) {
		for (i=1; i<numCols; i++) {
			uflow[i-1+(j-1) * (numCols-1)] = u[j][i];
			vflow[i-1+(j-1) * (numCols-1)] = v[j][i];
		}
	}

	//// Flow file writing
	//WriteFlowFile(img, ".\\\out.flo");

	//ofstream out(".\\salida.txt");
	//out << numCols << endl;
	//out.close();

	// Memory leaks are bad -----------------------------------------------------------
	free_dmatrix(u, 1, numRows-1, 1, numCols-1);
	free_dmatrix(v, 1, numRows-1, 1, numCols-1);

	for (k=1; k<JROWS; k++) {
		for (l=1; l<JCOLS; l++) {
			free_dmatrix(J[k][l], 1, numRows-1, 1, numCols-1);
		}
	}
	free_dmatrix(prevFrame, 1, numRows-1, 1, numCols-1);
	free_dmatrix(currFrame, 1, numRows-1, 1, numCols-1);

	return 1;
}

//Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
//u[1..n][1..n], using the right-hand side function rhs[1..n][1..n].
// Red-black, Gauss-Seidel relaxation, Neumann forced, pointwise coupled
void relax(double **u, double **v, double **J[JROWS][JCOLS], int rows, int cols, double alpha) {

	int i, ipass, isw, j, jsw = 1;
	double hx, hy, hx2, hy2, discr, detU, detV;

	hx = 1.0 / (cols-1);
	hy = hx;  //
	hx2 = hx * hx;
	hy2 = hy * hy;

	for (ipass=1; ipass<= 2; ipass++, jsw = 3-jsw) {	//Red and black sweeps.
		isw = jsw;
		for (j=1; j<=cols; j++, isw=3-isw)
			for (i=isw; i<=rows; i+=2) {		//Gauss-Seidel formula.
				if (i>1 && i<rows) {
					if (j>1 && j<cols) {
						discr = lin2by2det(alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j],
								J[1][2][i][j],
								J[1][2][i][j],
								alpha/hx2*2.0 + alpha/hy2*2.0 + J[2][2][i][j]);

						if (abs((double)discr)>EPS) {
							detU = lin2by2det(alpha/hx2*(u[i-1][j]+u[i+1][j]) + alpha/hy2*(u[i][j-1]+u[i][j+1]) - J[1][3][i][j],
									J[1][2][i][j],
									alpha/hx2*(v[i-1][j]+v[i+1][j]) + alpha/hy2*(v[i][j-1]+v[i][j+1]) - J[2][3][i][j],
									alpha/hx2*2.0 + alpha/hy2*2.0 + J[2][2][i][j]);
							u[i][j]= detU / discr;

							detV = lin2by2det(alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j],
									alpha/hx2 * (u[i-1][j]+u[i+1][j]) + alpha/hy2*(u[i][j-1]+u[i][j+1]) - J[1][3][i][j],
									J[1][2][i][j],
									alpha/hx2 * (v[i-1][j]+v[i+1][j]) + alpha/hy2*(v[i][j-1]+v[i][j+1]) - J[2][3][i][j]);
							v[i][j] = detV / discr;

						} else {
							u[i][j] = (alpha/hx2*(u[i-1][j]+u[i+1][j]) + alpha/hy2*(u[i][j-1]+u[i][j+1]) - (J[1][2][i][j]*v[i][j]+J[1][3][i][j])) / (alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j]);
							v[i][j] = (alpha/hx2*(v[i-1][j]+v[i+1][j]) + alpha/hy2*(v[i][j-1]+v[i][j+1]) - (J[1][2][i][j]*u[i][j]+J[2][3][i][j])) / (alpha/hx2*2.0 + alpha/hy2*2.0 + J[2][2][i][j]);
						}
					} else if (j==1) {
						u[i][j] = u[i][j+1];
						v[i][j] = v[i][j+1];
					} else if (j==cols) {
						u[i][j] = u[i][j-1];
						v[i][j] = v[i][j-1];
					} else {
						printf("relax: Bug: Index error(2).\n");
						return;
					}
				} else if (i==1) {
					if(j>1 && j<cols) {
						u[i][j] = u[i+1][j];
						v[i][j] = v[i+1][j];
					} else if (j==1) {
						u[i][j] = u[i+1][j+1];
						v[i][j] = v[i+1][j+1];
					} else if (j==cols) {
						u[i][j] = u[i+1][j-1];
						v[i][j] = v[i+1][j-1];
					} else {
						printf("relax: Bug: Index error.\n");
						return;
					}
				} else if (i==rows) {
					if(j>1 && j<cols) {
						u[i][j] = u[i-1][j];
						v[i][j] = v[i-1][j];
					} else if (j==1) {
						u[i][j] = u[i-1][j+1];
						v[i][j] = v[i-1][j+1];
					} else if (j==cols) {
						u[i][j] = u[i-1][j-1];
						v[i][j] = v[i-1][j-1];
					} else {
						printf("relax: Bug: Index error(3).\n");
						return;
					}
				}
			}
	}
}

void nrMatrixSmooth(double **nrMatrix, int rows, int cols, int filterSide, double filterSigma) {

	int i=1, j=1, k=0, l=0;
	double sum = 0.0;
	double **aux;

	fprintf(stderr, "allocating aux. matrix\n");
	// --- Matrix backup --------------------------------------------------------------------------
	aux = (double **) malloc((size_t) rows * sizeof(double));
	for (i=1; i<rows; i++) {
		aux[i] = malloc((size_t) cols * sizeof(double));
		for (j=1; j<cols; j++) {
			aux[i][j] = nrMatrix[i][j];
		}
	}

	fprintf(stderr, "gaussian filter 1D\n");
	// --- Gaussian filter 1D ---------------------------------------------------------------------
	int centralCoef = filterSide / 2;  // integer division, assuming odd filterSide
	double *filter = (double *) malloc((size_t) filterSide * sizeof(double));
	for (k=0; k<filterSide; k++) {
		filter[k] = exp(-(k-centralCoef) * (k-centralCoef)/2.0/filterSigma/filterSigma); //filter = exp(-(x-N-1).^2/2/sigma^2);
		sum += filter[k];
	}

	// --- Convolution -----------------------------------------------------------------------------
	// Comment: Removing some operations for the central pixels reduced the total derivation time
	// to one half.
	fprintf(stderr, "convolution\n");
	double momentSum = 0.0;
	sum = 0.0;
	int pos = 0;
	fprintf(stderr, " convolution, 1st dim\n");
	// --- First dimension ---
	for (i=1+centralCoef; i<rows-centralCoef; i++)
		for (j=1; j<cols; j++) {  // central pixels
			nrMatrix[i][j] = 0.0;
			for (k=0; k<filterSide; k++) {
				// central pixels
				nrMatrix[i][j] += aux[i+k-centralCoef][j] * filter[k];
			}
		}

	for (i=1; i<centralCoef+1; i++) {
		for (j=1; j<cols; j++) {  // upper boundary
			sum = 0.0;
			momentSum = 0.0;
			for (k=0; k<filterSide; k++) {
				pos = i+k-centralCoef;
				if (pos>=1 && pos<rows) {
					momentSum += aux[pos][j] * filter[k];
					sum += filter[k];
				}
			}
			nrMatrix[i][j] = momentSum / sum;
		}
	}

	for (i=rows-centralCoef; i<rows; i++) {
		for (j=1; j<cols; j++) {  // lower boundary
			sum = 0.0;
			momentSum = 0.0;
			for (k=0; k<filterSide; k++) {
				pos = i+k-centralCoef;
				if (pos>=1 && pos<rows) {
					momentSum += aux[pos][j] * filter[k];
					sum += filter[k];
				}
			}
			nrMatrix[i][j] = momentSum / sum;
		}
	}

	// Matrix backup
	fprintf(stderr, "  1st dim, matrix backup\n");
	for (i=1; i<rows; i++) {
		for (j=1; j<cols; j++) {
			aux[i][j] = nrMatrix[i][j];
		}
	}
	fprintf(stderr, " convolution, 2nd dim\n");
	// --- Second dimension ---
	for (i=1; i<rows; i++)
		for (j=1+centralCoef; j<cols-centralCoef; j++) {  // central pixels
			nrMatrix[i][j] = 0.0;
			for (k=0; k<filterSide; k++) {
				// central pixels
				nrMatrix[i][j] += aux[i][j+k-centralCoef] * filter[k];
			}
		}

	for (i=1; i<rows; i++) {
		for (j=1; j<centralCoef+1; j++) {  // left boundary
			sum = 0.0;
			momentSum = 0.0;
			for (k=0; k<filterSide; k++) {
				pos = j+k-centralCoef;
				if (pos>=1 && pos<cols) {
					momentSum += aux[i][pos] * filter[k];
					sum += filter[k];
				}
			}
			nrMatrix[i][j] = momentSum / sum;
		}
	}

	for (i=1; i<rows; i++) {
		for (j=cols-centralCoef; j<cols; j++) {  // right boundary
			sum = 0.0;
			momentSum = 0.0;
			for (k=0; k<filterSide; k++) {
				pos = j+k-centralCoef;
				if (pos>=1 && pos<cols) {
					momentSum += aux[i][pos] * filter[k];
					sum += filter[k];
				}
			}
			nrMatrix[i][j] = momentSum / sum;
		}
	}

	// --- Freeing memory --------------------------------------------------------
	for (i=1; i<rows; i++)
		free(aux[i]);
	free(filter);
	free(aux);
}

/*
 * Documentation here... soon :)
 */
int main(int argc, char *argv[]) {
	if (argc != 8) {
		printf("Usage: %s input_image_1 input_image_2 alpha rho sigma numIt output\n", argv[0]);
	} else {
		int w, h, pixeldim;

		// TODO hard coded test...
		// w = 5;
		// h = 5;
		// pixeldim = 1;
		//float im1[25] = {0.0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 255, 255, 0, 0,   0, 0, 0, 0, 0};
		//float im2[25] = {0.0, 0, 0, 0, 0,   0, 255, 255, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0, 0};

		// Load specified images
		float *im1 = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
		float *im2 = iio_read_image_float_vec(argv[2], &w, &h, &pixeldim);
		fprintf(stderr, "Two images loaded:\n\tim,: %dx%d image with %d channel(s)\n\tim2, %dx%d image with %d channel(s)\n", w, h, pixeldim, w, h, pixeldim);
                
                float *im1_gray = malloc(w*h*sizeof(float));
                float *im2_gray = malloc(w*h*sizeof(float));
                if (pixeldim==3){
                    int z;
                    for(z=0;z<w*h;z++){
                        im1_gray[z] =  (6968*im1[3*z] + 23434*im1[3*z + 1] + 2366*im1[3*z + 2])/32768;
                        im2_gray[z] =  (6968*im2[3*z] + 23434*im2[3*z + 1] + 2366*im2[3*z + 2])/32768;
                    }
                } else {
                    im1_gray = im1;
                    im2_gray = im2;
                }

		// parameters for CLG optical flow calculation
		// TODO read from console
		double alpha = (double)atof(argv[3]);//1.0;	//1150.0//2000.0//800.0// //pow(75.0,2)/3.0	// Alpha for yosemite seq. in Bruhn2005 = 800
		double rho = (double)atof(argv[4]); //0.0;	//5.5//2.25//1.0//0.0//6.75//0.0//1.8//
		double sigma = (double)atof(argv[5]); //0.0;	//2.0//1.50//2.54//1.0//1.95// 2.6 // 1.5
		int numIt = atoi(argv[6]); //50;

		// TODO make validation in order to check compatible kernel sizes with image size
		int filterSideSigma = (int) (2.0 * (int)(2.5 * sigma) + 1.0);
		int filterSideRho = 2 * ((int) (2.5 * rho)) + 1;

		// CLG optical flow computation
		// TODO jfd check this please ;)
		// (u, v) are initialized as the difference between images 1 and 2
		fprintf(stderr, "allocating memory for data structures\n");
		double* u;
		double* v;
		double* i1;
		double* i2;
		u = (double *) malloc((size_t) h * w * sizeof(double));
		v = (double *) malloc((size_t) h * w * sizeof(double));
		i1 = (double *) malloc((size_t) h * w * sizeof(double));
		i2 = (double *) malloc((size_t) h * w * sizeof(double));

		fprintf(stderr, "copying and casting data for double precision\n");
		int i = 0;
		int j = 0;
		for (i=0; i<w; i++) {
			for (j=0; j<h; j++) {
				u[w*j+i] = ((double) im2_gray[w*j+i]) - ((double) im1_gray[w*j+i]);
				v[w*j+i] = ((double) im2_gray[w*j+i]) - ((double) im1_gray[w*j+i]);
				i1[w*j+i] = ((double) im1_gray[w*j+i]);
				i2[w*j+i] = ((double) im2_gray[w*j+i]);
			}
		}

		fprintf(stderr, "call to calc_clg\n");
		int res = calc_clg(i2, i1, u, v, w, h, numIt, alpha, rho, sigma);

		// here we have flow already calculated and saved in two arrays
		// "u" and "v". now we save it in a two channels tiff file:

		// TODO debug OF fprinting
                /*
		for (i=0; i<w; i++) {
			for (j=0; j<h; j++) {
				fprintf(stderr, "u,v[%i,%i] = %g,\t%g\n", i, j, u[w*j+1], v[w*j+i]);
			}
		}
                 * */

		// Saving image
 
                float *uf = malloc(w * h * sizeof * uf);
                float *vf = malloc(w * h * sizeof * vf);
		float *x = malloc(w * h * 2 * sizeof * x);
		for (i=0; i<w; i++) {
			for (j=0; j<h; j++) {
                            uf[w*j+i] = (float) u[w*j+i];
                            vf[w*j+i] = (float) v[w*j+i];
                            x[2*(w*j+i)] = uf[w*j+i];
                            x[2*(w*j+i)+1] = vf[w*j+i];
			}
		}
		iio_save_image_float_vec(argv[7], x, w, h, 2);
                //iio_save_image_float_vec("../images/u.png", uf, w, h, 1);
                //iio_save_image_float_vec("../images/v.png", vf, w, h, 1);

		//free memory
		free(im1);
		free(im2);
                free(im1_gray);
                free(im2_gray);
		free(x);
		free(u);
		free(v);
                free(uf);
                free(vf);
        	fprintf(stderr, "clg-optical flow computation done.\n");
		return 0;
	}
}
