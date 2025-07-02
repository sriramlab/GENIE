#ifndef STATSFN_H

#define STATSFN_H

#include <iostream>

using namespace std;

namespace statsfn {
	#ifdef USE_DOUBLE
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
	#else
		typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
	#endif
	
   // Compute jackknife SEs
    // jack: matrix with statistics along rows and jackknife estimates along columns
    MatrixXdr jack_se (MatrixXdr& jack){
    	int nrows = jack.rows();
	    int ncols = jack.cols();
    	MatrixXdr sum_row = jack.rowwise().mean();
	    MatrixXdr SEjack;
    	SEjack = MatrixXdr::Zero(nrows,1);
	    double temp_val = 0;
    	for (int i = 0 ; i < nrows ; i++){
	    	for (int j = 0 ; j < ncols ; j++){
		    	temp_val = jack(i,j)-sum_row(i);
			    temp_val= temp_val* temp_val;
    			SEjack(i,0)+=temp_val;
	    	}
    		SEjack(i,0) = SEjack(i,0) * (Njack - 1) / Njack;
	    	SEjack(i,0) = sqrt(SEjack(i,0));
    	}

	    return SEjack;
    } 
}

#endif
