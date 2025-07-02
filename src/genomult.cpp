#include "genomult.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//
// Compute y^T X X^T y : output is a scalar
// X = X0[1:Nindv,index:(index+num_snp-1)]
// X0 : genotype matrix of Nindv X Nsnp
// vec (matrix of dimension Nindv X 1) (usually phenotype)
double compute_yXXy (int num_snp, MatrixXdr vec){
	MatrixXdr res = MatrixXdr::Zero (num_snp, 1);
	
	if (verbose >= 3){
		cout << "***In compute_yXXy***" << endl;
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;
		cout << "means = " << means.rows() << "," << means.cols () << "\t" << means.sum()<<endl;
		cout << "stds = " << stds.rows() << "," << stds.cols () << "\t" << stds.sum()<<endl;
	}

	mm.multiply_y_pre (vec, 1, res, false);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid = means.cwiseProduct(stds);
	resid = resid * vec.sum();
	
	if (verbose >= 4)
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;

	MatrixXdr Xy(num_snp,1);
	Xy = res - resid;

	if (verbose >= 4)
		cout << "Xy = " << Xy.rows() << "," << Xy.cols () << "\t" << Xy.sum()<<endl;

	double yXXy = (Xy.array() * Xy.array()).sum();
	return yXXy;
}


double compute_yVXXVy(int num_snp){
	MatrixXdr new_pheno_sum = new_pheno.colwise().sum();
	MatrixXdr res(num_snp, 1);

	mm.multiply_y_pre(new_pheno, 1, res, false);

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid = means.cwiseProduct(stds);
	resid = resid *new_pheno_sum;
	MatrixXdr Xy(num_snp,1);
	Xy = res - resid;
	double ytVXXVy = (Xy.array()* Xy.array()).sum();
	return ytVXXVy;
}

// Compute X X^T Z : Nindv X Nz matrix
// X = X0[1:Nindv,index:(index+num_snp-1)]
// X0 : genotype matrix of Nindv X Nsnp
// Z  : Zvec (matrix of dimension Nindv X Nz) (usually random vectors)
MatrixXdr  compute_XXz (int num_snp, MatrixXdr Zvec){
	MatrixXdr res = MatrixXdr::Zero (num_snp, Nz);

	if (verbose >= 3) {		
		cout << "***In compute_XXz***" << endl;
		cout << "res [" << res.rows() << "," << res.cols () << "]\t" << res.sum()<<endl;
		cout << "Zvec [" << Zvec.rows() << "," << Zvec.cols () << "]\t" << Zvec.sum()<<endl;
	}

	mm.multiply_y_pre(Zvec, Nz, res, false);

	if (verbose >= 4) {
		cout << "res [ " << res.rows() << ", " << res.cols () << "]\t" << res.sum()<<endl;
		cout << "means [ " << means.rows() << "," << means.cols () << "]\t" << means.sum()<<endl;
		cout << "stds [ " << stds.rows() << "," << stds.cols () << "]\t" << stds.sum()<<endl;
	}

	MatrixXdr zb_sum = Zvec.colwise().sum();

	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	if (verbose >= 4)
		cout << "res [" << res.rows() << "," << res.cols () << "]\t" << res.sum()<<endl;

	MatrixXdr resid(num_snp, Nz);
	MatrixXdr inter = means.cwiseProduct(stds);
	resid = inter * zb_sum;
	MatrixXdr inter_zb = res - resid;

	if (verbose >= 4) {
		cout << "resid [" << resid.rows() << "," << resid.cols () << "]\t" << resid.sum()<<endl;
		cout << "inter_zb [" << inter_zb.rows() << "," << inter_zb.cols () << "]\t" << inter_zb.sum()<<endl;
	}

	for(int k = 0; k < Nz; k++)
		for(int j = 0; j < num_snp ; j++)
			inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
	MatrixXdr new_zb = inter_zb.transpose();
	MatrixXdr new_res(Nz, Nindv);

	mm.multiply_y_post(new_zb, Nz, new_res, false);
	if (verbose >= 4)
		cout << "new_res = " << new_res.rows() << "," << new_res.cols () << "\t" << new_res.sum()<<endl;

	MatrixXdr new_resid(Nz, num_snp);
	MatrixXdr zb_scale_sum = new_zb * means;
	new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);
		/// new zb 
	MatrixXdr temp = new_res - new_resid;

	if (verbose >= 4)
		cout << "temp = " << temp.rows() << "," << temp.cols () << "\t" << temp.sum()<<endl;

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			temp(i,j) = temp(i,j) * mask(j,0);

	if (verbose >= 3)
		cout << "temp = " << temp.rows() << "," << temp.cols () << "\t" << temp.sum()<<endl;

	return temp.transpose();
}

MatrixXdr  compute_XXUz (int num_snp){
	res.resize(num_snp, Nz);

	mm.multiply_y_pre(all_Uzb,Nz,res, false);

	MatrixXdr zb_sum = all_Uzb.colwise().sum();


	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	MatrixXdr resid(num_snp, Nz);
	MatrixXdr inter = means.cwiseProduct(stds);
	resid = inter * zb_sum;
	MatrixXdr inter_zb = res - resid;


	for(int k = 0; k < Nz; k++)
		for(int j = 0; j < num_snp ; j++)
			inter_zb(j,k) =inter_zb(j,k) * stds(j,0);
	MatrixXdr new_zb = inter_zb.transpose();
	MatrixXdr new_res(Nz, Nindv);

	mm.multiply_y_post(new_zb, Nz, new_res, false);

	MatrixXdr new_resid(Nz, num_snp);
	MatrixXdr zb_scale_sum = new_zb * means;
	new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);

	MatrixXdr temp = new_res - new_resid;

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			temp(i,j) = temp(i,j) * mask(j,0);

	return temp.transpose();
}


MatrixXdr compute_yXXy_multi (int num_snp, MatrixXdr vec, int cur_pheno_count){

	MatrixXdr pheno_sum=vec.colwise().sum();
	MatrixXdr res = MatrixXdr::Zero (num_snp, cur_pheno_count);
	
	if (verbose >= 3){
		cout << "***In compute_yXXy_multi***" << endl;
	}

	mm.multiply_y_pre (vec, cur_pheno_count, res, false);

	for(int j=0; j<num_snp; j++)
		for(int k=0; k<cur_pheno_count;k++)
				res(j,k) = res(j,k)*stds(j,0);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	// res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, cur_pheno_count);
	// resid = means.cwiseProduct(stds);
	for (int j=0; j<num_snp; j++)
		for (int k=0; k<cur_pheno_count;k++)
			resid(j, k) = means(j, k)*stds(j, k);
	resid = resid * vec.sum();
	
	if (verbose >= 4)
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;

	MatrixXdr Xy(num_snp,cur_pheno_count);
	Xy = res - resid;

	if (verbose >= 4)
		cout << "Xy = " << Xy.rows() << "," << Xy.cols () << "\t" << Xy.sum()<<endl;

	// double yXXy = (Xy.array()* Xy.array()).sum();
	Xy=Xy.array()*Xy.array();
    MatrixXdr out_temp=Xy.colwise().sum();
	return out_temp;
}


MatrixXdr compute_yVXXVy_multi(int num_snp, MatrixXdr vec, int cur_pheno_count){
	MatrixXdr new_pheno_sum = vec.colwise().sum();
	MatrixXdr res(num_snp, cur_pheno_count);

	mm.multiply_y_pre(new_pheno, cur_pheno_count, res, false);

	// res = res.cwiseProduct(stds);

	for(int j=0; j<num_snp; j++)
		for(int k=0; k<cur_pheno_count;k++)
				res(j,k) = res(j,k)*stds(j,0);

	MatrixXdr resid(num_snp, cur_pheno_count);
	// resid = means.cwiseProduct(stds);

	for (int j=0; j<num_snp; j++)
		for (int k=0; k<cur_pheno_count;k++)
			resid(j, k) = means(j, k)*stds(j, k);

	resid = resid *new_pheno_sum;
	MatrixXdr Xy(num_snp, cur_pheno_count);
	
	Xy = res - resid;
	Xy=Xy.array()*Xy.array();
    MatrixXdr out_temp=Xy.colwise().sum();
    return out_temp;
}
