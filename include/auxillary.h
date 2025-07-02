#ifndef AUXILLARY_H
#define AUXILLARY_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <vector>

#ifdef USE_DOUBLE
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
#else
	typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
#endif

using namespace std;

extern int verbose;
extern int Nenv;
extern int Nbin; 
extern bool Annot_x_E;
extern vector<vector<bool> > annot_bool;
extern vector<int> len;
extern int Nsnp;
extern int Nsnp_annot;
extern vector<vector<int> > jack_bin;
extern vector<int>  snp_to_jackblock;
extern int Njack;
extern int step_size;
extern int step_size_rem;

extern vector<vector<int> > read_bin;
extern int Nreadblocks; 
extern vector<int> snp_to_read_block;

extern vector<int> SNP_annot;

extern bool gen_by_env;
extern bool cov_add_intercept;
extern bool add_env_to_cov;
extern MatrixXdr Enviro;
extern MatrixXdr pheno;
extern MatrixXdr mask;
extern MatrixXdr covariate;  

extern MatrixXdr new_pheno;
extern int phenocount;

int read_env (int Nind, std::string filename);
int read_cov (int Nind, std::string filename) ;

void read_pheno(int Nind, std::string filename);
int count_pheno(std::string filename);

void read_annot_1col (std::string filename);
void read_annot (std::string filename) ;

#endif
