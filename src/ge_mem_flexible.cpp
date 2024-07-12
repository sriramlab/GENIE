#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector> 
#include <random>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include "time.h"

#include "genotype.h"
#include "arguments.h"
#include "storage.h"
#include "matmult.h"
#include "io.h"
#include "std.h"
#include "functions.h"
#include "vectorfn.h"

// #if SSE_SUPPORT == 1
// 	#define fastmultiply fastmultiply_sse
// 	#define fastmultiply_pre fastmultiply_pre_sse
// #else
// 	#define fastmultiply fastmultiply_normal
// 	#define fastmultiply_pre fastmultiply_pre_normal
// #endif

using namespace Eigen;
using namespace std;

// Storing in RowMajor Form
#ifdef USE_DOUBLE
int use_double = 1;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;
#else
int use_double = 0;
typedef Matrix<float, Dynamic, Dynamic, RowMajor> MatrixXdr;
#endif

class data {
	public:
		MatrixXdr gen;
		int index;
};

// When true, writes to a file after pass 1 
// Reads from file in pass 2
bool opt1 = true;

// Reads SNPs in block that can be controlled based on memory constraints
int opt2 = 0 ;
int mem_Nsnp = -1;
int mem_alloc = -1;

//ENV
// Partition the GxE component with respect to the annotations in the annotation file
bool Annot_x_E = false;

// Add environment variables to covariates
bool add_env_to_cov = true;

int nongen_Nbin = 0;
int Nenv = 0;
MatrixXdr Enviro;
MatrixXdr point_est_adj_gxe;
MatrixXdr jack_adj_gxe;


//Intermediate Variables
int blocksize;
int hsegsize;
double *partialsums;
double *sum_op;		
double *yint_e;
double *yint_m;
double **y_e;
double **y_m;


struct timespec t0;

clock_t total_begin = clock();
MatrixXdr pheno;
MatrixXdr mask;
MatrixXdr covariate;  
MatrixXdr Q;
MatrixXdr v1; //W^ty
MatrixXdr v2;            //QW^ty
MatrixXdr v3;    //WQW^ty
MatrixXdr new_pheno;

genotype g;
genotype g1;
genotype g2;
MatrixXdr geno_matrix; //(p,n)
genotype* Geno;
int MAX_ITER;
int k,p,n;
int k_orig;

double score;

MatrixXdr c; //(p,k)
MatrixXdr x; //(k,n)
MatrixXdr v; //(p,k)
MatrixXdr means; //(p,1)
MatrixXdr stds; //(p,1)
MatrixXdr sum2;
MatrixXdr sum;  

Eigen::RowVectorXd means_na;
Eigen::RowVectorXd stds_na;
////////
//related to phenotype	
double y_sum; 
double y_mean;

options command_line_opts;

bool debug = false;
bool check_accuracy = false;
bool var_normalize = false;
double convergence_limit;
bool memory_efficient = false;
bool missing = false;
bool fast_mode = true;
bool text_version = false;
bool use_cov = false; 


//// jackknife index wich are computed based on annotation file
MatrixXdr dic_index;

// Not used
MatrixXdr jack_bin_size;

// Not used
vector<int> Annot;


int jack_scheme;
int Njack = 1000;
int jack_size;
bool use_ysum;
bool keep_xsum;

int Nreadblocks; 

int Nbin = 160;
// Total number of bins
// Number of VCs = Number of bins  + 1 (sigma_e)
// (1 for G model with single annotation)
int T_Nbin ;
int Nz = 10;
int Ncov;
int Nindv_mask;
int NC;
///////
MatrixXdr jack;
MatrixXdr point_est;
MatrixXdr enrich_jack;
MatrixXdr enrich_point_est;


//define random vector z's
MatrixXdr  all_zb;
MatrixXdr  all_Uzb;
MatrixXdr res;
MatrixXdr XXz;
MatrixXdr Xy;
MatrixXdr yXXy;
MatrixXdr jack_yXXy;


///
//Matrix<int, Dynamic, Dynamic, RowMajor> gen;
MatrixXdr gen;
bool read_header;
//read variables
unsigned char mask2;
int wordsize;
unsigned int unitsperword;
int unitsize;
int nrow, ncol;
unsigned char *gtype;
int Nsnp;

// Sum of number of SNPs assigned to each annotation 
// Can be greater than or less than or equal to the number of SNPs in the bim file
int Nsnp_annot = 0;

int Nindv;
bool **bin_annot;
int step_size;
int step_size_rem;
vector<vector<bool> > annot_bool;

// Number of SNPs in each bin
vector<int> len;

// Number of SNPs in each bin in a jackknife block
// Njack X (Total number of annotations (depends on the specific model used)
vector<vector<int> > jack_bin;
vector<int>  jack_block_size;
vector<int>  snp_to_jackblock;


vector<vector<int> > read_bin;
vector<int> jack_to_read_block;
vector<int> snp_to_read_block;

vector<vector<int> > mem_bin;
vector<int>  mem_block_size;

vector <data> allgen;
vector <genotype> allgen_mail;
int global_snp_index;
bool use_mailman = true;

///reading single col annot
vector <int> SNP_annot;
bool use_1col_annot = false;


///Variables for reg out cov on both side of LM
bool both_side_cov = true;
MatrixXdr UXXz;
MatrixXdr XXUz;
MatrixXdr Xz;
MatrixXdr trVK;

//CHANGE(10/20)
bool hetero_noise;
//CHANGE (2/27)
bool gen_by_env;
bool cov_add_intercept;
int verbose;
bool trace;
int nthreads;
MatMult mm;
std::ofstream outfile;
std::ofstream trace_file;
std::ofstream meta_file;

bool use_summary_genotypes = false;
string wgxsumfile;
string xsumfilepath;
string xsum_path;
string wgxsum_path;
string ysum_path;
std::ofstream xsum_ofs;
std::ofstream wgxsum_ofs;
std::ifstream wgxsum_ifs;
std::ifstream xsum_ifs;
std::ofstream ysum_ofs;
std::ifstream ysum_ifs;
string prefix;

int seed;
std::mt19937 seedr;
std::uniform_real_distribution<> udis (0,1);

std::istream& newline(std::istream& in)
{
	if ((in >> std::ws).peek() != std::char_traits<char>::to_int_type('\n')) {
		in.setstate(std::ios_base::failbit);
	}
	return in.ignore();
}


// Read environmental variable file
// Inputs: Number of individuals, filename
// Return number of environments
int read_env (int Nind, std::string filename){
	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	std::string line;
	std::istringstream in;
	int covIndex = 0;
	std::getline(ifs,line);
	in.str(line);
	string b;
	vector<vector<int> > missing;
	int covNum = 0;

	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
			missing.push_back(vector<int>()); //push an empty row  
			covNum++;
		}
	}
	vector<double> cov_sum(covNum, 0);
	if (gen_by_env == true) {
		Enviro = MatrixXdr::Zero(Nind, covNum);
		cout<< "Reading in "<< covNum << " environmental variables ..." << endl;
	}

	int j = 0;
	while(std::getline(ifs, line)){
		in.clear();
		in.str(line);
		string temp;
		in>>temp; in>>temp; //FID IID 
		for(int k = 0; k < covNum; k++){

			in>>temp;
			if(temp=="NA")
			{
				mask(j,0) = 0;
				continue;
			}
			double cur = atof(temp.c_str());
			if(cur==-9)
			{
				// Need to mask this case? 
				continue;
			} else {
				if (gen_by_env == true) {
					cov_sum[k]= cov_sum[k]+ cur;
					Enviro(j,k) = cur;
				}
			}
		}
		j++;
	}
	return covNum;
}

// Read covariate file. 
// Adds environmental variables to covariates if add_env_to_cov == True
// Adds intercept to covariates if cov_add_intercept == True
// Inputs: Number of individuals, filename
// Return number of covariates
int read_cov (int Nind, std::string filename) {
	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	std::string line;
	std::istringstream in;
	int covIndex = 0;
	std::getline(ifs,line);
	in.str(line);
	string b;
	vector<vector<int> > missing;
	int covNum = 0;
	vector<double> cov_sum;
	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
			missing.push_back(vector<int>()); //push an empty row  
			covNum++;
		}
	}

	if(covNum != 0)
		cov_sum.resize(covNum , 0);

	if (cov_add_intercept == true) {
		if(add_env_to_cov == true)
			covariate.resize(Nind, covNum + Nenv + 1);
		else
			covariate.resize(Nind, covNum + 1); 
	} else {
		if(add_env_to_cov == true)
			covariate.resize(Nind, covNum + Nenv);
		else
			covariate.resize(Nind, covNum);  
	}
	cout<< "Reading in "<< covNum << " covariates ..." << endl;

	int j = 0;
	while(std::getline(ifs, line)){
		in.clear();
		in.str(line);
		string temp;
		in>>temp; in>>temp; //FID IID 
		for(int k = 0; k < covNum; k++){

			in>>temp;
			if (temp == "NA")
			{
				missing[k].push_back(j);
				continue;
			}
			double cur = atof(temp.c_str());
			if (cur == -9)
			{
				missing[k].push_back(j);
				continue;
			}else{
				cov_sum[k]= cov_sum[k]+ cur;
				covariate(j,k) = cur;
			}
		}
		j++;
	}

	//compute cov mean and impute 
	for (int a = 0; a < covNum ; a++)
	{
		int missing_num = missing[a].size();
		cov_sum[a] = cov_sum[a] / (Nind - missing_num);

		for(int b = 0; b < missing_num; b++)
		{
			int index = missing[a][b];
			covariate(index, a) = cov_sum[a];
		}
	}

	if (gen_by_env == true) {
		if (verbose >= 1) {
			cout << "Shape of Env = " << Enviro.rows() << " " << Enviro.cols() << endl;
			cout << "Shape of Cov = " << covariate.rows() << " " << covariate.cols() << endl;
			cout << "Number of covariates = " << covNum<< ", number of environments = " << Nenv << endl;
		}
		if(add_env_to_cov == true){
			for(int i = 0 ; i < Nenv ; i++){
				covariate.col(covNum + i) = Enviro.col(i);
			}
		}
	} else {
		if (cov_add_intercept == true) {
			//adding col of all ones to covariates
			for (int i = 0 ; i < Nind ; i++)
				covariate (i,covNum + Nenv) = 1;
			return covNum + 1;
		} 
	}
	if (cov_add_intercept == true) {
		//adding col of all ones to covariates
		for (int i = 0 ; i < Nind ; i++)
			covariate(i,covNum + Nenv) = 1;
	}

	if (cov_add_intercept == true) {
		if(add_env_to_cov == true)
			return covNum + Nenv + 1;
		else
			return covNum + 1; 
	} else {
		if(add_env_to_cov == true)
			return covNum + Nenv;
		else
			return covNum;  
	}
}


void initial_var(){
	p = g.Nsnp;
	n = g.Nindv;


	c.resize(p,k);
	x.resize(k,n);
	v.resize(p,k);
	sum2.resize(p,1);
	sum.resize(p,1);


	if(!fast_mode && !memory_efficient){
		geno_matrix.resize(p,n);
		g.generate_eigen_geno(geno_matrix,var_normalize);
	}

	//TODO: Initialization of c with gaussian distribution
	c = MatrixXdr::Random(p,k);


	// Initial intermediate data structures
	blocksize = k;
	hsegsize = g.segment_size_hori;        // = log_3(n)
	int hsize = pow(3,hsegsize);
	int vsegsize = g.segment_size_ver;              // = log_3(p)
	int vsize = pow(3,vsegsize);

	partialsums = new double [blocksize];
	sum_op = new double[blocksize];
	yint_e = new double [hsize * blocksize];
	yint_m = new double [hsize * blocksize];
	memset (yint_m, 0, hsize * blocksize * sizeof(double));
	memset (yint_e, 0, hsize * blocksize * sizeof(double));

	y_e  = new double*[g.Nindv];
	for (int i = 0 ; i < g.Nindv ; i++) {
		y_e[i] = new double[blocksize];
		memset (y_e[i], 0, blocksize * sizeof(double));
	}

	y_m = new double*[hsegsize];
	for (int i = 0 ; i < hsegsize ; i++)
		y_m[i] = new double[blocksize];
}


void read_pheno2(int Nind, std::string filename){
	ifstream ifs(filename.c_str(), ios::in); 

	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}

	std::string line;
	std::istringstream in;  
	int phenocount = 0; 
	//read header
	std::getline(ifs,line); 
	in.str(line); 
	string b; 
	while(in>>b)
	{
		if(b!="FID" && b !="IID")
			phenocount++; 
	}
	pheno.resize(Nind, phenocount);
	mask.resize(Nind, phenocount);
	int i = 0;  
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line); 
		string temp;
		//fid,iid
		//todo: fid iid mapping; 
		//todo: handle missing phenotype
		in>>temp; in>>temp; 
		for(int j = 0; j < phenocount ; j++) {
			in>>temp;
			double cur = atof(temp.c_str());
			if(temp=="NA" || cur==-9){
				pheno(i,j) = 0;
				mask(i,j) = 0;
			}
			else{
				pheno(i,j) = atof(temp.c_str());
				mask(i,j) = 1;
			}
		}
		i++;
	}
}

// Compute y^T X X^T y : output is a scalar
// X = X0[1:Nindv,index:(index+num_snp-1)]
// X0 : genotype matrix of Nindv X Nsnp
// vec (matrix of dimension Nindv X 1) (usually phenotype)
double compute_yXXy(int num_snp, MatrixXdr vec){
	MatrixXdr res = MatrixXdr::Zero (num_snp, 1);
	
	if (verbose >= 3){
		cout << "***In compute_yXXy***" << endl;
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;
		cout << "means = " << means.rows() << "," << means.cols () << "\t" << means.sum()<<endl;
		cout << "stds = " << stds.rows() << "," << stds.cols () << "\t" << stds.sum()<<endl;
	}
	// if(use_mailman == true)
	//         multiply_y_pre_fast(vec,1,res,false);
	// else
	//          res = gen * vec;
	mm.multiply_y_pre(vec, 1, res, false);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid = means.cwiseProduct(stds);
	resid = resid *vec.sum();
	
	if (verbose >= 4)
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;

	MatrixXdr Xy(num_snp,1);
	Xy = res - resid;

	if (verbose >= 4)
		cout << "Xy = " << Xy.rows() << "," << Xy.cols () << "\t" << Xy.sum()<<endl;

	double yXXy = (Xy.array()* Xy.array()).sum();
	return yXXy;
}

double compute_yVXXVy(int num_snp){
	MatrixXdr new_pheno_sum = new_pheno.colwise().sum();
	MatrixXdr res(num_snp, 1);

		// if(use_mailman == true)
		//         multiply_y_pre_fast(new_pheno,1,res,false);
		// else
		//          res = gen * new_pheno;
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
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;
		cout << "Zvec = " << Zvec.rows() << "," << Zvec.cols () << "\t" << Zvec.sum()<<endl;
	}
		// if(use_mailman == true)
		//         multiply_y_pre_fast(Zvec,Nz,res, false);
		// else
		//         res = gen * Zvec;
	mm.multiply_y_pre(Zvec, Nz, res, false);

	if (verbose >= 4) {
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;
		cout << "means = " << means.rows() << "," << means.cols () << "\t" << means.sum()<<endl;
		cout << "stds = " << stds.rows() << "," << stds.cols () << "\t" << stds.sum()<<endl;
	}

	MatrixXdr zb_sum = Zvec.colwise().sum();

	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	MatrixXdr resid(num_snp, Nz);
	MatrixXdr inter = means.cwiseProduct(stds);
	resid = inter * zb_sum;
	MatrixXdr inter_zb = res - resid;

	if (verbose >= 4) {
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;
		cout << "inter_zb = " << inter_zb.rows() << "," << inter_zb.cols () << "\t" << inter_zb.sum()<<endl;
	}

	for(int k = 0; k < Nz; k++)
		for(int j = 0; j < num_snp ; j++)
			inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
	MatrixXdr new_zb = inter_zb.transpose();
	MatrixXdr new_res(Nz, Nindv);


		// if(use_mailman == true)
		//         multiply_y_post_fast(new_zb, Nz, new_res, false);
		// else
		//         new_res = new_zb * gen;
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

	// if(use_mailman == true)
	//         multiply_y_pre_fast(all_Uzb,Nz,res, false);
	// else
	//         res = gen * all_Uzb;
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


	// if(use_mailman == true)
	//         multiply_y_post_fast(new_zb, Nz, new_res, false);
	// else
	//         new_res = new_zb * gen;
	mm.multiply_y_post(new_zb, Nz, new_res, false);

	MatrixXdr new_resid(Nz, num_snp);
	MatrixXdr zb_scale_sum = new_zb * means;
	new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


	/// new zb 
	MatrixXdr temp = new_res - new_resid;

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			temp(i,j) = temp(i,j) * mask(j,0);


	return temp.transpose();
}


MatrixXdr  compute_Xz (int num_snp){
	MatrixXdr new_zb= MatrixXdr::Random(Nz,num_snp);
	new_zb = new_zb * sqrt(3);

	MatrixXdr new_res(Nz, Nindv);         

	// if(use_mailman == true)
	//         multiply_y_post_fast(new_zb,Nz,new_res, false);
	// else
	//         new_res = new_zb * gen;
	mm.multiply_y_post(new_zb,Nz,new_res, false);

	MatrixXdr new_resid(Nz, num_snp);
	MatrixXdr zb_scale_sum = new_zb * means;
	new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);

	/// new zb 
	MatrixXdr temp = new_res - new_resid;

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			temp(i,j) = temp(i,j) * mask(j,0);

	return temp.transpose();
}


void setup_read_blocks ()  {

	int snpindex = 0;
	int blockindex = 0 ;

	jack_to_read_block.resize (Njack);		
	snp_to_read_block.resize (Nsnp);

	for (int i = 0; i < Njack; i++) { 
		int jack_Nsnp = jack_block_size[i];
		int read_max_Nsnp = jack_Nsnp > mem_Nsnp ? mem_Nsnp: jack_Nsnp;
		int num_blocks = jack_Nsnp/read_max_Nsnp;

		jack_to_read_block[i] = num_blocks;
		for (int j = 0 ; j < num_blocks * read_max_Nsnp; j++) {
			snp_to_read_block[snpindex] = blockindex + j/read_max_Nsnp;
			snpindex ++ ;
		}
		for (int j = num_blocks * read_max_Nsnp; j < jack_Nsnp; j++){
			snp_to_read_block[snpindex] = blockindex + num_blocks - 1;
			snpindex ++ ;
		}	
		blockindex += num_blocks;
	}
	Nreadblocks = blockindex; 

	if (verbose >= 2) { 
		cout << "Number of read blocks = " << Nreadblocks << endl;

		if (verbose >= 4){
			vectorfn::printvector (jack_to_read_block); cout << endl;
			vectorfn::printvector (snp_to_read_block); cout << endl;
		}
	}
}

void read_annot (string filename){
	vector<bool> snp_annot;

	ifstream inp(filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int j = 0 ;
	int linenum = 0 ;
	int num_parti;
	stringstream check1(line);
	string intermediate;
	vector <string> tokens;
	while(std::getline (inp, line)){
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		if (line.empty())
			continue;
		j++;

		stringstream check1(line);
		string intermediate;
		vector <string> tokens;
				// Tokenizing w.r.t. space ' ' 
		while(getline(check1, intermediate, ' '))
		{
			tokens.push_back(intermediate);
		}
		if(linenum == 0){
			num_parti = tokens.size();
			Nbin = num_parti;
			if(Annot_x_E == false)
				snp_annot.resize(Nbin + Nenv,0);
			else
				snp_annot.resize(Nbin + (Nenv * Nbin),0);

			len.resize(num_parti,0);
		}
		int index_annot = 0;
		for(int i = 0; i < tokens.size(); i++){
			snp_annot[i]=0;
			if (tokens[i]=="1"){
				len[i]++;
				snp_annot[i]=1;
			}
		}
		if(Annot_x_E == false){
			for(int i = 0 ; i < Nenv ; i++)
				snp_annot[Nbin + i]=1;   /// need to modify for excluding snps
		}else{
			for(int i = 0 ; i < Nenv ; i++)
				for(int j = 0 ; j < Nbin ; j++)
					snp_annot[Nbin + (i * Nbin) + j]=snp_annot[j];
		}
		annot_bool.push_back(snp_annot);
		linenum++;
	}

	if(Nsnp != linenum){
		exitWithError ("Number of rows in bim file and annotation file does not match");
		exit (1);
	}

	Nsnp = linenum;
	//cout << "Total number of SNPs : " << Nsnp << endl;
	cout << "Number of annotations in annotation file = " << num_parti << endl;
	int selected_snps = 0;
	for (int i = 0 ; i < num_parti ; i++){
		cout << "Number of SNPs in annotation " << i << " = " << len[i] <<endl;
		selected_snps += len[i];
	}

	cout << "Number of SNPs selected according to annotation file = " <<selected_snps << endl;
	Nsnp_annot = selected_snps;


/*
	step_size = Nsnp / Njack;
	step_size_rem = Nsnp%Njack;
	cout << "Number of SNPs per jackknife block : " << step_size << endl;   
*/

	int Total_Nbin;
	if(Annot_x_E == false){
		jack_bin.resize(Njack, vector<int>(Nbin + Nenv + Nenv,0));
		read_bin.resize(Nreadblocks, vector<int>(Nbin + Nenv + Nenv,0));
		Total_Nbin = Nbin + Nenv;
	}else{
		jack_bin.resize(Njack, vector<int>(Nbin + (Nenv * Nbin) + Nenv,0));
		read_bin.resize(Nreadblocks, vector<int>(Nbin + (Nenv * Nbin) + Nenv,0));
		Total_Nbin = Nbin + (Nenv * Nbin);
	}
	for (int i = 0 ; i < Nsnp ; i++)
		for(int j = 0 ; j < Total_Nbin ; j++)
			if (annot_bool[i][j]==1){
				/*
				temp = i / step_size;
				if (temp>=Njack)
					temp = Njack - 1;
				jack_bin[temp][j]++;
				*/

				int temp = snp_to_jackblock[i];
				jack_bin[temp][j]++;

				temp = snp_to_read_block[i];
				read_bin[temp][j]++;
			}
}

void read_annot_1col (string filename){
	ifstream ifs(filename.c_str(), ios::in);

	std::string line;
	std::istringstream in;

	len.resize(Nbin,0);
	step_size = Nsnp / Njack;
	step_size_rem = Nsnp%Njack;
	cout << "Number of SNPs per block : " << step_size << endl;
	jack_bin.resize(Njack, vector<int>(Nbin,0));
	int i = 0;
	while(std::getline(ifs, line)){
		in.clear();
		in.str(line);
		string temp;

		in>>temp;        
		int cur = atoi(temp.c_str());
		SNP_annot.push_back(cur);
		len[cur - 1]++;

		int jack_val = i / step_size;
		if (jack_val == Njack)
			jack_val--;
		jack_bin[jack_val][SNP_annot[i]-1]++;
		i++;
	}

	if(Nsnp != i){
		exitWithError ("Number of rows in bim file and annotation file does not match");
		exit (1);
	}
	cout << "Total number of SNPs : " << Nsnp << endl;
	for (int i = 0 ; i < Nbin ; i++){
		cout << len[i]<<" SNPs in " << i<<"-th bin" << endl;
		Nsnp_annot += len[i];
	}
}


// Read bim file to get number of SNPs 
// filename: name of .bim file
// Return number of SNPs
// This is needed for defining jackknife blocks based on the total number (jackknife scheme == 1 or 2)
int get_number_of_snps (string bimfilename){
	ifstream inp(bimfilename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< bimfilename <<endl;
		exit(1);
	}
	string line;

	int j = 0 ;
	int linenum = 0 ;

	while(std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		if (line.empty())
			continue;

		j++;
	}
	int Nsnp = j;
	inp.close();
	if (verbose >= 1) { 
		cout << "Number of SNPs in bim file = "<< Nsnp << endl;
	}

	return Nsnp;
}

// Read bim file
// filename: name of .bim file
// Return number of SNPs
int read_bim (string filename){
	ifstream inp(filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;

	int j = 0 ;
	int jchr = 0 ;
	int jlastpos = 0;
	int jlastphyspos = 0;
	int jblockchr = 0 ;

	int jack_index = 0;
	int linenum = 0 ;
	string prevchrname = "";


	while(std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		if (line.empty())
			continue;

		vector<string> toks;
		functions::tokenize (line, toks, "\t ");
			
		if (toks.size () <6) {
			exitWithError ("Bad .bim file:");
			exit (1);	
		}

		bool newchr = false;
		string chrname = toks[0];
		string snpid = toks[1];
		int physpos = atoi(toks[3].c_str());

		if (jack_scheme == 1) {
			if ( j > 0 && j % step_size == 0) {
				jack_block_size.push_back (j - jlastpos);
				jack_index ++;
				jlastpos = j;
			}
		} else if (jack_scheme == 2) {
			if (j > 0 && chrname != prevchrname){ 
				jack_block_size.push_back (j - jlastpos);
				jack_index ++;
				jblockchr ++;
				jlastpos = j;
				jlastphyspos = physpos;
			} else if (jchr > 0 && jchr % step_size == 0) {
				jack_block_size.push_back (j - jlastpos);
				jack_index ++;
				jblockchr ++;
				jlastpos = j;
				jlastphyspos = physpos;		
			}

		} else if (jack_scheme == 3) { 
			if (j > 0 && chrname != prevchrname){ 
				jack_block_size.push_back (j - jlastpos);
				jack_index ++;
				jblockchr ++;
				jlastpos = j;
				jlastphyspos = physpos;
			} else if ( physpos - jlastphyspos >= jack_size) {
				jack_block_size.push_back (j - jlastpos);
				jack_index ++;
				jblockchr ++;
				jlastpos = j;
				jlastphyspos = physpos;		
			}
		}	

		if (chrname != prevchrname) {
			jchr = 0;
			jlastphyspos = physpos;
			jblockchr = 0;
		}
		snp_to_jackblock.push_back(jack_index);

		j++;
		jchr ++;

		prevchrname = chrname;
	}
	int Nsnp = j;
	if (jack_scheme == 1) {
		if ( Nsnp > 0 ) {
			if ( Nsnp % Njack == 0)
				jack_block_size.push_back (j - jlastpos);
			else {
				jack_block_size[jack_index-1] += (Nsnp%Njack);
				for (int k  = jlastpos; k < j; k++)
					snp_to_jackblock[k] = jack_index - 1;
			}
		}	
	} else if (jack_scheme == 2) {
		jack_block_size.push_back (j - jlastpos);
	} else if (jack_scheme == 3) { 
		jack_block_size.push_back (j - jlastpos);
	}	
		
	if (jack_scheme > 0)
		Njack = jack_block_size.size();

	inp.close();
	if (verbose >= 1) { 
		cout << "Number of SNPs in bim file = "<< Nsnp << endl;
		if (jack_scheme > 0) { 
			Njack = jack_block_size.size();
			cout << "Number of jackknife blocks = " << Njack << endl;
			if (jack_scheme == 1|| jack_scheme == 2)
				cout << "Number of SNPs per jackknife block = " << step_size << endl;   
		}

		if (verbose >= 5) {
			vectorfn::printvector(jack_block_size); cout << endl;
			vectorfn::printvector(snp_to_jackblock); cout << endl;
		}
	}

	return Nsnp;
}

// Count number of individuals in .pheno file
// filename: name of .pheno file
int count_pheno(std::string filename){
	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}

	std::string line;
	int i = 0;
	while(std::getline(ifs, line)){
		i++;
	}
	int Nindv = i - 1;
	return Nindv;
}

// Count number of individuals in .fam file. 
// filename: name of .fam file
int count_fam(std::string filename){
	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}

	std::string line;
	int i = 0;
	while(std::getline(ifs, line)){
		i++;
	}
	return i;
}

//// functions related to reading without mailman
template<typename T>
static std::istream & binary_read(std::istream& stream, T& value){
	return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

void set_metadata() {
	wordsize = sizeof(char) * 8;
	unitsize = 2;
	unitsperword = wordsize / unitsize;
	mask2 = 0;
	for (int i = 0 ; i < unitsize; i++)
			mask2 = mask2 |(0x1 << i);
	nrow = Nsnp;
	ncol = ceil(1.0 * Nindv / unitsperword);
}


int simulate2_geno_from_random(float p_j){
	float rval = udis (seedr);
	float dist_pj[3] = { (1 - p_j)*(1 - p_j), 2*p_j*(1 - p_j), p_j*p_j };
	if(rval < dist_pj[0] )
		return 0;
	else if( rval >= dist_pj[0] && rval < (dist_pj[0]+dist_pj[1]))
		return 1;
	else
		return 2;
}


float get_observed_pj(const unsigned char* line){
	int y[4];
	int observed_sum = 0;
	int observed_ct = 0;

	for (int k = 0 ;k < ncol ; k++) {
		unsigned char c = line [k];
		y[0] = (c)&mask2;
		y[1] = (c>>2)&mask2;
		y[2] = (c>>4)&mask2;
		y[3] = (c>>6)&mask2;
		int j0 = k * unitsperword;
		int lmax = 4;
		if (k == ncol - 1)  {
			lmax = Nindv%4;
			lmax = (lmax == 0)?4:lmax;
		}
		for ( int l = 0 ; l < lmax; l++){
			int j = j0 + l ;
						// Extract  PLINK coded genotype and convert into 0 / 1/2
						// // PLINK coding: 
						// // 00->0
						// // 01->missing
						// // 10->1
						// // 11->2
			int val = y[l];
			val-- ;
			if(val != 0){
				val =  (val < 0 ) ? 0 :val ;
				observed_sum += val;
				observed_ct ++;
			}
		}
	}
	return observed_sum * 0.5 / observed_ct;
}


void read_bed2 (std::istream& ifs, bool allow_missing, int num_snp)  {
	//ifstream ifs (filename.c_str(), ios::in|ios::binary);
	char magic[3];
	set_metadata ();

	gtype =  new unsigned char[ncol];

	if(read_header)
		binary_read(ifs,magic);

	int sum = 0;

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	int y[4];

	int bin_pointer;

	// for a given snp, the set of annotations (bins) that the SNP belongs to
	vector<int> pointer_bins;

	for(int i = 0 ; i < num_snp ; i++){
		global_snp_index++;
		ifs.read (reinterpret_cast<char*>(gtype), ncol * sizeof(unsigned char));
		float p_j = get_observed_pj(gtype);


		// for this snp, add all annotations (bins) it belongs to
		pointer_bins.clear();      
		for(int bin_index = 0 ; bin_index < Nbin ; bin_index++)
			if(annot_bool[global_snp_index][bin_index] == 1)
				pointer_bins.push_back(bin_index);
		//bin_pointer = bin_index;

		for (int k = 0 ;k < ncol ; k++) {
			unsigned char c = gtype [k];
			// Extract PLINK genotypes
			y[0] = (c)&mask2;
			y[1] = (c>>2)&mask2;
			y[2] = (c>>4)&mask2;
			y[3] = (c>>6)&mask2;
			int j0 = k * unitsperword;
			// Handle number of individuals not being a multiple of 4
			int lmax = 4;
			if (k == ncol - 1)  {
				lmax = Nindv%4;
				lmax = (lmax == 0)?4:lmax;
			}
			for ( int l = 0 ; l < lmax; l++){
				// j: index over individuals
				int j = j0 + l ;
				// Extract  PLINK coded genotype and convert into 0/1/2
				// PLINK coding: 
				// 00->0
				// 01->missing
				// 10->1
				// 11->2
				// val:  genotype
				int val = y[l];
				if(val == 1 && !allow_missing){
					// Impute by sampling from the distribution of genotype frequencies
					val = simulate2_geno_from_random(p_j);
					val++;
					val = (val == 1) ? 0 : val;
					//val = 0;
				}
				val-- ;
				val =  (val < 0 ) ? 0 :val ;
				sum += val;

				// Loop over all bins (=annotations) this SNP belong to
				for(int bin_index = 0 ; bin_index < pointer_bins.size();bin_index++){
					bin_pointer = pointer_bins[bin_index];

					// snp_index: index over SNPs that depends on the bin being considered
					int snp_index;
					if(use_mailman == true){
						snp_index = allgen_mail[bin_pointer].index;
						int horiz_seg_no = snp_index / allgen_mail[bin_pointer].segment_size_hori;
						allgen_mail[bin_pointer].p[horiz_seg_no][j] = 3 *allgen_mail[bin_pointer].p[horiz_seg_no][j]  + val;
						// computing sum for every snp to compute mean
						allgen_mail[bin_pointer].columnsum[snp_index]+=val;

					} else {
						snp_index = allgen[bin_pointer].index;
						allgen[bin_pointer].gen(snp_index,j) = val;
					}

				}

			}
		}

		// Update number of SNPs in each bin
		for(int bin_index = 0 ; bin_index < pointer_bins.size();bin_index++){
			bin_pointer = pointer_bins[bin_index];
			if(use_mailman == true)
				allgen_mail[bin_pointer].index++;
			else
				allgen[bin_pointer].index++;
		}
	}

	sum = 0 ;
	delete[] gtype;
}


void read_bed_1colannot (std::istream& ifs,bool allow_missing,int num_snp)  {
	//ifstream ifs (filename.c_str(), ios::in|ios::binary);
	char magic[3];
	set_metadata ();

	gtype =  new unsigned char[ncol];

	if(read_header)
		binary_read(ifs,magic);

	int sum = 0;

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	int y[4];

	int bin_pointer;


	for(int i = 0 ; i < num_snp ; i++){
		global_snp_index++;
		ifs.read (reinterpret_cast<char*>(gtype), ncol * sizeof(unsigned char));
		float p_j = get_observed_pj (gtype);

		for (int k = 0 ;k < ncol ; k++) {
			unsigned char c = gtype [k];
			// Extract PLINK genotypes
			y[0] = (c)&mask2;
			y[1] = (c>>2)&mask2;
			y[2] = (c>>4)&mask2;
			y[3] = (c>>6)&mask2;
			int j0 = k * unitsperword;
			// Handle number of individuals not being a multiple of 4
			int lmax = 4;
			if (k == ncol - 1)  {
				lmax = Nindv%4;
				lmax = (lmax == 0)?4:lmax;
			}
			for ( int l = 0 ; l < lmax; l++){
				int j = j0 + l ;
				// Extract  PLINK coded genotype and convert into 0 / 1/2
				// PLINK coding: 
				// 00->0
				// 01->missing
				// 10->1
				// 11->2
				int val = y[l];
				if(val == 1 && !allow_missing){
					val = simulate2_geno_from_random(p_j);
					val++;
					val = (val == 1) ? 0 : val;
					//val = 0;
				}
				val-- ;
				val =  (val < 0 ) ? 0 :val ;
				sum += val;

				bin_pointer = SNP_annot[global_snp_index]-1;

				int snp_index;
				if(use_mailman == true){
					snp_index = allgen_mail[bin_pointer].index;
					int horiz_seg_no = snp_index / allgen_mail[bin_pointer].segment_size_hori;
					allgen_mail[bin_pointer].p[horiz_seg_no][j] = 3 *allgen_mail[bin_pointer].p[horiz_seg_no][j]  + val;
					// computing sum for every snp to compute mean
					allgen_mail[bin_pointer].columnsum[snp_index]+=val;

				}else{
					snp_index = allgen[bin_pointer].index;
					allgen[bin_pointer].gen(snp_index,j) = val;
				}
			}
		}

		bin_pointer = SNP_annot[global_snp_index]-1;
		if(use_mailman == true)
			allgen_mail[bin_pointer].index++;
		else
			allgen[bin_pointer].index++;
	}

	sum = 0 ;
	delete[] gtype;

}

// Compute jackknife SEs
// jack: matrix with statistics along rows and jackknife estimates along columns
MatrixXdr jack_se(MatrixXdr jack){
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


void genotype_stream_pass_mem_efficient (string name){
	MatrixXdr output_yXXy;
	MatrixXdr output_XXz;
	MatrixXdr tmpoutput_XXz;		
	MatrixXdr output_XXUz;
	MatrixXdr output_env;
	MatrixXdr output;
	double temp_yXXy;


	if (hetero_noise == true) {
		T_Nbin = Nbin + nongen_Nbin + Nenv;

		output_XXz = MatrixXdr::Zero(Nindv,T_Nbin * Nz);
		tmpoutput_XXz = MatrixXdr::Zero(Nindv, Nz);

		if(both_side_cov == true){
			output_XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);
		}
		output_yXXy = MatrixXdr::Zero(T_Nbin, 1);
		jack_yXXy = MatrixXdr::Zero(T_Nbin, Njack);

	} else {
		T_Nbin = Nbin + nongen_Nbin;

		output_XXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);
		tmpoutput_XXz = MatrixXdr::Zero(Nindv, Nz);

		if(both_side_cov == true){
			output_XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);
		}
		output_yXXy = MatrixXdr::Zero(T_Nbin, 1);
		jack_yXXy = MatrixXdr::Zero(T_Nbin, Njack);

	}

	if (verbose >= 3) { 		
		cout << "both_side_cov = " << both_side_cov << endl;
		cout << "cov_add_intercept = " << cov_add_intercept << endl;
		cout << "T_Nbin = " << T_Nbin << " Nbin = " << Nbin << " nongen_Nbin = " << nongen_Nbin << endl;
		cout << "Nz = " << Nz << endl;
		cout << "tmpoutput(" << tmpoutput_XXz.rows() <<"," << tmpoutput_XXz.cols() << ") "<< tmpoutput_XXz.sum() << endl;
		cout << "output(" << output_XXz.rows() <<"," << output_XXz.cols() << ") "<< output_XXz.sum() << endl;
	}

	if (opt1){ 
		string prefix = command_line_opts.OUTPUT_FILE_PATH;

		if (!use_summary_genotypes) {
			xsum_path = prefix + ".xsum";
			xsum_ofs.open(xsum_path.c_str(), std::ios_base::out);
		}
		if (use_ysum) {
			ysum_path = prefix + ".ysum";
			ysum_ofs.open(ysum_path.c_str(), std::ios_base::out);
		}
	}

	ifstream ifs (name.c_str(), ios::in|ios::binary);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< name <<endl;
		exit(1);
	}

	read_header = true;
	global_snp_index=-1;

	MatrixXdr vec1;
	MatrixXdr w1;
	MatrixXdr w2;
	MatrixXdr w3;

	MatrixXdr  A_trs(T_Nbin,T_Nbin);
	MatrixXdr b_trk(T_Nbin,1);
	MatrixXdr c_yky(T_Nbin,1);

	MatrixXdr X_l(T_Nbin + 1,T_Nbin + 1);
	MatrixXdr Y_r(T_Nbin + 1,1);
	MatrixXdr B1;
	MatrixXdr B2;
	MatrixXdr C1;
	MatrixXdr C2;
	double trkij;
	double yy = (pheno.array() * pheno.array()).sum();

	if(both_side_cov == true){
		yy = (new_pheno.array() * new_pheno.array()).sum();
	}

	Nindv_mask = mask.sum();
	if(both_side_cov == true)
		NC = Nindv_mask - Ncov;
	else
		NC = Nindv_mask;

	MatrixXdr herit;

	// Matrix of jacknife estimates
	// Rows: statistics (=variance components including sigma_e). Hence number of rows = T_Nbin + 1
	// Columns: jackknife estimates
	jack.resize(T_Nbin + 1,Njack);

	// Matrix of variance components
	// Includes sigma_e. Hence number of entries = T_Nbin + 1
	point_est.resize(T_Nbin + 1,1);

	point_est_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,1);
	jack_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,Njack);

	enrich_jack.resize(T_Nbin,Njack);
	enrich_point_est.resize(T_Nbin,1);


	MatrixXdr h1;
	MatrixXdr h2;
	MatrixXdr h3;

	double trkij_res1;
	double trkij_res2;
	double trkij_res3;
	double tk_res;

	int global_block_index = 0;
	for (int jack_index = 0 ; jack_index < Njack ; jack_index++){

		if (hetero_noise == true) {
			output_XXz = MatrixXdr::Zero(Nindv,T_Nbin * Nz);

			if(both_side_cov == true){
				output_XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);
			}
			output_yXXy = MatrixXdr::Zero(T_Nbin,1);
		} else {
			T_Nbin = Nbin + nongen_Nbin;

			output_XXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);

			if(both_side_cov == true){
				output_XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz);
			}
			output_yXXy = MatrixXdr::Zero(T_Nbin, 1);
		}	

		int jack_Nsnp = jack_block_size [jack_index];	
		int read_Nsnp = jack_Nsnp > mem_Nsnp ? mem_Nsnp: jack_Nsnp;
		int num_blocks = jack_to_read_block [jack_index];
//		int read_Nsnp = jack_Nsnp/num_blocks;
		int rem = jack_Nsnp - num_blocks * read_Nsnp;

		if (verbose >= 1) 
			cout << "************Pass 1: Reading jackknife block " << jack_index << " ************" <<endl;

		cout << "Pass 1: Reading jackknife block " << jack_index << endl;
		cout << "Dividing jackknife block " << jack_index << " into " << num_blocks << " read blocks each with " << read_Nsnp << " SNPs ";
		if (rem > 0)
			cout << " (last block has " << read_Nsnp + rem << " SNPs)";
		cout << endl;

		for (int block_index = 0; block_index < num_blocks; block_index ++, global_block_index++ ){
			read_Nsnp += (block_index < (num_blocks - 1))? 0 : rem;

			cout << "Pass 1: Reading SNP block " << block_index << " in jackknife block " << jack_index << ", global SNP block index " << global_block_index << endl;

			if (verbose >= 1)
				cout << "read_Nsnp = " << read_Nsnp << endl;

			if(use_mailman == true){
				if (verbose >= 1) { 
					cout << "Number of SNPs in each bin in jackknife block "<< jack_index << endl;
					vectorfn::printvector (jack_bin[jack_index]); cout << endl;
				}
				for (int i = 0 ; i < Nbin ; i++){
					allgen_mail[i].segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
					allgen_mail[i].Nsegments_hori = ceil(read_bin[global_block_index][i]*1.0 / (allgen_mail[i].segment_size_hori * 1.0));
					allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,vector<int>(Nindv));
					//allgen_mail[i].not_O_i.resize(jack_bin[jack_index][i]);
					//allgen_mail[i].not_O_j.resize(Nindv);
					
					// number of SNPs in this annotation (bin) in this jackknife block
					allgen_mail[i].index = 0;
					allgen_mail[i].Nsnp = read_bin[global_block_index][i];
					allgen_mail[i].Nindv = Nindv;

					allgen_mail[i].columnsum.resize(read_bin[global_block_index][i],1);
					for (int index_temp = 0 ; index_temp < read_bin[global_block_index][i];index_temp++)
						allgen_mail[i].columnsum[index_temp]=0;
				}
			} else {
				for (int k = 0 ; k < Nbin ; k++){
					allgen[k].gen.resize(read_bin[global_block_index][k],Nindv);
					// number of SNPs in this annotation (bin) in this jackknife block
					allgen[k].index = 0;
				}
			}

			if(use_1col_annot == true)
				read_bed_1colannot(ifs, missing, read_Nsnp);
			else
				read_bed2(ifs, missing, read_Nsnp);
			read_header = false;

			for (int bin_index = 0 ; bin_index < Nbin ; bin_index++){
				int num_snp;
				if (use_mailman == true)
					num_snp = allgen_mail[bin_index].index;
				else
					num_snp = allgen[bin_index].index;

				if (verbose >= 2)
					cout << "Number of SNPs in bin " << bin_index << "  = " << num_snp << endl; 

				// The case where the read block for this bin is empty
				// Can skip this bin unless this is the last read block inside jackknife block
				// In that case, we need to update the whole-genome statistics computed so far 
				// We also need to write the relevant statistics to the file
				if (num_snp == 0)  {
					if (block_index == num_blocks - 1){

						if (!use_summary_genotypes){
						for (int z_index = 0 ; z_index < Nz ; z_index++){
							XXz.col((bin_index * 2 * Nz) + Nz + z_index) += output_XXz.col(bin_index * Nz + z_index);   /// save whole sample contribution

							if(both_side_cov == true) {
								vec1 = output_XXz.col(bin_index * Nz + z_index);
								w1 = covariate.transpose() * vec1;
								w2 = Q * w1;
								w3 = covariate * w2;
								UXXz.col((bin_index * 2 * Nz) + Nz + z_index) += w3;
								UXXz.col((bin_index * 2 * Nz) + z_index) += w3;
							}
						}
						}

						yXXy(bin_index,1) += output_yXXy(bin_index, 0);


						if (opt1) {
							// Check the number of SNPs in this bin within the jackknife block
							int tmpsnp = jack_bin[jack_index][bin_index];
							// If the bin is also empty for the jackknife block 
							// Only need to write the number of SNPs (=0)
							// So that this bin within the jackknnife block can be skipped in pass number 2
							if (use_summary_genotypes){
								double temp_yXXy = output_yXXy(bin_index, 0);
								
								if (use_ysum)
									ysum_ofs.write((char *) (&temp_yXXy), sizeof(double));
								else
									jack_yXXy(bin_index, jack_index) = temp_yXXy;
							} else {
							xsum_ofs.write((char *) (&tmpsnp), sizeof(int));
							if (tmpsnp != 0){

								// Handle non-empty bin within jackknife block
								// First write XXz, XXUz.
								MatrixXdr tmpoutput = output_XXz.block (0, bin_index * Nz, output_XXz.rows(), Nz);	

								write_matrix (xsum_ofs, tmpoutput);
								if (both_side_cov == true ) {
									MatrixXdr tmpoutput = output_XXUz.block (0, bin_index * Nz, output_XXUz.rows(), Nz);	
									write_matrix (xsum_ofs, tmpoutput);
								}

								double temp_yXXy = output_yXXy(bin_index, 0);
								if (use_ysum)
									ysum_ofs.write((char *) (&temp_yXXy), sizeof(double));
								else
									jack_yXXy(bin_index, jack_index) = temp_yXXy;
							}
							}
						}
					}
					continue; 
				} // End of code to handle empty blocks

				stds.resize(num_snp,1);
				means.resize(num_snp,1);

				if(use_mailman == true){
					for (int i = 0 ; i < num_snp ; i++)
						means(i,0) = (double)allgen_mail[bin_index].columnsum[i]/Nindv;
				} else	  
					means = allgen[bin_index].gen.rowwise().mean();


				for (int i = 0 ; i < num_snp ; i++)
					stds(i,0) = 1 / sqrt((means(i,0) * (1-(0.5 * means(i,0)))));

				if (use_mailman == true){
					g = allgen_mail[bin_index];
					g.segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
					g.Nsegments_hori = ceil(read_bin[global_block_index][bin_index]*1.0 / (g.segment_size_hori * 1.0));
					g.p.resize(g.Nsegments_hori,vector<int>(Nindv));
					//g.not_O_i.resize(jack_bin[jack_index][bin_index]);
					//g.not_O_j.resize(Nindv);
					initial_var();

				} else {
					gen = allgen[bin_index].gen;
				}

				// Setup data structures for mailman
				mm = MatMult(g, gen, debug, var_normalize, memory_efficient, missing, use_mailman, nthreads, Nz);

				if (!use_summary_genotypes){
				tmpoutput_XXz = compute_XXz(num_snp, all_zb);

				for (int z_index = 0 ; z_index < Nz ; z_index++)
					output_XXz.col(bin_index * Nz + z_index) = output_XXz.col(bin_index * Nz + z_index) + tmpoutput_XXz.col(z_index);

				if (verbose >= 3) {
					cout << "all_zb = " << all_zb.rows() << "," << all_zb.cols() << "\t" << all_zb.sum () << endl;
					cout << "Pass 1: " << jack_index << " " << block_index << " " << bin_index << "\ttmpoutput(" << tmpoutput_XXz.rows() <<"," << tmpoutput_XXz.cols() << ") "<< tmpoutput_XXz.sum() << endl;
					cout << "Pass 1: " << jack_index << " " << block_index << " " << bin_index << "\toutput(" << output_XXz.rows() <<"," << output_XXz.cols() << ") "<< output_XXz.sum() << endl;
				}

				if (opt1 && block_index == (num_blocks - 1)) {
					int tmpsnp = jack_bin[jack_index][bin_index];
					MatrixXdr tmpoutput = output_XXz.block (0, bin_index * Nz, output_XXz.rows(), Nz);	
					if (verbose >= 3) {
						cout << "output(" << tmpoutput.rows() <<"," << tmpoutput.cols() << ") "<< tmpoutput.sum() << endl;
					}

					xsum_ofs.write((char *) (&tmpsnp), sizeof(int));
					write_matrix (xsum_ofs, tmpoutput);
				}

				// begin gxe computations
				// This code block has not been tested for the memory optimized setting
				MatrixXdr scaled_pheno;
				if (gen_by_env == true) {
				// This code block not tested
					for (int env_index = 0 ; env_index < Nenv ; env_index++){
						MatrixXdr env_all_zb = all_zb.array().colwise() * Enviro.col(env_index).array();
						output_env = compute_XXz(num_snp,env_all_zb);
						output_env = output_env.array().colwise() * Enviro.col(env_index).array();

						int gxe_bin_index;
						if(Annot_x_E == true)
							gxe_bin_index = Nbin + (env_index * Nbin) + bin_index;
						else
							gxe_bin_index = Nbin + env_index;

						for (int z_index = 0 ; z_index < Nz ; z_index++){
							XXz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=output_env.col(z_index);   /// save whole sample
							XXz.col((gxe_bin_index * 2*Nz) + z_index)+=output_env.col(z_index);

							if(both_side_cov == true) {
								vec1 = output_env.col(z_index);
								w1 = covariate.transpose() * vec1;
								w2 = Q * w1;
								w3 = covariate * w2;

								UXXz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=w3;
								UXXz.col((gxe_bin_index * 2*Nz) + z_index)+=w3;
							}
						}

						if (both_side_cov == true){
							MatrixXdr env_all_Uzb = all_Uzb.array().colwise() * Enviro.col(env_index).array();
							output_env = compute_XXz(num_snp,env_all_Uzb);
							output_env = output_env.array().colwise() * Enviro.col(env_index).array();

							for (int z_index = 0 ; z_index < Nz ; z_index++){
								XXUz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=output_env.col(z_index);   /// save whole sample
								XXUz.col((gxe_bin_index * 2*Nz) + z_index)+=output_env.col(z_index); 
							}
						}
						if(both_side_cov == true)
							scaled_pheno = new_pheno.array() * Enviro.col(env_index).array();
						else
							scaled_pheno= pheno.array() * Enviro.col(env_index).array();


						double temp_e_yxxy;
						temp_e_yxxy= compute_yXXy(num_snp,scaled_pheno);

						yXXy(gxe_bin_index,1)+=temp_e_yxxy;
						yXXy(gxe_bin_index,0)+=temp_e_yxxy;

					}
				}//end gxe computations

				if (block_index == num_blocks - 1){
					for (int z_index = 0 ; z_index < Nz ; z_index++){
						XXz.col((bin_index * 2 * Nz) + Nz + z_index) += output_XXz.col(bin_index * Nz + z_index);   /// save whole sample contribution

						if(both_side_cov == true) {
							vec1 = output_XXz.col(bin_index * Nz + z_index);
							w1 = covariate.transpose() * vec1;
							w2 = Q * w1;
							w3 = covariate * w2;
							UXXz.col((bin_index * 2 * Nz) + Nz + z_index) += w3;
							UXXz.col((bin_index * 2 * Nz) + z_index) += w3;
						}
					}
				}

				if (both_side_cov == true){
					output = compute_XXUz(num_snp); 
					for (int z_index = 0 ; z_index < Nz ; z_index++)
						output_XXUz.col(bin_index * Nz + z_index) = output_XXUz.col(bin_index * Nz + z_index) + output.col(z_index);

					if (opt1 && block_index == num_blocks - 1) { 
						MatrixXdr tmpoutput = output_XXUz.block (0, bin_index * Nz, output_XXUz.rows(), Nz);	

						int rows = static_cast<int>(tmpoutput.rows());
					    int cols = static_cast<int>(tmpoutput.cols());
						write_matrix (xsum_ofs, tmpoutput);
					}
					if (block_index == num_blocks - 1){
						for (int z_index = 0 ; z_index < Nz ; z_index++){
							XXUz.col((bin_index * 2 * Nz) + Nz + z_index) += output_XXUz.col(bin_index * Nz + z_index);   /// save whole sample
						}	 
					}
				}

				if (verbose >= 2)
					cout << "pheno(" << pheno.rows() << "," << pheno.cols () << ") " << pheno.sum() << endl;

				}

				if(both_side_cov == false) {
					temp_yXXy = compute_yXXy(num_snp, pheno);
					output_yXXy(bin_index, 0) += temp_yXXy;
				} else {
					temp_yXXy = compute_yVXXVy(num_snp);
					output_yXXy(bin_index, 0) += temp_yXXy;
				}				
				if (opt1 && block_index == num_blocks - 1) {
					temp_yXXy = output_yXXy(bin_index, 0);
					if (use_ysum)
						ysum_ofs.write((char *) (&temp_yXXy), sizeof(double));
					else
						jack_yXXy(bin_index, jack_index) = temp_yXXy;
				}
				if (block_index == num_blocks - 1) 
					yXXy(bin_index,1) += output_yXXy(bin_index, 0);


				if (verbose >= 2) { 
					cout << "Pass 1: " << jack_index << " " << block_index << " " << bin_index << "\tXXz(" << XXz.rows() <<"," << XXz.cols() << ") "<< XXz.sum() << endl;
					cout << "Pass 1: " << jack_index << " " << block_index << " " << bin_index << "\ttemp_yXXy : " << temp_yXXy << endl;
					cout << "Pass 1: " << jack_index <<" " << block_index << " " << bin_index << "\tyXXy(" << yXXy.rows() <<"," << yXXy.cols() << ") "<< yXXy.sum() << endl;
					if (yXXy.rows () > 1) {
						cout << "Pass 1: " << jack_index <<" " << block_index << " " << bin_index << "\toutput_yXXy " << output_yXXy(0,0) <<", " << output_yXXy(1,0) << endl;
						cout << "Pass 1: " << jack_index <<" " << block_index << " " << bin_index << "\tyXXy " << yXXy(0,1) <<", " << yXXy(1,1) << endl;
					}
				}

				mm.clean_up();
				if(use_mailman == true){
					delete[] sum_op;
					delete[] partialsums;
					delete[] yint_e;
					delete[] yint_m;
					for (int i  = 0 ; i < hsegsize; i++)
						delete[] y_m [i];
					delete[] y_m;

					for (int i  = 0 ; i < g.Nindv; i++)
						delete[] y_e[i];
					delete[] y_e;

					vector< vector<int> >().swap(g.p);
					vector< vector<int> >().swap(allgen_mail[bin_index].p);
					g.columnsum.clear();
					g.columnsum2.clear();
					g.columnmeans.clear();
					g.columnmeans2.clear();
					allgen_mail[bin_index].columnsum.clear();
					allgen_mail[bin_index].columnsum2.clear();
					allgen_mail[bin_index].columnmeans.clear();
					allgen_mail[bin_index].columnmeans2.clear();
				}
			} // loop over bins

		} // loop over read blocks

	}//end loop over jackknife blocks
	cout << "Finished reading and computing over all blocks" << endl;
	cout << endl;
	cout << endl;

	if (hetero_noise == true) {
		MatrixXdr hetro_all_Uzb;
		for (int env_index = 0 ; env_index < Nenv ; env_index++){
		/// add hetero env noise
			MatrixXdr hetro_all_zb = all_zb.array().colwise() * Enviro.col(env_index).array();
			hetro_all_zb = hetro_all_zb.array().colwise() * Enviro.col(env_index).array();

			if(both_side_cov == true){
				hetro_all_Uzb = all_Uzb.array().colwise() * Enviro.col(env_index).array();
				hetro_all_Uzb = hetro_all_Uzb.array().colwise() * Enviro.col(env_index).array();
			}

			int hetro_index;
			if(Annot_x_E == true)
				hetro_index = Nbin + (Nenv * Nbin) + env_index;
			else
				hetro_index = Nbin + Nenv + env_index;
			for (int z_index = 0 ; z_index < Nz ; z_index++){
				XXz.col(((hetro_index) * 2*Nz) + Nz + z_index) = hetro_all_zb.col(z_index);
				if(both_side_cov == true){
					vec1 = hetro_all_zb.col(z_index);
					w1 = covariate.transpose() * vec1;
					w2 = Q * w1;
					w3 = covariate * w2;
					UXXz.col(((hetro_index) * 2*Nz) + Nz + z_index) = w3;
					XXUz.col(((hetro_index) * 2*Nz) + Nz + z_index) = hetro_all_Uzb.col(z_index);
				}
			}

			MatrixXdr scaled_pheno;
			if(both_side_cov == true)
				scaled_pheno = new_pheno.array() * Enviro.col(env_index).array();
			else
				scaled_pheno= pheno.array() * Enviro.col(env_index).array();

			yXXy(hetro_index,1) = (scaled_pheno.array() * scaled_pheno.array()).sum();
			len.push_back(1);
		}
	}

	if (verbose == 1) {
		cout << "Size of bins :" << endl;
		if (hetero_noise == true) {
			for(int i = 0 ; i < Nbin + nongen_Nbin + Nenv ; i++)
				cout << "bin " << i<<" : " << len[i]<<endl;
		} else {
			for(int i = 0 ; i < Nbin + nongen_Nbin ; i++)
				cout << "bin " << i<<" : " << len[i]<<endl;
		}
	}

	cout << "Number of individuals without missing phenotype and enviroment = " << mask.sum() << endl;
	cout << endl;
	cout << endl;

	//handle when jackknife block does not include any SNPs from a bin//refill
	for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
		for (int z_index = 0 ; z_index < Nz ; z_index++){
			XXz.col((bin_index * 2 * Nz) + z_index) = XXz.col((bin_index * 2 * Nz) + Nz + z_index);
			if(both_side_cov == true){
				UXXz.col((bin_index * 2*Nz) + z_index) = UXXz.col((bin_index * 2 * Nz) + Nz + z_index);
				XXUz.col((bin_index * 2*Nz) + z_index) = XXUz.col((bin_index * 2 * Nz) + Nz + z_index);  
			}
		}
		yXXy(bin_index,0)= yXXy(bin_index,1);
	}

	if (opt1)  {
		if (!use_summary_genotypes){
			string prefix = command_line_opts.OUTPUT_FILE_PATH;
			string wgxsum_path = prefix + ".wgxsum";
			wgxsum_ofs.open(wgxsum_path.c_str(), std::ios_base::out);
			for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
				MatrixXdr tmpoutput = XXz.block (0, bin_index * 2 * Nz + Nz, XXz.rows(), Nz);
				write_matrix (wgxsum_ofs, tmpoutput);
				
				if (both_side_cov) {
					tmpoutput = UXXz.block (0, bin_index * 2 * Nz + Nz, UXXz.rows(), Nz);
					write_matrix (wgxsum_ofs, tmpoutput);
					tmpoutput = XXUz.block (0, bin_index * 2 * Nz + Nz, XXUz.rows(), Nz);
					write_matrix (wgxsum_ofs, tmpoutput);
				}
			}
			wgxsum_ofs.close ();
		}

		if (!use_summary_genotypes)
			xsum_ofs.close();
		if (use_ysum)
			ysum_ofs.close();
	}
}


void genotype_stream_pass (string name, int pass_num){
	if (verbose >= 3) {
			cout << "both_side_cov = " << both_side_cov << endl;
			cout << "cov_add_intercept = " << cov_add_intercept << endl;
			cout << "T_Nbin = " << T_Nbin << " Nbin = " << Nbin << " nongen_Nbin = " << nongen_Nbin << endl;
			cout << "Nz = " << Nz << endl;
	}	
	if (opt1){ 
		string prefix = command_line_opts.OUTPUT_FILE_PATH;
		if (use_summary_genotypes) {
			xsum_path = xsumfilepath + ".xsum";
			wgxsum_path = xsumfilepath + ".wgxsum";
		} else
			xsum_path = prefix + ".xsum";
		ysum_path = prefix + ".ysum";
		if (pass_num==1) {
			if (!use_summary_genotypes)  
				xsum_ofs.open(xsum_path.c_str(), std::ios_base::out);
			if (use_ysum)	
				ysum_ofs.open(ysum_path.c_str(), std::ios_base::out);
		} else {
			xsum_ifs.open(xsum_path.c_str(), std::ios_base::in);
			if (!xsum_ifs.is_open()) {	
				cerr << "Error reading file "<< xsum_path <<endl;
				exit(1);
			}
			if (use_ysum)	
				ysum_ifs.open(ysum_path.c_str(), std::ios_base::in);
		}
	}


	ifstream ifs (name.c_str(), ios::in|ios::binary);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< name <<endl;
		exit(1);
	}

	read_header = true;
	global_snp_index=-1;

	MatrixXdr output;
	MatrixXdr output_env;

	MatrixXdr vec1;
	MatrixXdr w1;
	MatrixXdr w2;
	MatrixXdr w3;

	MatrixXdr  A_trs(T_Nbin,T_Nbin);
	MatrixXdr b_trk(T_Nbin,1);
	MatrixXdr c_yky(T_Nbin,1);

	MatrixXdr X_l(T_Nbin + 1,T_Nbin + 1);
	MatrixXdr Y_r(T_Nbin + 1,1);
	MatrixXdr B1;
	MatrixXdr B2;
	MatrixXdr C1;
	MatrixXdr C2;
	double trkij;
	double yy = (pheno.array() * pheno.array()).sum();

	if(both_side_cov == true){
		yy = (new_pheno.array() * new_pheno.array()).sum();
	}

	Nindv_mask = mask.sum();
	if(both_side_cov == true)
		NC = Nindv_mask - Ncov;
	else
		NC = Nindv_mask;

	MatrixXdr herit;

	// Matrix of jacknife estimates
	// Rows: statistics (=variance components including sigma_e). Hence number of rows = T_Nbin + 1
	// Columns: jackknife estimates
	jack.resize(T_Nbin + 1,Njack);

	// Matrix of variance components
	// Includes sigma_e. Hence number of entries = T_Nbin + 1
	point_est.resize(T_Nbin + 1,1);

	point_est_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,1);
	jack_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,Njack);

	enrich_jack.resize(T_Nbin,Njack);
	enrich_point_est.resize(T_Nbin,1);


	MatrixXdr h1;
	MatrixXdr h2;
	MatrixXdr h3;

	double trkij_res1;
	double trkij_res2;
	double trkij_res3;
	double tk_res;


	if (pass_num == 2) {
		if (use_summary_genotypes){

			wgxsum_ifs.open(wgxsum_path.c_str(), std::ios_base::in);
			if (!wgxsum_ifs.is_open()) {	
				cerr << "Error reading file "<< wgxsum_path <<endl;
				exit(1);
			}

			for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
				read_matrix (wgxsum_ifs, output);
				XXz.block (0, bin_index * 2 * Nz + Nz, Nindv, Nz) = output;
				if (both_side_cov) {
					read_matrix (wgxsum_ifs, output);
					UXXz.block (0, bin_index * 2 * Nz + Nz, Nindv, Nz) = output;
					read_matrix (wgxsum_ifs, output);
					XXUz.block (0, bin_index * 2 * Nz + Nz, Nindv, Nz) = output;
				}
			}


			// Handle when a jackknife block does not include any SNPs from a bin
			// Only need to do this for genotype matrices (XXz, UXXz, XXUz) since
			// the phenotype summaries (yXXy) are being generated

			for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
				for (int z_index = 0 ; z_index < Nz ; z_index++){
					XXz.col((bin_index * 2*Nz) + z_index) = XXz.col((bin_index * 2*Nz) + Nz + z_index);
					if(both_side_cov == true){
						UXXz.col((bin_index * 2*Nz) + z_index) = UXXz.col((bin_index * 2*Nz) + Nz + z_index);
						XXUz.col((bin_index * 2*Nz) + z_index) = XXUz.col((bin_index * 2*Nz) + Nz + z_index);
					}
				}
//				yXXy(bin_index,0)= yXXy(bin_index,1);
			}

			wgxsum_ifs.close ();
		}
	}

	for (int jack_index = 0 ; jack_index < Njack ; jack_index++){

		int read_Nsnp = jack_block_size[jack_index];	
		cout << "Pass "<< pass_num << ": Reading jackknife block " << jack_index << endl;
		
		if (verbose >= 1)  {
			cout << "************Pass " << pass_num << ": Reading jackknife block " << jack_index << " ************" <<endl;
			if (verbose >= 2)
			cout << "read_Nsnp = " << read_Nsnp << endl;
		}

		if(use_mailman == true){
			for (int i = 0 ; i < Nbin ; i++){
				allgen_mail[i].segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
				allgen_mail[i].Nsegments_hori = ceil(jack_bin[jack_index][i]*1.0 / (allgen_mail[i].segment_size_hori * 1.0));
				allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,vector<int>(Nindv));
				//allgen_mail[i].not_O_i.resize(jack_bin[jack_index][i]);
				//allgen_mail[i].not_O_j.resize(Nindv);
				
				// number of SNPs in this annotation (bin) in this jackknife block
				allgen_mail[i].index = 0;
				allgen_mail[i].Nsnp = jack_bin[jack_index][i];
				allgen_mail[i].Nindv = Nindv;

				allgen_mail[i].columnsum.resize(jack_bin[jack_index][i],1);
				for (int index_temp = 0 ; index_temp < jack_bin[jack_index][i];index_temp++)
					allgen_mail[i].columnsum[index_temp]=0;
			}
		} else {
			for (int k = 0 ; k < Nbin ; k++){
				allgen[k].gen.resize(jack_bin[jack_index][k],Nindv);
				// number of SNPs in this annotation (bin) in this jackknife block
				allgen[k].index = 0;
			}
		}

		if (opt1  && pass_num == 2) {}
		else { 
			if(use_1col_annot == true)
				read_bed_1colannot(ifs, missing, read_Nsnp);
			else
				read_bed2(ifs, missing, read_Nsnp);
			read_header = false;
		}

		if (verbose >= 1) { 
			cout << "Number of SNPs in each bin in jackknife block "<< jack_index << endl;
			vectorfn::printvector (jack_bin[jack_index]); cout << endl;
		}
		for (int bin_index = 0 ; bin_index < Nbin ; bin_index++){
			int num_snp;
			if (opt1 && pass_num == 2) { 
				xsum_ifs.read((char *) (&num_snp), sizeof(int)); 
			} else { 
				if (use_mailman == true)
					num_snp = allgen_mail[bin_index].index;
				else
					num_snp = allgen[bin_index].index;
			}

			if (verbose >= 2)
				cout << "Number of SNPs in bin " << bin_index << "  = " << num_snp << endl; 
	
			// Skip empty bin
			if (num_snp == 0)  {
				if (opt1 && pass_num == 1) 
					xsum_ofs.write((char *) (&num_snp), sizeof(int));
				continue;
			}

			if (opt1 && pass_num==2) {
				read_matrix (xsum_ifs, output);
				
				if (verbose >= 2) {
					cout << "Pass "<< pass_num << ": " << jack_index << " " << bin_index << "\toutput(" << output.rows() <<"," << output.cols() << ") "<< output.sum() << endl;
				}
			} else {
				stds.resize(num_snp,1);
				means.resize(num_snp,1);

				if(use_mailman == true){
					for (int i = 0 ; i < num_snp ; i++)
						means(i,0) = (double)allgen_mail[bin_index].columnsum[i]/Nindv;
				} else	  
					means = allgen[bin_index].gen.rowwise().mean();


				for (int i = 0 ; i < num_snp ; i++)
					stds(i,0) = 1 / sqrt((means(i,0) * (1-(0.5 * means(i,0)))));

				if (use_mailman == true){
					g = allgen_mail[bin_index];
					g.segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
					g.Nsegments_hori = ceil(jack_bin[jack_index][bin_index]*1.0 / (g.segment_size_hori * 1.0));
					g.p.resize(g.Nsegments_hori,vector<int>(Nindv));
					//g.not_O_i.resize(jack_bin[jack_index][bin_index]);
					//g.not_O_j.resize(Nindv);
					initial_var();

				} else {
					gen = allgen[bin_index].gen;
				}

				mm = MatMult(g, gen, debug, var_normalize, memory_efficient, missing, use_mailman, nthreads, Nz);
				output = compute_XXz(num_snp, all_zb);

				if (verbose >= 2) 
					cout << "Pass "<< pass_num << ": " << jack_index << " " << bin_index << "\toutput(" << output.rows() <<"," << output.cols() << ") " << output.sum() << endl;

				if (opt1) {
					xsum_ofs.write((char *) (&num_snp), sizeof(int));
					write_matrix (xsum_ofs, output);
				}
			}

			// begin gxe computations
			// This code block has not been tested for the memory optimized setting
			MatrixXdr scaled_pheno;
			if (gen_by_env == true) { 
				// Not tested
				for (int env_index = 0 ; env_index < Nenv ; env_index++){
					MatrixXdr env_all_zb = all_zb.array().colwise() * Enviro.col(env_index).array();
					output_env = compute_XXz(num_snp,env_all_zb);
					output_env = output_env.array().colwise() * Enviro.col(env_index).array();

					int gxe_bin_index;
					if(Annot_x_E == true)
						gxe_bin_index = Nbin + (env_index * Nbin) + bin_index;
					else
						gxe_bin_index = Nbin + env_index;

					for (int z_index = 0 ; z_index < Nz ; z_index++){
						if(pass_num == 1){
							XXz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=output_env.col(z_index);   /// save whole sample
							XXz.col((gxe_bin_index * 2*Nz) + z_index)+=output_env.col(z_index);
						}
						else if(num_snp != len[bin_index])
							XXz.col((gxe_bin_index * 2*Nz) + z_index) = XXz.col((gxe_bin_index * 2*Nz) + z_index)-output_env.col(z_index);   /// save corresponding jack contrib

						if(both_side_cov == true) {
							vec1 = output_env.col(z_index);
							w1 = covariate.transpose() * vec1;
							w2 = Q * w1;
							w3 = covariate * w2;

							if(pass_num == 1){
								UXXz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=w3;
								UXXz.col((gxe_bin_index * 2*Nz) + z_index)+=w3;
							}
							else //if(num_snp != len[bin_index])
								UXXz.col((gxe_bin_index * 2*Nz) + z_index) = UXXz.col((gxe_bin_index * 2*Nz) + z_index)-w3;
						}
					}

					if (both_side_cov == true){
						MatrixXdr env_all_Uzb = all_Uzb.array().colwise() * Enviro.col(env_index).array();
						output_env = compute_XXz(num_snp,env_all_Uzb);
						output_env = output_env.array().colwise() * Enviro.col(env_index).array();

						for (int z_index = 0 ; z_index < Nz ; z_index++){
							if(pass_num == 1){
								XXUz.col((gxe_bin_index * 2*Nz) + Nz + z_index)+=output_env.col(z_index);   /// save whole sample
								XXUz.col((gxe_bin_index * 2*Nz) + z_index)+=output_env.col(z_index); 
							}
							else
								XXUz.col((gxe_bin_index * 2*Nz) + z_index) = XXUz.col((gxe_bin_index * 2*Nz) + z_index)-output_env.col(z_index);
						}
					}
					if(both_side_cov == true)
						scaled_pheno = new_pheno.array() * Enviro.col(env_index).array();
					else
						scaled_pheno= pheno.array() * Enviro.col(env_index).array();


					double temp_e_yxxy;
					temp_e_yxxy= compute_yXXy(num_snp,scaled_pheno);

					if(pass_num == 1){
						yXXy(gxe_bin_index,1)+=temp_e_yxxy;
						yXXy(gxe_bin_index,0)+=temp_e_yxxy;
					}
					else 
						yXXy(gxe_bin_index,0)= yXXy(gxe_bin_index,0)-temp_e_yxxy;      

				}
			}//end gxe computations

			// Number of columns = 2 * Number of annotations * Number of random vectors
			// For each annotation, the second half of XXz contains the whole genome statistic
			// The first half contains the jackknife statistic for the current jackknife block
			// These jackknife statistics are used to compute the jackknife variance components before moving onto the next jackknife block
			for (int z_index = 0 ; z_index < Nz ; z_index++){
				if(pass_num == 1){
					XXz.col((bin_index * 2 * Nz) + Nz + z_index) += output.col(z_index);   /// save whole genome contribution
				}
				else //if(num_snp != len[bin_index])
					XXz.col((bin_index * 2 * Nz) + z_index) = XXz.col((bin_index * 2 * Nz) + Nz + z_index)-output.col(z_index);   /// save corresponding jack contribution

				if(both_side_cov == true) {
					vec1 = output.col(z_index);
					w1 = covariate.transpose() * vec1;
					w2 = Q * w1;
					w3 = covariate * w2;
					if(pass_num == 1){
						UXXz.col((bin_index * 2*Nz) + Nz + z_index)+=w3;
						UXXz.col((bin_index * 2*Nz) + z_index)+=w3;
					}
					else if(num_snp != len[bin_index])
						UXXz.col((bin_index * 2*Nz) + z_index) = UXXz.col((bin_index * 2*Nz) + Nz + z_index)-w3;
				}
			}

			if (both_side_cov == true){
				if (opt1 && pass_num == 2) { 
					read_matrix (xsum_ifs, output);
				} else {
					output = compute_XXUz(num_snp); 
					if (opt1) {
						write_matrix (xsum_ofs, output);
					}
				}
				for (int z_index = 0 ; z_index < Nz ; z_index++){
					if(pass_num == 1){
						XXUz.col((bin_index * 2 * Nz) + Nz + z_index) += output.col(z_index);   /// save whole sample
					}
					else //if(num_snp != len[bin_index])
						XXUz.col((bin_index * 2 * Nz) + z_index) = XXUz.col((bin_index * 2*Nz) + Nz + z_index)-output.col(z_index);
				}	 
			}

			//compute yXXy
			double temp_yXXy;
			if (opt1 && pass_num == 2) { 
				if (use_ysum) 
					ysum_ifs.read((char *) (&temp_yXXy), sizeof(double)); 
				else
					temp_yXXy = jack_yXXy(bin_index, jack_index);
			} else {
				if (verbose >= 2)
					cout << "pheno(" << pheno.rows() << "," << pheno.cols () << ") " << pheno.sum() << endl;
				if(both_side_cov == false)
					temp_yXXy = compute_yXXy(num_snp, pheno);
				else
					temp_yXXy = compute_yVXXVy(num_snp);
					
				if (opt1) {
					jack_yXXy(bin_index, jack_index) = temp_yXXy;
					ysum_ofs.write((char *) (&temp_yXXy), sizeof(double));
				}	
			}

			// In the first pass, compute the whole-sample statistics: yXXy(,1)
			// In the second pass, compute the jackknife statistic for the current block
			// These jackknife statistics are used to compute the jackknife variance components before moving onto the next jackknife block
			if(pass_num == 1){
				yXXy(bin_index,1) += temp_yXXy;
			}
			else 
				yXXy(bin_index,0)= yXXy(bin_index , 1)-temp_yXXy;


			if (verbose >= 2) {
				cout << "Pass "<< pass_num << ": " << jack_index << " " << bin_index << "\tXXz(" << XXz.rows() <<"," << XXz.cols() << ") "<< XXz.sum() << endl;
				cout << "Pass "<< pass_num << ": " << jack_index << " " <<  bin_index << "\ttemp_yXXy: " << temp_yXXy << endl;
				cout << "Pass "<< pass_num << ": " << jack_index << " " << bin_index << "\tyXXy(" << yXXy.rows() <<"," << yXXy.cols() << ") "<< yXXy.sum() << endl;
				if (yXXy.rows () > 1)
					cout << "Pass " << pass_num << ": " << jack_index <<" " << bin_index << "\tyXXy " << yXXy(0,1) <<"," << yXXy(1,1) << endl;
			}

			if (opt1 && pass_num == 2){ 
			} else {
				mm.clean_up();
				if(use_mailman == true){

					delete[] sum_op;
					delete[] partialsums;
					delete[] yint_e;
					delete[] yint_m;
					for (int i  = 0 ; i < hsegsize; i++)
						delete[] y_m [i];
					delete[] y_m;

					for (int i  = 0 ; i < g.Nindv; i++)
						delete[] y_e[i];
					delete[] y_e;

					vector< vector<int> >().swap(g.p);
					vector< vector<int> >().swap(allgen_mail[bin_index].p);
					g.columnsum.clear();
					g.columnsum2.clear();
					g.columnmeans.clear();
					g.columnmeans2.clear();
					allgen_mail[bin_index].columnsum.clear();
					allgen_mail[bin_index].columnsum2.clear();
					allgen_mail[bin_index].columnmeans.clear();
					allgen_mail[bin_index].columnmeans2.clear();
				}
			}
		} // loop over bins

		if(pass_num == 2){
			// 	
			// Compute variance components for each jackknife subsample
			//
			for(int l = 0 ; l < T_Nbin ; l++)
				if( len[l]==jack_bin[jack_index][l])
					jack_bin[jack_index][l]=0;

			for (int i = 0 ; i < T_Nbin ; i++){
				b_trk(i,0) = Nindv_mask;

				if(i >= (T_Nbin-(nongen_Nbin + Nenv)) ){
					B1 = XXz.block(0, (i * 2 * Nz), Nindv, Nz);
					B1 =all_zb.array() * B1.array();
					b_trk(i,0) = B1.sum() / (len[i]-jack_bin[jack_index][i]) / Nz;
				}

				if (jack_index == Njack - 1){
					if (verbose >= 2)
						cout << "yXXy(bin,0) = " << yXXy(i,0) << ", Number of SNPs in bin[" <<i<<"] " << len[i] << " Number of SNPs in jackknife block[" << jack_index << "], bin[" <<i<<"] " << jack_bin[jack_index][i] << endl;
				}

				c_yky(i,0) = yXXy(i,0) / (len[i]-jack_bin[jack_index][i]);

				if(both_side_cov == true){
					B1 = XXz.block(0,(i * 2*Nz),Nindv,Nz);
					C1 = B1.array() * all_Uzb.array();
					C2 = C1.colwise().sum();	
					tk_res = C2.sum();  

					tk_res = tk_res / (len[i]-jack_bin[jack_index][i]) / Nz;

					b_trk(i,0) = b_trk(i,0)-tk_res;
				}

				for (int j = i ; j < T_Nbin ; j++){
					B1 = XXz.block(0,(i * 2*Nz),Nindv,Nz);
					B2 = XXz.block(0,(j * 2*Nz),Nindv,Nz);
					C1 = B1.array() * B2.array();
					C2 = C1.colwise().sum();
					trkij = C2.sum();


					if(both_side_cov == true){
						h1 = covariate.transpose() * B1;
						h2 = Q * h1;
						h3 = covariate * h2;
						C1 = h3.array() * B2.array();
						C2 = C1.colwise().sum();
						trkij_res1 = C2.sum();


						B1 = XXUz.block(0,(i * 2*Nz),Nindv,Nz);
						B2 = UXXz.block(0,(j * 2*Nz),Nindv,Nz);
						C1 = B1.array() * B2.array();
						C2 = C1.colwise().sum();
						trkij_res3 = C2.sum();


						trkij += trkij_res3 - trkij_res1 - trkij_res1 ;

					}

					trkij = trkij / (len[i]-jack_bin[jack_index][i]) / (len[j]-jack_bin[jack_index][j]) / Nz;
					A_trs(i,j) = trkij;
					A_trs(j,i) = trkij;

				}
			}


			X_l << A_trs,b_trk,b_trk.transpose(),NC;
			Y_r << c_yky,yy;
			herit = X_l.colPivHouseholderQr().solve(Y_r);


			if(jack_index == 0){
				outfile << "Number of individuals after filtering = " << Nindv_mask << endl;
				outfile << "Number of covariates = " << Ncov << endl;
				outfile << "Number of environments = " << Nenv << endl;

				if (verbose >= 2) {
					cout << "Jackknife block = " << jack_index << endl;
					cout << "Xl[" << jack_index << "] = " << X_l << endl;
					cout << "Yr[" << jack_index << "] = " << Y_r << endl;
					double relative_error = (X_l * herit - Y_r).norm() / Y_r.norm(); // norm() is L2 norm
					cout << "The relative error is: " << relative_error << endl;
					
					#ifdef USE_DOUBLE
						JacobiSVD<MatrixXd> svd(X_l);
					#else
						JacobiSVD<MatrixXf> svd(X_l);
					#endif
					double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
					cout << "Condition number:  "<< cond << endl;

					outfile << "Jackknife block = " << jack_index << endl;
					outfile << "Xl[" << jack_index << "] = " << X_l << endl;
					outfile << "Yr[" << jack_index << "] = " << Y_r << endl;
					outfile << "Normal Equations info:" << endl;
					outfile << "The relative error is: " << relative_error << endl;
					outfile << "Max sing.val: " << svd.singularValues()(0) << endl;
					outfile << "Min sing.val: " << svd.singularValues()(svd.singularValues().size()-1) << endl;
					outfile << "Condition number: " << cond << endl;
					cout << endl;
					cout << endl;
					outfile << endl;
					outfile << endl;
				}
			}

			if (jack_index == 2){
				if (verbose >= 2){
					cout << "Jackknife block = " << jack_index << endl;
					cout << "Xl[" << jack_index << "] = " << X_l << endl;
					cout << "Yr[" << jack_index << "] = " << Y_r << endl;
				}
			}

			// Fill in the jackknife statistics
			for(int i = 0 ; i<(T_Nbin + 1);i++)
				jack(i,jack_index) = herit(i,0);


			for(int i = 0 ; i<(T_Nbin + 1);i++)
				if(i == T_Nbin)
					jack_adj_gxe(i,jack_index) = jack(i,jack_index) * NC;
				else    
					jack_adj_gxe(i,jack_index) = jack(i,jack_index) * b_trk(i,0);


			//handle when a jackknife block does not include any SNPs from a bin

			for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
				for (int z_index = 0 ; z_index < Nz ; z_index++){
					XXz.col((bin_index * 2*Nz) + z_index) = XXz.col((bin_index * 2*Nz) + Nz + z_index);
					if(both_side_cov == true){
						UXXz.col((bin_index * 2*Nz) + z_index) = UXXz.col((bin_index * 2*Nz) + Nz + z_index);
						XXUz.col((bin_index * 2*Nz) + z_index) = XXUz.col((bin_index * 2*Nz) + Nz + z_index);
					}
				}
				yXXy(bin_index,0)= yXXy(bin_index,1);
			}

			if (trace){
				trace_file << X_l.block(0, 0, Nbin, Nbin);
				for (int j = 0; j< Nbin; j++) // TODO: is this the right way to count SNPs in the jn block for partitioned heritability
					trace_file << "," <<  Nsnp_annot - jack_bin[jack_index][j];
				trace_file << endl;
			}

		} //end if pass_num = 2

	}//end loop over jackknife blocks
	cout << "Finished reading and computing over all blocks" << endl;
	cout << endl;
	cout << endl;

	if(pass_num == 1){
		if (hetero_noise == true) {
			MatrixXdr hetro_all_Uzb;
			for (int env_index = 0 ; env_index < Nenv ; env_index++){
			/// add hetero env noise
				MatrixXdr hetro_all_zb = all_zb.array().colwise() * Enviro.col(env_index).array();
				hetro_all_zb = hetro_all_zb.array().colwise() * Enviro.col(env_index).array();

				if(both_side_cov == true){
					hetro_all_Uzb = all_Uzb.array().colwise() * Enviro.col(env_index).array();
					hetro_all_Uzb = hetro_all_Uzb.array().colwise() * Enviro.col(env_index).array();
				}

				int hetro_index;
				if(Annot_x_E == true)
					hetro_index = Nbin + (Nenv * Nbin) + env_index;
				else
					hetro_index = Nbin + Nenv + env_index;
				for (int z_index = 0 ; z_index < Nz ; z_index++){
					XXz.col(((hetro_index) * 2*Nz) + Nz + z_index) = hetro_all_zb.col(z_index);
					if(both_side_cov == true){
						vec1 = hetro_all_zb.col(z_index);
						w1 = covariate.transpose() * vec1;
						w2 = Q * w1;
						w3 = covariate * w2;
						UXXz.col(((hetro_index) * 2*Nz) + Nz + z_index) = w3;
						XXUz.col(((hetro_index) * 2*Nz) + Nz + z_index) = hetro_all_Uzb.col(z_index);
					}
				}

				MatrixXdr scaled_pheno;
				if(both_side_cov == true)
					scaled_pheno = new_pheno.array() * Enviro.col(env_index).array();
				else
					scaled_pheno= pheno.array() * Enviro.col(env_index).array();

				yXXy(hetro_index,1) = (scaled_pheno.array() * scaled_pheno.array()).sum();
				len.push_back(1);
			}
		}

		if (verbose >= 1) { 
			cout << "Size of bins :" << endl;
			if (hetero_noise == true) {
				for(int i = 0 ; i < Nbin + nongen_Nbin + Nenv ; i++)
					cout << "bin " << i<<" : " << len[i]<<endl;
			} else {
				for(int i = 0 ; i < Nbin + nongen_Nbin ; i++)
					cout << "bin " << i<<" : " << len[i]<<endl;
			}
		}

		cout << "Number of individuals without missing phenotype and enviroment: " << mask.sum() << endl;
		cout << endl;
		cout << endl;



		//handle when jackknife block does not include any SNPs from a bin//refill
		for (int bin_index = 0 ; bin_index < T_Nbin ; bin_index++){
			for (int z_index = 0 ; z_index < Nz ; z_index++){
				XXz.col((bin_index * 2*Nz) + z_index) = XXz.col((bin_index * 2*Nz) + Nz + z_index);
				if(both_side_cov == true){
					UXXz.col((bin_index * 2*Nz) + z_index) = UXXz.col((bin_index * 2*Nz) + Nz + z_index);
					XXUz.col((bin_index * 2*Nz) + z_index) = XXUz.col((bin_index * 2*Nz) + Nz + z_index);  
				}
			}
			yXXy(bin_index,0)= yXXy(bin_index,1);
		}

	} // pass 1

	if(pass_num == 2){
		//
		// Compute variance components for the full sample
		//
		for (int i = 0 ; i < T_Nbin ; i++){
			c_yky(i,0) = yXXy(i,1) / len[i];
			//if(both_side_cov == false)
			b_trk(i,0) = Nindv_mask;

			if(i>=(T_Nbin-(nongen_Nbin + Nenv)) ){
				B1 = XXz.block(0,(i * 2*Nz) + Nz,Nindv,Nz);
				B1 =all_zb.array() * B1.array();
				b_trk(i,0) = B1.sum() / len[i]/Nz;
			}

			if(both_side_cov == true){
				B1 = XXz.block(0,(i * 2*Nz) + Nz,Nindv,Nz);
				C1 = B1.array() * all_Uzb.array();
				C2 = C1.colwise().sum();
				tk_res = C2.sum();
				tk_res = tk_res / len[i]/Nz;

				b_trk(i,0) = b_trk(i,0)-tk_res;
			}

			for (int j = i ; j < T_Nbin ; j++){
				B1 = XXz.block(0,(i * 2*Nz) + Nz,Nindv,Nz);
				B2 = XXz.block(0,(j * 2*Nz) + Nz,Nindv,Nz);
				C1 = B1.array() * B2.array();
				C2 = C1.colwise().sum();
				trkij = C2.sum();

				if(both_side_cov == true){
					h1 = covariate.transpose() * B1;
					h2 = Q * h1;
					h3 = covariate * h2;
					C1 = h3.array() * B2.array();
					C2 = C1.colwise().sum();
					trkij_res1 = C2.sum();
					B1 = XXUz.block(0,(i * 2*Nz) + Nz,Nindv,Nz);
					B2 = UXXz.block(0,(j * 2*Nz) + Nz,Nindv,Nz);
					C1 = B1.array() * B2.array();
					C2 = C1.colwise().sum();
					trkij_res3 = C2.sum();
					trkij += trkij_res3 - trkij_res1 - trkij_res1 ;

				}
				trkij = trkij / len[i]/len[j]/Nz;
				A_trs(i,j) = trkij;
				A_trs(j,i) = trkij;
			}
		} // loop over bins     


		X_l << A_trs,b_trk,b_trk.transpose(),NC;
		Y_r << c_yky,yy;
		if (trace){
			trace_file << X_l.block(0, 0, Nbin, Nbin);
			for (int j = 0; j< Nbin; j++)
				trace_file << "," << len[j];
			trace_file << endl;
		}

		herit = X_l.colPivHouseholderQr().solve(Y_r);

		if (verbose >= 2 ) { 
			cout << "Whole-genome normal equations" << endl;
			cout << "Xl" << endl << X_l << endl;
			cout << "Yr" << endl << Y_r << endl;
		}
		for(int i = 0 ; i<(T_Nbin + 1);i++)
			point_est(i,0) = herit(i,0);

		//point_est_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,1);
		//jack_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,Njack);

		for(int i = 0 ; i<(T_Nbin + 1);i++)
			if(i == T_Nbin)
				point_est_adj_gxe(i,0) = point_est(i,0) * NC;
			else
				point_est_adj_gxe(i,0) = point_est(i,0) * b_trk(i,0);
	} // pass 2


	if (opt1){ 
		if (pass_num==1)  {
			if (!use_summary_genotypes) 
				xsum_ofs.close();
			if (use_ysum) 
				ysum_ofs.close();
		} else {
			xsum_ifs.close();
			if (!keep_xsum)
				std::remove (xsum_path.c_str());
			if (use_ysum)  {
				ysum_ifs.close();
				std::remove (ysum_path.c_str());
			}
		}
	}
}

// Regress covariates from phenotypes
//  
void regress_covariates (){
	bool normalize_proj_pheno = command_line_opts.normalize_proj_pheno;
	if(use_cov == true){
		MatrixXdr mat_mask = mask.replicate(1,Ncov);
		covariate = covariate.cwiseProduct(mat_mask);

		MatrixXdr WtW= covariate.transpose() * covariate;
		Q = WtW.inverse(); // Q = (W^tW)^-1

		if (both_side_cov == false){
			MatrixXdr v1 = covariate.transpose() * pheno; //W^ty
			MatrixXdr v2 = Q * v1;            //QW^ty
			MatrixXdr v3 = covariate * v2;    //WQW^ty
			new_pheno = pheno - v3;
			pheno = new_pheno.cwiseProduct(mask);

			y_sum = pheno.sum();
			y_mean = y_sum / mask.sum();

			for(int i = 0; i < Nindv; i++){
				if(mask(i,0)!=0)
					pheno(i,0) =pheno(i,0) - y_mean; //center phenotype
			}
			y_sum = pheno.sum();
		}

		if (both_side_cov == true){
			y_sum = pheno.sum();
			y_mean = y_sum / mask.sum();
			double phen_sd = 0;
			for(int i = 0; i < Nindv; i++){
				phen_sd+=(pheno(i,0)-y_mean) * (pheno(i,0)-y_mean);
				if(mask(i,0)!=0)
					pheno(i,0) =pheno(i,0) - y_mean; //center phenotype
			}
			phen_sd = sqrt(phen_sd / (mask.sum()-1));
			pheno = pheno / phen_sd;
			y_sum = pheno.sum();

			v1 = covariate.transpose() * pheno; //W^ty
			v2 = Q * v1;            //QW^ty
			v3 = covariate * v2;    //WQW^ty
			new_pheno = pheno - v3;
			new_pheno = new_pheno.cwiseProduct(mask);

			//normalize the projected phenotype
			if (normalize_proj_pheno == true) {
				y_sum = new_pheno.sum();
				y_mean = y_sum / mask.sum();
				phen_sd = 0;
				for(int i = 0; i < Nindv; i++){
					phen_sd+=(new_pheno(i,0)-y_mean) * (new_pheno(i,0)-y_mean);
					if(mask(i,0)!=0)
						new_pheno(i,0) =new_pheno(i,0) - y_mean; //center phenotype
				}
				phen_sd = sqrt(phen_sd / (mask.sum()-1));
				new_pheno = new_pheno / phen_sd;
			}

		}
	} else {
		y_sum = pheno.sum();
		y_mean = y_sum / mask.sum();
		for(int i = 0; i < Nindv; i++){
			if(mask(i,0)!=0)
				pheno(i,0) = pheno(i,0) - y_mean; //center phenotype
		}
		y_sum = pheno.sum();

	}
}

// Read .bim, .pheno, .env, .fam, .annot, .cov file.
void read_auxillary_files () { 

	// Read bim file to count number of SNPs
	string geno_name = command_line_opts.GENOTYPE_FILE_PATH;
	std::stringstream f1;
	f1 << geno_name << ".bim";

	if (jack_scheme == 1|| jack_scheme == 2){ 
		Nsnp = get_number_of_snps (f1.str());
		step_size = Nsnp / Njack;
		step_size_rem = Nsnp%Njack;
		read_bim (f1.str());
	} else {
		Nsnp = read_bim (f1.str());
	}

	setup_read_blocks();

	// Read phenotype and save the number of indvs
	string filename = command_line_opts.PHENOTYPE_FILE_PATH;
	Nindv = count_pheno (filename);
	read_pheno2 (Nindv, filename);
	cout << "Number of individuals = "<< Nindv << endl;
	y_sum = pheno.sum();


	// Read .env file depending on the model
	// Model to fit
	gen_by_env = command_line_opts.gen_by_env;
	hetero_noise = command_line_opts.hetero_noise;
	if (gen_by_env == false) {
		Nenv = 0;
		hetero_noise = false;
	} else {
		//Read environment file 
		std::string envfile = command_line_opts.ENV_FILE_PATH;
		Nenv = read_env (Nindv, envfile);
	}

	// Read annotation files
	filename = command_line_opts.Annot_PATH;

	if(use_1col_annot == true){
		read_annot_1col(filename);
	}else{
		read_annot(filename);
	}

	// Read .fam file
	std::stringstream f0;
	f0 << geno_name << ".fam";
	string name_fam = f0.str();
	int fam_lines = count_fam(name_fam);

	if (fam_lines != Nindv) {
		exitWithError ("Number of individuals in fam file and pheno file does not match ");
		exit (1);
	}

	// Read covariate file
	std::string covfile = command_line_opts.COVARIATE_FILE_PATH;
	std::string covname = "";
	if(covfile != "" ){
		use_cov = true;
		Ncov = read_cov (Nindv, covfile);
	} else if (covfile == ""){
		cout << "No covariate file specified" << endl;

		if (cov_add_intercept == true) {
			covariate.resize(Nindv,1);
			for(int i = 0 ; i < Nindv ; i++)
				covariate(i,0) = 1;
			Ncov = 1;
			use_cov = true;
			cout << "Intercept included" << endl;
		} else {
			both_side_cov = false;
			use_cov = false;
			cout << "No intercept included" << endl;
		}
	}


	xsumfilepath = command_line_opts.XSUM_FILE_PATH;
	if(xsumfilepath != "" ){
		use_summary_genotypes = true;	
		string xsum_path = xsumfilepath + ".xsum";
		string wgxsum_path = xsumfilepath + ".wgxsum";
		xsum_ifs.open(xsum_path.c_str(), std::ios_base::in);
		if (!xsum_ifs.is_open()) {	
			cerr << "Error reading file "<< xsum_path <<endl;
			exit(1);
		}
		xsum_ifs.close();

		wgxsum_ifs.open(wgxsum_path.c_str(), std::ios_base::in);
		if (!wgxsum_ifs.is_open()) {	
			cerr << "Error reading file "<< wgxsum_path <<endl;
			exit(1);
		}
		wgxsum_ifs.close();

	}	

}


void print_results () {
	int T_Nbin;
	if (gen_by_env == false) {
		T_Nbin = Nbin;
	} else if (hetero_noise == true) {
		T_Nbin = Nbin + nongen_Nbin + Nenv;
	} else {
		T_Nbin = Nbin + nongen_Nbin;
	}

	MatrixXdr point_se;
	point_se = jack_se(jack);

	cout << "*****" << endl;
	outfile << "*****" << endl;
	for (int i = 0 ; i < T_Nbin ; i++){
		cout << "Number of SNPs in bin " << i<<" = " << len[i]<<endl;
		outfile << "Number of SNPs in bin " << i<<" = " << len[i]<<endl;
	}
	cout << "*****" << endl;
	outfile << "*****" << endl;
	cout << "Number of G variance components = " << Nbin << endl;
	cout << "Number of GxE variance components = " << nongen_Nbin << endl;
	outfile << "Number of G variance components = " << Nbin << endl;
	outfile << "Number of GxE variance components = " << nongen_Nbin << endl;
	if (hetero_noise == true) {
		cout << "Number of NxE variance components = " << Nenv << endl;
		outfile << "Number of NxE variance components = " << Nenv << endl;
	} else {
		cout << "Number of NxE variance components = 0" << endl;
		outfile << "Number of NxE variance components = 0" << endl;
	}

	cout << "*****" << endl;
	outfile << "*****" << endl;


	cout << endl << "OUTPUT: " << endl << "Variance components: " << endl;
	outfile << "OUTPUT: " << endl << "Variance components: " << endl;
	for (int j = 0 ; j < T_Nbin ; j++){
		if(j < Nbin){
			cout << "Sigma^2_g[" << j<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
			outfile << "Sigma^2_g[" << j<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
		}
		else if (j <(Nbin + nongen_Nbin)){
			int k = j - Nbin;
			cout << "Sigma^2_gxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
			outfile << "Sigma^2_gxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
		}
		else if (j<(T_Nbin)){
			int k = j - Nbin - nongen_Nbin;
			cout << "Sigma^2_nxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
			outfile << "Sigma^2_nxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
		}
	}
	cout << "Sigma^2_e : " << point_est(T_Nbin,0) << "  SE : " << point_se(T_Nbin,0) << endl;
	outfile << "Sigma^2_e : " << point_est(T_Nbin,0) << "  SE : " << point_se(T_Nbin,0) << endl;

	double temp_sig = 0;
	double temp_sum = point_est.sum();
	for (int j = 0 ; j < T_Nbin ; j++){
		point_est(j,0) = point_est(j,0) / temp_sum;
		temp_sig += point_est(j,0);
	}
	point_est(T_Nbin,0) = temp_sig;


	for (int i = 0 ; i < Njack ; i++){
		temp_sig = 0;
		temp_sum = jack.col(i).sum();
		for (int j = 0 ; j < T_Nbin ; j++){
			jack(j,i) = jack(j,i) / temp_sum;
			temp_sig += jack(j,i);
		}
		jack(T_Nbin,i) = temp_sig;
	}

	////adj for GXE
	temp_sum = point_est_adj_gxe.sum();
	temp_sig = 0;
	for (int j = 0 ; j < T_Nbin ; j++){
		point_est_adj_gxe(j,0) = point_est_adj_gxe(j,0) / temp_sum;
		temp_sig += point_est_adj_gxe(j,0);
	}
	point_est_adj_gxe(T_Nbin,0) = temp_sig;
	temp_sig = 0;
	if(Annot_x_E == true){
		for(int k = Nbin ; k<(2 * Nbin);k++)
			temp_sig += point_est_adj_gxe(k,0);
		point_est_adj_gxe(T_Nbin + 1,0) = temp_sig;// total GxE
	}else{
			point_est_adj_gxe(T_Nbin + 1,0) = point_est_adj_gxe(Nbin,0);

	}
	temp_sig = 0;
	for(int k = 0 ; k < Nbin ; k++)
		temp_sig += point_est_adj_gxe(k,0);
	point_est_adj_gxe(T_Nbin + 2,0) = temp_sig;//total G

	double temp2_sig = 0;
	for (int i = 0 ; i < Njack ; i++){
		temp_sig = 0;
		temp_sum = jack_adj_gxe.col(i).sum();
		for (int j = 0 ; j < T_Nbin ; j++){
			jack_adj_gxe(j,i) = jack_adj_gxe(j,i) / temp_sum;
			temp_sig += jack_adj_gxe(j,i);
		}
		jack_adj_gxe(T_Nbin,i) = temp_sig;
		temp2_sig = 0;
		if(Annot_x_E == true){
			for(int k = Nbin ; k<(2 * Nbin);k++)
				temp2_sig += jack_adj_gxe(k,i);
			jack_adj_gxe(T_Nbin + 1,i) = temp2_sig;
		}else{
			jack_adj_gxe(T_Nbin + 1,i) = jack_adj_gxe(Nbin,i);
		}
		temp2_sig = 0;
		for(int k = 0 ; k < Nbin ; k++)
			temp2_sig += jack_adj_gxe(k,i);
		jack_adj_gxe(T_Nbin + 2,i) = temp2_sig;
	}

	MatrixXdr SEjack_adj_gxe = jack_se(jack_adj_gxe);
	cout << "*****" << endl;
	outfile << "*****" << endl;
	cout << "Heritabilities: " << endl;
	outfile << "Heritabilities: " << endl;
	for (int j = 0 ; j < T_Nbin ; j++){
		if(j < Nbin){
			cout << "h2_g[" << j<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
			outfile << "h2_g[" << j<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
		}
		else if (j<(Nbin + nongen_Nbin)){
			int k = j - Nbin;
			cout << "h2_gxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
			outfile << "h2_gxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
		}
		else if (j<(T_Nbin)){
			int k = j - Nbin - nongen_Nbin;
			cout << "h2_nxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
			outfile << "h2_nxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
		}
	}
	cout << "Total h2 : " << point_est_adj_gxe(T_Nbin,0) << " SE: " << SEjack_adj_gxe(T_Nbin,0) << endl;
	outfile << "Total h2 : " << point_est_adj_gxe(T_Nbin,0) << " SE: " << SEjack_adj_gxe(T_Nbin,0) << endl;
	cout << "Total h2_g : " << point_est_adj_gxe(T_Nbin + 2,0) << " SE: " << SEjack_adj_gxe(T_Nbin + 2,0) << endl;
	outfile << "Total h2_g : " << point_est_adj_gxe(T_Nbin + 2,0) << " SE: " << SEjack_adj_gxe(T_Nbin + 2,0) << endl;
	if (gen_by_env == true) {
		cout << "Total h2_gxe : " << point_est_adj_gxe(T_Nbin + 1,0) << " SE: " << SEjack_adj_gxe(T_Nbin + 1,0) << endl;
		outfile << "Total h2_gxe : " << point_est_adj_gxe(T_Nbin + 1,0) << " SE: " << SEjack_adj_gxe(T_Nbin + 1,0) << endl;
	}
	cout << "*****" << endl;
	outfile << "*****" << endl;

	MatrixXdr enrich_g;
	MatrixXdr jack_enrich_g;
	jack_enrich_g.resize(Nbin,Njack);
	enrich_g.resize(Nbin,1);
	double total_g_h2 = 0;
	int total_snp = 0;
	for(int i = 0 ; i < Nbin ; i++){
		total_g_h2 += point_est_adj_gxe(i,0);
		total_snp += len[i];
	}
	double numi;
	double denom = (double)total_g_h2 / total_snp;
	for(int i = 0 ; i < Nbin ; i++){
		numi = point_est_adj_gxe(i,0) / len[i];
		enrich_g(i,0) = numi / denom;
	}
	for(int j = 0 ; j < Njack ; j++){
		total_g_h2 = 0;
		total_snp = 0;
		for (int k = 0 ; k < Nbin ; k++){
			total_snp += len[k]-jack_bin[j][k];
			total_g_h2 += jack_adj_gxe(k,j);
		}
		denom = (double)total_g_h2 / total_snp;
		for(int k = 0 ; k < Nbin ; k++){
			numi = (double)jack_adj_gxe(k,j) / (len[k]-jack_bin[j][k]);
			jack_enrich_g(k,j) = (double)numi / denom;
			}
	}
	MatrixXdr enrich_g_se;
	enrich_g_se = MatrixXdr::Zero(Nbin,1);
	enrich_g_se = jack_se(jack_enrich_g);


	cout << "Enrichments:" << endl;
	outfile << "Enrichments:" << endl;
	cout << "G enrichment" << endl;
	outfile << "G enrichment" << endl;
	for(int i = 0 ; i < Nbin ; i++){
		cout << "Enrichment g[" << i<<"] : " << enrich_g(i,0) << " SE : " << enrich_g_se(i,0) << endl;
		outfile << "Enrichment g[" << i<<"] : " << enrich_g(i,0) << " SE : " << enrich_g_se(i,0) << endl;
	}
	// compute enrich GxE
	if(Annot_x_E == true){
		if (gen_by_env == true) {
			MatrixXdr enrich_gxe;
			MatrixXdr jack_enrich_gxe;
			jack_enrich_gxe.resize(Nbin,Njack);
			enrich_gxe.resize(Nbin,1);
			double total_gxe_h2 = 0;
			int total_snp = 0;
			for(int i = 0 ; i < Nbin ; i++){
				total_gxe_h2 += point_est_adj_gxe(Nbin + i,0);
				total_snp += len[Nbin + i];
			}
			double numi;
			double denom = (double)total_gxe_h2 / total_snp;
			for(int i = 0 ; i < Nbin ; i++){
				numi = point_est_adj_gxe(Nbin + i,0) / len[Nbin + i];
				enrich_gxe(i,0) = numi / denom;
			}
			for(int j = 0 ; j < Njack ; j++){
				total_gxe_h2 = 0;
				total_snp = 0;
				for (int k = 0 ; k < Nbin ; k++){
					total_snp += len[Nbin + k]-jack_bin[j][Nbin + k];
					total_gxe_h2 += jack_adj_gxe(Nbin + k,j);
				}
				denom = (double)total_gxe_h2 / total_snp;
				for(int k = 0 ; k < Nbin ; k++){
					numi = jack_adj_gxe(Nbin + k,j) / (len[Nbin + k]-jack_bin[j][Nbin + k]);
					jack_enrich_gxe(k,j) = numi / denom;
				}
			}
			MatrixXdr enrich_gxe_se;
			enrich_gxe_se = MatrixXdr::Zero(Nbin,1);
			enrich_gxe_se = jack_se(jack_enrich_gxe);

			cout << "GxE enrichment" << endl;
			outfile << "GxE enrichment" << endl;
			for(int i = 0 ; i < Nbin ; i++){
				cout << "Enrichment gxe[" << i<<"] : " << enrich_gxe(i,0) << " SE : " << enrich_gxe_se(i,0) << endl;
				outfile << "Enrichment gxe[" << i<<"] : " << enrich_gxe(i,0) << " SE : " << enrich_gxe_se(i,0) << endl;
			}
		}
	}
	cout << "*****" << endl;
	outfile << "*****" << endl;

	///compute parameters for overlapping annotations based  on s-ldsc definition :

	MatrixXdr her_per_snp;
	MatrixXdr her_cat_ldsc;
	MatrixXdr point_her_cat_ldsc;;
	her_cat_ldsc = MatrixXdr::Zero(T_Nbin,Njack);
	MatrixXdr her_per_snp_inbin(T_Nbin,1);
	point_her_cat_ldsc = MatrixXdr::Zero(T_Nbin,1);

	for (int k = 0 ; k<=Njack ; k++){
		if(k == Njack){
			for(int i = 0 ; i < T_Nbin ; i++){
				her_per_snp_inbin(i,0) = (double)point_est_adj_gxe(i,0) / len[i];
			}
		}else{
			for(int i = 0 ; i < T_Nbin ; i++){
				her_per_snp_inbin(i,0) = (double)jack_adj_gxe(i,k) / (len[i]-jack_bin[k][i]);
			}

		}
		her_per_snp = MatrixXdr::Zero(Nsnp,2);


		for(int i = 0 ; i < Nsnp ; i++){
			for(int j = 0 ; j < Nbin ; j++){
				if(annot_bool[i][j]==1)
					her_per_snp(i,0)+=her_per_snp_inbin(j,0);
				if((annot_bool[i][Nbin + j]==1) && (Annot_x_E == true))
					her_per_snp(i,1)+=her_per_snp_inbin(Nbin + j,0);
			}        
			if(k == Njack){
				for(int j = 0 ; j < Nbin ; j++){
					if(annot_bool[i][j]==1)
						point_her_cat_ldsc(j,0)+=her_per_snp(i,0);
					if((annot_bool[i][Nbin + j]==1) && (Annot_x_E == true))
						point_her_cat_ldsc(Nbin + j,0)+=her_per_snp(i,1);
				}                
			}else{
				int temp = i / step_size;
				if(temp>=Njack)
					temp = Njack - 1;
				for(int j = 0 ; j < Nbin ; j++){
					if(annot_bool[i][j]==1 && temp != k)
						her_cat_ldsc(j,k)+=her_per_snp(i,0);
					if((annot_bool[i][Nbin + j]==1) && (temp != k) && (Annot_x_E == true))
						her_cat_ldsc(Nbin + j,k)+=her_per_snp(i,1);
				}                   
			}     
		}
	}

	MatrixXdr se_her_cat_ldsc = jack_se(her_cat_ldsc);

	cout << "*****" << endl;
	outfile << "*****" << endl;
	cout << "Heritabilities and enrichments computed based on overlapping setting" << endl;
	outfile << "Heritabilities and enrichments computed based on overlapping setting" << endl;

	cout << "Heritabilities: " << endl;
	outfile << "Heritabilities: " << endl;
	for (int j = 0 ; j < Nbin ; j++){
		cout << "h2_g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
		outfile << "h2_g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
	}
	if ((gen_by_env == true) && (Annot_x_E == true)) {
		for (int j = 0 ; j < Nbin ; j++){
			int k = j + Nbin;
			cout << "h2_gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;
			outfile << "h2_gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;

		}
	}

	int all_snp = 0;
	for(int i = 0 ; i < Nbin ; i++){
		all_snp += len[i];
	}
	double snp_por;
	for (int i = 0 ; i < Nbin ; i++){
		point_her_cat_ldsc(i,0) = (double)point_her_cat_ldsc(i,0) / point_est_adj_gxe(T_Nbin + 2,0);  //additive
		if((gen_by_env == true) && (Annot_x_E == true))
			point_her_cat_ldsc(i + Nbin,0) = (double)point_her_cat_ldsc(i + Nbin,0) / point_est_adj_gxe(T_Nbin + 1,0);  //GxE

		snp_por = (double)len[i]/all_snp;
		point_her_cat_ldsc(i,0) = (double)point_her_cat_ldsc(i,0) / snp_por;
		if((gen_by_env == true) && (Annot_x_E == true))
			point_her_cat_ldsc(i + Nbin,0) = (double)point_her_cat_ldsc(i + Nbin,0) / snp_por;
	}

	double temp_size;
	for(int i = 0 ; i < Njack ; i++){
		temp_size = all_snp;
		for(int k = 0 ; k < Nbin ; k++)
			temp_size = temp_size - jack_bin[i][k];
		for(int j = 0 ; j < Nbin ; j++){
			her_cat_ldsc(j,i) = (double)her_cat_ldsc(j,i) / jack_adj_gxe(T_Nbin + 2,i);
			if((gen_by_env == true) && (Annot_x_E == true))
				her_cat_ldsc(j + Nbin,i) = (double)her_cat_ldsc(j + Nbin,i) / jack_adj_gxe(T_Nbin + 1,i);

			snp_por = (double)(len[j]-jack_bin[i][j]) / temp_size;
			her_cat_ldsc(j,i) = (double)her_cat_ldsc(j,i) / snp_por;
			if((gen_by_env == true) && (Annot_x_E == true))
				her_cat_ldsc(j + Nbin,i) = (double)her_cat_ldsc(j + Nbin,i) / snp_por;

		}
	}

	se_her_cat_ldsc = jack_se(her_cat_ldsc);

	cout << "Enrichments (overlapping def): " << endl;
	outfile << "Enrichments (overlapping def): " << endl;
	for (int j = 0 ; j < Nbin ; j++){
		cout << "Enrichment g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
		outfile << "Enrichment g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
	}
	if ((gen_by_env == true) && (Annot_x_E == true)) {
		for (int j = 0 ; j < Nbin ; j++){
			int k = j + Nbin;
			cout << "Enrichment gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;
			outfile << "Enrichment gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE: " << se_her_cat_ldsc(k,0) << endl;
		}
	}
}

void print_trace () { 
	string prefix = command_line_opts.OUTPUT_FILE_PATH;
	string trpath = prefix + ".trace";
	string mnpath = prefix + ".MN";
	trace_file.open(trpath.c_str(), std::ios_base::out);
	trace_file << "TRACE,NSNPS_JACKKNIFE" << endl;
	meta_file.open(mnpath.c_str(), std::ios_base::out);
	meta_file << "NSAMPLE,NSNP,NBLK,NBIN" << endl << Nindv << "," << Nsnp << "," << Njack << "," << Nbin;
	meta_file.close();
}

void print_input_parameters(){ 
	string outpath = command_line_opts.OUTPUT_FILE_PATH;
	outfile.open(outpath.c_str(), std::ios_base::out);
	
	outfile << "##################################" << endl;
	outfile << "#                                #" << endl;
	outfile << "#          GENIE (v1.1.0)        #" << endl;
	outfile << "#                                #" << endl;
	outfile << "##################################" << endl;
	outfile << endl;
	outfile << endl;

	outfile << "Active essential options: " << endl;
	if (command_line_opts.GENOTYPE_FILE_PATH != "")
		outfile << "\t - g (genotype) " << command_line_opts.GENOTYPE_FILE_PATH << endl;
	if (command_line_opts.Annot_PATH != "")
		outfile << "\t - annot (annotation) " << command_line_opts.Annot_PATH << endl;
	if (command_line_opts.PHENOTYPE_FILE_PATH != "")
		outfile << "\t - p (phenotype) " << command_line_opts.PHENOTYPE_FILE_PATH << endl;
	if (command_line_opts.COVARIATE_FILE_PATH != "")
		outfile << "\t - c (covariates) " << command_line_opts.COVARIATE_FILE_PATH << endl;
	if (command_line_opts.OUTPUT_FILE_PATH != "")
		outfile << "\t - o (output) " << command_line_opts.OUTPUT_FILE_PATH << endl;
	if (command_line_opts.ENV_FILE_PATH != "")
		outfile << "\t - e (environment) " << command_line_opts.ENV_FILE_PATH << endl;
	if (command_line_opts.model != "")
		outfile << "\t - m (model) " << command_line_opts.model << endl;
	if (command_line_opts.num_of_evec > 0) 
		outfile << "\t - k (# random vectors) " << std::to_string(command_line_opts.num_of_evec) << endl;
	if (command_line_opts.jack_number > 0)
		outfile << "\t - jn (# jackknife blocks) " << std::to_string(command_line_opts.jack_number) << endl;
	if (command_line_opts.nthreads > 0)
		outfile << "\t - t (# threads) " << std::to_string(command_line_opts.nthreads) << endl;
	if (command_line_opts.seed != -1)
		outfile << "\t - s (seed) " << std::to_string(command_line_opts.seed) << endl;
	if (command_line_opts.exannot == true)
		outfile << "\t - eXannt (paritioned GxE)" << endl;
	if (verbose >= 1) {
		outfile << "Other options: " << endl;
		outfile << "\t - norm_proj_pheno (normalize pheno after projection on covariates) " << std::to_string(command_line_opts.normalize_proj_pheno) << endl;
		outfile << "\t - cov_add_intercept (intercept term added to covariates) " << std::to_string(command_line_opts.cov_add_intercept) << endl;
		outfile << "\t - v (verbose) " << std::to_string(command_line_opts.verbose) << endl;
	}
	if (trace)
		outfile << "\t - tr (trace summaries) " << endl;

	outfile << endl;
	outfile << endl;
	if ( verbose >= 3 ) {
		cout << "use_double = " << use_double << endl;
	}
}

void init_params () {

	verbose = command_line_opts.verbose;
	debug = command_line_opts.debugmode ||  verbose>= 3 ;

	trace = command_line_opts.print_trace;
	k_orig = command_line_opts.num_of_evec ;
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = false;
	k = k_orig + command_line_opts.l;
	k = (int)ceil(k / 10.0) * 10;
	command_line_opts.l = k - k_orig;
	srand((unsigned int) time(0));
	Nz = command_line_opts.num_of_evec;
	k = Nz;

	jack_scheme = command_line_opts.jack_scheme;
	Njack = command_line_opts.jack_number;
	jack_size = command_line_opts.jack_size;
	jack_size *= 1e6; 

	opt1 = command_line_opts.opt1;
	opt2 = command_line_opts.opt2;
	mem_Nsnp = command_line_opts.mem_Nsnp;
	use_mailman = command_line_opts.fast_mode;

	seed = command_line_opts.seed;
	if (command_line_opts.exannot == true)
		Annot_x_E = true;
	nthreads = command_line_opts.nthreads;
	
	cov_add_intercept = command_line_opts.cov_add_intercept;
	use_ysum = command_line_opts.use_ysum;
	keep_xsum = command_line_opts.keep_xsum;
}

int main(int argc, char const *argv[]){
    struct timeval now;
    gettimeofday(&now, NULL);
    long starttime = now.tv_sec * UMILLION + now.tv_usec;

	parse_args (argc,argv);
	init_params ();

	read_auxillary_files ();
	regress_covariates ();	

	if (gen_by_env == true) {    
		///mask out indv with missingness from Enviro
		Enviro = Enviro.array().colwise() * mask.col(0).array();
	}

	//define random vector z's
	all_zb= MatrixXdr::Random(Nindv,Nz);

	if (seed == -1) {
		std::random_device rd;
		seedr.seed(rd());
	} else {
		seedr.seed(seed);
	}

	std::normal_distribution<> dist(0,1);

	auto z_vec = std::bind(dist, seedr);

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			all_zb(j,i) = z_vec();

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			all_zb(j,i) = all_zb(j,i) * mask(j,0);

	if(both_side_cov == true){
		all_Uzb.resize(Nindv , Nz);
		for (int j = 0 ; j < Nz ; j++){
			MatrixXdr w1 = covariate.transpose() * all_zb.col(j);
			MatrixXdr w2 = Q * w1;
			MatrixXdr w3 = covariate * w2;
			all_Uzb.col(j) = w3;
		}
	}

	if(Annot_x_E == true){
		nongen_Nbin = Nenv * Nbin;
		for (int i = 0 ; i < Nenv ; i++)
			for (int j = 0 ; j < Nbin ; j++)
				len.push_back(len[j]);
	}else{
		nongen_Nbin = Nenv;
		for (int i = 0 ; i < Nenv ; i++)
			len.push_back(Nsnp);
	}

	// XXz : Nindv X (Nbins * Nz * 2): 
	// First half of columns has the computation from the current jackknife block.
	// Last half of columns has the computation from the whole genome.
	if(hetero_noise == true) {
		T_Nbin = Nbin + nongen_Nbin + Nenv;
	} else {
		T_Nbin = Nbin + nongen_Nbin;
	}

	XXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);

	if(both_side_cov == true){
		UXXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);
		XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);
	}
	yXXy = MatrixXdr::Zero(T_Nbin,2);

	if (!use_ysum)
		jack_yXXy = MatrixXdr::Zero(T_Nbin, Njack);

	if(use_mailman == true) 
		allgen_mail.resize(Nbin);
	else
		allgen.resize(Nbin);

	int bin_index = 0;
	///// code for handling overlapping annotations
	string geno_name = command_line_opts.GENOTYPE_FILE_PATH;
	std::stringstream f3;
	f3 << geno_name << ".bed";
	string name = f3.str();
	ifstream ifs (name.c_str(), ios::in|ios::binary);
	read_header = true;
	global_snp_index=-1;

	if (!ifs.is_open()){
		cerr << "Error reading file "<< name  <<endl;
		exit(1);
	}

	cout << endl;
	cout << "Reading genotypes ..." << endl;

	print_input_parameters ();
	
	//CHANGE(03 / 05): add trace summary files. input to -o is now just the prefix (all output file endings are fixed to .log)
	if (trace){
		print_trace ();
	}

	// This is where most of the computation happens
	if (opt2){
		genotype_stream_pass_mem_efficient (name);
	} else {
		genotype_stream_pass (name,1);
	}
	genotype_stream_pass (name,2);
	print_results ();

    gettimeofday(&now, NULL);
    long endtime = now.tv_sec * UMILLION + now.tv_usec;
    double elapsed = endtime - starttime;
    elapsed/=1.e6;
    cout << "GENIE ran successfully. Time elapsed = " << elapsed << " seconds " << endl;
	
	return 0;
}
