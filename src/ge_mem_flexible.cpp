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
// #include "mailman.h"
#include "arguments.h"
//#include "helper.h"
#include "storage.h"
#include "matmult.h"
#include "io.h"
#include "std.h"

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
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;
//typedef Matrix<int, Dynamic, Dynamic, RowMajor> MatrixXdrInt;

class data {
	public:
		MatrixXdr gen;
		int index;
};

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
MatrixXdr jack_bin_size;
vector<int> len;
vector<int> Annot;
int Njack = 100;
int Nbin = 160;
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
int Nsnp_annot = 0;
int Nindv;
bool **bin_annot;
int step_size;
int step_size_rem;
std::vector<std::vector<bool> > annot_bool;
std::vector<std::vector<int> > jack_bin;
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
bool verbose;
bool trace;
int nthreads;
MatMult mm;
std::ofstream outfile;
std::ofstream trace_file;
std::ofstream meta_file;
string prefix;

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
		if (verbose == true) {
			cout << "The shape of Env: " << Enviro.rows() << " " << Enviro.cols() << endl;
			cout << "The shape of Cov: " << covariate.rows() << " " << covariate.cols() << endl;
			cout << "Number of covariates: " << covNum<< ", the number of environments: " << Nenv << endl;
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

double compute_yXXy(int num_snp,MatrixXdr vec){
	MatrixXdr res(num_snp, 1);
	// if(use_mailman == true)
	//         multiply_y_pre_fast(vec,1,res,false);
	// else
	//          res = gen * vec;
	mm.multiply_y_pre(vec,1,res,false);

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid = means.cwiseProduct(stds);
	resid = resid *vec.sum();
	MatrixXdr Xy(num_snp,1);
	Xy = res - resid;

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
	mm.multiply_y_pre(new_pheno,1,res,false);

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid = means.cwiseProduct(stds);
	resid = resid *new_pheno_sum;
	MatrixXdr Xy(num_snp,1);
	Xy = res - resid;
	double ytVXXVy = (Xy.array()* Xy.array()).sum();
	return ytVXXVy;
}


MatrixXdr  compute_XXz (int num_snp,MatrixXdr Zvec){
	res.resize(num_snp, Nz);

		// if(use_mailman == true)
		//         multiply_y_pre_fast(Zvec,Nz,res, false);
		// else
		//         res = gen * Zvec;
	mm.multiply_y_pre(Zvec,Nz,res, false);


	MatrixXdr zb_sum = Zvec.colwise().sum();


	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	MatrixXdr resid(num_snp, Nz);
	MatrixXdr inter = means.cwiseProduct(stds);
	resid = inter * zb_sum;
	MatrixXdr inter_zb = res - resid;

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
			inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
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
	int selected_snps = 0;
	for (int i = 0 ; i < num_parti ; i++){
		cout << len[i]<<" SNPs in bin " << i << endl;
		selected_snps += len[i];
	}

	cout << "Number of SNPs selected according to annotation file : " <<selected_snps << endl;
	Nsnp_annot = selected_snps;


	step_size = Nsnp / Njack;
	step_size_rem = Nsnp%Njack;
	cout << "Number of SNPs per jackknife block : " << step_size << endl;   

	int Total_Nbin;
	if(Annot_x_E == false){
		jack_bin.resize(Njack, vector<int>(Nbin + Nenv + Nenv,0));
		Total_Nbin = Nbin + Nenv;
	}else{
		jack_bin.resize(Njack, vector<int>(Nbin + (Nenv * Nbin) + Nenv,0));
		Total_Nbin = Nbin + (Nenv * Nbin);
	}
	int temp;
	for (int i = 0 ; i < Nsnp ; i++)
		for(int j = 0 ; j < Total_Nbin ; j++)
			if (annot_bool[i][j]==1){
				temp = i / step_size;
				if (temp>=Njack)
					temp = Njack - 1;
				jack_bin[temp][j]++;
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
	if (verbose == true)
		cout << "Number of SNPs in bim file: "<< Nsnp <<endl;
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
	float rval = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float dist_pj[3] = { (1 - p_j) * (1 - p_j), 2 * p_j * (1 - p_j), p_j * p_j };
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


void read_bed2 (std::istream& ifs,bool allow_missing,int num_snp)  {
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

	vector<int> pointer_bins;

	for(int i = 0 ; i < num_snp ; i++){
		global_snp_index++;
		ifs.read (reinterpret_cast<char*>(gtype), ncol * sizeof(unsigned char));
		float p_j = get_observed_pj(gtype);


		pointer_bins.clear();      
		for(int bin_index = 0 ; bin_index < Nbin ; bin_index++)
			if(annot_bool[global_snp_index][bin_index]==1)
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

				for(int bin_index = 0 ; bin_index < pointer_bins.size();bin_index++){
					bin_pointer = pointer_bins[bin_index];
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
		float p_j = get_observed_pj(gtype);

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

void genotype_stream_pass(string name, int pass_num){
	int T_Nbin;
	if (hetero_noise == true)
		T_Nbin = Nbin + nongen_Nbin + Nenv;
	else
		T_Nbin = Nbin + nongen_Nbin;

	ifstream ifs (name.c_str(), ios::in|ios::binary);
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

	jack.resize(T_Nbin + 1,Njack);
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

	for (int jack_index = 0 ; jack_index < Njack ; jack_index++){

		int read_Nsnp = (jack_index<(Njack - 1)) ? (step_size) : (step_size + step_size_rem);
		cout << "Reading block " << jack_index << " in pass " << pass_num << endl;

		if(use_mailman == true){
			for (int i = 0 ; i < Nbin ; i++){
				allgen_mail[i].segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
				allgen_mail[i].Nsegments_hori = ceil(jack_bin[jack_index][i]*1.0 / (allgen_mail[i].segment_size_hori * 1.0));
				allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,std::vector<int>(Nindv));
				//allgen_mail[i].not_O_i.resize(jack_bin[jack_index][i]);
				//allgen_mail[i].not_O_j.resize(Nindv);
				allgen_mail[i].index = 0;
				allgen_mail[i].Nsnp = jack_bin[jack_index][i];
				allgen_mail[i].Nindv = Nindv;

				allgen_mail[i].columnsum.resize(jack_bin[jack_index][i],1);
				for (int index_temp = 0 ; index_temp < jack_bin[jack_index][i];index_temp++)
					allgen_mail[i].columnsum[index_temp]=0;
			}
		}else{
			for (int k = 0 ; k < Nbin ; k++){
				allgen[k].gen.resize(jack_bin[jack_index][k],Nindv);
				allgen[k].index = 0;
			}
		}

		if(use_1col_annot == true)
			read_bed_1colannot(ifs,missing,read_Nsnp);
		else
			read_bed2(ifs,missing,read_Nsnp);
		read_header = false;

		for (int bin_index = 0 ; bin_index < Nbin ; bin_index++){
			int num_snp;
			if (use_mailman == true)
				num_snp = allgen_mail[bin_index].index;
			else
				num_snp = allgen[bin_index].index;

			if(num_snp != 0){
				stds.resize(num_snp,1);
				means.resize(num_snp,1);

				if(use_mailman == true){
					for (int i = 0 ; i < num_snp ; i++)
						means(i,0) = (double)allgen_mail[bin_index].columnsum[i]/Nindv;
				}			
				else	  
					means = allgen[bin_index].gen.rowwise().mean();


				for (int i = 0 ; i < num_snp ; i++)
					stds(i,0) = 1 / sqrt((means(i,0) * (1-(0.5 * means(i,0)))));

				if (use_mailman == true){
					g = allgen_mail[bin_index];
					g.segment_size_hori = floor(log(Nindv) / log(3)) - 2 ;
					g.Nsegments_hori = ceil(jack_bin[jack_index][bin_index]*1.0 / (g.segment_size_hori * 1.0));
					g.p.resize(g.Nsegments_hori,std::vector<int>(Nindv));
					//g.not_O_i.resize(jack_bin[jack_index][bin_index]);
					//g.not_O_j.resize(Nindv);
					initial_var();

				}else{
					gen = allgen[bin_index].gen;
				}

				mm = MatMult(g, gen, false, var_normalize, memory_efficient, missing, use_mailman, nthreads, Nz);
				output = compute_XXz(num_snp,all_zb);

				MatrixXdr scaled_pheno;
				if (gen_by_env == true) { 
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

					}//end gxe computations
				}

				for (int z_index = 0 ; z_index < Nz ; z_index++){
					if(pass_num == 1){
						XXz.col((bin_index * 2*Nz) + Nz + z_index)+=output.col(z_index);   /// save whole sample
					}
					else //if(num_snp != len[bin_index])
						XXz.col((bin_index * 2*Nz) + z_index) = XXz.col((bin_index * 2*Nz) + Nz + z_index)-output.col(z_index);   /// save corresponding jack contrib

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
					output = compute_XXUz(num_snp); 
					for (int z_index = 0 ; z_index < Nz ; z_index++){
						if(pass_num == 1){
							XXUz.col((bin_index * 2*Nz) + Nz + z_index)+=output.col(z_index);   /// save whole sample
						}
						else //if(num_snp != len[bin_index])
							XXUz.col((bin_index * 2*Nz) + z_index) = XXUz.col((bin_index * 2*Nz) + Nz + z_index)-output.col(z_index);
						}	 
				}

				//compute yXXy
				double temp_yxxy;
				if(both_side_cov == false)
					temp_yxxy= compute_yXXy(num_snp,pheno);
				else
					temp_yxxy= compute_yVXXVy(num_snp);

				if(pass_num == 1){
					yXXy(bin_index,1)+=temp_yxxy;
				}
				else 
					yXXy(bin_index,0)= yXXy(bin_index,1)-temp_yxxy;

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

					std::vector< std::vector<int> >().swap(g.p);
					std::vector< std::vector<int> >().swap(allgen_mail[bin_index].p);
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
		}

		if(pass_num == 2){
			for(int l = 0 ; l < T_Nbin ; l++)
				if( len[l]==jack_bin[jack_index][l])
			jack_bin[jack_index][l]=0;

			for (int i = 0 ; i < T_Nbin ; i++){
				b_trk(i,0) = Nindv_mask;

				if(i>=(T_Nbin-(nongen_Nbin + Nenv)) ){
					B1 = XXz.block(0,(i * 2*Nz),Nindv,Nz);
					B1 =all_zb.array() * B1.array();
					b_trk(i,0) = B1.sum() / (len[i]-jack_bin[jack_index][i]) / Nz;
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
				outfile << "Number of individuals after filtering: " << Nindv_mask << endl;
				outfile << "Number of covariates: " << Ncov << endl;
				outfile << "Number of environments: " << Nenv << endl;

				if (verbose == true) {
					cout << "LHS of Normal Eq" << endl << X_l << endl;
					cout << "RHS of Normal Eq" << endl << Y_r << endl;
					double relative_error = (X_l * herit - Y_r).norm() / Y_r.norm(); // norm() is L2 norm
					cout << "The relative error is: " << relative_error << endl;
					JacobiSVD<MatrixXd> svd(X_l);
					double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
					cout << "condition number:  "<< cond << endl;

					outfile << "LHS of Normal Eq" << endl << X_l << endl;
					outfile << "RHS of Normal Eq" << endl << Y_r << endl;
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

	}//end loop over jack
	cout << "Finished reading and computing over all blocks" << endl;
	cout << endl;
	cout << endl;

	if(pass_num == 1){
		cout << "pass = " << pass_num << endl;
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
		cout << "Size of bins :" << endl;
		if (hetero_noise == true) {
			for(int i = 0 ; i < Nbin + nongen_Nbin + Nenv ; i++)
				cout << "bin " << i<<" : " << len[i]<<endl;
		} else {
			for(int i = 0 ; i < Nbin + nongen_Nbin ; i++)
				cout << "bin " << i<<" : " << len[i]<<endl;
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

	}

	if(pass_num == 2){
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
		}     


		X_l << A_trs,b_trk,b_trk.transpose(),NC;
		Y_r << c_yky,yy;
		if (trace){
			trace_file << X_l.block(0, 0, Nbin, Nbin);
			for (int j = 0; j< Nbin; j++)
				trace_file << "," << len[j];
			trace_file << endl;
		}

		herit = X_l.colPivHouseholderQr().solve(Y_r);

		if (verbose == true ) { 
			cout << "Xl" << endl << X_l << endl;
			cout << "Yl" << endl << Y_r << endl;
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
	Nsnp = read_bim (f1.str());

	// Read phenotype and save the number of indvs
	string filename = command_line_opts.PHENOTYPE_FILE_PATH;
	Nindv = count_pheno (filename);
	read_pheno2 (Nindv, filename);
	cout << "Number of individuals : "<< Nindv << endl;
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
	std::string covname="";
	if(covfile!=""){
		use_cov = true;
		Ncov = read_cov(Nindv, covfile);
	} else if (covfile==""){
		cout << "No Covariate File Specified" << endl;

		if (cov_add_intercept == true) {
			covariate.resize(Nindv,1);
			for(int i = 0 ; i < Nindv ; i++)
				covariate(i,0) = 1;
				Ncov = 1;
		} else {
			both_side_cov = false;
			use_cov = false;
		}
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
		cout << "Number of features in bin " << i<<" : " << len[i]<<endl;
		outfile << "Number of features in bin " << i<<" : " << len[i]<<endl;
	}
	cout << "*****" << endl;
	outfile << "*****" << endl;
	cout << "Number of G variance components : " << Nbin << endl;
	cout << "Number of GxE variance components : " << nongen_Nbin << endl;
	outfile << "Number of G variance components : " << Nbin << endl;
	outfile << "Number of GxE variance components : " << nongen_Nbin << endl;
	if (hetero_noise == true) {
		cout << "Number of NxE variance components : " << Nenv << endl;
		outfile << "Number of NxE variance components : " << Nenv << endl;
	} else {
		cout << "Number of NxE variance components : 0" << endl;
		outfile << "Number of NxE variance components : 0" << endl;
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
	outfile << "#          GENIE (v1.0.0)        #" << endl;
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
	if (verbose) {
		outfile << "Other options: " << endl;
		outfile << "\t - norm_proj_pheno (normalize pheno after projection on covariates) " << std::to_string(command_line_opts.normalize_proj_pheno) << endl;
		outfile << "\t - cov_add_intercept (intercept term added to covariates) " << std::to_string(command_line_opts.cov_add_intercept) << endl;
		outfile << "\t - v (verbose) " << std::to_string(command_line_opts.verbose) << endl;
	}
	if (trace)
		outfile << "\t - tr (trace summaries) " << endl;

	outfile << endl;
	outfile << endl;
}


int main(int argc, char const *argv[]){

	parse_args (argc,argv);
	verbose = command_line_opts.verbose;
	trace = command_line_opts.print_trace;

	////////////////////////////////////////////
	///////////////////////////////////////////
	//int B = command_line_opts.batchNum;
	k_orig = command_line_opts.num_of_evec ;
	debug = command_line_opts.debugmode ;
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = false;
	k = k_orig + command_line_opts.l;
	k = (int)ceil(k / 10.0) * 10;
	command_line_opts.l = k - k_orig;
	srand((unsigned int) time(0));
	Nz = command_line_opts.num_of_evec;
	k = Nz;
	Njack = command_line_opts.jack_number;
	int seed = command_line_opts.seed;
	if (command_line_opts.exannot == true)
		Annot_x_E = true;
	nthreads = command_line_opts.nthreads;
	////
	//////////////////////////// Read multi genotypes
	string line;
	int num_files = 0;
	cov_add_intercept = command_line_opts.cov_add_intercept;

	read_auxillary_files ();
	regress_covariates ();	

	if (gen_by_env == true) {    
		///mask out indv with missingness from Enviro
		Enviro = Enviro.array().colwise() * mask.col(0).array();
	}

	//define random vector z's
	all_zb= MatrixXdr::Random(Nindv,Nz);

	std::mt19937 seedr;
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
		all_Uzb.resize(Nindv,Nz);
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

	if(hetero_noise == true) {
		XXz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin + Nenv) * Nz * 2);

		if(both_side_cov == true){
			UXXz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin + Nenv) * Nz * 2);
			XXUz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin + Nenv) * Nz * 2);
		}
		yXXy = MatrixXdr::Zero(Nbin + nongen_Nbin + Nenv,2);
	} else {
		XXz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin) * Nz * 2);

		if(both_side_cov == true){
			UXXz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin) * Nz * 2);
			XXUz = MatrixXdr::Zero(Nindv,(Nbin + nongen_Nbin) * Nz * 2);
		}
		yXXy = MatrixXdr::Zero(Nbin + nongen_Nbin,2);
	}

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
	genotype_stream_pass(name,1);
	genotype_stream_pass(name,2);

	print_results ();
	////////////////////////////////////////////////////////
	return 0;
}
