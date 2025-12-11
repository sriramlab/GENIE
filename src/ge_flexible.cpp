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

#include "auxillary.h"
#include "genotype.h"
#include "genomult.h"
#include "arguments.h"
#include "storage.h"
#include "matmult.h"
#include "io.h"
#include "std.h"
#include "functions.h"
#include "vectorfn.h"
#include "statsfn.h"
#include "sample_matching.h"

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

// If memeff = true, uses the memory efficient version
bool memeff = false;

// If opt1 = true, writes to a file after pass 1 
// Reads from file in pass 2
bool opt1 = true;

// If opt2 = true, reads groups of SNPs in block (termed read block)
// Number of SNPs in a read block can be controlled based on memory constraints by setting mem_Nsnp
bool opt2 = true ;
int mem_Nsnp = 10;
int mem_alloc = -1;

// Number of read blocks
// Approximately: Nsnp/mem_Nsnp
// Will  vary depending on the structure of jackknife blocks
int Nreadblocks; 

//ENV
// Partition the GxE component with respect to the annotations in the annotation file
bool Annot_x_E = false;

// Add environment variables to covariates
bool add_env_to_cov = true;

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
MatrixXdr geno_matrix; //(p,n)
int k,p,n;

MatrixXdr means; //(p,1)
MatrixXdr stds; //(p,1)
MatrixXdr sum2;
MatrixXdr sum;  

////////
//related to phenotype	
double y_sum; 
double y_mean;

options command_line_opts;

bool debug = false;
bool var_normalize = false;
bool memory_efficient = false;
bool missing = false;
bool fast_mode = true;
bool use_cov = false; 

// Jackknife-related variables
// Different approaches for defining the jackknife blocks
// jack_scheme  =
// 1: each block has constant number of SNPs (except the last). a block might span distinct chromosomes
// 2: each block has constant number of SNPs with each block contained within a single chromosome
// 3: each block has constant physical length with each block contained within a single chromosome. 
// For 3: physical length is specified in Mbs.
// Default : 1
int jack_scheme = 1;

// Number of jackknife blocks
// When jack_scheme = 2 or 3, this number will be set to 
// actual number based on the number of SNPs
int Njack = 1000;

// Size of jackknife block (in 1 Mb)
// Relevant for jack_scheme = 3
int jack_size;

int step_size;
int step_size_rem;

bool use_ysum;
bool keep_xsum;

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
int Nindv;
int Nsnp;
//
// Number of environmental variables
int Nenv = 0;

// Number of genetic bins (annotations)
// Set in read_annot
int Nbin = 1;

// Number of non-genetic bins
int nongen_Nbin = 0;

// Number of genetic bins
int gen_Nbin = 0;

// Total number of bins (related to the number of annotations)
// Number of VCs = T_Nbin  + 1 (sigma_e)
// (1 for G model with single annotation)
// For heterogeneous noise: T_Nbin = Nbin + nongen_Nbin + Nenv;
// For no heterogeneous noise: T_Nbin = Nbin + nongen_Nbin;
int T_Nbin ;

// Number of covariates
int Ncov;
int Nindv_mask;
int NC;

// Number of random vectors
int Nz = 10;


// Sum of number of SNPs assigned to each annotation 
// Can be greater than or less than or equal to the number of SNPs in the bim file (for overlapping annotations)
int Nsnp_annot = 0;

// Annotation matrix for binary annotations (Number of SNPs x Number of annotations) 
vector<vector<bool> > annot_bool;

// Number of SNPs in each bin
vector<int> len;

// Number of SNPs in each bin in a jackknife block
// Njack X (Total number of annotations (depends on the specific model used)
vector<vector<int> > jack_bin;
vector<int>  jack_block_size;
vector<int>  snp_to_jackblock;

vector<vector<int> > read_bin;

// Assigned in setup_read_blocks
vector<int> jack_to_read_block;
vector<int> snp_to_read_block;

vector<vector<int> > mem_bin;
vector<int>  mem_block_size;

vector <data> allgen;
vector <genotype> allgen_mail;
int global_snp_index;
bool use_mailman = true;

///reading single col annot
vector<int> SNP_annot;
bool use_1col_annot = false;


///Variables for reg out cov on both side of LM
bool both_side_cov = true;
MatrixXdr UXXz;
MatrixXdr XXUz;
MatrixXdr Xz;
MatrixXdr trVK;

//CHANGE (2/27)
bool gen_by_env;
//CHANGE(10/20)
bool hetero_noise;
bool cov_add_intercept;
int verbose;
bool trace;
bool use_dummy = false;
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

int phenocount = 0;

std::istream& newline(std::istream& in)
{
	if ((in >> std::ws).peek() != std::char_traits<char>::to_int_type('\n')) {
		in.setstate(std::ios_base::failbit);
	}
	return in.ignore();
}

void initial_var(){
	p = g.Nsnp;
	n = g.Nindv;

	sum2.resize(p,1);
	sum.resize(p,1);

	if(!fast_mode && !memory_efficient){
		geno_matrix.resize (p,n);
		g.generate_eigen_geno (geno_matrix, var_normalize);
	}

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

// Needed for the memory efficient vesion (opt1 = true && opt2 = true).
// This is only used when memeff == 1
// Needed to handle the cases where jackknife blocks do not all have equal number of SNPs
// Sets up data structures that map SNPs and jackknife blocks to read blocks
void setup_read_blocks ()  {
	if (verbose >= 3) {
		cout << "In setup_read_blocks" << endl;
	}

	int snpindex = 0;
	int blockindex = 0 ;

    // How many read blocks in each jackknife block
	jack_to_read_block.resize (Njack);		

    // Which read block each SNP belongs to
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
	if (verbose >= 3) {
		cout << "Finished setting up read blocks" << endl;
	}
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

	MatrixXdr A_trs(T_Nbin,T_Nbin);
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
	jack.resize(T_Nbin + 1, Njack);

	// Matrix of variance components
	// Includes sigma_e. Hence number of entries = T_Nbin + 1
	point_est.resize(T_Nbin + 1, 1);

	point_est_adj_gxe = MatrixXdr::Zero(T_Nbin + 3, 1);
	jack_adj_gxe = MatrixXdr::Zero(T_Nbin + 3, Njack);

	enrich_jack.resize(T_Nbin, Njack);
	enrich_point_est.resize(T_Nbin, 1);


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

	cout << "Number of individuals without missing phenotype and enviromental variables = " << mask.sum() << endl;
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

	MatrixXdr A_trs(T_Nbin,T_Nbin);
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
	jack.resize(T_Nbin + 1, Njack);

	// Matrix of variance components
	// Includes sigma_e. Hence number of entries = T_Nbin + 1
	point_est.resize(T_Nbin + 1, 1);

	point_est_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,1);
	jack_adj_gxe = MatrixXdr::Zero(T_Nbin + 3,Njack);

	enrich_jack.resize(T_Nbin, Njack);
	enrich_point_est.resize(T_Nbin, 1);


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
			}

			wgxsum_ifs.close ();
		}
	}

	cout << "************ Making pass number "<< pass_num << " over genotypes ************" << endl;
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
				if (verbose >= 3) {
					cout << "In gxe" << endl;
				}
				for (int env_index = 0 ; env_index < Nenv ; env_index++){


					if (opt1  && pass_num == 2) {
						read_matrix (xsum_ifs, output_env);
					} else {
						MatrixXdr env_all_zb = all_zb.array().colwise() * Enviro.col(env_index).array();
						output_env = compute_XXz(num_snp, env_all_zb);
						output_env = output_env.array().colwise() * Enviro.col(env_index).array();
						if (opt1) {
							write_matrix (xsum_ofs, output_env);
						}	
					}
					

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
						else 
						// if(num_snp != len[bin_index])
							XXz.col((gxe_bin_index * 2 * Nz) + z_index) = XXz.col((gxe_bin_index * 2*Nz) + z_index) - output_env.col (z_index);   /// save corresponding jack contrib

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

						if (opt1  && pass_num == 2) {
							read_matrix (xsum_ifs, output_env);
						} else {
							MatrixXdr env_all_Uzb = all_Uzb.array().colwise() * Enviro.col(env_index).array();
							output_env = compute_XXz(num_snp, env_all_Uzb);
							output_env = output_env.array().colwise() * Enviro.col(env_index).array();
							if (opt1) {
								write_matrix (xsum_ofs, output_env);
							}	
						}

						for (int z_index = 0 ; z_index < Nz ; z_index++){
							if(pass_num == 1){
								XXUz.col((gxe_bin_index * 2 * Nz) + Nz + z_index) += output_env.col(z_index);   /// save whole sample
								XXUz.col((gxe_bin_index * 2 * Nz) + z_index) += output_env.col(z_index); 
							}
							else
								XXUz.col((gxe_bin_index * 2 * Nz) + z_index) = XXUz.col((gxe_bin_index * 2*Nz) + z_index)-output_env.col(z_index);
						}
					}
					if(both_side_cov == true)
						scaled_pheno = new_pheno.array() * Enviro.col(env_index).array();
					else
						scaled_pheno= pheno.array() * Enviro.col(env_index).array();


					double temp_e_yxxy;
					if (opt1 && pass_num == 2)  {
						if (use_ysum) 
							ysum_ifs.read((char *) (&temp_e_yxxy), sizeof(double)); 
						else
							temp_e_yxxy = jack_yXXy(gxe_bin_index, jack_index);
					} else {

						temp_e_yxxy= compute_yXXy(num_snp, scaled_pheno);

						if (opt1) {
							jack_yXXy(gxe_bin_index, jack_index) = temp_e_yxxy;
							ysum_ofs.write((char *) (&temp_e_yxxy), sizeof(double));
						}
					}

					if(pass_num == 1){
						yXXy(gxe_bin_index, 1) += temp_e_yxxy;
						yXXy(gxe_bin_index, 0) += temp_e_yxxy;
					}
					else 
						yXXy(gxe_bin_index, 0) = yXXy(gxe_bin_index,0) - temp_e_yxxy;      

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
					else // if(num_snp != len[bin_index])
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
				yXXy(bin_index, 1) += temp_yXXy;
			}
			else 
				yXXy(bin_index, 0)= yXXy(bin_index , 1) - temp_yXXy;


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
				for (int i=0; i< Nbin; i++){
						for (int j=0; j < Nbin; j++){
								trace_file << (X_l(i, j) - Nindv_mask) * (len[i] - jack_bin[jack_index][i])*(len[j] - jack_bin[jack_index][j])/pow(Nindv_mask, 2)
								<< ",";
						}
						trace_file << len[i] - jack_bin[jack_index][i] << endl;
				}
			}

		} //end if pass_num = 2

	}//end loop over jackknife blocks
	cout << "Finished reading and computing over all blocks" << endl;
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

		cout << "Number of individuals without missing phenotype and environmental variables: " << mask.sum() << endl;
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
			for (int i=0; i< Nbin; i++){
				for (int j=0; j < Nbin; j++){
					trace_file << (X_l(i, j) - Nindv_mask) * len[i]*len[j]/pow(Nindv_mask,2) << ",";
				}
				trace_file << len[i] << endl;
			}
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
	gen_Nbin = Nbin;
}


/* This is the key code in ge_hetero_flexible.cpp*/
void genotype_stream_single_pass (string name) {
	if (verbose >= 3) {
		cout << "both_side_cov = " << both_side_cov << endl;
		cout << "cov_add_intercept = " << cov_add_intercept << endl;
		cout << "T_Nbin = " << T_Nbin << " Nbin = " << Nbin << " nongen_Nbin = " << nongen_Nbin << endl;
		cout << "Nz = " << Nz << endl;
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

	for (int jack_index = 0; jack_index < Njack; jack_index ++){
		int read_Nsnp = jack_block_size[jack_index];	
		cout << "Reading jackknife block " << jack_index << endl;
		if (verbose >= 1)  {
			cout << "************Reading jackknife block " << jack_index << " ************" <<endl;
			if (verbose >= 2)
				cout << "read_Nsnp = " << read_Nsnp << endl;
		}


		//		int read_Nsnp = (jack_index<(Njack-1)) ? (step_size) : (step_size+step_size_rem);
		if(use_mailman == true){
			for (int i = 0; i < Nbin; i++){
				allgen_mail[i].segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
				allgen_mail[i].Nsegments_hori = ceil(jack_bin[jack_index][i]*1.0/(allgen_mail[i].segment_size_hori*1.0));
				allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,std::vector<int>(Nindv));
				allgen_mail[i].not_O_i.resize(jack_bin[jack_index][i]);
				allgen_mail[i].not_O_j.resize(Nindv);
				allgen_mail[i].index = 0;
				allgen_mail[i].Nsnp = jack_bin[jack_index][i];
				allgen_mail[i].Nindv = Nindv;

				allgen_mail[i].columnsum.resize(jack_bin[jack_index][i],1);
				for (int index_temp = 0; index_temp < jack_bin[jack_index][i]; index_temp++)
					allgen_mail[i].columnsum[index_temp] = 0;
			}
		}
		else{
			for (int bin_index = 0;bin_index < Nbin; bin_index++){
				allgen[bin_index].gen.resize(jack_bin[jack_index][bin_index],Nindv);
				allgen[bin_index].index = 0;
			}
		}

		if(use_1col_annot == true)
			read_bed_1colannot (ifs, missing, read_Nsnp);
		else
			read_bed2 (ifs, missing, read_Nsnp);
		read_header = false;

		for (int bin_index = 0; bin_index < Nbin; bin_index++){
			int num_snp;
			if (use_mailman == true)
				num_snp = allgen_mail[bin_index].index;
			else
				num_snp = allgen[bin_index].index;

			if(num_snp != 0){
				stds.resize(num_snp, 1);
				means.resize(num_snp, 1);

				if(use_mailman == true){
					for (int i = 0; i < num_snp; i++)
						means(i,0) = (double)allgen_mail[bin_index].columnsum[i]/Nindv;
				}			
				else	  
					means = allgen[bin_index].gen.rowwise().mean();


				for (int i = 0; i < num_snp; i++)
					stds(i,0) = 1/sqrt((means(i,0)*(1-(0.5*means(i,0)))));

				if (use_mailman == true){
					g = allgen_mail[bin_index];
					g.segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
					g.Nsegments_hori = ceil(jack_bin[jack_index][bin_index]*1.0/(g.segment_size_hori*1.0));
					g.p.resize(g.Nsegments_hori,std::vector<int>(Nindv));
					g.not_O_i.resize(jack_bin[jack_index][bin_index]);
					g.not_O_j.resize(Nindv);
					initial_var();

				}else{
					gen = allgen[bin_index].gen;

				} 
				mm = MatMult(g, gen, debug, var_normalize, memory_efficient, missing, use_mailman, nthreads, Nz);
				// cout << "here1" << endl;
				output = compute_XXz(num_snp,all_zb);
				// cout << "here2" << endl;

				///gxe computations
				MatrixXdr scaled_pheno;
				if (gen_by_env == true) {
					for (int env_index = 0; env_index < Nenv; env_index++){
						MatrixXdr env_all_zb = all_zb.array().colwise()*Enviro.col(env_index).array();                                   
						output_env = compute_XXz (num_snp, env_all_zb);                                   
						output_env = output_env.array().colwise()*Enviro.col(env_index).array();

						int gxe_bin_index;
						if(Annot_x_E == true)  
							gxe_bin_index = Nbin+(env_index*Nbin)+bin_index;
						else
							gxe_bin_index = Nbin+env_index;			

						for (int z_index = 0; z_index < Nz; z_index++){
							XXz.col(((gxe_bin_index)*(Njack+1)*Nz)+(jack_index*Nz)+z_index) += output_env.col(z_index);	
							XXz.col(((gxe_bin_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index) += output_env.col(z_index);

							if(both_side_cov == true) {
								vec1 = output_env.col(z_index);
								w1 = covariate.transpose() * vec1;
								w2 = Q * w1;
								w3 = covariate * w2;
								//if(num_snp!=len[bin_index])
								UXXz.col(((gxe_bin_index)*(Njack+1)*Nz)+(jack_index*Nz)+z_index) += w3;
								UXXz.col(((gxe_bin_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index) += w3;
							}

						}


						if (both_side_cov==true){
							MatrixXdr env_all_Uzb = all_Uzb.array().colwise()*Enviro.col(env_index).array();
							output_env = compute_XXz(num_snp,env_all_Uzb);
							output_env = output_env.array().colwise()*Enviro.col(env_index).array();

							for (int z_index = 0; z_index < Nz; z_index++){
								XXUz.col(((gxe_bin_index)*(Njack+1)*Nz)+(jack_index*Nz)+z_index) += output_env.col(z_index);
								XXUz.col(((gxe_bin_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index) += output_env.col(z_index);   /// save whole sample
							}
						}
						if(both_side_cov==true)
							scaled_pheno = new_pheno.array()*Enviro.col(env_index).array();
						else
							scaled_pheno = pheno.array()*Enviro.col(env_index).array();

						double temp = compute_yXXy(num_snp,scaled_pheno);
						yXXy(gxe_bin_index,jack_index) += temp;
						yXXy(gxe_bin_index,Njack) += temp;
					}
				}
				////end gxe computation

				for (int z_index = 0; z_index < Nz; z_index++){
					// per-block contribution (overwrite)
					XXz.col((bin_index*(Njack+1)*Nz) + (jack_index*Nz) + z_index) = output.col(z_index);
					// whole-sample accumulator (+=)
					XXz.col((bin_index*(Njack+1)*Nz) + (Njack*Nz)   + z_index) += output.col(z_index);
					if(both_side_cov == true) {
						vec1 = output.col(z_index);
						w1 = covariate.transpose()*vec1;
						w2 = Q * w1;
						w3 = covariate * w2;
						// if(num_snp != len[bin_index])
							UXXz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index) = w3;
						UXXz.col((bin_index*(Njack+1)*Nz)+(Njack*Nz)+z_index) += w3;
					}

				}

				if (both_side_cov == true){
					output=compute_XXUz(num_snp); 
					for (int z_index = 0; z_index < Nz; z_index++){
						// if(num_snp != len[bin_index])
							XXUz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index) = output.col(z_index);
						XXUz.col((bin_index*(Njack+1)*Nz)+(Njack*Nz)+z_index) += output.col(z_index);   /// save whole sample
					}
				}

				//compute yXXy

				if(both_side_cov == false)
					yXXy(bin_index, jack_index) = compute_yXXy(num_snp,pheno);
				else
					yXXy(bin_index, jack_index) = compute_yVXXVy(num_snp);

				yXXy(bin_index, Njack) += yXXy(bin_index,jack_index);

				if(num_snp == len[bin_index])
					yXXy(bin_index, jack_index) = 0;



				if (verbose >= 2) {
					cout << jack_index << " " << bin_index << "\tXXz(" << XXz.rows() <<"," << XXz.cols() << ") "<< XXz.sum() << endl;
					cout << jack_index << " " << bin_index << "\tyXXy(" << yXXy.rows() <<"," << yXXy.cols() << ") "<< yXXy.sum() << endl;
					if (yXXy.rows () > 1)
						cout << jack_index <<" " << bin_index << "\tyXXy " << yXXy(0,1) <<"," << yXXy(1,1) << endl;
				}

				//compute Xz

				if (verbose >= 2) {
					cout << num_snp << " SNPs in bin "<< bin_index<< " of jackknife block  " << jack_index << endl;   
					cout << "Reading and computing bin " << bin_index << "  of jackknife block " << jack_index << " completed" << endl;
				}
				mm.clean_up();
				if(use_mailman==true){
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
					std::vector< std::vector<int> >().swap(g.not_O_j);
					std::vector< std::vector<int> >().swap(g.not_O_i);
					std::vector< std::vector<int> >().swap(allgen_mail[bin_index].p);
					std::vector< std::vector<int> >().swap(allgen_mail[bin_index].not_O_j);
					std::vector< std::vector<int> >().swap(allgen_mail[bin_index].not_O_i);
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
	}//end loop over jackknife blocks
	cout << "Finished reading and computing over all blocks" << endl;
	cout << endl;


	if (hetero_noise == true) {
		MatrixXdr hetro_all_Uzb;
		for (int env_index=0;env_index<Nenv;env_index++){
			/// add hetero env noise
			MatrixXdr hetro_all_zb=all_zb.array().colwise()*Enviro.col(env_index).array();
			hetro_all_zb=hetro_all_zb.array().colwise()*Enviro.col(env_index).array();

			if(both_side_cov==true){
				hetro_all_Uzb=all_Uzb.array().colwise()*Enviro.col(env_index).array();
				hetro_all_Uzb=hetro_all_Uzb.array().colwise()*Enviro.col(env_index).array();
			}

			int hetro_index;
			if(Annot_x_E==true)
				hetro_index=Nbin+(Nenv*Nbin)+env_index;
			else
				hetro_index=Nbin+Nenv+env_index;
			for (int z_index=0;z_index<Nz;z_index++){
				XXz.col(((hetro_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index)=hetro_all_zb.col(z_index);		 
				if(both_side_cov==true){
					vec1=hetro_all_zb.col(z_index);
					w1=covariate.transpose()*vec1;
					w2=Q*w1;
					w3=covariate*w2;
					UXXz.col(((hetro_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index)=w3;
					XXUz.col(((hetro_index)*(Njack+1)*Nz)+(Njack*Nz)+z_index)=hetro_all_Uzb.col(z_index);
				}

			}

			MatrixXdr scaled_pheno;
			if(both_side_cov==true)
				scaled_pheno=new_pheno.array()*Enviro.col(env_index).array();
			else
				scaled_pheno= pheno.array()*Enviro.col(env_index).array();

			yXXy(hetro_index,Njack)=(scaled_pheno.array()*scaled_pheno.array()).sum();
			len.push_back(1);
		} 
	}

	cout<<"Size of bins :"<<endl;
	//CHANGE(10/20)
	if (hetero_noise == true) {
		for(int i = 0; i < Nbin + nongen_Nbin + Nenv; i++)
			cout << "Bin " << i << " : " << len[i] << endl;
	} else {
		for(int i = 0; i < Nbin + nongen_Nbin; i++)
			cout << "Bin " << i << " : " << len[i] << endl;
	}

	cout << "Number of individuals without missing phenotype and enviromental variables: " << mask.sum() << endl;
	cout << endl;

	gen_Nbin = Nbin;

	//CHANGE(10/20)
	if (hetero_noise == true)
		Nbin = Nbin + nongen_Nbin + Nenv;
	else
		Nbin = Nbin + nongen_Nbin;

	for(int bin_index = 0; bin_index < Nbin; bin_index++){
		for(int jack_index = 0; jack_index < Njack; jack_index++){
			for (int z_index = 0; z_index < Nz; z_index++){
				MatrixXdr v1 = XXz.col((bin_index*(Njack+1)*Nz)+(Njack*Nz)+z_index);
				MatrixXdr v2 = XXz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index);
				XXz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index) = v1-v2;
				if(both_side_cov == true){
					v1 = XXUz.col((bin_index*(Njack+1)*Nz)+(Njack*Nz)+z_index);
					v2 = XXUz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index);
					XXUz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index) = v1-v2;

					v1 = UXXz.col((bin_index*(Njack+1)*Nz)+(Njack*Nz)+z_index);
					v2 = UXXz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index);
					UXXz.col((bin_index*(Njack+1)*Nz)+(jack_index*Nz)+z_index) = v1-v2;
				}    
			}
			yXXy(bin_index,jack_index) = yXXy(bin_index,Njack)-yXXy(bin_index,jack_index);
		}
	}
	//// all XXy and yXXy and contributions of every jackknife subsamples  were computed till this line.

	/// normal equations LHS
	MatrixXdr  A_trs(Nbin,Nbin);
	MatrixXdr b_trk(Nbin,1);
	MatrixXdr c_yky(Nbin,1);

	MatrixXdr X_l(Nbin+1,Nbin+1);
	MatrixXdr Y_r(Nbin+1,1);

	int jack_index=Njack;
	MatrixXdr B1;
	MatrixXdr B2;
	MatrixXdr C1;
	MatrixXdr C2;
	double trkij;
	double yy = (pheno.array() * pheno.array()).sum();

	if(both_side_cov == true)
		yy = (new_pheno.array()*new_pheno.array()).sum();

	int Nindv_mask = mask.sum();
	if(both_side_cov == true)
		NC = Nindv_mask - Ncov;
	else
		NC = Nindv_mask;

	point_est_adj_gxe=MatrixXdr::Zero(Nbin+3,1);
	jack_adj_gxe=MatrixXdr::Zero(Nbin+3,Njack);


	jack.resize(Nbin+1,Njack);
	point_est.resize(Nbin+1,1);

	enrich_jack.resize(Nbin,Njack);
	enrich_point_est.resize(Nbin,1);

	MatrixXdr h1;
	MatrixXdr h2;
	MatrixXdr h3;

	double trkij_res1;
	double trkij_res2;
	double trkij_res3;
	double tk_res;


	for (jack_index = 0; jack_index <= Njack ; jack_index ++){
		for(int k = 0; k < Nbin; k++)
			if( jack_index < Njack && len[k] == jack_bin[jack_index][k])
				jack_bin[jack_index][k] = 0;

		for (int i = 0;i < Nbin; i++){
			b_trk(i,0) = Nindv_mask;
			// CHANGE (2/17)
			int sum_num_nongen_bin = 0;
			if (hetero_noise == true) 
				sum_num_nongen_bin = nongen_Nbin + Nenv;
			else 
				sum_num_nongen_bin = nongen_Nbin;

			if(i>=(Nbin-sum_num_nongen_bin) ){

				B1=XXz.block(0,(i*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
				B1 =all_zb.array()*B1.array();

				if(jack_index==Njack)
					b_trk(i,0)=B1.sum()/len[i]/Nz;
				else
					b_trk(i,0)=B1.sum()/(len[i]-jack_bin[jack_index][i])/Nz;

			}

			if(jack_index==Njack)
				c_yky(i,0)=yXXy(i,jack_index)/len[i];
			else
				c_yky(i,0)=yXXy(i,jack_index)/(len[i]-jack_bin[jack_index][i]);


			if(both_side_cov==true){
				B1=XXz.block(0,(i*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
				C1=B1.array()*all_Uzb.array();
				C2=C1.colwise().sum();	
				tk_res=C2.sum();  

				if(jack_index==Njack)
					tk_res=tk_res/len[i]/Nz;
				else
					tk_res=tk_res/(len[i]-jack_bin[jack_index][i])/Nz;

				b_trk(i,0)=b_trk(i,0)-tk_res;
			}

			for (int j=i;j<Nbin;j++){
				B1=XXz.block(0,(i*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
				B2=XXz.block(0,(j*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
				C1=B1.array()*B2.array();
				C2=C1.colwise().sum();
				trkij=C2.sum();


				if(both_side_cov==true){

					h1=covariate.transpose()*B1;
					h2=Q*h1;
					h3=covariate*h2;
					C1=h3.array()*B2.array();
					C2=C1.colwise().sum();
					trkij_res1=C2.sum();

					B1=XXUz.block(0,(i*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
					B2=UXXz.block(0,(j*(Njack+1)*Nz)+(jack_index*Nz),Nindv,Nz);
					C1=B1.array()*B2.array();
					C2=C1.colwise().sum();
					trkij_res3=C2.sum();


					trkij+=trkij_res3-trkij_res1-trkij_res1 ;
				}

				if(jack_index==Njack)
					trkij=trkij/len[i]/len[j]/Nz;
				else
					trkij=trkij/(len[i]-jack_bin[jack_index][i])/(len[j]-jack_bin[jack_index][j])/Nz;
				A_trs(i,j)=trkij;
				A_trs(j,i)=trkij;

			}
		}


		X_l<<A_trs,b_trk,b_trk.transpose(),NC;
		Y_r<<c_yky,yy;

		if (trace){
			if (jack_index < Njack){
				for (int i=0; i< Nbin; i++){
					for (int j=0; j < Nbin; j++){
						trace_file << (X_l(i, j) - Nindv_mask) * (len[i] - jack_bin[jack_index][i])*(len[j] - jack_bin[jack_index][j])/pow(Nindv_mask, 2)
							<< ",";
					}
					trace_file << len[i] - jack_bin[jack_index][i] << endl;
				}
			}
			else{
				for (int i=0; i< Nbin; i++){
					for (int j=0; j < Nbin; j++){
						trace_file << (X_l(i, j) - Nindv_mask) * len[i]*len[j]/pow(Nindv_mask,2) << ",";
					}
					trace_file << len[i] << endl;
				}

			}
		}

		MatrixXdr herit = X_l.fullPivHouseholderQr().solve(Y_r);


		if(jack_index == Njack){
			outfile << "Number of individuals after filtering: " << Nindv_mask << endl;
			outfile << "Number of covariates: " << Ncov << endl;
			outfile << "Number of environments: " << Nenv << endl;
			if (verbose == true) {
				cout<<"LHS of Normal Eq"<<endl<<X_l<<endl;
				cout<<"RHS of Normal Eq"<<endl<<Y_r<<endl;
				double relative_error = (X_l*herit - Y_r).norm() / Y_r.norm(); // norm() is L2 norm
				cout << "The relative error is: " << relative_error << endl;
#ifdef USE_DOUBLE
				JacobiSVD<MatrixXd> svd(X_l);
#else
				JacobiSVD<MatrixXf> svd(X_l);
#endif

				double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
				cout<<"condition number:  "<< cond<<endl;

				outfile<<"LHS of Normal Eq"<<endl<<X_l<<endl;
				outfile<<"RHS of Normal Eq"<<endl<<Y_r<<endl;
				outfile<<"Normal Equations info:"<<endl;
				outfile<<"The relative error is: " << relative_error << endl;
				outfile<<"Max sing.val: "<<svd.singularValues()(0)<<endl;
				outfile<<"Min sing.val: "<<svd.singularValues()(svd.singularValues().size()-1)<<endl;
				outfile<<"Condition number: "<<cond<<endl;
				cout << endl;
				cout << endl;
				outfile << endl;
				outfile << endl;
			}
			for(int i = 0; i < (Nbin + 1); i++)
				point_est(i,0) = herit(i,0);

			//adj gxe

			for(int i=0;i<(Nbin+1);i++){
				if(i==Nbin)
					point_est_adj_gxe(i,0)=point_est(i,0)*NC;
				else
					point_est_adj_gxe(i,0)=point_est(i,0)*b_trk(i,0);
			} 
		}else{
			for(int i=0;i<(Nbin+1);i++)
				jack(i,jack_index)=herit(i,0);

			//adj gxe
			for(int i=0;i<(Nbin+1);i++){
				if(i==Nbin)
					jack_adj_gxe(i,jack_index)=jack(i,jack_index)*NC;
				else
					jack_adj_gxe(i,jack_index)=jack(i,jack_index)*b_trk(i,0);
			}
		}
	}//end of loop over jack

}



// Regress covariates from phenotype
// If no covariates are provided, center phenotype
//  
void regress_covariates () {
	if (verbose >= 3) {
		cout << "In regress_covariates" << endl;
	}

	bool normalize_proj_pheno = command_line_opts.normalize_proj_pheno;
	if(use_cov == true){
		if (verbose >= 3) {
			cout << "Regressing covariates" << endl;
		}
		MatrixXdr mat_mask = mask.replicate(1, Ncov);
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
		if (verbose >= 3) {
			cout << "Centering phenotype" << endl;
		}

		y_sum = pheno.sum();
		y_mean = y_sum / mask.sum();
		for(int i = 0; i < Nindv; i++){
			if(mask(i,0)!=0)
				pheno(i,0) = pheno(i,0) - y_mean; //center phenotype
		}
		y_sum = pheno.sum();
	}
	if (verbose >= 3) {
		cout << "Finished regressing covariates" << endl;
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
	cout << "Number of SNPs = "<< Nsnp << endl;

	setup_read_blocks();

	// Get file paths
	use_dummy = command_line_opts.use_dummy_pheno;
	std::stringstream f0;
	f0 << geno_name << ".fam";
	string fam_path = f0.str();
	string pheno_path = command_line_opts.PHENOTYPE_FILE_PATH;
	string cov_path = command_line_opts.COVARIATE_FILE_PATH;
	string env_path = command_line_opts.ENV_FILE_PATH;
	string annot_path = command_line_opts.Annot_PATH;

	// Model to fit
	gen_by_env = command_line_opts.gen_by_env;
	hetero_noise = command_line_opts.hetero_noise;
	if (gen_by_env == false) {
		Nenv = 0;
		hetero_noise = false;
	}

	if (command_line_opts.no_match_ids) {
		// ============================================================
		// LEGACY MODE: Position-based matching (original behavior)
		// ============================================================
		cerr << "Warning: Using legacy position-based sample matching (--no-match-ids)." << endl;
		cerr << "  Results may be incorrect if files are not pre-aligned." << endl;

		int fam_lines = count_fam(fam_path);

		if (use_dummy) {
			Nindv = fam_lines;
			// Initialize dummy phenotype and mask
			pheno.resize(Nindv, 1);
			mask.resize(Nindv, 1);
			new_pheno.resize(Nindv, 1);
			pheno.setZero();
			mask.setOnes();
			new_pheno.setZero();
			phenocount = 1;
		} else {
			Nindv = count_pheno(pheno_path);
			read_pheno(Nindv, pheno_path);
		}

		if (fam_lines != Nindv) {
			exitWithError("Number of individuals in fam file and pheno file does not match ");
			exit(1);
		}

		cout << "Number of individuals = " << Nindv << endl;
		y_sum = pheno.sum();

		// Read environment file (if needed)
		if (gen_by_env) {
			Nenv = read_env(Nindv, env_path);
		}

		// Read annotation files
		if (use_1col_annot == true) {
			read_annot_1col(annot_path);
		} else {
			read_annot(annot_path);
		}

		// Read covariate file
		if (cov_path != "") {
			use_cov = true;
			Ncov = read_cov(Nindv, cov_path);
		} else {
			cout << "No covariate file specified" << endl;
			if ((cov_add_intercept == true) && (!use_dummy)) {
				covariate.resize(Nindv, 1);
				for (int i = 0; i < Nindv; i++)
					covariate(i, 0) = 1;
				Ncov = 1;
				use_cov = true;
				cout << "Intercept included" << endl;
			} else {
				both_side_cov = false;
				use_cov = false;
				cout << "No intercept included" << endl;
			}
		}

	} else {
		// ============================================================
		// NEW MODE: ID-based sample matching (default)
		// BUG-003 fix: Match samples by FID/IID across all input files
		// ============================================================

		// Step 1: Read FAM IDs (master sample list)
		FamIndex fam_idx = read_fam_ids(fam_path);
		cout << "Read " << fam_idx.count << " samples from FAM file" << endl;

		// Step 2: Initialize tracking
		SampleMatchResult match_result;
		match_result.fam_count = fam_idx.count;
		match_result.all_aligned = true;
		match_result.covar_count = 0;
		match_result.env_count = 0;
		match_result.fam_not_in_covar = 0;
		match_result.covar_not_in_fam = 0;
		match_result.fam_not_in_env = 0;
		match_result.env_not_in_fam = 0;

		// Step 3: Read environment file FIRST (if needed, before covariates)
		MatchStats env_stats = {0, 0, 0, true};
		if (gen_by_env) {
			Nenv = read_env_matched(fam_idx, env_path, env_stats);
			match_result.env_count = env_stats.total_in_file;
			match_result.fam_not_in_env = fam_idx.count - env_stats.matched;
			match_result.env_not_in_fam = env_stats.not_in_fam;
			if (!env_stats.already_aligned) match_result.all_aligned = false;
		}

		// Step 4: Read phenotype file
		MatchStats pheno_stats = {0, 0, 0, true};
		if (use_dummy) {
			Nindv = fam_idx.count;
			// Initialize dummy phenotype and mask
			pheno.resize(Nindv, 1);
			mask.resize(Nindv, 1);
			new_pheno.resize(Nindv, 1);
			pheno.setZero();
			mask.setOnes();
			new_pheno.setZero();
			phenocount = 1;
			pheno_stats = {fam_idx.count, fam_idx.count, 0, true};
		} else {
			pheno_stats = read_pheno_matched(fam_idx, pheno_path, phenocount);
			Nindv = fam_idx.count;  // Matrix sized to FAM count
		}
		match_result.pheno_count = pheno_stats.total_in_file;
		match_result.fam_not_in_pheno = fam_idx.count - pheno_stats.matched;
		match_result.pheno_not_in_fam = pheno_stats.not_in_fam;
		if (!pheno_stats.already_aligned) match_result.all_aligned = false;

		cout << "Number of individuals (FAM) = " << Nindv << endl;
		y_sum = pheno.sum();

		// Step 5: Read annotation file (unchanged - SNP-level, no sample IDs)
		if (use_1col_annot == true) {
			read_annot_1col(annot_path);
		} else {
			read_annot(annot_path);
		}

		// Step 6: Read covariate file
		MatchStats cov_stats = {0, 0, 0, true};
		if (cov_path != "") {
			use_cov = true;
			Ncov = read_cov_matched(fam_idx, cov_path, cov_stats);
			match_result.covar_count = cov_stats.total_in_file;
			match_result.fam_not_in_covar = fam_idx.count - cov_stats.matched;
			match_result.covar_not_in_fam = cov_stats.not_in_fam;
			if (!cov_stats.already_aligned) match_result.all_aligned = false;
		} else {
			cout << "No covariate file specified" << endl;
			if ((cov_add_intercept == true) && (!use_dummy)) {
				covariate.resize(Nindv, 1);
				for (int i = 0; i < Nindv; i++) {
					covariate(i, 0) = 1;
				}
				Ncov = 1;
				use_cov = true;
				cout << "Intercept included" << endl;
			} else {
				both_side_cov = false;
				use_cov = false;
				cout << "No intercept included" << endl;
			}
		}

		// Step 7: Compute intersection count and report
		match_result.intersection_count = compute_intersection_count();
		report_sample_matching(match_result);

		// Step 8: Validate intersection is not empty
		if (match_result.intersection_count == 0) {
			cerr << "ERROR: No samples in common across all input files." << endl;
			cerr << "  Hint: Check that FID/IID columns match across files." << endl;
			cerr << "        Sample IDs are case-sensitive." << endl;
			exit(1);
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
		T_Nbin = gen_Nbin;
	} else if (hetero_noise == true) {
		T_Nbin = gen_Nbin + nongen_Nbin + Nenv;
	} else {
		T_Nbin = gen_Nbin + nongen_Nbin;
	}

	MatrixXdr point_se;
	point_se = statsfn::jack_se(jack);

	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;

	for (int i = 0 ; i < T_Nbin ; i++){
		cout << "Number of SNPs in bin " << i<<" = " << len[i]<<endl;
		outfile << "Number of SNPs in bin " << i<<" = " << len[i]<<endl;
	}
	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;
	cout << "Number of G variance components = " << gen_Nbin << endl;
	cout << "Number of GxE variance components = " << nongen_Nbin << endl;
	outfile << "Number of G variance components = " << gen_Nbin << endl;
	outfile << "Number of GxE variance components = " << nongen_Nbin << endl;
	if (hetero_noise == true) {
		cout << "Number of NxE variance components = " << Nenv << endl;
		outfile << "Number of NxE variance components = " << Nenv << endl;
	} else {
		cout << "Number of NxE variance components = 0" << endl;
		outfile << "Number of NxE variance components = 0" << endl;
	}
	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;

	if (use_dummy){
		cout << "!!! A dummy phenotype is used for this GENIE run. The heritability estimates are NOT meaningful (please use the trace summaries only) !!!" << endl;
	    cout << "**************************************************" << endl;
		outfile << "!!! A dummy phenotype is used for this GENIE run. The heritability estimates are NOT meaningful (please use the trace summaries only) !!!" << endl;
	    outfile << "**************************************************" << endl;
	}


	cout << "OUTPUT: " << endl << "Variance components: " << endl;
	outfile << "OUTPUT: " << endl << "Variance components: " << endl;
	if (verbose >= 3) { 
		cout << "gen_Nbin = " << gen_Nbin << "\tT_Nbin = " << T_Nbin << endl;
	}		
	for (int j = 0 ; j < T_Nbin ; j++){
		if(j < gen_Nbin){
			cout << "Sigma^2_g[" << j<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
			outfile << "Sigma^2_g[" << j<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
		}
		else if (j <(gen_Nbin + nongen_Nbin)){
			int k = j - gen_Nbin;
			cout << "Sigma^2_gxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
			outfile << "Sigma^2_gxe[" << k<<"] : " << point_est(j,0) << "  SE : " << point_se(j,0) << endl;
		}
		else if (j<(T_Nbin)){
			int k = j - gen_Nbin - nongen_Nbin;
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
		for(int k = gen_Nbin ; k<(2 * gen_Nbin);k++)
			temp_sig += point_est_adj_gxe(k,0);
		point_est_adj_gxe(T_Nbin + 1,0) = temp_sig;// total GxE
	}else{
			point_est_adj_gxe(T_Nbin + 1,0) = point_est_adj_gxe(gen_Nbin,0);

	}
	temp_sig = 0;
	for(int k = 0 ; k < gen_Nbin ; k++)
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
			for(int k = gen_Nbin ; k<(2 * gen_Nbin);k++)
				temp2_sig += jack_adj_gxe(k,i);
			jack_adj_gxe(T_Nbin + 1,i) = temp2_sig;
		}else{
			jack_adj_gxe(T_Nbin + 1,i) = jack_adj_gxe(gen_Nbin,i);
		}
		temp2_sig = 0;
		for(int k = 0 ; k < gen_Nbin ; k++)
			temp2_sig += jack_adj_gxe(k,i);
		jack_adj_gxe(T_Nbin + 2,i) = temp2_sig;
	}

	MatrixXdr SEjack_adj_gxe = statsfn::jack_se(jack_adj_gxe);
	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;
	cout << "Heritabilities: " << endl;
	outfile << "Heritabilities: " << endl;
	for (int j = 0 ; j < T_Nbin ; j++){
		if(j < gen_Nbin){
			cout << "h2_g[" << j<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
			outfile << "h2_g[" << j<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
		}
		else if (j<(gen_Nbin + nongen_Nbin)){
			int k = j - gen_Nbin;
			cout << "h2_gxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
			outfile << "h2_gxe[" << k<<"] : " << point_est_adj_gxe(j,0) << " SE : " << SEjack_adj_gxe(j,0) << endl;
		}
		else if (j<(T_Nbin)){
			int k = j - gen_Nbin - nongen_Nbin;
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
	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;

	MatrixXdr enrich_g;
	MatrixXdr jack_enrich_g;
	jack_enrich_g.resize(gen_Nbin,Njack);
	enrich_g.resize(gen_Nbin,1);
	double total_g_h2 = 0;
	int total_snp = 0;
	for(int i = 0 ; i < gen_Nbin ; i++){
		total_g_h2 += point_est_adj_gxe(i,0);
		total_snp += len[i];
	}
	double numi;
	double denom = (double)total_g_h2 / total_snp;
	for(int i = 0 ; i < gen_Nbin ; i++){
		numi = point_est_adj_gxe(i,0) / len[i];
		enrich_g(i,0) = numi / denom;
	}
	for(int j = 0 ; j < Njack ; j++){
		total_g_h2 = 0;
		total_snp = 0;
		for (int k = 0 ; k < gen_Nbin ; k++){
			total_snp += len[k]-jack_bin[j][k];
			total_g_h2 += jack_adj_gxe(k,j);
		}
		denom = (double)total_g_h2 / total_snp;
		for(int k = 0 ; k < gen_Nbin ; k++){
			numi = (double)jack_adj_gxe(k,j) / (len[k]-jack_bin[j][k]);
			jack_enrich_g(k,j) = (double)numi / denom;
			}
	}
	MatrixXdr enrich_g_se;
	enrich_g_se = MatrixXdr::Zero(gen_Nbin,1);
	enrich_g_se = statsfn::jack_se(jack_enrich_g);


	cout << "Enrichments:" << endl;
	outfile << "Enrichments:" << endl;
	cout << "G enrichment" << endl;
	outfile << "G enrichment" << endl;
	for(int i = 0 ; i < gen_Nbin ; i++){
		cout << "Enrichment g[" << i<<"] : " << enrich_g(i,0) << " SE : " << enrich_g_se(i,0) << endl;
		outfile << "Enrichment g[" << i<<"] : " << enrich_g(i,0) << " SE : " << enrich_g_se(i,0) << endl;
	}
	// compute enrich GxE
	if(Annot_x_E == true){
		if (gen_by_env == true) {
			MatrixXdr enrich_gxe;
			MatrixXdr jack_enrich_gxe;
			jack_enrich_gxe.resize(gen_Nbin,Njack);
			enrich_gxe.resize(gen_Nbin,1);
			double total_gxe_h2 = 0;
			int total_snp = 0;
			for(int i = 0 ; i < gen_Nbin ; i++){
				total_gxe_h2 += point_est_adj_gxe(gen_Nbin + i,0);
				total_snp += len[gen_Nbin + i];
			}
			double numi;
			double denom = (double)total_gxe_h2 / total_snp;
			for(int i = 0 ; i < gen_Nbin ; i++){
				numi = point_est_adj_gxe(gen_Nbin + i,0) / len[gen_Nbin + i];
				enrich_gxe(i,0) = numi / denom;
			}
			for(int j = 0 ; j < Njack ; j++){
				total_gxe_h2 = 0;
				total_snp = 0;
				for (int k = 0 ; k < gen_Nbin ; k++){
					total_snp += len[gen_Nbin + k]-jack_bin[j][gen_Nbin + k];
					total_gxe_h2 += jack_adj_gxe(gen_Nbin + k,j);
				}
				denom = (double)total_gxe_h2 / total_snp;
				for(int k = 0 ; k < gen_Nbin ; k++){
					numi = jack_adj_gxe(gen_Nbin + k,j) / (len[gen_Nbin + k]-jack_bin[j][gen_Nbin + k]);
					jack_enrich_gxe(k,j) = numi / denom;
				}
			}
			MatrixXdr enrich_gxe_se;
			enrich_gxe_se = MatrixXdr::Zero(gen_Nbin,1);
			enrich_gxe_se = statsfn::jack_se(jack_enrich_gxe);

			cout << "GxE enrichment" << endl;
			outfile << "GxE enrichment" << endl;
			for(int i = 0 ; i < gen_Nbin ; i++){
				cout << "Enrichment gxe[" << i<<"] : " << enrich_gxe(i,0) << " SE : " << enrich_gxe_se(i,0) << endl;
				outfile << "Enrichment gxe[" << i<<"] : " << enrich_gxe(i,0) << " SE : " << enrich_gxe_se(i,0) << endl;
			}
		}
	}
	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;

	///compute parameters for overlapping annotations based on s-ldsc definition :

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
			for(int j = 0 ; j < gen_Nbin ; j++){
				if(annot_bool[i][j]==1)
					her_per_snp(i,0)+=her_per_snp_inbin(j,0);
				if((annot_bool[i][gen_Nbin + j]==1) && (Annot_x_E == true))
					her_per_snp(i,1)+=her_per_snp_inbin(gen_Nbin + j,0);
			}        
			if(k == Njack){
				for(int j = 0 ; j < gen_Nbin ; j++){
					if(annot_bool[i][j]==1)
						point_her_cat_ldsc(j,0)+=her_per_snp(i,0);
					if((annot_bool[i][gen_Nbin + j]==1) && (Annot_x_E == true))
						point_her_cat_ldsc(gen_Nbin + j,0)+=her_per_snp(i,1);
				}                
			}else{
				int temp = i / step_size;
				if(temp>=Njack)
					temp = Njack - 1;
				for(int j = 0 ; j < gen_Nbin ; j++){
					if(annot_bool[i][j]==1 && temp != k)
						her_cat_ldsc(j,k)+=her_per_snp(i,0);
					if((annot_bool[i][gen_Nbin + j]==1) && (temp != k) && (Annot_x_E == true))
						her_cat_ldsc(gen_Nbin + j,k)+=her_per_snp(i,1);
				}                   
			}     
		}
	}

	MatrixXdr se_her_cat_ldsc = statsfn::jack_se(her_cat_ldsc);

	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;
	cout << "Heritabilities and enrichments computed (overlapping def)" << endl;
	outfile << "Heritabilities and enrichments computed (overlapping def)" << endl;

	cout << "Heritabilities: " << endl;
	outfile << "Heritabilities: " << endl;
	for (int j = 0 ; j < gen_Nbin ; j++){
		cout << "h2_g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
		outfile << "h2_g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
	}
	if ((gen_by_env == true) && (Annot_x_E == true)) {
		for (int j = 0 ; j < gen_Nbin ; j++){
			int k = j + gen_Nbin;
			cout << "h2_gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;
			outfile << "h2_gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;

		}
	}

	int all_snp = 0;
	for(int i = 0 ; i < gen_Nbin ; i++){
		all_snp += len[i];
	}
	double snp_por;
	for (int i = 0 ; i < gen_Nbin ; i++){
		point_her_cat_ldsc(i,0) = (double)point_her_cat_ldsc(i,0) / point_est_adj_gxe(T_Nbin + 2,0);  //additive
		if((gen_by_env == true) && (Annot_x_E == true))
			point_her_cat_ldsc(i + gen_Nbin,0) = (double)point_her_cat_ldsc(i + gen_Nbin,0) / point_est_adj_gxe(T_Nbin + 1,0);  //GxE

		snp_por = (double)len[i]/all_snp;
		point_her_cat_ldsc(i,0) = (double)point_her_cat_ldsc(i,0) / snp_por;
		if((gen_by_env == true) && (Annot_x_E == true))
			point_her_cat_ldsc(i + gen_Nbin,0) = (double)point_her_cat_ldsc(i + gen_Nbin,0) / snp_por;
	}

	double temp_size;
	for(int i = 0 ; i < Njack ; i++){
		temp_size = all_snp;
		for(int k = 0 ; k < gen_Nbin ; k++)
			temp_size = temp_size - jack_bin[i][k];
		for(int j = 0 ; j < gen_Nbin ; j++){
			her_cat_ldsc(j,i) = (double)her_cat_ldsc(j,i) / jack_adj_gxe(T_Nbin + 2,i);
			if((gen_by_env == true) && (Annot_x_E == true))
				her_cat_ldsc(j + gen_Nbin,i) = (double)her_cat_ldsc(j + gen_Nbin,i) / jack_adj_gxe(T_Nbin + 1,i);

			snp_por = (double)(len[j]-jack_bin[i][j]) / temp_size;
			her_cat_ldsc(j,i) = (double)her_cat_ldsc(j,i) / snp_por;
			if((gen_by_env == true) && (Annot_x_E == true))
				her_cat_ldsc(j + gen_Nbin,i) = (double)her_cat_ldsc(j + gen_Nbin,i) / snp_por;

		}
	}

	se_her_cat_ldsc = statsfn::jack_se(her_cat_ldsc);

	cout << "**************************************************" << endl;
	outfile << "**************************************************" << endl;
	cout << "Enrichments (overlapping def): " << endl;
	outfile << "Enrichments (overlapping def): " << endl;
	for (int j = 0 ; j < gen_Nbin ; j++){
		cout << "Enrichment g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
		outfile << "Enrichment g[" << j<<"] : " << point_her_cat_ldsc(j,0) << " SE : " << se_her_cat_ldsc(j,0) << endl;
	}
	if ((gen_by_env == true) && (Annot_x_E == true)) {
		for (int j = 0 ; j < gen_Nbin ; j++){
			int k = j + gen_Nbin;
			cout << "Enrichment gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE : " << se_her_cat_ldsc(k,0) << endl;
			outfile << "Enrichment gxe[" << j<<"] : " << point_her_cat_ldsc(k,0) << " SE: " << se_her_cat_ldsc(k,0) << endl;
		}
	}
}

void print_trace () {
	string prefix = command_line_opts.OUTPUT_FILE_PATH;
	string trpath=prefix + ".tr";
	string mnpath=prefix + ".MN";
	trace_file.open(trpath.c_str(), std::ios_base::out);
	stringstream ss;
	for (int i=0; i<Nbin; i++){
		ss << "LD_SUM_" << i << ",";
	}
	ss << "NSNPS_JACKKNIFE";
	trace_file << ss.str() << endl;
	meta_file.open(mnpath.c_str(), std::ios_base::out);
	meta_file << "NSAMPLE,NSNPS,NBLKS,NBINS,K" << endl << Nindv << "," << Nsnp << "," << Njack << "," << Nbin << "," << Nz;
	meta_file.close();
	trace_file.close();
}


void print_input_parameters() {
	string outpath = command_line_opts.OUTPUT_FILE_PATH;
	outfile.open(outpath.c_str(), std::ios_base::out);
	if (!outfile.is_open()){
		cerr << "Error writing to output file : "<< outpath <<endl;
		exit(1);
	}
	outfile << command_line_opts.version_string << endl;
    outfile << "Seed = " << seed << endl;

}

void init_params () {

	verbose = command_line_opts.verbose;
	debug = command_line_opts.debugmode ||  verbose >= 4 ;

	trace = command_line_opts.print_trace;
	srand((unsigned int) time(0));
	Nz = command_line_opts.num_of_vec;
	k = Nz;

	jack_scheme = command_line_opts.jack_scheme;
	Njack = command_line_opts.jack_number;
	jack_size = command_line_opts.jack_size;
	jack_size *= 1e6; 

	memeff = command_line_opts.memeff;
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

template <typename Func>
void dummy_pheno(int Nind, Func& func){
	// fill in dummy phenotype
	cout << "Filling in dummy" << endl;
	mask.resize(Nindv, 1);
	pheno.resize(Nind, 1);
	for (int i=0; i < Nind; i++){
			pheno(i,0) = func();
			mask(i,0) = 1;

	}
}

int main(int argc, char const *argv[]){
    struct timeval now;
    gettimeofday(&now, NULL);
    long starttime = now.tv_sec * UMILLION + now.tv_usec;

	parse_args (argc,argv);
	init_params ();

	read_auxillary_files ();

	if (gen_by_env == true) {    
		///mask out indv with missingness from Enviro
		Enviro = Enviro.array().colwise() * mask.col(0).array();
	}
	//define random vector z's
	all_zb = MatrixXdr::Random (Nindv, Nz);

	if (seed == -1) {
		std::random_device rd;
        seed =  rd ();
		seedr.seed(seed);
	} else {
		seedr.seed(seed);
	}
    cout << "Seed = " << seed << endl;

	std::normal_distribution<> dist(0,1);
	auto z_vec = std::bind(dist, seedr);

	for (int i = 0 ; i < Nz ; i++)
		for(int j = 0 ; j < Nindv ; j++)
			all_zb(j,i) = z_vec();
	
	// fill in dummy (for trace summaries)
	if (use_dummy)
		dummy_pheno(Nindv, z_vec);

	cout << "Regressing covariates" << endl;
	regress_covariates ();	

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

	if (memeff) { 
		XXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);

		if(both_side_cov == true){
			UXXz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);
			XXUz = MatrixXdr::Zero(Nindv, T_Nbin * Nz * 2);
		}
		yXXy = MatrixXdr::Zero(T_Nbin, 2);

	} else {
		XXz=MatrixXdr::Zero(Nindv, T_Nbin * (Njack+1) * Nz);

		if(both_side_cov==true){
			UXXz=MatrixXdr::Zero(Nindv, T_Nbin * (Njack+1) * Nz);
			XXUz=MatrixXdr::Zero(Nindv, T_Nbin * (Njack+1) * Nz);
		}
		yXXy=MatrixXdr::Zero(T_Nbin, Njack + 1);
	}

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
	read_header = true;
	global_snp_index=-1;

	ifstream ifs (name.c_str(), ios::in|ios::binary);
	if (!ifs.is_open()){
		cerr << "Error reading file "<< name  <<endl;
		exit(1);
	}

	print_input_parameters ();

	cout << endl;
	cout << "Reading genotypes ..." << endl;

	
	//CHANGE(03/05): add trace summary files. input to -o is now just the prefix (all output file endings are fixed to .log)
	if (trace){
		print_trace ();
	}

	// This is where most of the computation happens
    // memeff = 1: Memory-efficient two-pass version
	if (memeff) {

        // In memory-efficient version, pass 2 is common
        // Setting opt2 = 1: makes Pass 1 even more efficient 
        //
		if (opt2){
			genotype_stream_pass_mem_efficient (name); // flexible read block size
		} else {
			genotype_stream_pass (name, 1);
		}

		genotype_stream_pass (name, 2);
	} else {
		genotype_stream_single_pass (name);
	}

	print_results ();

    gettimeofday(&now, NULL);
    long endtime = now.tv_sec * UMILLION + now.tv_usec;
    double elapsed = endtime - starttime;
    elapsed /= 1.e6;
    cout << "GENIE ran successfully. Time elapsed = " << elapsed << " seconds " << endl;
	
    outfile.close();
	return 0;
}
