#include "genotype.h"
//#include "arguments.h"
#include "storage.h"
#include "functions.h"
#include "vectorfn.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void genotype::init_means(bool is_missing){
	columnmeans.resize (Nsnp);
	for(int i=0;i< Nsnp;i++){
		double sum = columnsum[i];
		if(is_missing)
			columnmeans[i] = sum*1.0/(Nindv-not_O_i[i].size());
		else
			columnmeans[i] = sum*1.0/Nindv;
	}
}

void genotype::read_txt_naive (std::string filename,bool allow_missing){
	
	ifstream ifs (filename.c_str(), ios::in);                                       
	
	std::string line;
	std::getline(ifs, line);
    std::istringstream iss(line);
    if (!(iss >> Nsnp >> Nindv)) { 
		cout<<"ERROR: Header with number of SNPs and individuals not present"<<endl; 
		exit(-1);
	}
	
	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}

	int i=0;

	vector <bool> m;
	vector <bool> l;
	while(std::getline(ifs,line)){
		int sum=0;
		for(int j=0;j<line.size();j++){
			int val = int(line[j]-'0');	
			if(val==0){
				l.push_back(false);
				m.push_back(false);
			}
			else if(val==1){
				sum+=1;
				l.push_back(true);
				m.push_back(false);
			}
			else if(val==2){
				sum+=2;
				l.push_back(false);
				m.push_back(true);
			}
			else if(val==9 && allow_missing){
				not_O_i[i].push_back(j);
				not_O_j[j].push_back(i);
				l.push_back(false);
				m.push_back(false);
			}
			else{
				cout<<"Invalid entry in Genotype Matrix"<<endl;
				cout<<"If there is Missing data, run with -miss flag"<<endl;
				exit(-1);				
			}
		}
		i++;
		columnsum.push_back(sum);
		msb.push_back(m);
		lsb.push_back(l);
		m.clear();
		l.clear();
	}
	init_means(allow_missing);
}

void genotype::read_txt_mailman (std::string filename,bool allow_missing){
   	ifstream ifs (filename.c_str(), ios::in);                                       
	
	// Calculating the sizes and other stuff for genotype matrix
	std::string line;
	std::getline(ifs, line);
    std::istringstream iss(line);
    if (!(iss >> Nsnp >> Nindv)) { 
		cout<<"ERROR: Header with number of SNPs and individuals not present"<<endl; 
		exit(-1);
	}
	segment_size_hori = ceil(log(Nindv)/log(3));
	segment_size_ver = ceil(log(Nsnp)/log(3));
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	Nsegments_ver = ceil(Nindv*1.0/(segment_size_ver*1.0));
	p.resize(Nsegments_hori,std::vector<int>(Nindv));

	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}

	int i=0;

	while(std::getline(ifs,line)){
		int horiz_seg_no = i/segment_size_hori ;
		int sum=0;
		for(int j=0;j<line.size();j++){
			int val = int(line[j]-'0');
			if(val==0 || val==1 || val==2){
				sum+=val;
				p[horiz_seg_no][j] = (3 * p[horiz_seg_no][j]) + val ;
			}	
			else if(val==9 && allow_missing){
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j] ;
				not_O_i[i].push_back(j);
				not_O_j[j].push_back(i);
			}
			else{
				cout<<"ERROR: Invalid character in genotype file"<<endl;
				cout<<"If there is Missing data, run with -miss flag"<<endl;				
				exit(-1);
			}			
		}
		i++;
		columnsum.push_back(sum);
	}
	init_means(allow_missing);	
}


template<typename T>
static std::istream & binary_read(std::istream& stream, T& value){
	return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template <class T>
inline void printvector(vector<T> &t, string delim = " ", bool newline = false){
		for (int i = 0; i < t.size(); i++)
				cout << t[i] << delim;
        if (newline)
            cout << endl;
}

template <class T>
inline void printvectornl(vector<T> &t, string delim = " "){
    printvector (t, delim, true);
}


void genotype::read_bim (string filename){
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
	Nsnp = j;
	inp.close();
}

void genotype::read_fam (string filename){
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
	Nindv = j;
	inp.close();
}

// Accessor Functions

double genotype::get_geno(int snpindex,int indvindex,bool var_normalize=false){
	double m = msb[snpindex][indvindex];
	double l = lsb[snpindex][indvindex];
	double geno = (m*2.0+l) - get_col_mean(snpindex);
	if(var_normalize)
		return geno/get_col_std(snpindex);
	else
		return geno;
}

double genotype::get_col_mean(int snpindex){
	double temp = columnmeans[snpindex];
	return temp;
}

double genotype::get_col_sum(int snpindex){
	double temp = columnsum[snpindex];
	return temp;
}


double genotype::get_col_sum2(int snpindex){
	double temp=columnsum2[snpindex]; 
	return temp; 
}
double genotype::get_col_std(int snpindex){
	double p_i = get_col_mean(snpindex);
	if(p_i == 0 || p_i == 2)
		return 1.0;
	double temp = sqrt(p_i*(1-(0.5*p_i))) ; 
	return temp;
}

void genotype::generate_eigen_geno(MatrixXdr &geno_matrix,bool var_normalize){
	for(int i=0;i<Nsnp;i++){
		for(int j=0;j<Nindv;j++){
			double m = msb[i][j];
			double l = lsb[i][j];
			double geno = (m*2.0+l) - get_col_mean(i);
			if(var_normalize)
				geno_matrix(i,j) = geno/get_col_std(i);
			else
				geno_matrix(i,j) = geno;
		}
	}
}

// Modifier Functions

void genotype::update_col_mean(int snpindex,double value){
	columnmeans[snpindex] = value;
}

//
// Read bim file to get number of SNPs 
// filename: name of .bim file
// Return number of SNPs
// This is needed for defining jackknife blocks based on the total number (jackknife scheme == 1 or 2)
int get_number_of_snps (string bimfilename) {
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
int read_bim (string filename) {
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
			
		if (toks.size () < 6) {
//			exitWithError ("Bad .bim file:");
			cerr << "Bad .bim file:" << endl;
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

// Compute MAF from genotypes at a given SNP. 
// Genotypes are encoded in bed format in line. 
float get_observed_pj(const unsigned char* line){
	int y[4];
	int observed_sum = 0;
	int observed_ct = 0;

	for (int k = 0 ; k < ncol ; k++) {
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

    if (observed_ct == 0 ) {
		cerr <<"Cannot compute frequency at monomorphic SNP" << endl;
		exit (1);
    }

	return observed_sum * 0.5 / observed_ct;
}



// Draw genotype given allele frequency
int simulate_geno_from_random (float p_j, std::mt19937 &seedr){
    std::uniform_real_distribution<> udis (0,1);
	float rval = udis (seedr);
	float dist_pj[3] = { (1 - p_j)*(1 - p_j), 2*p_j*(1 - p_j), p_j*p_j };
	if(rval < dist_pj[0] )
		return 0;
	else if( rval >= dist_pj[0] && rval < (dist_pj[0]+dist_pj[1]))
		return 1;
	else
		return 2;
}


// Read specified numer of SNPs from bed file.
// ifs: stream object associated with bed file.
// num_snp: number of SNPs to read
void read_bed2 (std::istream& ifs, bool allow_missing, int num_snp)  {
	char magic[3];
    set_metadata ();

	gtype =  new unsigned char[ncol];

	if(read_header)
		binary_read(ifs, magic);

	int sum = 0;

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	int y[4];

//	int bin_pointer;

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
					val = simulate_geno_from_random(p_j, seedr);
					val++;
					val = (val == 1) ? 0 : val;
				}
				val-- ;
				val =  (val < 0 ) ? 0 :val ;
				sum += val;

				// Loop over all bins (=annotations) this SNP belong to
				for(int bin_index = 0 ; bin_index < pointer_bins.size();bin_index++){
				    int bin_pointer = pointer_bins[bin_index];

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
			int bin_pointer = pointer_bins[bin_index];
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
					val = simulate_geno_from_random(p_j, seedr);
					val++;
					val = (val == 1) ? 0 : val;
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
