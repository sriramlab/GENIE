#include "auxillary.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

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
	int xNenv = 0;

	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
			missing.push_back(vector<int>()); //push an empty row  
			xNenv++;
		}
	}
	vector<double> cov_sum(xNenv, 0);
	if (gen_by_env == true) {
		Enviro = MatrixXdr::Zero(Nind, xNenv);
		cout<< "Reading in "<< xNenv << " environmental variables ..." << endl;
	}

	int j = 0;
	while(std::getline(ifs, line)){
		in.clear();
		in.str(line);
		string temp;
		in>>temp; in>>temp; //FID IID 
		for(int k = 0; k < xNenv ; k++){

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
	return xNenv;
}

// Read covariate file. 
// Read environmental variables (and Enviro is initialized) before calling this function.
// Adds environmental variables to covariates if add_env_to_cov == True
// Adds intercept to covariates if cov_add_intercept == True
// Inputs: Number of individuals, filename
// Return total number of covariates (covariates specified in covariate file [ + environmental variables] [ + intercept])
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
    int Nenv = 0;
	if(add_env_to_cov == true)
        Nenv = Enviro.cols ();

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
				covariate (i, covNum + Nenv) = 1;
			return covNum + 1;
		} 
	}
	if (cov_add_intercept == true) {
		//adding col of all ones to covariates
		for (int i = 0 ; i < Nind ; i++)
			covariate(i, covNum + Nenv) = 1;
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


// Read pheno file
// Inputs: Number of individuals, filename
void read_pheno(int Nind, std::string filename){
	ifstream ifs(filename.c_str(), ios::in); 

	if (!ifs.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}

	std::string line;
	std::istringstream in;  
	phenocount = 0; 
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
	new_pheno.resize(Nind, phenocount);
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

//
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

void read_annot (string filename) {
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
	vector<bool> snp_annot;

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
		cerr << "Number of rows in bim file and annotation file does not match\n" << endl;
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

	int Total_Nbin;
	if(Annot_x_E == false){
		jack_bin.resize(Njack, vector<int>(Nbin + Nenv + Nenv,0));
		read_bin.resize(Nreadblocks, vector<int>(Nbin + Nenv + Nenv,0));
		Total_Nbin = Nbin + Nenv;
	} else {
		jack_bin.resize(Njack, vector<int>(Nbin + (Nenv * Nbin) + Nenv,0));
		read_bin.resize(Nreadblocks, vector<int>(Nbin + (Nenv * Nbin) + Nenv,0));
		Total_Nbin = Nbin + (Nenv * Nbin);
	}

	for (int i = 0 ; i < Nsnp ; i++)
		for(int j = 0 ; j < Total_Nbin ; j++)
			if (annot_bool[i][j]==1){
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
	step_size_rem = Nsnp % Njack;
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
        cerr << "Number of rows in bim file and annotation file does not match" << endl;
		exit (1);
	}
	cout << "Total number of SNPs : " << Nsnp << endl;
	for (int i = 0 ; i < Nbin ; i++){
		cout << len[i]<<" SNPs in " << i<<"-th bin" << endl;
		Nsnp_annot += len[i];
	}
}
