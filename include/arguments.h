#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <bits/stdc++.h>
#include "io.h"
#include "std.h"

using namespace std;
int io::debug = 0;

struct options{

	std::string GENOTYPE_FILE_PATH;
	std::string OUTPUT_FILE_PATH;
	std::string PHENOTYPE_FILE_PATH; 
	std::string COVARIATE_FILE_PATH; 
	std::string COVARIATE_NAME; 
	std::string ENV_FILE_PATH;
	 std::string Annot_PATH;
	

	std::string simul_par_PATH;
	  std::string maf_ld_PATH;


	//int batchNum; 
	int num_of_evec ;
	int jack_number;

	bool getaccuracy ;
	bool debugmode;
	int accelerated_em;
	int l;
	bool memory_efficient;
	bool fast_mode;
	bool exannot;
	bool missing;
	bool text_version;
	double beta;
	//CHANGE(10/20)
	bool hetero_noise;
	bool gen_by_env;
	std::string model;
	int seed;
	bool normalize_proj_pheno;
	bool cov_add_intercept;
	int nthreads;
	bool verbose;
	//CHANGE(11/12)
	// bool perm_E_in_GxE;
    //CHANGE(03/04)
    bool print_trace;
};

/*template<typename T, typename U>
struct is_same {
    static const bool value = false; 
};

template<typename T>
struct is_same<T,T> { 
   static const bool value = true; 
};

*/
extern options command_line_opts;


void exitWithError(const std::string &error) {
	std::cout << error;
	std::cout.flush();

	exit(1);
}


std::string usage ( ) {
	std::ostringstream ostr;

	ostr <<"Usage: GENIE/GENIE_mem " << "[input flags]" << endl;
	ostr <<"\t-g [genotype] -annot [annotation] -p [phenotype] -c [covariates] -o [output]" << endl;
	ostr  << "\t[-e [environment] -m [G|G+GxE|G+GxE+NxE] -k [# random vectors] -jn [# jackknife subsamples]	-t [# threads]]" << endl;
	ostr <<  "\t[-s [random seed] -eXannot -norm_proj_pheno [0|1] -cov_add_intercept [0|1] -v [0|1]]"<<endl;

	return ostr.str ();
}


class Convert{
public:
	template <typename T>
	static std::string T_to_string(T const &val){
		std::ostringstream ostr;
		ostr << val;
		return ostr.str();
	}
		
	template <typename T>
	static T string_to_T(std::string const &val){
		std::istringstream istr(val);
		T returnVal;
        int num;
        //CHANGE(03/04): allow 0, 1 boolean
		if(is_same<T,bool>::value){
			if (!(istr >> std::boolalpha >> returnVal)){
                std::istringstream istr(val);
                if (istr >> num){
                    if (num == 0 || num == 1){
                        returnVal = num;
                        return returnVal;
                    }
                }
                exitWithError("CFG: Not a valid bool received!\n");
            }
			return returnVal;
		}
		else{
			if (!(istr >> returnVal))
				exitWithError("CFG: Not a valid " + (std::string)typeid(T).name() + " received!\n");
			return returnVal;
		}
	}

	static std::string string_to_T(std::string const &val){
		return val;
	}
};


class ConfigFile{
private:
	std::map<std::string, std::string> contents;
	std::string fName;

	void removeComment(std::string &line) const{
		if (line.find('#') != line.npos)
			line.erase(line.find('#'));
	}

	bool onlyWhitespace(const std::string &line) const{
		return (line.find_first_not_of(' ') == line.npos);
	}
	bool validLine(const std::string &line) const{
		std::string temp = line;
		temp.erase(0, temp.find_first_not_of("\t "));
		if (temp[0] == '=')
			return false;

		for (size_t i = temp.find('=') + 1; i < temp.length(); i++)
			if (temp[i] != ' ')
				return true;

		return false;
	}

	void extractKey(std::string &key, size_t const &sepPos, const std::string &line) const{
		key = line.substr(0, sepPos);
		if (key.find('\t') != line.npos || key.find(' ') != line.npos)
			key.erase(key.find_first_of("\t "));
	}
	void extractValue(std::string &value, size_t const &sepPos, const std::string &line) const{
		value = line.substr(sepPos + 1);
		value.erase(0, value.find_first_not_of("\t "));
		value.erase(value.find_last_not_of("\t ") + 1);
	}

	void extractContents(const std::string &line){
		std::string temp = line;
		temp.erase(0, temp.find_first_not_of("\t "));
		size_t sepPos = temp.find('=');

		std::string key, value;
		extractKey(key, sepPos, temp);
		extractValue(value, sepPos, temp);

		if (!keyExists(key))
			contents.insert(std::pair<std::string, std::string>(key, value));
		else
			exitWithError("CFG: Can only have unique key names!\n");
	}

	void parseLine(const std::string &line, size_t const lineNo){
        //TODO: Allow for boolean flags without "=". might have to do some input filtering
		if (line.find('=') == line.npos)
			exitWithError("CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n");

		if (!validLine(line))
			exitWithError("CFG: Bad format for line: " + Convert::T_to_string(lineNo) + "\n");

		extractContents(line);
	}

	void ExtractKeys(){
		std::ifstream file;
		file.open(fName.c_str());
		if (!file)
			exitWithError("CFG: File " + fName + " couldn't be found!\n");

		std::string line;
		size_t lineNo = 0;
		while (std::getline(file, line)){
			lineNo++;
			std::string temp = line;

			if (temp.empty())
				continue;

			removeComment(temp);
			if (onlyWhitespace(temp))
				continue;

			parseLine(temp, lineNo);
		}

		file.close();
	}
public:
	ConfigFile(const std::string &fName){
		this->fName = fName;
		ExtractKeys();
	}

	ConfigFile () {}

	bool keyExists(const std::string &key) const{
		return contents.find(key) != contents.end();
	}

	template <typename ValueType>
	ValueType getValueOfKey(const std::string &key, ValueType const &defaultValue = ValueType()) const{
		if (!keyExists(key))
			return defaultValue;

		return Convert::string_to_T<ValueType>(contents.find(key)->second);
	}

	void insertKey (const std::string &key, const std::string &value) {
		contents[key] = value;		
	}

	void printVersion () {
		io::println ("##################################",0);
		io::println ("#                                #",0);
		io::println ("#          GENIE (v1.0.0)        #",0);
		io::println ("#                                #",0);
		io::println ("##################################",0);
		io::println ("",0);
	}

	void printParameters ()  {

	    io::println ("###########Parameters#############",0);                                            
    	map<string,string>::iterator i;                                                                  
    	for (i = contents.begin(); i!=contents.end(); i++){                                            
        	string s = "#" + i->first +"\t" + i->second ;                                                
			io::println (s,0); 
    	}   
		io::println ("#################################",0); 
	}
};

void parse_args(int argc, char const *argv[]){
	
	// Setting Default Values
	command_line_opts.num_of_evec=10;
	command_line_opts.debugmode=false;
	command_line_opts.OUTPUT_FILE_PATH = "";
	bool got_genotype_file=false;
	bool got_output_file = false;
	bool got_phenotype_file = false;
	bool got_annot_file = false;

	command_line_opts.memory_efficient=false;
	command_line_opts.fast_mode=true;
	command_line_opts.missing=false;
	command_line_opts.text_version = false;
	command_line_opts.exannot=false;
	//CHANGE(10/20)
	command_line_opts.hetero_noise=true;
	command_line_opts.gen_by_env=true;
	command_line_opts.jack_number = 100;
	command_line_opts.model = "G+GxE+NxE";
	command_line_opts.seed = -1;
	command_line_opts.normalize_proj_pheno = true;
	command_line_opts.cov_add_intercept = true;
	command_line_opts.nthreads = 1;
	command_line_opts.verbose = false;
	//CHANGE(11/12)
	// command_line_opts.perm_E_in_GxE=false;
    // CHANGE(03/04)
    command_line_opts.print_trace=false;

	if(argc<2){
		// cout<<"Correct Usage is "<<argv[0]<<" -p <parameter file>"<<endl;
		
		exitWithError (usage ());
//		cout<<"Correct Usage is "<<argv[0] << ":" << endl;
//		cout <<"\t-g [genotype] -annot [annotation] -p [phenotype] -c [covariates] -o [output]" << endl;
//		cout  << "\t[-e [environment] -m [G|G+GxE|G+GxE+NxE] -k [# random vectors] -jn [# jackknife subsamples]	-t [# threads]]" << endl;
//		cout <<  "\t[-s [random seed] -eXannot -norm_proj_pheno [0|1] -cov_add_intercept [0|1] -v [0|1]]"<<endl;
//		exit(-1);
	}
    // using a config file instead of cmd-line args. TODO: have all the current options as config version. remove deprecated/redundant options
	if(strcmp(argv[1],"--config")==0){
		std::string cfg_filename = std::string(argv[2]);
		ConfigFile cfg(cfg_filename);
		cfg.printVersion();
		got_genotype_file=cfg.keyExists("genotype");
	    command_line_opts.jack_number = cfg.getValueOfKey<int>("num_jack", 10);
		//command_line_opts.batchNum = cfg.getValueOfKey<int> ("batchNum",10); 
		command_line_opts.num_of_evec=cfg.getValueOfKey<int>("num_vec",2);
		command_line_opts.getaccuracy=cfg.getValueOfKey<bool>("accuracy",false);
		//command_line_opts.getaccuracy=cfg.keyExists("accuracy"); // config file parsing doesn't allow boolean flags... needs to be fixed
		//command_line_opts.debugmode=cfg.getValueOfKey<bool>("debug",false);
		command_line_opts.debugmode=cfg.getValueOfKey<bool>("debug",false);
		command_line_opts.l=cfg.getValueOfKey<int>("l",0);
		command_line_opts.OUTPUT_FILE_PATH = cfg.getValueOfKey<string>("output",string(""));
		command_line_opts.GENOTYPE_FILE_PATH = cfg.getValueOfKey<string>("genotype",string(""));
		command_line_opts.PHENOTYPE_FILE_PATH= cfg.getValueOfKey<string>("phenotype", string("")); 
		command_line_opts.COVARIATE_FILE_PATH= cfg.getValueOfKey<string>("covariate", string(""));
        command_line_opts.ENV_FILE_PATH = cfg.getValueOfKey<string>("environment", string(""));
        command_line_opts.Annot_PATH = cfg.getValueOfKey<string>("annotation", string(""));
        command_line_opts.nthreads = cfg.getValueOfKey<int>("nthreads", 1); 
		command_line_opts.COVARIATE_NAME=cfg.getValueOfKey<string>("covariateName", string(""));  
		command_line_opts.accelerated_em = cfg.getValueOfKey<int>("accelerated_em",0);
		command_line_opts.memory_efficient = cfg.getValueOfKey<bool>("memory_efficient",false);	
		command_line_opts.fast_mode = cfg.getValueOfKey<bool>("fast_mode",true);
		command_line_opts.missing = cfg.getValueOfKey<bool>("missing",false);	
		command_line_opts.text_version = cfg.getValueOfKey<bool>("text_version",false);						
        command_line_opts.print_trace = cfg.getValueOfKey<bool>("trace", false);
        command_line_opts.verbose = cfg.getValueOfKey<bool>("verbose", false);
        //command_line_opts.print_trace = cfg.keyExists("print_trace");
        //command_line_opts.verbose = cfg.keyExists("verbose");
        command_line_opts.model = cfg.getValueOfKey<string>("model", string(""));
        if (cfg.keyExists("model")){
            const char* model_arg = command_line_opts.model.c_str();
            if (strcmp(model_arg, "G")==0) {
                command_line_opts.hetero_noise=false;
                command_line_opts.gen_by_env=false;
                command_line_opts.normalize_proj_pheno = false;
                cout << "Estimation of G heritability" << endl;
            } else if (strcmp(model_arg, "G+GxE")==0) {
                command_line_opts.hetero_noise=false;
                command_line_opts.gen_by_env=true;
                cout << "Estimation of G and GxE heritability (no heterogeneous noise)" << endl;
            } else if (strcmp(model_arg, "G+GxE+NxE")==0) {
                command_line_opts.hetero_noise=true;
                command_line_opts.gen_by_env=true;
                cout << "Estimation of G and GxE heritability (with heterogeneous noise)" << endl;
            } else {
                cout << "model can only be G, G+GxE, or G+GxE+NxE. Using default G+GxE+NxE model" << endl;
                command_line_opts.hetero_noise=true;
                command_line_opts.gen_by_env=true;
                cout << "Estimation of G and GxE heritability (with heterogeneous noise)" << endl;
            }
	    }
		cfg.printParameters ();
    }
        //currently, trace is only available for "G" model. (TODO for SUMRHE to have other options)
	else {
		//Add other default parameters here
		ConfigFile cfg;
		cfg.printVersion();
		cfg.insertKey ("num_vec", "10");
		cfg.insertKey ("num_block", "100");
		cfg.insertKey ("debug", "0");
		cfg.insertKey ("seed",  Convert::T_to_string (command_line_opts.seed));

		for (int i = 1; i < argc; i++) { 
			if (i + 1 != argc){
				if(strcmp(argv[i],"-g")==0){
					command_line_opts.GENOTYPE_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("genotype", argv[i+1]);
					got_genotype_file=true;
					i++;
				}
				else if(strcmp(argv[i],"-o")==0){
					command_line_opts.OUTPUT_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("output", argv[i+1]);
					got_output_file=true;
					i++;
				}
				else if(strcmp(argv[i],"-annot")==0){
					command_line_opts.Annot_PATH= string(argv[i+1]);
					cfg.insertKey ("annotation", argv[i+1]);
            	    got_annot_file=true;
					i++;
				}
				else if(strcmp(argv[i],"-simul_par")==0){
					command_line_opts.simul_par_PATH= string(argv[i+1]);
					cfg.insertKey ("simul_par", argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-maf_ld")==0){
					command_line_opts.maf_ld_PATH= string(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i], "-p")==0){
					command_line_opts.PHENOTYPE_FILE_PATH =string(argv[i+1]); 
					cfg.insertKey ("phenotype", argv[i+1]);
                	got_phenotype_file=true;
					i++; 
				}
				else if (strcmp(argv[i], "-e")==0){
					command_line_opts.ENV_FILE_PATH=string(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-c")==0){
					command_line_opts.COVARIATE_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("covariate", argv[i+1]);
					i++; 
				}
				else if(strcmp(argv[i],"-cn")==0){
					command_line_opts.COVARIATE_NAME = string(argv[i+1]);
					cfg.insertKey ("covariate name", argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-k")==0){
					command_line_opts.num_of_evec = atoi(argv[i+1]);
					cfg.insertKey ("num_vec", argv[i+1]);
					i++;
				}
                //CHANGE(03/04): batchNum is deprecated
                /*
				else if(strcmp(argv[i],"-b")==0){
					command_line_opts.batchNum=atoi(argv[i+1]); 
					i++; 
				}
                */
				else if(strcmp(argv[i],"-l")==0){
					command_line_opts.l = atoi(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-jn")==0){
					command_line_opts.jack_number= atoi(argv[i+1]);
					cfg.insertKey ("num_jackknife_blocks", argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-h2")==0){
					command_line_opts.beta= atof(argv[i+1]);
					i++;
				}
				//CHANGE(2/27)
				else if(strcmp(argv[i],"-m")==0) {
					if (strcmp(argv[i+1], "G")==0) {
						command_line_opts.model = "G";
						command_line_opts.hetero_noise=false;
						command_line_opts.gen_by_env=false;
						command_line_opts.normalize_proj_pheno = false;
						cout << "Estimation of G heritability" << endl;
				} else if (strcmp(argv[i+1], "G+GxE")==0) {
						command_line_opts.model = "G+GxE";
						command_line_opts.hetero_noise=false;
						command_line_opts.gen_by_env=true;
						cout << "Estimation of G and GxE heritability (no heterogeneous noise)" << endl;
				} else if (strcmp(argv[i+1], "G+GxE+NxE")==0) {
						command_line_opts.hetero_noise=true;
						command_line_opts.gen_by_env=true;
						cout << "Estimation of G and GxE heritability (with heterogeneous noise)" << endl;
				} else {
						cout << "model can only be G, G+GxE, or G+GxE+NxE. Using default G+GxE+NxE model" << endl;
						command_line_opts.hetero_noise=true;
						command_line_opts.gen_by_env=true;
						cout << "Estimation of G and GxE heritability (with heterogeneous noise)" << endl;
					}
					i++;
				} else if (strcmp(argv[i], "-s")==0) {
					command_line_opts.seed = atoi(argv[i+1]);
					cfg.insertKey ("seed", argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-norm_proj_pheno")==0) {
					int flag = atoi(argv[i+1]);
					if (flag == 0) {
						command_line_opts.normalize_proj_pheno = false;
						cout << "The phenotype vector after projection on covariates will be be normalized" << endl;
					} else {
						command_line_opts.normalize_proj_pheno = true;
					}
					i++;
				} else if (strcmp(argv[i], "-cov_add_intercept")==0) {
					int flag = atoi(argv[i+1]);
					if (flag == 0) {
						command_line_opts.cov_add_intercept = false;
						cout << "The one vector (intercept term) will not be added to the covariates" << endl;
					} else {
						command_line_opts.cov_add_intercept = true;
					}
					i++;
				} else if (strcmp(argv[i], "-t")==0) {
					command_line_opts.nthreads = atoi(argv[i+1]);
					i++;
                   // CHANGE(03/04): makes more sense to have flags as a boolean true. TODO: also, why aren't we using getopt_long?
				} else if ((strcmp(argv[i], "-v")==0) || (strcmp(argv[i], "--verbose")==0)) {
					command_line_opts.verbose = true;
				} else if ((strcmp(argv[i], "-tr")==0) || (strcmp(argv[i], "--trace")==0)){
                    command_line_opts.print_trace=true;
                }
							//CHANGE(11/12)
							// else if (strcmp(argv[i], "-perm_E_in_GxE")==0) {
							// 	command_line_opts.perm_E_in_GxE=true;
							// 	cout << "permute E in GxE flag set" << endl;
							// }
							// else if(strcmp(argv[i],"-v")==0)
							// 	command_line_opts.debugmode=true;
				else if(strcmp(argv[i],"-mem")==0)
					command_line_opts.memory_efficient=true;
				else if(strcmp(argv[i],"-miss")==0)
					command_line_opts.missing=true;
				else if(strcmp(argv[i],"-nfm")==0)
					command_line_opts.fast_mode=false;
				else if(strcmp(argv[i],"-eXannot")==0)		
					command_line_opts.exannot=true;	
				else{
					cout<<"Not Enough or Invalid arguments"<<endl;
									// cout<<"Correct Usage is "<<argv[0]<<" -g <genotype file> -k <num_of_vec> -b <num_of_zb/10>  -v (for debugmode) -a (for getting accuracy)"<<endl;
					exitWithError (usage());
/*
					cout<<"Correct Usage is "<<argv[0] << ":" << endl;
					cout <<"\t-g [genotype] -annot [annotation] -p [phenotype] -c [covariates] -o [output]" << endl;
					cout  << "\t[-e [environment] -m [G|G+GxE|G+GxE+NxE] -k [# random vectors] -jn [# jackknife subsamples]	-t [# threads]]" << endl;
					cout <<  "\t[-s [random seed] -eXannot -norm_proj_pheno [0|1] -cov_add_intercept [0|1] -v [0|1]]"<<endl;
					exit(-1);
*/
				}
			}
			else if(strcmp(argv[i],"-v")==0) {
				command_line_opts.debugmode=true;
				cfg.insertKey ("debug", "1");
			}
			else if(strcmp(argv[i],"-a")==0)
				command_line_opts.getaccuracy=true;
					//	else if(strcmp(argv[i],"-vn")==0)
					//			command_line_opts.var_normalize=true;
			else if(strcmp(argv[i],"-mem")==0)
				command_line_opts.memory_efficient=true;
			else if(strcmp(argv[i],"-nfm")==0)
				command_line_opts.fast_mode=false;
			else if(strcmp(argv[i],"-miss")==0)
				command_line_opts.missing=true;
			else if(strcmp(argv[i],"-txt")==0)
				command_line_opts.text_version=true;
		}
		cfg.printParameters ();
	}
	
	if(got_genotype_file==false){
		cout<<"Genotype file missing"<<endl;
		exitWithError (usage ());
//		cout<<"Correct Usage is "<<argv[0]<<" -g <genotype file> -k <num_of_vec> -m <max_iterations> -v (for debugmode) -a (for getting accuracy)"<<endl;
		exit(-1);
	}

}

#endif
