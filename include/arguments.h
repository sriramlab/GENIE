#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
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
	std::string XSUM_FILE_PATH;
	
	int num_of_vec ;

	bool debugmode;
	bool fast_mode;
	bool exannot;
	bool missing;

	bool hetero_noise;
	bool gen_by_env;
	std::string model;
	int seed;
	bool normalize_proj_pheno;
	bool cov_add_intercept;
	int nthreads;
	int verbose;

	// If memeff = true, uses the memory efficient version
	bool memeff;

	// If opt1 = true, writes to a file after pass 1 
	// Reads from file in pass 2
	bool opt1;

	// If opt2 = true, reads groups of SNPs in block (termed read block)
	// Number of SNPs in a read block can be controlled based on memory constraints by setting mem_Nsnp
	bool opt2;

	int mem_Nsnp;
	int mem_alloc;

	bool use_ysum;
	bool keep_xsum;
	

	// How to do jackknife
	// 0: no jackknife
	// 1: constant number of SNPs
	// 2: same as 0 (adjusting for chromosome boundaries)
	// 3: constant physical length (adjusting for chromosome boundaries)
	// Default: 1
	int jack_scheme;
	int jack_number;

	// Length of jackknife block in MB
	int jack_size;

    bool print_trace;
	bool use_dummy_pheno; // in case no phenotype is provided; use dummy pheno
	std::string TRACE_FILE_PATH;
	string version_string;
};

extern options command_line_opts;

void exitWithError(const std::string &error) {
	std::cerr << error;
	std::cerr.flush();

	exit(1);
}

void noTrace(options& opt){  // trace summaries are only supported for G only model
	if (opt.print_trace == true){
		cerr << "Trace summary is only supported for G only model! Unsetting trace summary option..." << endl;
		opt.print_trace = false;
	}	
}

// CHANGE (07/16): IMO we should be using getopt_long, but for the moment we can use this shortcut for argument processing
int isarg(const char* arg, const char* shortarg, const char* longarg){
	return ((strcmp(shortarg, arg) == 0) || (strcmp(longarg, arg)==0));
}

std::string usage ( ) {
	std::ostringstream ostr;

	ostr <<"Usage: GENIE " << "[--config parameter file|command-line options]" << endl;
	ostr << "\nOptions:\n"
		 << "  -h, --help 					Show this message and exit\n"
		 << "  -g, --genotype				The path of PLINK BED genotype file\n"
		 << "  -p, --phenotype				The path of phenotype file\n"
		 << "  -c, --covariate				The path of covariate file\n"
		 << "  -a, --annot					The path of input annotation for partitioned heritability\n"
		 << "  -o, --output					Output file path (prefix)\n"
		 << "  -e, --environment				The path of environment file\n"
		 << "  -m, --model					Specification of the model. Currently there are 3 options:\n"
		 << "						1. additive genetic (G) effects only (arg: G)\n"
		 << "							The model reduces to RHE-mc (Pazokitoroudi et al. Nat Commun (2020). https://doi.org/10.1038/s41467-020-17576-9).\n"
		 << "						2. additive genetic (G) and gene-environment (GxE) effects (arg: G+GxE)\n"
		 << "							The model treats noise/environment effects as homogeneous.\n"
		 << "						3. additive genetic (G), gene-environment (GxE) and heterogeneous noise (NxE) effects (arg: G+GxE+NxE)\n"
		 << "							The model treats noise/environment effects as heterogeneous.\n"
		 << "  -k, --num-vec					The number of random vectors (10 is recommended).\n"
		 << "  -jn, --num-jack				The number of jackknife blocks (100 is recommended).\n"
		 << "  -t, --nthreads				The number of threads for multithreading\n"
		 << "  -tr, --trace					Flag for printing trace summary files (.tr) and metadata (.MN), compatible with SUM-RHE.\n"
		 << "						If this flag is set, a phenotype file is no longer required (a dummy phenotype without covariates will be used). Currently supports G-model only.\n"
		 << "  -v, --verbose					Verbose mode; Output extra information (Normal equation, number of samples, etc.).\n"
		 << "  -eXa, --eXannot				By default, GENIE fits a single GxE variance component. To partition the GxE component w.r.t the annotation file, add '-eXannot' flag.\n"
		 << "  -np, --norm-proj-pheno			By default, the phenotype vector is standardized after regressing covariates. Turn this off by setting '--norm-proj-pheno 0'.\n"
		 << "  -i, --cov-add-intercept			By default, a vector of ones is appended to the covariates (the intercept term). Turn this off by setting '--cov-add-intercept 0'.\n"
		 << endl;
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
				cout << val << endl;
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
	std::map<std::string, int> levels;
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
		int level = 0;

		contents[key] = value;
		levels[key] = 0;
	}

	void parseLine(const std::string &line, size_t const lineNo){
        //TODO: Allow for boolean flags without "=". might have to do some input filtering
		if (line.find('=') == line.npos)
			exitWithError("ERROR in parameter file: Could not find separator on line: " + Convert::T_to_string(lineNo) + "\n");

		if (!validLine(line))
			exitWithError("ERROR in parameter file: Bad format for line: " + Convert::T_to_string(lineNo) + "\n");

		extractContents(line);
	}

	void ExtractKeys(){
		std::ifstream file;
		file.open(fName.c_str());
		if (!file)
			exitWithError("ERROR: File " + fName + " not found!\n");

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
    string modelstring;
    string outputstring; 
	ConfigFile(const std::string &fName){
		this->fName = fName;
		ExtractKeys();

        modelstring = "";
        outputstring  = "";
	}

	ConfigFile () {
        modelstring = "";
        outputstring  = "";
	}

	void readFile (const std::string &fName) {
		this->fName = fName;
		ExtractKeys();
	}

	bool keyExists(const std::string &key) const{
		return contents.find(key) != contents.end();
	}

	template <typename ValueType>
	ValueType getValueOfKey(const std::string &key, ValueType const &defaultValue = ValueType()) const{
		if (!keyExists(key))
			return defaultValue;

		return Convert::string_to_T<ValueType>(contents.find(key)->second);
	}

	void insertKey (const std::string &key, const std::string &value, int level = 1 ) {
		contents[key] = value;		
		levels[key] = level; 
	}

	std::string versionstring () {
		std::ostringstream ostr;
		ostr << "############################################################" << endl;
		ostr << "#                                                          #" << endl;
		ostr << "#                       GENIE (v1.2.1)                     #" << endl;
		ostr << "#                                                          #" << endl;
		ostr << "############################################################" << endl;
		return ostr.str ();
	}

	std::string paramstring (int verbose = 0 )  {
		std::ostringstream ostr;
	    ostr << "####################### Parameters #########################" << endl; 

    	map<string,string>::iterator i;
		for (i = contents.begin(); i!=contents.end(); i++){
			string s = "#" + i->first +"\t" + i->second;
			int level = levels [i->first];
			if (level <= verbose)
				ostr << s << endl;
    	}
		ostr << "############################################################" << endl;
		ostr << modelstring << endl;
		ostr << outputstring << endl;
		return ostr.str();
	}
};


void parse_args(int argc, char const *argv[]){
	
	// Setting Default Values
	command_line_opts.fast_mode = true;
	command_line_opts.missing = false;

	bool got_genotype_file = false;
	bool got_output_file = false;
	bool got_phenotype_file = false;
	bool got_covariate_file = false;
	bool got_environment_file = false;
	bool got_annot_file = false;


	command_line_opts.num_of_vec = 10;
	command_line_opts.debugmode = false;
	command_line_opts.OUTPUT_FILE_PATH = "";

	command_line_opts.exannot = false;
	command_line_opts.hetero_noise = true;
	command_line_opts.gen_by_env = true;

	command_line_opts.jack_scheme = 1;
	command_line_opts.jack_number = 100;
	command_line_opts.jack_size = 10;

	command_line_opts.use_ysum = false;
	command_line_opts.keep_xsum = true;

	command_line_opts.model = "G+GxE+NxE";
	command_line_opts.seed = -1;
	command_line_opts.normalize_proj_pheno = true;
	command_line_opts.cov_add_intercept = true;
	command_line_opts.nthreads = 1;
	command_line_opts.verbose = 0;
    command_line_opts.print_trace = false;
	command_line_opts.use_dummy_pheno = false;

	command_line_opts.memeff = false;
	command_line_opts.opt1 = true;
	command_line_opts.opt2 = true;
	command_line_opts.mem_Nsnp = 10; 
	command_line_opts.mem_alloc = -1; 

	ConfigFile cfg;
	string param_string; 
	string version_string = cfg.versionstring();
	command_line_opts.version_string = version_string;

	if(argc<2 || (argc == 2 && (isarg(argv[1], "-h", "--help")))){
		cerr << version_string << endl;
		exitWithError (usage ());
	}
    // using a config file instead of cmd-line args. TODO: have all the current options as config version. remove deprecated/redundant options

	cfg.insertKey ("num_vec", "10");

	cfg.insertKey ("jack_scheme", "1");
	cfg.insertKey ("jack_size", "10");
	cfg.insertKey ("num_jack", "100");

	cfg.insertKey ("nthreads", "1");
	cfg.insertKey ("model", "G+GxE+NxE");
	cfg.insertKey ("norm-proj-pheno", "1");
	cfg.insertKey ("cov-add-intercept", "1");
	cfg.insertKey ("verbose", "0");
	cfg.insertKey ("memeff", "0");
	cfg.insertKey ("opt1", "1");
	cfg.insertKey ("opt2", "1");
	cfg.insertKey ("eXannot", "0");
	cfg.insertKey ("hetero_noise", "1");
	cfg.insertKey ("gene_by_env", "1");

	if (strcmp(argv[1],"--config")==0) {
		std::string cfg_filename = std::string(argv[2]);
		cfg.readFile (cfg_filename);
		
	    command_line_opts.jack_scheme = cfg.getValueOfKey<int>("jack_scheme", 1);
	    command_line_opts.jack_number = cfg.getValueOfKey<int>("num_jack", 100);
	    command_line_opts.jack_size = cfg.getValueOfKey<int>("jack_size", 10);

		command_line_opts.num_of_vec = cfg.getValueOfKey<int>("num_vec", 10);
		command_line_opts.debugmode = cfg.getValueOfKey<bool>("debug", false);

		got_output_file = cfg.keyExists("output");
		command_line_opts.OUTPUT_FILE_PATH = cfg.getValueOfKey<string>("output",string(""));

		got_genotype_file = cfg.keyExists("genotype");
		command_line_opts.GENOTYPE_FILE_PATH = cfg.getValueOfKey<string>("genotype",string(""));

		command_line_opts.PHENOTYPE_FILE_PATH= cfg.getValueOfKey<string>("phenotype", string(""));
		if (command_line_opts.PHENOTYPE_FILE_PATH != "")
			got_phenotype_file = true;

		command_line_opts.COVARIATE_FILE_PATH= cfg.getValueOfKey<string>("covariate", string(""));
		if (command_line_opts.COVARIATE_FILE_PATH != "")
			got_covariate_file = true;

        command_line_opts.ENV_FILE_PATH = cfg.getValueOfKey<string>("environment", string(""));
		if (command_line_opts.ENV_FILE_PATH != "")
			got_environment_file = true;

        command_line_opts.Annot_PATH = cfg.getValueOfKey<string>("annotation", string(""));
		if (command_line_opts.Annot_PATH != "")
			got_annot_file = true;

		command_line_opts.XSUM_FILE_PATH= cfg.getValueOfKey<string>("summary_genotype", string(""));
		command_line_opts.TRACE_FILE_PATH= cfg.getValueOfKey<string>("trace_input", string("")); 
		command_line_opts.COVARIATE_NAME=cfg.getValueOfKey<string>("covariateName", string(""));  
		command_line_opts.seed = cfg.getValueOfKey<int>("seed", -1);
		command_line_opts.fast_mode = cfg.getValueOfKey<bool>("fast_mode", true);
		command_line_opts.missing = cfg.getValueOfKey<bool>("missing", false);	
        command_line_opts.print_trace = cfg.getValueOfKey<bool>("trace", false);
        command_line_opts.verbose = cfg.getValueOfKey<int>("verbose", 0);
        command_line_opts.nthreads = cfg.getValueOfKey<int>("nthreads", 1); 

		command_line_opts.memeff  = cfg.getValueOfKey<bool>("memeff", false);
		command_line_opts.opt1 = cfg.getValueOfKey<bool>("opt1", true);
        command_line_opts.opt2 = cfg.getValueOfKey<bool>("opt2", true); 
        command_line_opts.mem_Nsnp = cfg.getValueOfKey<int>("mem_Nsnp", 10); 
        command_line_opts.mem_alloc = cfg.getValueOfKey<int>("mem_alloc", -1); 

        command_line_opts.use_ysum = cfg.getValueOfKey<bool>("use_phenosum", false);
        command_line_opts.keep_xsum = cfg.getValueOfKey<bool>("output_genosum", true);

        command_line_opts.normalize_proj_pheno = cfg.getValueOfKey<bool> ("norm_proj_pheno", true);
		command_line_opts.cov_add_intercept = cfg.getValueOfKey<bool>("cov_add_intercept", true);
		command_line_opts.exannot = cfg.getValueOfKey<bool>("eXannot", false);

        command_line_opts.model = cfg.getValueOfKey<string>("model", string("G+GxE+NxE"));
        if (cfg.keyExists("model")){
            const char* model_arg = command_line_opts.model.c_str();
            if (strcmp(model_arg, "G")==0) {
                command_line_opts.hetero_noise = false;
                command_line_opts.gen_by_env = false;
                command_line_opts.normalize_proj_pheno = false;
                cfg.modelstring = "# Estimating G heritability";	
            } else if (strcmp(model_arg, "G+GxE")==0) {
                command_line_opts.hetero_noise = false;
                command_line_opts.gen_by_env = true;
                cfg.modelstring = "# Estimating G and GxE heritability (no heterogeneous noise)";
				// noTrace(command_line_opts);
            } else if (strcmp(model_arg, "G+GxE+NxE")==0) {
                command_line_opts.hetero_noise = true;
                command_line_opts.gen_by_env = true;
                cfg.modelstring = "# Estimating G and GxE heritability (with heterogeneous noise)";
				// noTrace(command_line_opts);			
			
            } else {
                command_line_opts.hetero_noise = true;
                command_line_opts.gen_by_env = true;
                cfg.modelstring = "# Choice of models must be one of G, G+GxE, or G+GxE+NxE. Using default G+GxE+NxE model\n# Estimating G and GxE heritability (with heterogeneous noise)";
				// noTrace(command_line_opts);
            }
	    }

    //currently, trace is only available for "G" model. (TODO for SUMRHE to have other options)
	} else {

		for (int i = 1; i < argc; i++) {
			if (i + 1 != argc){
				if(isarg(argv[i], "-g", "--genotype")){
					command_line_opts.GENOTYPE_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("genotype", argv[i+1], 0);
					got_genotype_file = true;
					i++;
				} else if(isarg(argv[i], "-o", "--output")){
					command_line_opts.OUTPUT_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("output", argv[i+1], 0);
					got_output_file = true;
					i++;
				} else if(isarg(argv[i], "-a", "--annot")){
					command_line_opts.Annot_PATH= string(argv[i+1]);
					cfg.insertKey ("annotation", argv[i+1], 0);
            	    got_annot_file = true;
					i++;
				} else if(isarg(argv[i], "-p", "--phenotype")){
					command_line_opts.PHENOTYPE_FILE_PATH =string(argv[i+1]); 
					cfg.insertKey ("phenotype", argv[i+1], 0);
                	got_phenotype_file = true;
					i++; 
				} else if (isarg(argv[i], "-e", "--environment")){
					command_line_opts.ENV_FILE_PATH=string(argv[i+1]);
					cfg.insertKey ("environment", argv[i+1], 0);
                	got_environment_file = true;
					i++;
				} else if(isarg(argv[i], "-c", "--covariate")){
					command_line_opts.COVARIATE_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("covariate", argv[i+1], 0);
                	got_covariate_file = true;
					i++; 
				} else if(strcmp(argv[i],"-sg")==0){
					command_line_opts.XSUM_FILE_PATH = string(argv[i+1]);
					i++; 
				} else if(isarg(argv[i], "-cn", "--covariate-name")){
					command_line_opts.COVARIATE_NAME = string(argv[i+1]);
					cfg.insertKey ("covariate name", argv[i+1], 0);
					i++;
				} else if(isarg(argv[i], "-k", "--num-vec")){
					command_line_opts.num_of_vec = atoi(argv[i+1]);
					cfg.insertKey ("num_vec", argv[i+1], 0 );
					i++;
				} else if(isarg(argv[i], "-jn", "--num-jack")){
					command_line_opts.jack_number = atoi(argv[i+1]);
					cfg.insertKey ("num_jack", argv[i+1], 0);
					i++;
				} else if(isarg(argv[i], "-js", "--jack-scheme")){
					command_line_opts.jack_scheme = atoi(argv[i+1]);
					cfg.insertKey ("jack_scheme", argv[i+1], 0);
					i++;
				} else if(isarg(argv[i], "-m", "--model")) {
                    cfg.insertKey ("model", argv[i+1], 0);
					if (strcmp(argv[i+1], "G")==0) {
						command_line_opts.model = "G";
						command_line_opts.hetero_noise = false;
						command_line_opts.gen_by_env = false;
						command_line_opts.normalize_proj_pheno = false;
                        cfg.modelstring = "# Estimating G heritability";	
				    } else if (strcmp(argv[i+1], "G+GxE")==0) {
						command_line_opts.model = "G+GxE";
						command_line_opts.hetero_noise = false;
						command_line_opts.gen_by_env = true;
                        cfg.modelstring = "# Estimating G and GxE heritability (no heterogeneous noise)";
						// noTrace(command_line_opts);
				    } else if (strcmp(argv[i+1], "G+GxE+NxE")==0) {
						command_line_opts.hetero_noise = true;
						command_line_opts.gen_by_env = true;
                        cfg.modelstring = "# Estimating and GxE heritability (with heterogeneous noise)";
						// noTrace(command_line_opts);
				    } else {
						command_line_opts.hetero_noise = true;
						command_line_opts.gen_by_env = true;
                        cfg.modelstring = "# Choice of models must be one of G, G+GxE, or G+GxE+NxE. Using default G+GxE+NxE model\n# Estimating G and GxE heritability (with heterogeneous noise)";
					}
					i++;
				} else if (isarg(argv[i], "-s", "--seed")) {
					command_line_opts.seed = atoi(argv[i+1]);
					cfg.insertKey ("seed", argv[i+1], 0);
					i++;
				} else if (isarg(argv[i], "-np", "--norm-proj-pheno")) {
					int flag = atoi(argv[i+1]);
					if (flag == 0) 
						command_line_opts.normalize_proj_pheno = false;
					else 
						command_line_opts.normalize_proj_pheno = true;
					cfg.insertKey ("norm-proj-pheno", argv[i+1], 0);
					i++;
				} else if (isarg(argv[i], "-i", "--cov-add-intercept")) {
					int flag = atoi(argv[i+1]);
					if (flag == 0) {
						command_line_opts.cov_add_intercept = false;
					} else {
						command_line_opts.cov_add_intercept = true;
					}
					cfg.insertKey ("cov-add-intercept", argv[i+1], 0);
					i++;
				} else if (isarg(argv[i], "-t", "--nthreads")) {
					command_line_opts.nthreads = atoi(argv[i+1]);
                    cfg.insertKey ("nthreads", argv[i+1], 0);
					i++;
				} else if (isarg(argv[i], "-v", "--verbose")) {
					command_line_opts.verbose = atoi(argv[i+1]);
                    cfg.insertKey ("verbose", argv[i+1], 0);
					i++;
				} else if (isarg(argv[i], "-tr", "--trace")){
                    command_line_opts.print_trace = true;
                    cfg.insertKey ("trace", "1", 0);
					cfg.outputstring += "# Printing trace summaries";
				} else if (strcmp(argv[i], "-tr_input")==0) {
					command_line_opts.TRACE_FILE_PATH = string(argv[i+1]);
					cfg.insertKey ("trace_input", argv[i+1], 0);
					i++;
				}
				else if(isarg(argv[i], "-miss", "--missing"))
					command_line_opts.missing = true;
				else if(isarg(argv[i], "-nfm", "--no-fast-mode"))
					command_line_opts.fast_mode = false;
				else if(isarg(argv[i],"-eXa", "--eXannot")) {
					command_line_opts.exannot = true;
                    cfg.insertKey ("eXannot", "1", 0);
                } else if (isarg(argv[i],"-mem", "--memeff") ) {
                    command_line_opts.memeff = true;        
                    cfg.insertKey ("memeff", "1", 0);
                } else if (isarg(argv[i],"-opt2", "--opt2") ) {
                    command_line_opts.opt2 = (atoi (argv[i+1])>0);        
                    cfg.insertKey ("opt2",Convert::T_to_string (command_line_opts.opt2), 0);
                } else if (isarg(argv[i],"-opt1", "--opt1") ) {
                    command_line_opts.opt1 = (atoi (argv[i+1])>0);        
                    cfg.insertKey ("opt1",Convert::T_to_string (command_line_opts.opt1), 0);
                } else if (isarg(argv[i],"-mem_Nsnp", "--mem_Nsnp") ) {
                    command_line_opts.mem_Nsnp = atoi(argv[i+1]);     
                    cfg.insertKey ("mem_Nsnp",Convert::T_to_string (command_line_opts.mem_Nsnp), 0);
                } else {
					cerr << version_string << endl;
					cerr << "ERROR: Not Enough or Invalid arguments: '"<< argv[i] << "'" << endl;
					exitWithError (usage());
				}
			} else if (isarg(argv[i], "-tr", "--trace")){
				command_line_opts.print_trace = true;
				cfg.insertKey ("trace", "1", 0);
				cfg.outputstring += "# Printing trace summaries";
            } else if (isarg(argv[i],"-mem", "--memeff") ) {
                command_line_opts.memeff = true;
                cfg.insertKey ("memeff", "1", 0);
			} else if(isarg(argv[i], "-nfm", "--no-fast-mode"))
				command_line_opts.fast_mode = false;
			else if(isarg(argv[i], "-miss", "--missing"))
				command_line_opts.missing = true;
		}
	}
	

	if (command_line_opts.cov_add_intercept)
		cfg.outputstring += "# The intercept term *will* be added to any fixed-effect covariates";
	else
		cfg.outputstring += "# The intercept term *will not* be added to any fixed-effect covariates";

	if (command_line_opts.normalize_proj_pheno)
		cfg.outputstring += "# After projecting out the covariates, the phenotype vector *will* be standardized\n";
	else
		cfg.outputstring += "# After projecting out the covariates, the phenotype vector *will not* be standardized\n";
	param_string = cfg.paramstring (command_line_opts.verbose);
	
	cout << version_string << endl;

	if (!got_genotype_file){
		cerr << "ERROR: Missing genotype file" << endl;
		exitWithError (usage());
	}
	// if ((command_line_opts.print_trace) && (command_line_opts.model != "G")){
	// 	cerr << "Trace summary is only supported for G only model. Disabling --trace option." << endl;
	// 	command_line_opts.print_trace = false;
	// }

	if (!got_phenotype_file){
		if (command_line_opts.print_trace){
			cout << "Estimating trace summary without phenotype input (will be using dummy phenotype)" << endl;
			command_line_opts.use_dummy_pheno = true;
			if (got_covariate_file){
				cerr << "ERROR: Trace summary estimation with dummy phenotype does not support covariates at the moment." << endl;
				exitWithError (usage());
			}
		}
		else {
			cerr << "ERROR: Missing phenotype file" << endl;
			exitWithError (usage());
		}
	}
	const char* model_arg = command_line_opts.model.c_str();
	if (strcmp(model_arg, "G+GxE")==0 || strcmp (model_arg, "G+GxE+NxE")==0) {
		if (!got_environment_file){
			cerr << "ERROR: Missing file with environmental variables" << endl;
			exitWithError (usage());
		}
	}

	cout << param_string << endl;
}

#endif
