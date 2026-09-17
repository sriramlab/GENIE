#include "sample_matching.h"
#include "auxillary.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//=============================================================================
// read_fam_ids: Read FAM file and build sample index
//=============================================================================
FamIndex read_fam_ids(const std::string& fam_path) {
	FamIndex idx;
	idx.count = 0;

	ifstream ifs(fam_path.c_str(), ios::in);
	if (!ifs.is_open()) {
		cerr << "ERROR: Cannot open FAM file: " << fam_path << endl;
		exit(1);
	}

	string line;
	while (getline(ifs, line)) {
		if (line.empty()) continue;

		istringstream iss(line);
		string fid, iid;
		iss >> fid >> iid;

		string id = fid + "_" + iid;

		// Check for duplicate IDs
		if (idx.id_to_idx.count(id) > 0) {
			cerr << "ERROR: Duplicate sample ID in FAM file: " << id << endl;
			cerr << "  First occurrence: row " << (idx.id_to_idx[id] + 1) << endl;
			cerr << "  Second occurrence: row " << (idx.count + 1) << endl;
			cerr << "  Hint: Sample IDs must be unique." << endl;
			exit(1);
		}

		idx.ids.push_back(id);
		idx.id_to_idx[id] = idx.count;
		idx.count++;
	}

	ifs.close();
	return idx;
}

//=============================================================================
// read_pheno_matched: Read phenotype with ID-based matching
//=============================================================================
MatchStats read_pheno_matched(
	const FamIndex& fam_idx,
	const std::string& filename,
	int& phenocount_out
) {
	MatchStats stats = {0, 0, 0, true};

	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()) {
		cerr << "ERROR: Cannot open phenotype file: " << filename << endl;
		exit(1);
	}

	// Parse header to count phenotypes (skip FID, IID columns)
	string line;
	getline(ifs, line);
	istringstream header_iss(line);
	string col;
	phenocount_out = 0;
	while (header_iss >> col) {
		if (col != "FID" && col != "IID") {
			phenocount_out++;
		}
	}

	if (phenocount_out > 1) {
		cerr << "ERROR: The phenotype file can only have a single phenotype" << endl;
		exit(1);
	}

	// Initialize matrices sized to FAM count
	int Nindv = fam_idx.count;
	pheno.resize(Nindv, phenocount_out);
	new_pheno.resize(Nindv, phenocount_out);
	mask.resize(Nindv, phenocount_out);

	pheno.setZero();
	new_pheno.setZero();
	mask.setZero();  // Default: all masked out

	// Track order for alignment check
	int last_matched_fam_row = -1;

	// Read data lines
	while (getline(ifs, line)) {
		if (line.empty()) continue;

		istringstream iss(line);
		string fid, iid;
		iss >> fid >> iid;
		string id = fid + "_" + iid;

		stats.total_in_file++;

		// Lookup in FAM index
		auto it = fam_idx.id_to_idx.find(id);
		if (it == fam_idx.id_to_idx.end()) {
			// Sample not in FAM - skip
			stats.not_in_fam++;
			continue;
		}

		int fam_row = it->second;
		stats.matched++;

		// Check alignment: rows must appear in ascending FAM order
		if (fam_row <= last_matched_fam_row) {
			stats.already_aligned = false;
		}
		last_matched_fam_row = fam_row;

		// Read phenotype values into FAM-order position
		for (int j = 0; j < phenocount_out; j++) {
			string temp;
			iss >> temp;
			double val = atof(temp.c_str());

			if (temp == "NA" || val == -9) {
				pheno(fam_row, j) = 0;
				mask(fam_row, j) = 0;  // Missing phenotype
			} else {
				pheno(fam_row, j) = val;
				mask(fam_row, j) = 1;  // Valid
			}
		}
	}

	ifs.close();

	// Check alignment: must have matched all FAM samples in order with no extras
	if (stats.matched != fam_idx.count || stats.not_in_fam > 0) {
		stats.already_aligned = false;
	}

	return stats;
}

//=============================================================================
// read_env_matched: Read environment with ID-based matching
//=============================================================================
int read_env_matched(
	const FamIndex& fam_idx,
	const std::string& filename,
	MatchStats& stats_out
) {
	stats_out = {0, 0, 0, true};

	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()) {
		cerr << "ERROR: Cannot open environment file: " << filename << endl;
		exit(1);
	}

	// Parse header
	string line;
	getline(ifs, line);
	istringstream header_iss(line);
	string col;
	int envNum = 0;
	while (header_iss >> col) {
		if (col != "FID" && col != "IID") {
			envNum++;
		}
	}

	// Initialize environment matrix
	int Nindv = fam_idx.count;

	if (gen_by_env) {
		Enviro.resize(Nindv, envNum);
		Enviro.setZero();
		cout << "Reading in " << envNum << " environmental variables ..." << endl;
	}

	int last_matched_fam_row = -1;

	// Read data
	while (getline(ifs, line)) {
		if (line.empty()) continue;

		istringstream iss(line);
		string fid, iid;
		iss >> fid >> iid;
		string id = fid + "_" + iid;

		stats_out.total_in_file++;

		auto it = fam_idx.id_to_idx.find(id);
		if (it == fam_idx.id_to_idx.end()) {
			stats_out.not_in_fam++;
			continue;
		}

		int fam_row = it->second;
		stats_out.matched++;

		if (fam_row <= last_matched_fam_row) {
			stats_out.already_aligned = false;
		}
		last_matched_fam_row = fam_row;

		// Read environment values
		for (int k = 0; k < envNum; k++) {
			string temp;
			iss >> temp;

			if (temp == "NA") {
				mask(fam_row, 0) = 0;  // Exclude sample
				continue;
			}
			double val = atof(temp.c_str());
			if (val == -9) {
				mask(fam_row, 0) = 0;  // Exclude sample
				continue;
			}

			if (gen_by_env) {
				Enviro(fam_row, k) = val;
			}
		}
	}

	ifs.close();

	if (stats_out.matched != fam_idx.count || stats_out.not_in_fam > 0) {
		stats_out.already_aligned = false;
	}

	return envNum;
}

//=============================================================================
// read_cov_matched: Read covariate with ID-based matching
//=============================================================================
int read_cov_matched(
	const FamIndex& fam_idx,
	const std::string& filename,
	MatchStats& stats_out
) {
	stats_out = {0, 0, 0, true};

	ifstream ifs(filename.c_str(), ios::in);
	if (!ifs.is_open()) {
		cerr << "ERROR: Cannot open covariate file: " << filename << endl;
		exit(1);
	}

	// Parse header
	string line;
	getline(ifs, line);
	istringstream header_iss(line);
	string col;
	int covNum = 0;
	while (header_iss >> col) {
		if (col != "FID" && col != "IID") {
			covNum++;
		}
	}

	// Determine final covariate count (with env and intercept if applicable)
	int Nenv_cols = 0;
	if (add_env_to_cov && gen_by_env && Enviro.cols() > 0) {
		Nenv_cols = Enviro.cols();
	}

	int total_cov_cols = covNum;
	if (add_env_to_cov) total_cov_cols += Nenv_cols;
	if (cov_add_intercept) total_cov_cols += 1;

	// Initialize covariate matrix
	int Nindv = fam_idx.count;
	covariate.resize(Nindv, total_cov_cols);
	covariate.setZero();

	// Track for missing value imputation
	vector<vector<int> > missing(covNum);
	vector<double> cov_sum(covNum, 0.0);
	vector<int> cov_count(covNum, 0);

	// Track which FAM samples have covariate data
	vector<bool> has_cov(Nindv, false);

	cout << "Reading in " << covNum << " covariates ..." << endl;

	int last_matched_fam_row = -1;

	// Read data
	while (getline(ifs, line)) {
		if (line.empty()) continue;

		istringstream iss(line);
		string fid, iid;
		iss >> fid >> iid;
		string id = fid + "_" + iid;

		stats_out.total_in_file++;

		auto it = fam_idx.id_to_idx.find(id);
		if (it == fam_idx.id_to_idx.end()) {
			stats_out.not_in_fam++;
			continue;
		}

		int fam_row = it->second;
		stats_out.matched++;
		has_cov[fam_row] = true;

		if (fam_row <= last_matched_fam_row) {
			stats_out.already_aligned = false;
		}
		last_matched_fam_row = fam_row;

		// Read covariate values
		for (int k = 0; k < covNum; k++) {
			string temp;
			iss >> temp;

			if (temp == "NA") {
				missing[k].push_back(fam_row);
				continue;
			}
			double val = atof(temp.c_str());
			if (val == -9) {
				missing[k].push_back(fam_row);
				continue;
			}

			cov_sum[k] += val;
			cov_count[k]++;
			covariate(fam_row, k) = val;
		}
	}

	ifs.close();

	// Impute missing values with mean (only for samples that were matched)
	for (int k = 0; k < covNum; k++) {
		if (cov_count[k] > 0) {
			double mean_val = cov_sum[k] / cov_count[k];
			for (size_t i = 0; i < missing[k].size(); i++) {
				int idx = missing[k][i];
				covariate(idx, k) = mean_val;
			}
		}
	}

	// Exclude samples not in covariate file from analysis
	for (int i = 0; i < Nindv; i++) {
		if (!has_cov[i]) {
			mask(i, 0) = 0;
		}
	}

	// Add environment to covariates if requested
	if (add_env_to_cov && Nenv_cols > 0) {
		for (int i = 0; i < Nenv_cols; i++) {
			covariate.col(covNum + i) = Enviro.col(i);
		}
	}

	// Add intercept if requested
	if (cov_add_intercept) {
		int intercept_col = covNum + (add_env_to_cov ? Nenv_cols : 0);
		for (int i = 0; i < Nindv; i++) {
			covariate(i, intercept_col) = 1.0;
		}
	}

	if (stats_out.matched != fam_idx.count || stats_out.not_in_fam > 0) {
		stats_out.already_aligned = false;
	}

	return total_cov_cols;
}

//=============================================================================
// report_sample_matching: Print PLINK2-style report
//=============================================================================
void report_sample_matching(const SampleMatchResult& result) {
	cout << endl;
	cout << "=== Sample Matching Report ===" << endl;
	cout << result.fam_count << " samples in genotype file (FAM)" << endl;
	cout << result.pheno_count << " samples in phenotype file" << endl;

	if (result.covar_count > 0) {
		cout << result.covar_count << " samples in covariate file" << endl;
	}
	if (result.env_count > 0) {
		cout << result.env_count << " samples in environment file" << endl;
	}

	cout << result.intersection_count << " samples remaining after intersection" << endl;

	// Report drops
	bool any_drops = false;
	if (result.fam_not_in_pheno > 0) {
		cout << "  (" << result.fam_not_in_pheno << " genotype samples not in phenotype)" << endl;
		any_drops = true;
	}
	if (result.pheno_not_in_fam > 0) {
		cout << "  (" << result.pheno_not_in_fam << " phenotype samples not in genotype)" << endl;
		any_drops = true;
	}
	if (result.covar_count > 0 && result.fam_not_in_covar > 0) {
		cout << "  (" << result.fam_not_in_covar << " genotype samples not in covariate)" << endl;
		any_drops = true;
	}
	if (result.covar_count > 0 && result.covar_not_in_fam > 0) {
		cout << "  (" << result.covar_not_in_fam << " covariate samples not in genotype)" << endl;
		any_drops = true;
	}
	if (result.env_count > 0 && result.fam_not_in_env > 0) {
		cout << "  (" << result.fam_not_in_env << " genotype samples not in environment)" << endl;
		any_drops = true;
	}
	if (result.env_count > 0 && result.env_not_in_fam > 0) {
		cout << "  (" << result.env_not_in_fam << " environment samples not in genotype)" << endl;
		any_drops = true;
	}

	if (result.all_aligned && !any_drops) {
		cout << "  (all files already aligned - no reordering needed)" << endl;
	}

	cout << "==============================" << endl;
	cout << endl;
}

//=============================================================================
// compute_intersection_count: Count samples with mask=1
//=============================================================================
int compute_intersection_count() {
	return static_cast<int>(mask.sum());
}
