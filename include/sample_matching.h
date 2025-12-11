#ifndef SAMPLE_MATCHING_H
#define SAMPLE_MATCHING_H

// BUG-003: Sample ID matching across input files
// This module implements PLINK2-style sample matching by FID/IID
// instead of assuming files are pre-aligned by row position.

#include <string>
#include <vector>
#include <unordered_map>

// Structure to hold FAM file sample index
// This is the "master" sample list - genotype defines the sample universe
struct FamIndex {
    std::vector<std::string> ids;                    // FID_IID for each FAM row
    std::unordered_map<std::string, int> id_to_idx;  // Lookup: ID -> FAM row index
    int count;                                        // Total samples in FAM
};

// Structure to track matching statistics for a single file
struct MatchStats {
    int total_in_file;       // Total samples in the file
    int matched;             // Samples found in FAM
    int not_in_fam;          // Samples in file but not in FAM (dropped)
    bool already_aligned;    // True if file order matches FAM order exactly
};

// Structure to hold overall sample matching results
struct SampleMatchResult {
    int fam_count;           // Total samples in FAM (genotype)
    int pheno_count;         // Total samples in phenotype file
    int covar_count;         // Total samples in covariate file (0 if no file)
    int env_count;           // Total samples in environment file (0 if no file)

    int intersection_count;  // Samples in ALL required files
    int fam_not_in_pheno;    // FAM samples missing from phenotype
    int fam_not_in_covar;    // FAM samples missing from covariate
    int fam_not_in_env;      // FAM samples missing from environment
    int pheno_not_in_fam;    // Phenotype samples not in FAM
    int covar_not_in_fam;    // Covariate samples not in FAM
    int env_not_in_fam;      // Environment samples not in FAM

    bool all_aligned;        // True if all files were already aligned
};

// Read FAM file and build sample index
// Returns FamIndex with IDs and hash map for O(1) lookup
FamIndex read_fam_ids(const std::string& fam_path);

// Read phenotype file with ID-based matching
// Loads data directly into FAM-order positions
// Sets mask=0 for FAM samples not in phenotype (or with NA/-9 values)
// Returns MatchStats for this file
// phenocount_out: output parameter for number of phenotype columns
MatchStats read_pheno_matched(
    const FamIndex& fam_idx,
    const std::string& filename,
    int& phenocount_out
);

// Read covariate file with ID-based matching
// Returns number of covariate columns (same semantics as original read_cov)
// stats_out: output parameter for matching statistics
int read_cov_matched(
    const FamIndex& fam_idx,
    const std::string& filename,
    MatchStats& stats_out
);

// Read environment file with ID-based matching
// Returns number of environment columns (same semantics as original read_env)
// stats_out: output parameter for matching statistics
int read_env_matched(
    const FamIndex& fam_idx,
    const std::string& filename,
    MatchStats& stats_out
);

// Print PLINK2-style sample matching report
void report_sample_matching(const SampleMatchResult& result);

// Compute final intersection count from mask
// Returns number of samples with mask=1
int compute_intersection_count();

#endif // SAMPLE_MATCHING_H
