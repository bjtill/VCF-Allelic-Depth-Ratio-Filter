#ifndef VCF_FILTER_H
#define VCF_FILTER_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

// Enum for VCF format types
enum class VCFType {
    STANDARD,   // Standard VCF with AD field (comma-separated allele depths)
    VARSCAN     // VarScan2 VCF with RD, AD, and FREQ fields
};

// Structure to hold filter parameters
struct FilterParams {
    double homozygous_purity;  // Required purity for homozygous calls (e.g., 100.0 = 100%)
    double het_min_ratio;       // Minimum ratio for heterozygous calls (e.g., 40.0 = 40%)
    double het_max_ratio;       // Maximum ratio for heterozygous calls (e.g., 60.0 = 60%)
    VCFType vcf_type;          // Type of VCF format
    
    FilterParams() : homozygous_purity(100.0), het_min_ratio(40.0), het_max_ratio(60.0), 
                     vcf_type(VCFType::STANDARD) {}
};

// Structure to hold sample information
struct SampleInfo {
    std::string name;
    bool include;  // true if sample should be evaluated, false if ignored
    
    SampleInfo() : name(""), include(false) {}
    SampleInfo(const std::string& n, bool inc) : name(n), include(inc) {}
};

// Structure to hold filtering statistics
struct FilterStats {
    int total_variants;
    int retained_variants;
    int filtered_variants;
    std::map<std::string, int> filter_reasons;  // Reason -> count
    
    FilterStats() : total_variants(0), retained_variants(0), filtered_variants(0) {}
    
    void addFilterReason(const std::string& reason) {
        filtered_variants++;
        filter_reasons[reason]++;
    }
};

// Class to handle VCF filtering
class VCFFilter {
private:
    std::string vcf_file;
    std::string sample_list_file;
    std::string output_vcf_file;
    std::string log_file;
    FilterParams params;
    
    std::vector<std::string> vcf_header;
    std::vector<SampleInfo> samples;
    std::map<std::string, int> sample_indices;  // sample name -> column index in VCF
    FilterStats stats;
    
    // Parse the sample list file
    bool parseSampleList();
    
    // Parse VCF header and identify sample columns
    bool parseVCFHeader(std::ifstream& vcf_in);
    
    // Extract AD values from FORMAT field
    std::vector<int> extractAD(const std::string& format, const std::string& sample_data);
    
    // Extract frequency from VarScan FREQ field
    double extractFREQ(const std::string& format, const std::string& sample_data);
    
    // Extract RD and AD from VarScan format
    std::pair<int, int> extractVarScanDepths(const std::string& format, const std::string& sample_data);
    
    // Extract GT (genotype) from sample data
    std::string extractGT(const std::string& format, const std::string& sample_data);
    
    // Check if a variant passes filters for a specific sample
    bool passesFilter(const std::string& gt, const std::vector<int>& ad, std::string& reason);
    
    // Check if a variant passes filters for VarScan format (using FREQ)
    bool passesFilterVarScan(const std::string& gt, double freq, std::string& reason);
    
    // Check if variant line should be retained (checks all included samples)
    bool shouldRetainVariant(const std::vector<std::string>& fields, std::string& filter_reason);
    
    // Write log file with statistics
    void writeLogFile();
    
public:
    VCFFilter(const std::string& vcf, const std::string& samples, 
              const std::string& output, const std::string& log,
              const FilterParams& p = FilterParams());
    
    // Main filtering function
    bool filterVCF();
    
    // Setters for parameters
    void setHomozygousPurity(double purity) { params.homozygous_purity = purity; }
    void setHetRange(double min_ratio, double max_ratio) { 
        params.het_min_ratio = min_ratio; 
        params.het_max_ratio = max_ratio; 
    }
};

// Utility functions
std::vector<std::string> split(const std::string& str, char delimiter);
std::string trim(const std::string& str);

#endif // VCF_FILTER_H
