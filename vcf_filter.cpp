#include "vcf_filter.h"
#include <iomanip>
#include <ctime>

// Utility function to split strings
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Utility function to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

// Constructor
VCFFilter::VCFFilter(const std::string& vcf, const std::string& samples, 
                     const std::string& output, const std::string& log,
                     const FilterParams& p)
    : vcf_file(vcf), sample_list_file(samples), output_vcf_file(output), 
      log_file(log), params(p) {}

// Parse the sample list file
bool VCFFilter::parseSampleList() {
    std::ifstream infile(sample_list_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open sample list file: " << sample_list_file << std::endl;
        return false;
    }
    
    std::string line;
    int line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        line = trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 2) {
            std::cerr << "Warning: Line " << line_num << " in sample list has fewer than 2 columns, skipping." << std::endl;
            continue;
        }
        
        std::string sample_name = trim(fields[0]);
        std::string symbol = trim(fields[1]);
        
        bool include = (symbol == "+" || symbol == "1" || symbol == "include" || symbol == "INCLUDE");
        
        samples.push_back(SampleInfo(sample_name, include));
    }
    
    infile.close();
    
    std::cout << "Loaded " << samples.size() << " samples from sample list." << std::endl;
    int included = 0;
    for (const auto& s : samples) {
        if (s.include) included++;
    }
    std::cout << "  " << included << " samples will be evaluated for filtering." << std::endl;
    std::cout << "  " << (samples.size() - included) << " samples will be ignored." << std::endl;
    
    return true;
}

// Parse VCF header and identify sample columns
bool VCFFilter::parseVCFHeader(std::ifstream& vcf_in) {
    std::string line;
    while (std::getline(vcf_in, line)) {
        if (line.empty()) continue;
        
        // Store all header lines
        if (line[0] == '#') {
            vcf_header.push_back(line);
            
            // Check if this is the column header line
            if (line.substr(0, 6) == "#CHROM") {
                std::vector<std::string> fields = split(line, '\t');
                
                // Sample columns start at index 9
                for (size_t i = 9; i < fields.size(); i++) {
                    sample_indices[fields[i]] = i;
                }
                
                std::cout << "Found " << sample_indices.size() << " samples in VCF." << std::endl;
                return true;
            }
        } else {
            // Non-header line reached, but we haven't found the column header
            std::cerr << "Error: Could not find #CHROM header line in VCF." << std::endl;
            return false;
        }
    }
    
    std::cerr << "Error: VCF file appears to be empty or malformed." << std::endl;
    return false;
}

// Extract frequency from VarScan FREQ field
double VCFFilter::extractFREQ(const std::string& format, const std::string& sample_data) {
    std::vector<std::string> format_fields = split(format, ':');
    std::vector<std::string> data_fields = split(sample_data, ':');
    
    if (format_fields.size() != data_fields.size()) {
        return -1.0;  // Return -1 on mismatch
    }
    
    // Find FREQ field
    for (size_t i = 0; i < format_fields.size(); i++) {
        if (format_fields[i] == "FREQ") {
            std::string freq_str = data_fields[i];
            // Remove the '%' character if present
            if (!freq_str.empty() && freq_str.back() == '%') {
                freq_str.pop_back();
            }
            try {
                return std::stod(freq_str);
            } catch (...) {
                return -1.0;
            }
        }
    }
    
    return -1.0;  // FREQ not found
}

// Extract RD and AD from VarScan format (returns pair of <RD, AD>)
std::pair<int, int> VCFFilter::extractVarScanDepths(const std::string& format, const std::string& sample_data) {
    std::vector<std::string> format_fields = split(format, ':');
    std::vector<std::string> data_fields = split(sample_data, ':');
    
    if (format_fields.size() != data_fields.size()) {
        return {0, 0};  // Return zeros on mismatch
    }
    
    int rd = 0, ad = 0;
    bool found_rd = false, found_ad = false;
    
    // Find RD and AD fields
    for (size_t i = 0; i < format_fields.size(); i++) {
        if (format_fields[i] == "RD") {
            try {
                rd = std::stoi(data_fields[i]);
                found_rd = true;
            } catch (...) {
                rd = 0;
            }
        } else if (format_fields[i] == "AD") {
            try {
                ad = std::stoi(data_fields[i]);
                found_ad = true;
            } catch (...) {
                ad = 0;
            }
        }
        
        if (found_rd && found_ad) break;
    }
    
    return {rd, ad};
}

// Extract AD values from FORMAT field
std::vector<int> VCFFilter::extractAD(const std::string& format, const std::string& sample_data) {
    std::vector<int> ad_values;
    
    std::vector<std::string> format_fields = split(format, ':');
    std::vector<std::string> data_fields = split(sample_data, ':');
    
    if (format_fields.size() != data_fields.size()) {
        return ad_values;  // Return empty vector on mismatch
    }
    
    // Find AD field
    for (size_t i = 0; i < format_fields.size(); i++) {
        if (format_fields[i] == "AD") {
            std::vector<std::string> ad_strings = split(data_fields[i], ',');
            for (const auto& ad_str : ad_strings) {
                try {
                    if (ad_str == "." || ad_str.empty()) {
                        ad_values.push_back(0);
                    } else {
                        ad_values.push_back(std::stoi(ad_str));
                    }
                } catch (...) {
                    ad_values.push_back(0);
                }
            }
            break;
        }
    }
    
    return ad_values;
}

// Extract GT (genotype) from sample data
std::string VCFFilter::extractGT(const std::string& format, const std::string& sample_data) {
    std::vector<std::string> format_fields = split(format, ':');
    std::vector<std::string> data_fields = split(sample_data, ':');
    
    if (format_fields.size() != data_fields.size()) {
        return "./.";  // Return missing genotype on mismatch
    }
    
    // GT is typically the first field, but we'll search for it
    for (size_t i = 0; i < format_fields.size(); i++) {
        if (format_fields[i] == "GT") {
            return data_fields[i];
        }
    }
    
    return "./.";
}

// Check if a variant passes filters for a specific sample
bool VCFFilter::passesFilter(const std::string& gt, const std::vector<int>& ad, std::string& reason) {
    // Missing genotype or AD
    if (gt == "./." || gt == ".|." || ad.empty()) {
        reason = "Missing GT or AD";
        return true;  // Don't filter on missing data
    }
    
    // Calculate total depth
    int total_depth = 0;
    for (int depth : ad) {
        total_depth += depth;
    }
    
    if (total_depth == 0) {
        reason = "Zero total depth";
        return true;  // Don't filter on zero depth
    }
    
    // Parse genotype (assuming diploid: 0/0, 0/1, 1/1, etc.)
    std::string gt_clean = gt;
    std::replace(gt_clean.begin(), gt_clean.end(), '|', '/');
    std::vector<std::string> alleles = split(gt_clean, '/');
    
    if (alleles.size() != 2) {
        reason = "Non-diploid genotype";
        return true;  // Don't filter non-diploid calls
    }
    
    int allele1 = -1, allele2 = -1;
    try {
        allele1 = (alleles[0] == ".") ? -1 : std::stoi(alleles[0]);
        allele2 = (alleles[1] == ".") ? -1 : std::stoi(alleles[1]);
    } catch (...) {
        reason = "Invalid genotype format";
        return true;
    }
    
    if (allele1 < 0 || allele2 < 0) {
        reason = "Missing allele";
        return true;
    }
    
    // Check if we have AD values for these alleles
    if (allele1 >= (int)ad.size() || allele2 >= (int)ad.size()) {
        reason = "AD array too small for genotype";
        return true;
    }
    
    // Homozygous reference (0/0)
    if (allele1 == 0 && allele2 == 0) {
        double ref_pct = (ad[0] * 100.0) / total_depth;
        if (ref_pct < params.homozygous_purity) {
            std::stringstream ss;
            ss << "Hom ref purity " << std::fixed << std::setprecision(1) << ref_pct 
               << "% < " << params.homozygous_purity << "%";
            reason = ss.str();
            return false;
        }
        return true;
    }
    
    // Homozygous alternate (1/1, 2/2, etc.)
    if (allele1 == allele2 && allele1 > 0) {
        double alt_pct = (ad[allele1] * 100.0) / total_depth;
        if (alt_pct < params.homozygous_purity) {
            std::stringstream ss;
            ss << "Hom alt purity " << std::fixed << std::setprecision(1) << alt_pct 
               << "% < " << params.homozygous_purity << "%";
            reason = ss.str();
            return false;
        }
        return true;
    }
    
    // Heterozygous (0/1, 1/2, etc.)
    int depth1 = ad[allele1];
    int depth2 = ad[allele2];
    
    double ratio1 = (depth1 * 100.0) / total_depth;
    double ratio2 = (depth2 * 100.0) / total_depth;
    
    // Check if either allele falls within the acceptable range
    bool passes = (ratio1 >= params.het_min_ratio && ratio1 <= params.het_max_ratio) ||
                  (ratio2 >= params.het_min_ratio && ratio2 <= params.het_max_ratio);
    
    if (!passes) {
        std::stringstream ss;
        ss << "Het ratio " << std::fixed << std::setprecision(1) << ratio1 << "/" << ratio2 
           << " outside " << params.het_min_ratio << "-" << params.het_max_ratio << "%";
        reason = ss.str();
        return false;
    }
    
    return true;
}

// Check if a variant passes filters for VarScan format (using FREQ)
bool VCFFilter::passesFilterVarScan(const std::string& gt, double freq, std::string& reason) {
    // Missing genotype or frequency
    if (gt == "./." || gt == ".|." || freq < 0) {
        reason = "Missing GT or FREQ";
        return true;  // Don't filter on missing data
    }
    
    // Parse genotype (assuming diploid: 0/0, 0/1, 1/1, etc.)
    std::string gt_clean = gt;
    std::replace(gt_clean.begin(), gt_clean.end(), '|', '/');
    std::vector<std::string> alleles = split(gt_clean, '/');
    
    if (alleles.size() != 2) {
        reason = "Non-diploid genotype";
        return true;  // Don't filter non-diploid calls
    }
    
    int allele1 = -1, allele2 = -1;
    try {
        allele1 = (alleles[0] == ".") ? -1 : std::stoi(alleles[0]);
        allele2 = (alleles[1] == ".") ? -1 : std::stoi(alleles[1]);
    } catch (...) {
        reason = "Invalid genotype format";
        return true;
    }
    
    if (allele1 < 0 || allele2 < 0) {
        reason = "Missing allele";
        return true;
    }
    
    // Homozygous reference (0/0)
    // FREQ represents variant frequency, so for 0/0 we want FREQ to be near 0%
    if (allele1 == 0 && allele2 == 0) {
        double allowed_variant_freq = 100.0 - params.homozygous_purity;
        if (freq > allowed_variant_freq) {
            std::stringstream ss;
            ss << "Hom ref has variant freq " << std::fixed << std::setprecision(1) << freq 
               << "% > " << allowed_variant_freq << "%";
            reason = ss.str();
            return false;
        }
        return true;
    }
    
    // Homozygous alternate (1/1, 2/2, etc.)
    // FREQ represents variant frequency, so for 1/1 we want FREQ to be near 100%
    if (allele1 == allele2 && allele1 > 0) {
        if (freq < params.homozygous_purity) {
            std::stringstream ss;
            ss << "Hom alt variant freq " << std::fixed << std::setprecision(1) << freq 
               << "% < " << params.homozygous_purity << "%";
            reason = ss.str();
            return false;
        }
        return true;
    }
    
    // Heterozygous (0/1, 1/2, etc.)
    // For heterozygous, FREQ should be within the acceptable range
    if (freq < params.het_min_ratio || freq > params.het_max_ratio) {
        std::stringstream ss;
        ss << "Het variant freq " << std::fixed << std::setprecision(1) << freq 
           << "% outside " << params.het_min_ratio << "-" << params.het_max_ratio << "%";
        reason = ss.str();
        return false;
    }
    
    return true;
}

// Check if variant line should be retained (checks all included samples)
bool VCFFilter::shouldRetainVariant(const std::vector<std::string>& fields, std::string& filter_reason) {
    if (fields.size() < 10) {
        filter_reason = "Malformed VCF line";
        return false;
    }
    
    std::string format = fields[8];
    
    // Check each sample that should be included
    for (const auto& sample : samples) {
        if (!sample.include) continue;
        
        // Find this sample's column in the VCF
        auto it = sample_indices.find(sample.name);
        if (it == sample_indices.end()) {
            std::cerr << "Warning: Sample '" << sample.name << "' not found in VCF." << std::endl;
            continue;
        }
        
        int col_idx = it->second;
        if (col_idx >= (int)fields.size()) {
            filter_reason = "Sample column out of range";
            return false;
        }
        
        std::string sample_data = fields[col_idx];
        std::string gt = extractGT(format, sample_data);
        
        std::string reason;
        bool passes = false;
        
        if (params.vcf_type == VCFType::VARSCAN) {
            // Use FREQ field for VarScan
            double freq = extractFREQ(format, sample_data);
            passes = passesFilterVarScan(gt, freq, reason);
        } else {
            // Use AD field for standard VCF
            std::vector<int> ad = extractAD(format, sample_data);
            passes = passesFilter(gt, ad, reason);
        }
        
        if (!passes) {
            filter_reason = "Sample " + sample.name + ": " + reason;
            return false;
        }
    }
    
    return true;
}

// Write log file with statistics
void VCFFilter::writeLogFile() {
    std::ofstream log_out(log_file);
    if (!log_out.is_open()) {
        std::cerr << "Error: Cannot open log file: " << log_file << std::endl;
        return;
    }
    
    // Get current time
    std::time_t now = std::time(nullptr);
    log_out << "VCF Filter Log\n";
    log_out << "Generated: " << std::ctime(&now);
    log_out << "===============================================\n\n";
    
    log_out << "Input Files:\n";
    log_out << "  VCF: " << vcf_file << "\n";
    log_out << "  Sample List: " << sample_list_file << "\n\n";
    
    log_out << "Output Files:\n";
    log_out << "  Filtered VCF: " << output_vcf_file << "\n";
    log_out << "  Log: " << log_file << "\n\n";
    
    log_out << "Filter Parameters:\n";
    log_out << "  VCF Type: " << (params.vcf_type == VCFType::VARSCAN ? "VarScan2" : "Standard") << "\n";
    log_out << "  Homozygous purity: " << params.homozygous_purity << "%\n";
    log_out << "  Heterozygous ratio range: " << params.het_min_ratio << "% - " 
            << params.het_max_ratio << "%\n\n";
    
    log_out << "Samples Evaluated:\n";
    for (const auto& sample : samples) {
        if (sample.include) {
            log_out << "  - " << sample.name << "\n";
        }
    }
    log_out << "\n";
    
    log_out << "Summary Statistics:\n";
    log_out << "  Total variants: " << stats.total_variants << "\n";
    log_out << "  Retained variants: " << stats.retained_variants << " ("
            << std::fixed << std::setprecision(2) 
            << (stats.total_variants > 0 ? (stats.retained_variants * 100.0 / stats.total_variants) : 0)
            << "%)\n";
    log_out << "  Filtered variants: " << stats.filtered_variants << " ("
            << std::fixed << std::setprecision(2)
            << (stats.total_variants > 0 ? (stats.filtered_variants * 100.0 / stats.total_variants) : 0)
            << "%)\n\n";
    
    if (!stats.filter_reasons.empty()) {
        log_out << "Filter Reasons:\n";
        for (const auto& pair : stats.filter_reasons) {
            log_out << "  " << pair.first << ": " << pair.second << "\n";
        }
    }
    
    log_out.close();
    std::cout << "Log file written to: " << log_file << std::endl;
}

// Main filtering function
bool VCFFilter::filterVCF() {
    // Parse sample list
    if (!parseSampleList()) {
        return false;
    }
    
    // Open input VCF
    std::ifstream vcf_in(vcf_file);
    if (!vcf_in.is_open()) {
        std::cerr << "Error: Cannot open VCF file: " << vcf_file << std::endl;
        return false;
    }
    
    // Parse VCF header
    if (!parseVCFHeader(vcf_in)) {
        vcf_in.close();
        return false;
    }
    
    // Open output VCF
    std::ofstream vcf_out(output_vcf_file);
    if (!vcf_out.is_open()) {
        std::cerr << "Error: Cannot create output VCF file: " << output_vcf_file << std::endl;
        vcf_in.close();
        return false;
    }
    
    // Write header to output
    for (const auto& header_line : vcf_header) {
        vcf_out << header_line << "\n";
    }
    
    // Add filter info to header
    vcf_out << "##FILTER=<ID=AD_FILTER,Description=\"Filtered based on allelic depth ratios\">\n";
    
    // Process variants
    std::string line;
    std::cout << "Processing variants..." << std::endl;
    
    while (std::getline(vcf_in, line)) {
        if (line.empty()) continue;
        
        stats.total_variants++;
        
        std::vector<std::string> fields = split(line, '\t');
        std::string filter_reason;
        
        if (shouldRetainVariant(fields, filter_reason)) {
            vcf_out << line << "\n";
            stats.retained_variants++;
        } else {
            stats.addFilterReason(filter_reason);
        }
        
        // Progress indicator
        if (stats.total_variants % 10000 == 0) {
            std::cout << "  Processed " << stats.total_variants << " variants..." << std::endl;
        }
    }
    
    vcf_in.close();
    vcf_out.close();
    
    std::cout << "Filtering complete!" << std::endl;
    std::cout << "  Total variants: " << stats.total_variants << std::endl;
    std::cout << "  Retained: " << stats.retained_variants << std::endl;
    std::cout << "  Filtered: " << stats.filtered_variants << std::endl;
    
    // Write log file
    writeLogFile();
    
    return true;
}
