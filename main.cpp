#include "vcf_filter.h"
#include <iostream>
#include <cstring>

void printUsage(const char* program_name) {
    std::cout << "VCF Allelic Depth Filter\n";
    std::cout << "Usage: " << program_name << " [options]\n\n";
    std::cout << "Required arguments:\n";
    std::cout << "  -v, --vcf FILE          Input multi-sample VCF file\n";
    std::cout << "  -s, --samples FILE      Sample list file (tab-delimited: sample_name symbol)\n";
    std::cout << "  -o, --output FILE       Output filtered VCF file\n\n";
    std::cout << "Optional arguments:\n";
    std::cout << "  -t, --vcf-type TYPE     VCF format type: 'standard' or 'varscan' (default: standard)\n";
    std::cout << "  -l, --log FILE          Log file (default: output.log)\n";
    std::cout << "  -p, --hom-purity PCT    Homozygous purity percentage (default: 100.0)\n";
    std::cout << "  -m, --het-min PCT       Heterozygous minimum ratio (default: 40.0)\n";
    std::cout << "  -M, --het-max PCT       Heterozygous maximum ratio (default: 60.0)\n";
    std::cout << "  -h, --help              Display this help message\n\n";
    std::cout << "Sample List Format:\n";
    std::cout << "  Column 1: Sample name (must match VCF header)\n";
    std::cout << "  Column 2: Symbol (+, 1, include, INCLUDE = evaluate; anything else = ignore)\n\n";
    std::cout << "VCF Types:\n";
    std::cout << "  standard - Standard VCF with AD field (comma-separated allele depths)\n";
    std::cout << "           Example: GT:AD format with AD=100,50 for het call\n";
    std::cout << "  varscan  - VarScan2 VCF with RD, AD, and FREQ fields\n";
    std::cout << "           Example: GT:RD:AD:FREQ format with FREQ=33.33%\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << program_name << " -v input.vcf -s samples.txt -o filtered.vcf -p 95 -m 35 -M 65\n";
    std::cout << "  " << program_name << " -v varscan.vcf -s samples.txt -o filtered.vcf -t varscan\n\n";
    std::cout << "Filter Logic:\n";
    std::cout << "  - Homozygous calls (0/0, 1/1): Require specified purity percentage\n";
    std::cout << "    Standard: 0/0 with AD=100,0 passes; AD=99,1 fails at 100% purity\n";
    std::cout << "    VarScan: 0/0 with FREQ=0% passes; FREQ=5% fails at 100% purity\n";
    std::cout << "  - Heterozygous calls (0/1, etc.): Require ratio within min-max range\n";
    std::cout << "    Standard: 0/1 with AD=40,60 passes; AD=70,30 fails with 40-60% range\n";
    std::cout << "    VarScan: 0/1 with FREQ=45% passes; FREQ=70% fails with 40-60% range\n";
    std::cout << "  - Variants are retained only if ALL included samples pass filters\n";
}

int main(int argc, char* argv[]) {
    std::string vcf_file;
    std::string sample_list;
    std::string output_file;
    std::string log_file = "output.log";
    
    FilterParams params;
    
    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if ((arg == "-v" || arg == "--vcf") && i + 1 < argc) {
            vcf_file = argv[++i];
        } else if ((arg == "-s" || arg == "--samples") && i + 1 < argc) {
            sample_list = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_file = argv[++i];
        } else if ((arg == "-t" || arg == "--vcf-type") && i + 1 < argc) {
            std::string type_str = argv[++i];
            // Convert to lowercase for case-insensitive comparison
            std::transform(type_str.begin(), type_str.end(), type_str.begin(), ::tolower);
            if (type_str == "varscan" || type_str == "varscan2") {
                params.vcf_type = VCFType::VARSCAN;
            } else if (type_str == "standard" || type_str == "std") {
                params.vcf_type = VCFType::STANDARD;
            } else {
                std::cerr << "Error: Invalid VCF type. Use 'standard' or 'varscan'." << std::endl;
                return 1;
            }
        } else if ((arg == "-l" || arg == "--log") && i + 1 < argc) {
            log_file = argv[++i];
        } else if ((arg == "-p" || arg == "--hom-purity") && i + 1 < argc) {
            try {
                params.homozygous_purity = std::stod(argv[++i]);
                if (params.homozygous_purity < 0 || params.homozygous_purity > 100) {
                    std::cerr << "Error: Homozygous purity must be between 0 and 100." << std::endl;
                    return 1;
                }
            } catch (...) {
                std::cerr << "Error: Invalid homozygous purity value." << std::endl;
                return 1;
            }
        } else if ((arg == "-m" || arg == "--het-min") && i + 1 < argc) {
            try {
                params.het_min_ratio = std::stod(argv[++i]);
                if (params.het_min_ratio < 0 || params.het_min_ratio > 100) {
                    std::cerr << "Error: Heterozygous minimum ratio must be between 0 and 100." << std::endl;
                    return 1;
                }
            } catch (...) {
                std::cerr << "Error: Invalid heterozygous minimum ratio value." << std::endl;
                return 1;
            }
        } else if ((arg == "-M" || arg == "--het-max") && i + 1 < argc) {
            try {
                params.het_max_ratio = std::stod(argv[++i]);
                if (params.het_max_ratio < 0 || params.het_max_ratio > 100) {
                    std::cerr << "Error: Heterozygous maximum ratio must be between 0 and 100." << std::endl;
                    return 1;
                }
            } catch (...) {
                std::cerr << "Error: Invalid heterozygous maximum ratio value." << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate required arguments
    if (vcf_file.empty() || sample_list.empty() || output_file.empty()) {
        std::cerr << "Error: Missing required arguments." << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    // Validate filter parameters
    if (params.het_min_ratio > params.het_max_ratio) {
        std::cerr << "Error: Heterozygous minimum ratio cannot be greater than maximum ratio." << std::endl;
        return 1;
    }
    
    // Display configuration
    std::cout << "VCF Allelic Depth Filter\n";
    std::cout << "========================\n";
    std::cout << "Configuration:\n";
    std::cout << "  Input VCF: " << vcf_file << "\n";
    std::cout << "  Sample list: " << sample_list << "\n";
    std::cout << "  Output VCF: " << output_file << "\n";
    std::cout << "  Log file: " << log_file << "\n";
    std::cout << "  VCF Type: " << (params.vcf_type == VCFType::VARSCAN ? "VarScan2" : "Standard") << "\n";
    std::cout << "  Homozygous purity: " << params.homozygous_purity << "%\n";
    std::cout << "  Heterozygous range: " << params.het_min_ratio << "% - " 
              << params.het_max_ratio << "%\n\n";
    
    // Create filter and run
    VCFFilter filter(vcf_file, sample_list, output_file, log_file, params);
    
    if (!filter.filterVCF()) {
        std::cerr << "Error: Filtering failed." << std::endl;
        return 1;
    }
    
    std::cout << "\nFiltering completed successfully!" << std::endl;
    return 0;
}
