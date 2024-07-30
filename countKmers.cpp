#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cctype>

void check_args(int argc, char* argv[]) {
    if (argc != 3) {
        throw std::runtime_error(std::string("Usage: ") + argv[0] + " <input_file> <output_file>");
    }
    std::ifstream infile(argv[1]);
    if (!infile) {
        throw std::runtime_error(std::string("Error: File '") + argv[1] + "' not found");
    }
    infile.seekg(0, std::ios::end);
    if (infile.tellg() == 0) {
        throw std::runtime_error(std::string("Error: File '") + argv[1] + "' is empty");
    }
}

std::vector<std::string> read_file(const std::string& infile) {
    std::ifstream f(infile.c_str());
    if (!f) {
        throw std::runtime_error("Error: Unable to open input file");
    }
    std::vector<std::string> data;
    std::string line;
    while (std::getline(f, line)) {
        data.push_back(line);
    }
    if (f.bad()) {
        throw std::runtime_error("Error: An error occurred while reading the input file");
    }
    return data;
}

bool is_valid_dna_sequence(const std::string& sequence) {
    for (size_t i = 0; i < sequence.size(); ++i) {
        char c = sequence[i];
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            return false;
        }
    }
    return true;
}

std::map<std::string, int> count_kmers(const std::vector<std::string>& data) {
    std::map<std::string, int> kmer_dict;
    for (size_t i = 0; i < data.size(); ++i) {
        const std::string& line = data[i];
        if (line.empty() || line[0] == '>') {
            continue;
        }
        if (!is_valid_dna_sequence(line)) {
            std::cerr << "Warning: Invalid DNA sequence found: " << line << std::endl;
            continue;
        }
        if (line.size() < 4) {
            std::cerr << "Warning: Line too short to contain any k-mers: " << line << std::endl;
            continue;
        }
        for (size_t j = 0; j <= line.size() - 4; ++j) {
            std::string kmer = line.substr(j, 4);
            kmer_dict[kmer]++;
        }
    }
    return kmer_dict;
}

bool compare_kmer_counts(const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
    return a.second > b.second || (a.second == b.second && a.first < b.first);
}

std::vector<std::pair<std::string, int> > sort_kmers(const std::map<std::string, int>& kmer_dict) {
    std::vector<std::pair<std::string, int> > sorted_kmers(kmer_dict.begin(), kmer_dict.end());
    std::sort(sorted_kmers.begin(), sorted_kmers.end(), compare_kmer_counts);
    return sorted_kmers;
}

void write_output(const std::vector<std::pair<std::string, int> >& sorted_kmers, const std::string& outfile) {
    std::ofstream f(outfile.c_str());
    if (!f) {
        throw std::runtime_error("Error: Unable to open output file");
    }
    for (size_t i = 0; i < sorted_kmers.size(); ++i) {
        const std::pair<std::string, int>& kmer_count_pair = sorted_kmers[i];
        f << kmer_count_pair.first << "\t" << kmer_count_pair.second << "\n";
    }
    if (f.bad()) {
        throw std::runtime_error("Error: An error occurred while writing to the output file");
    }
}

int main(int argc, char* argv[]) {
    try {
        check_args(argc, argv);
        std::vector<std::string> data = read_file(argv[1]);
        std::map<std::string, int> kmer_dict = count_kmers(data);
        std::vector<std::pair<std::string, int> > sorted_kmers = sort_kmers(kmer_dict);
        write_output(sorted_kmers, argv[2]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
