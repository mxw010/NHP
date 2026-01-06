#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

struct SequenceAnalysis {
    std::string upstream8bp;
    std::string downstream8bp;
    int nnkLength;
    size_t nnkStartPos;
    size_t nnkEndPos;
    
    void print(const std::string& label) const {
        std::cout << "\n=== " << label << " ===" << std::endl;
        std::cout << "NNK start position: " << nnkStartPos << std::endl;
        std::cout << "NNK end position: " << nnkEndPos << std::endl;
        std::cout << "NNK length: " << nnkLength << " bp (" << nnkLength/3 << " codons)" << std::endl;
        std::cout << "Upstream: " << upstream8bp << std::endl;
        std::cout << "Downstream: " << downstream8bp << std::endl;
    }
};

class NNKAnalyzer {
public:
    static std::string toUpperCase(const std::string& seq) {
        std::string result = seq;
        std::transform(result.begin(), result.end(), result.begin(), ::toupper);
        return result;
    }
    
    static char complement(char base) {
        switch(std::toupper(base)) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            case 'N': return 'N';
            case 'K': return 'M';
            case 'M': return 'K';
            default: return base;
        }
    }
    
    static std::string reverseComplement(const std::string& seq) {
        std::string rc;
        rc.reserve(seq.length());
        
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            rc += complement(*it);
        }
        
        return rc;
    }
    
    static SequenceAnalysis analyzeSequence(const std::string& seq, int flankingLength = 8) {
        SequenceAnalysis result;
        
        // Find NNK repeat region
        size_t nnkStart = seq.find("NNK");
        if (nnkStart == std::string::npos) {
            result.nnkLength = 0;
            result.nnkStartPos = 0;
            result.nnkEndPos = 0;
            return result;
        }
        
        result.nnkStartPos = nnkStart;
        
        // Calculate length of NNK repeats
        size_t pos = nnkStart;
        int count = 0;
        
        while (pos + 2 < seq.length() && 
               seq[pos] == 'N' && 
               seq[pos + 1] == 'N' && 
               seq[pos + 2] == 'K') {
            count++;
            pos += 3;
        }
        
        result.nnkLength = count * 3;
        result.nnkEndPos = nnkStart + result.nnkLength - 1;
        
        // Extract upstream flanking sequence (closest to NNK start)
        if (nnkStart >= static_cast<size_t>(flankingLength)) {
            result.upstream8bp = seq.substr(nnkStart - flankingLength, flankingLength);
            result.upstream8bp = toUpperCase(result.upstream8bp);
        } else if (nnkStart > 0) {
            result.upstream8bp = seq.substr(0, nnkStart);
            result.upstream8bp = toUpperCase(result.upstream8bp);
        } else {
            result.upstream8bp = "";
        }
        
        // Extract downstream flanking sequence (closest to NNK end)
        size_t downstreamStart = nnkStart + result.nnkLength;
        if (downstreamStart + flankingLength <= seq.length()) {
            result.downstream8bp = seq.substr(downstreamStart, flankingLength);
            result.downstream8bp = toUpperCase(result.downstream8bp);
        } else if (downstreamStart < seq.length()) {
            result.downstream8bp = seq.substr(downstreamStart);
            result.downstream8bp = toUpperCase(result.downstream8bp);
        } else {
            result.downstream8bp = "";
        }
        
        return result;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <sequence> [flanking_length]" << std::endl;
        std::cerr << "  sequence: DNA sequence containing NNK repeats" << std::endl;
        std::cerr << "  flanking_length: number of base pairs to extract (default: 8)" << std::endl;
        std::cerr << "\nExample: " << argv[0] << " \"gtcctatggacaagtggccNNKNNKNNKgcacaggc\" 10" << std::endl;
        return 1;
    }
    
    std::string sequence = argv[1];
    int flankingLength = 8; // default
    
    if (argc == 3) {
        try {
            flankingLength = std::stoi(argv[2]);
            if (flankingLength <= 0) {
                std::cerr << "Error: flanking_length must be a positive integer" << std::endl;
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid flanking_length. Must be an integer." << std::endl;
            return 1;
        }
    }
    
    // Remove underscores if present
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '_'), sequence.end());
    
    std::cout << "Original sequence:" << std::endl;
    std::cout << sequence << std::endl;
    std::cout << "Flanking length: " << flankingLength << " bp" << std::endl;
    
    // Analyze forward sequence
    SequenceAnalysis forward = NNKAnalyzer::analyzeSequence(sequence, flankingLength);
    forward.print("Forward Sequence Analysis");
    
    // Create reverse complement analysis by reverse complementing the flanking sequences
    SequenceAnalysis reverse;
    reverse.nnkLength = forward.nnkLength;
    reverse.nnkStartPos = forward.nnkStartPos;
    reverse.nnkEndPos = forward.nnkEndPos;
    reverse.upstream8bp = NNKAnalyzer::reverseComplement(forward.upstream8bp);
    reverse.downstream8bp = NNKAnalyzer::reverseComplement(forward.downstream8bp);
    
    reverse.print("Reverse Complement of Flanking Sequences");
    
    // std::cout << "\nNote: The upstream RC is the reverse complement of the forward upstream." << std::endl;
    // std::cout << "      The downstream RC is the reverse complement of the forward downstream." << std::endl;
    
    return 0;
}