#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <zlib.h>
#include <cctype>
#include <map>

struct FastqRead {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
    std::vector<int> qualityScores;
};

struct SequenceAnalysis {
    std::string upstream8bp;
    std::string downstream8bp;
    int nnkLength;
    size_t nnkStartPos;
    size_t nnkEndPos;
};

struct FlankingSequences {
    std::string upstream;
    std::string downstream;
    std::string upstreamRC;
    std::string downstreamRC;
    int nnkLength;
};

// NNK Analyzer functionality
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
        
        size_t nnkStart = seq.find("NNK");
        if (nnkStart == std::string::npos) {
            result.nnkLength = 0;
            result.nnkStartPos = 0;
            result.nnkEndPos = 0;
            return result;
        }
        
        result.nnkStartPos = nnkStart;
        
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
        
        if (nnkStart >= static_cast<size_t>(flankingLength)) {
            result.upstream8bp = seq.substr(nnkStart - flankingLength, flankingLength);
            result.upstream8bp = toUpperCase(result.upstream8bp);
        } else if (nnkStart > 0) {
            result.upstream8bp = seq.substr(0, nnkStart);
            result.upstream8bp = toUpperCase(result.upstream8bp);
        } else {
            result.upstream8bp = "";
        }
        
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

// FASTQ processing
class FastqExtractor {
private:
    gzFile file;
    static const int BUFFER_SIZE = 65536;
    int phredOffset;
    
public:
    FastqExtractor(const std::string& filename, int offset = 33) : phredOffset(offset) {
        file = gzopen(filename.c_str(), "rb");
        if (!file) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
    }
    
    ~FastqExtractor() {
        if (file) {
            gzclose(file);
        }
    }
    
    bool getNextRead(FastqRead& read) {
        char buffer[BUFFER_SIZE];
        
        if (gzgets(file, buffer, BUFFER_SIZE) == nullptr) return false;
        read.header = std::string(buffer);
        if (read.header.back() == '\n') read.header.pop_back();
        
        if (gzgets(file, buffer, BUFFER_SIZE) == nullptr) return false;
        read.sequence = std::string(buffer);
        if (read.sequence.back() == '\n') read.sequence.pop_back();
        
        if (gzgets(file, buffer, BUFFER_SIZE) == nullptr) return false;
        read.plus = std::string(buffer);
        if (read.plus.back() == '\n') read.plus.pop_back();
        
        if (gzgets(file, buffer, BUFFER_SIZE) == nullptr) return false;
        read.quality = std::string(buffer);
        if (read.quality.back() == '\n') read.quality.pop_back();
        
        read.qualityScores = convertQualityToScores(read.quality);
        
        return true;
    }
    
    std::vector<int> convertQualityToScores(const std::string& qualityString) {
        std::vector<int> scores;
        scores.reserve(qualityString.length());
        
        for (char c : qualityString) {
            scores.push_back(static_cast<int>(c) - phredOffset);
        }
        
        return scores;
    }
    
    double getAverageQuality(const std::vector<int>& scores) {
        if (scores.empty()) return 0.0;
        double sum = 0.0;
        for (int score : scores) sum += score;
        return sum / scores.size();
    }
};

class SequenceProcessor {
public:
    static char translateCodon(const std::string& codon) {
        // Standard genetic code
        static const std::map<std::string, char> codonTable = {
            {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
            {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
            {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
            {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
            {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
            {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
            {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
            {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
            {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
            {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
            {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
            {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
            {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
            {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
            {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
            {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
        };
        
        auto it = codonTable.find(codon);
        if (it != codonTable.end()) {
            return it->second;
        }
        return 'X'; // Unknown amino acid
    }
    
    static std::string translateToPeptide(const std::string& dnaSeq) {
        if (dnaSeq.length() % 3 != 0) {
            return ""; // Invalid sequence length
        }
        
        std::string peptide;
        for (size_t i = 0; i < dnaSeq.length(); i += 3) {
            std::string codon = dnaSeq.substr(i, 3);
            char aa = translateCodon(codon);
            peptide += aa;
            // Don't break - keep translating even after stop codons
        }
        
        return peptide;
    }
    
    static bool isValidNNKSequence(const std::string& seq) {
        if (seq.length() % 3 != 0) return false;
        
        for (size_t i = 0; i < seq.length(); i += 3) {
            if (i + 2 >= seq.length()) return false;
            
            char k = std::toupper(seq[i + 2]);
            if (k != 'G' && k != 'T') {
                return false;
            }
        }
        return true;
    }

    static bool isSTOP(const std::string& pep)
    {
        return pep.find('*') != std::string::npos;
    }
    
    static size_t findSequence(const std::string& text, const std::string& pattern) {
        return text.find(pattern);
    }
    
    static std::pair<std::string, size_t> extractBetweenFlanking(const std::string& seq, 
                                              const std::string& upstream,
                                              const std::string& downstream,
                                              int expectedLength) {
        size_t upstreamPos = seq.find(upstream);
        if (upstreamPos == std::string::npos) return {"", 0};
        
        size_t extractStart = upstreamPos + upstream.length();
        size_t downstreamPos = seq.substr(extractStart).find(downstream);
        
        if (downstreamPos == std::string::npos) return {"", 0};
        
        std::string extracted = seq.substr(extractStart, downstreamPos);
        
        if (static_cast<int>(extracted.length()) != expectedLength) return {"", 0};
        
        return {extracted, extractStart};
    }
};

struct ExtractedSequence {
    std::string readHeader;
    std::string extractedSeq;
    std::string peptideSeq;
    std::string orientation;
    double avgQuality;
    size_t startPosition;
    bool isValidNNK;
    bool isStopCodon;
};

void processAndExtract(const std::string& fastqFile, 
                       const FlankingSequences& flanking,
                       double minQuality,
                       const std::string& outputFile) {
    FastqExtractor extractor(fastqFile);
    std::vector<ExtractedSequence> validSequences;
    
    FastqRead read;
    int totalReads = 0;
    int passedQuality = 0;
    int matchedFlanking = 0;
    int validNNK = 0;
    int invalidNNK = 0;
    int normalPeptides = 0;
    int StopPeptides = 0;
    
    while (extractor.getNextRead(read)) {
        totalReads++;
        
        double avgQual = extractor.getAverageQuality(read.qualityScores);
        if (avgQual < minQuality) continue;
        passedQuality++;
        
        std::string upperSeq = NNKAnalyzer::toUpperCase(read.sequence);
        std::string extracted;
        std::string orientation;
        size_t startPos = 0;
        
        std::pair<std::string, size_t> forwardResult = SequenceProcessor::extractBetweenFlanking(
            upperSeq, flanking.upstream, flanking.downstream, flanking.nnkLength);
        
        if (!forwardResult.first.empty()) {
            extracted = forwardResult.first;
            startPos = forwardResult.second;
            orientation = "forward";
        } else {
            // For reverse orientation, downstream RC comes first, then upstream RC
            std::pair<std::string, size_t> reverseResult = SequenceProcessor::extractBetweenFlanking(
                upperSeq, flanking.downstreamRC, flanking.upstreamRC, flanking.nnkLength);
            
            if (!reverseResult.first.empty()) {
                extracted = NNKAnalyzer::reverseComplement(reverseResult.first);
                startPos = reverseResult.second;
                orientation = "reverse";
            }
        }
        
        if (extracted.empty()) continue;
        matchedFlanking++;
        
        // counting valid and invalid nnks
        if (!SequenceProcessor::isValidNNKSequence(extracted)) {
            invalidNNK++;
        } else {
            validNNK++;
        }
        
        // Translate to peptide
        std::string peptide = SequenceProcessor::translateToPeptide(extracted);

        // counting STOP codons
        if (!SequenceProcessor::isSTOP(peptide)) {
            normalPeptides++;
        } else {
            StopPeptides++;
        }

        ExtractedSequence result;
        result.readHeader = read.header;
        result.extractedSeq = extracted;
        result.peptideSeq = peptide;
        result.orientation = orientation;
        result.avgQuality = avgQual;
        result.startPosition = startPos;
        result.isValidNNK = SequenceProcessor::isValidNNKSequence(extracted);
        result.isStopCodon = SequenceProcessor::isSTOP(peptide);

        validSequences.push_back(result);
    }
    
    std::ofstream outFile(outputFile);
    if (!outFile) {
        throw std::runtime_error("Failed to open output file: " + outputFile);
    }
    
    outFile << "# Total reads: " << totalReads << std::endl;
    outFile << "# Passed quality (>" << minQuality << "): " << passedQuality << std::endl;
    outFile << "# Matched flanking sequences: " << matchedFlanking << std::endl;
    outFile << "# Valid NNK sequences (K=G/T): " << validNNK << std::endl;
    outFile << "# Invalid NNK sequences (K!=G/T): " << invalidNNK << std::endl;
    outFile << "# Final extracted sequences: " << validSequences.size() << std::endl;
    outFile << std::endl;
    outFile << "ReadHeader\tStartPosition\tDNASequence\tPeptideSequence\tOrientation\tAvgQuality\tValidNNK\tStopCodon" << std::endl;
    
    for (const auto& seq : validSequences) {
        outFile << seq.readHeader << "\t" 
                << seq.startPosition << "\t"
                << seq.extractedSeq << "\t" 
                << seq.peptideSeq << "\t"
                << seq.orientation << "\t" 
                << seq.avgQuality << "\t"
                << (seq.isValidNNK ? "TRUE" : "FALSE") << "\t"
                << (seq.isStopCodon ? "TRUE" : "FALSE") <<std::endl;
    }
    
    outFile.close();
    
    std::cout << R"("total_reads": ")" << totalReads << R"(",)" << std::endl;
    std::cout << R"("passed_qc": ")" << passedQuality << R"(",)" << std::endl;
    std::cout << R"("matched_flanking": ")" << matchedFlanking << R"(",)" << std::endl;
    std::cout << R"("extracted": ")" << validSequences.size() << R"(",)" << std::endl;
    std::cout << R"("valid_NNK": ")" << validNNK << R"(",)" << std::endl;
    std::cout << R"("invalid_NNK": ")" << invalidNNK << R"(",)" << std::endl;
    std::cout << R"("valid_codons": ")" << normalPeptides << R"(",)" << std::endl;
    std::cout << R"("stop_codons": ")" << StopPeptides << R"(",)" << std::endl;
    std::cout << R"("output": ")" << outputFile << R"(")" << std::endl;
    std::cout << "}" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " <reference_sequence> <fastq.gz> <min_quality> <output.txt> [flanking_length]" << std::endl;
        std::cerr << "\nParameters:" << std::endl;
        std::cerr << "  reference_sequence : Reference sequence containing NNK repeats" << std::endl;
        std::cerr << "  fastq.gz          : Input FASTQ file (gzipped)" << std::endl;
        std::cerr << "  min_quality       : Minimum average quality score (e.g., 30)" << std::endl;
        std::cerr << "  output.txt        : Output file for extracted sequences" << std::endl;
        std::cerr << "  flanking_length   : Optional - bp to extract as flanking (default: 8)" << std::endl;
        std::cerr << "\nExample:" << std::endl;
        std::cerr << "  " << argv[0] << " \"gtcctatggacaagtggccNNKNNKNNKgcacaggc\" reads.fastq.gz 30 extracted.txt 8" << std::endl;
        return 1;
    }
    
    try {
        std::string refSequence = argv[1];
        std::string fastqFile = argv[2];
        double minQuality = std::stod(argv[3]);
        std::string outputFile = argv[4];
        int flankingLength = 8;
        
        if (argc == 6) {
            flankingLength = std::stoi(argv[5]);
            if (flankingLength <= 0) {
                std::cerr << "Error: flanking_length must be a positive integer" << std::endl;
                return 1;
            }
        }
        
        // Remove underscores from reference sequence
        refSequence.erase(std::remove(refSequence.begin(), refSequence.end(), '_'), refSequence.end());
        
        // Analyze reference sequence to get flanking sequences
        SequenceAnalysis analysis = NNKAnalyzer::analyzeSequence(refSequence, flankingLength);
        
        if (analysis.nnkLength == 0) {
            std::cerr << "Error: No NNK repeats found in reference sequence" << std::endl;
            return 1;
        }
        
        // Prepare flanking sequences
        FlankingSequences flanking;
        flanking.upstream = analysis.upstream8bp;
        flanking.downstream = analysis.downstream8bp;
        flanking.upstreamRC = NNKAnalyzer::reverseComplement(analysis.upstream8bp);
        flanking.downstreamRC = NNKAnalyzer::reverseComplement(analysis.downstream8bp);
        flanking.nnkLength = analysis.nnkLength;
        
        std::cout << "{" << std::endl;
        std::cout << R"("reference": ")" << refSequence << R"(",)" << std::endl;
        std::cout << R"("fastq": ")" << fastqFile << R"(",)" << std::endl;
        std::cout << R"("min_quality": ")" << std::to_string(minQuality) << R"(",)" << std::endl;
        std::cout << R"("start_position": ")" << analysis.nnkStartPos << R"(",)" << std::endl;
        std::cout << R"("NNK_length": ")" << analysis.nnkLength << R"(",)" << std::endl;
        std::cout << R"("flanking": ")" << flanking.upstream << R"(",)" << std::endl;
        std::cout << R"("upstream_rc": ")" << flanking.upstreamRC << R"(",)" << std::endl;
        std::cout << R"("downstream_rc": ")" << flanking.downstreamRC << R"(",)" << std::endl;

        processAndExtract(fastqFile, flanking, minQuality, outputFile);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}