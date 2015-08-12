//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_consensus.cpp -- entry point to consensus functions
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <queue>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "nanopolish_math.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_khmm_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "nanopolish_variants.h"
#include "nanopolish_haplotype.h"
#include "profiler.h"
#include "progress.h"
#include "stdaln.h"

// Macros
#define max3(x,y,z) std::max(std::max(x,y), z)

// Flags to turn on/off debugging information

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1
//#define DEBUG_TRANSITION 1
//#define DEBUG_PATH_SELECTION 1
//#define DEBUG_SINGLE_SEGMENT 1
//#define DEBUG_SHOW_TOP_TWO 1
//#define DEBUG_SEGMENT_ID 193
//#define DEBUG_BENCHMARK 1

//
// Getopt
//
#define SUBPROGRAM "consensus"

static const char *CONSENSUS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CONSENSUS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Compute a new consensus sequence for an assembly using a signal-level HMM\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -m, --min-read-evidence=N            require at least N reads to have a variant to try calling it\n"
"      --snps                           only call SNPs\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -e, --event-bam=FILE                 the events aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -o, --outfile=FILE                   write result to FILE [default: stdout]\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --vcf=FILE                       write called variants to vcf FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string event_bam_file;
    static std::string genome_file;
    static std::string output_file;
    static std::string output_vcf;
    static std::string window;
    static int snps_only = 0;
    static int min_read_evidence = 1;
    static int show_progress = 0;
    static int num_threads = 1;
}

static const char* shortopts = "r:b:g:t:w:o:e:m:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VCF, OPT_PROGRESS, OPT_SNPS_ONLY };

static const struct option longopts[] = {
    { "verbose",           no_argument,       NULL, 'v' },
    { "reads",             required_argument, NULL, 'r' },
    { "bam",               required_argument, NULL, 'b' },
    { "event-bam",         required_argument, NULL, 'e' },
    { "genome",            required_argument, NULL, 'g' },
    { "window",            required_argument, NULL, 'w' },
    { "outfile",           required_argument, NULL, 'o' },
    { "threads",           required_argument, NULL, 't' },
    { "min-read-evidence", required_argument, NULL, 'm' },
    { "vcf",               required_argument, NULL, OPT_VCF },
    { "snps",              no_argument,       NULL, OPT_SNPS_ONLY },
    { "progress",          no_argument,       NULL, OPT_PROGRESS },
    { "help",              no_argument,       NULL, OPT_HELP },
    { "version",           no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Handy wrappers for scoring/debugging functions
// The consensus algorithms call into these so we can switch
// scoring functions without writing a bunch of code
double score_sequence(const std::string& sequence, const HMMInputData& data)
{
    //return score_skip_merge(sequence, state);
    //return score_khmm_model_postmerge(sequence, state);
    //return khmm_score(sequence, state, AP_GLOBAL);
    return profile_hmm_score(sequence, data);
    //return score_emission_dp(sequence, state);
}

int get_contig_length(const std::string& contig)
{
    faidx_t *fai = fai_load(opt::genome_file.c_str());
    int len = faidx_seq_len(fai, contig.c_str());
    fai_destroy(fai);
    return len;
}

std::vector<AlignmentState> hmm_align(const std::string& sequence, const HMMInputData& data)
{
    return profile_hmm_align(sequence, data);
//    return khmm_posterior_decode(sequence, state);
}

void debug_sequence(const std::string& name, uint32_t seq_id, uint32_t read_id, const std::string& sequence, const HMMInputData& data)
{
    std::vector<AlignmentState> alignment = hmm_align(sequence, data);
    print_alignment(name, seq_id, read_id, sequence, data, alignment);
}

void update_training_with_segment(const std::string& sequence, const HMMInputData& data)
{
    profile_hmm_update_training(sequence, data);
}

//
// Outlier filtering
//
void filter_outlier_data(std::vector<HMMInputData>& input, const std::string& sequence)
{
    std::vector<HMMInputData> out_rs;
    for(uint32_t ri = 0; ri < input.size(); ++ri) {
        const HMMInputData& rs = input[ri];

        double curr = score_sequence(sequence, rs);
        double n_events = abs(rs.event_start_idx - rs.event_stop_idx) + 1.0f;
        double lp_per_event = curr / n_events;

        if(opt::verbose >= 3) {
            fprintf(stderr, "OUTLIER_FILTER %d %.2lf %.2lf %.2lf\n", ri, curr, n_events, lp_per_event);
        }

        if(fabs(lp_per_event) < 3.5f) {
            out_rs.push_back(rs);
        }
    }
    input.swap(out_rs);
}

void filter_outlier_data_by_kmers(std::vector<HMMInputData>& input, const std::string& sequence)
{
    std::vector<HMMInputData> out_rs;
    for(uint32_t ri = 0; ri < input.size(); ++ri) {
        const HMMInputData& rs = input[ri];

        double curr = score_sequence(sequence, rs);
        double n_kmers = sequence.size() - K + 1;
        double lp_per_kmer = curr / n_kmers;

        if(opt::verbose >= 3) {
            fprintf(stderr, "OUTLIER_FILTER %d %.2lf %.2lf %.2lf\n", ri, curr, n_kmers, lp_per_kmer);
        }

        if(fabs(lp_per_kmer) < 3.5f) {
            out_rs.push_back(rs);
        }
    }
    input.swap(out_rs);
}


// update the training data on the current segment
void train_segment(HMMRealignmentInput& window, uint32_t segment_id)
{
#if 0
    // Get the segments
    assert(segment_id + 2 < window.anchored_columns.size());
    HMMAnchoredColumn& start_column = window.anchored_columns[segment_id];
    HMMAnchoredColumn& middle_column = window.anchored_columns[segment_id + 1];
    HMMAnchoredColumn& end_column = window.anchored_columns[segment_id + 2];

    std::string s_m_base = start_column.base_sequence;
    std::string m_e_base = middle_column.base_sequence;

    std::string segment_sequence = join_sequences_at_kmer(s_m_base, m_e_base);

    // Set up the the input data for the HMM
    std::vector<HMMInputData> input = get_input_for_columns(window, start_column, end_column);
     
    for(uint32_t ri = 0; ri < input.size(); ++ri) {
        std::vector<AlignmentState> decodes = hmm_align(segment_sequence, input[ri]);
        update_training_with_segment(segment_sequence, input[ri]);
    }
#endif
}

void train(HMMRealignmentInput& window)
{
    // train on current consensus
    uint32_t num_segments = window.anchored_columns.size();
    for(uint32_t segment_id = 0; segment_id < num_segments - 2; ++segment_id) {
        train_segment(window, segment_id);
    }

    // Update model parameters
    for(uint32_t ri = 0; ri < window.reads.size(); ++ri) {
        window.reads[ri]->parameters[0].train();
        window.reads[ri]->parameters[1].train();
    }
}

std::vector<Variant> generate_all_snps(const std::string& reference)
{
    std::vector<Variant> out;
    for(size_t i = 0; i < reference.size(); ++i) {
        for(size_t bi = 0; bi < 4; ++bi) {
            char b = "ACGT"[bi];
            if(reference[i] != b) {
                Variant v;
                v.ref_name = "noctg";
                v.ref_position = i;
                v.ref_seq = reference[i];
                v.alt_seq = b;
                out.push_back(v);
            }
        }
    }
    return out;
}

std::vector<Variant> generate_variants_from_kLCS(const std::string& reference, const std::string& sequence, const int kmer)
{
    std::vector<Variant> out;
    kLCSResult result = kLCS(reference, sequence, kmer);

    uint32_t match_idx = 0;
    uint32_t last_idx = result.size() - 1;
    while(match_idx < result.size()) {

        // advance the match to the next point of divergence
        while(match_idx != last_idx && 
              result[match_idx].i == result[match_idx + 1].i - 1 &&
              result[match_idx].j == result[match_idx + 1].j - 1) {
            match_idx++;
        }

        // no more divergences to process
        if(match_idx == last_idx)
            break;

        uint32_t bl = result[match_idx + 1].i - result[match_idx].i + kmer;
        uint32_t rl = result[match_idx + 1].j - result[match_idx].j + kmer;

        std::string ref_subseq = reference.substr(result[match_idx].i, bl);
        std::string read_subseq = sequence.substr(result[match_idx].j, rl);
        
        // substrings must share a k-mer at the beginning/end
        assert(ref_subseq.substr(0, kmer) == read_subseq.substr(0, kmer));    
        assert(ref_subseq.substr(bl - kmer) == read_subseq.substr(rl - kmer));    
        assert(ref_subseq != read_subseq);

        // Find the left boundary of the difference
        int ref_s = 0;
        int read_s = 0;
        while(ref_subseq[ref_s] == read_subseq[read_s]) {
            ref_s += 1;
            read_s += 1;
        }

        // if its an indel or complex variant, include one matching base
        if(ref_subseq.length() != read_subseq.length() || ref_s != read_s) {
            ref_s -= 1;
            read_s -= 1;
        }
        
        std::string tmp_ref = ref_subseq.substr(ref_s);
        std::string tmp_read = read_subseq.substr(read_s);

        // trim unnecessary bases from the end
        while(tmp_ref.size() > 1 && tmp_read.size() > 1 && 
              tmp_ref.back() == tmp_read.back()) 
        {
            tmp_ref.pop_back();
            tmp_read.pop_back();
            assert(!tmp_ref.empty());
            assert(!tmp_read.empty());
        }

        Variant v;
        v.ref_name = "noctg";
        v.ref_position = ref_s + result[match_idx].i;
        v.ref_seq = tmp_ref;
        v.alt_seq = tmp_read;
        v.quality = 0.0f;
        out.push_back(v);

        match_idx += 1;
    }
    return out;
}

std::vector<Variant> generate_variants_from_reads(const std::string& reference, const std::vector<std::string>& reads)
{
    std::vector<Variant> out;

    // Generate alternatives
    for(uint32_t ri = 0; ri < reads.size(); ++ri) {

        if(reads[ri].size() < K)
            continue;

        std::vector<Variant> read_variants = generate_variants_from_kLCS(reference, reads[ri], K);
        out.insert(out.end(), read_variants.begin(), read_variants.end());
    }
    return out;
}

double parallel_score(const std::string& sequence, const std::vector<HMMInputData>& event_sequences)
{
    double score = 0.0f;

    #pragma omp parallel for
    for(size_t i = 0; i < event_sequences.size(); ++i) {
        double s = profile_hmm_score(sequence, event_sequences[i]);

        #pragma omp critical
        score += s;
    }
    return score;
}

struct BranchSequence
{
    std::string sequence;
    double score;

    static bool sortByScore(const BranchSequence& a, const BranchSequence& b) { return a.score > b.score; }
};

std::vector<std::string> search_for_extensions(const std::string& root, 
                                               const std::string& ref_subseq, 
                                               const std::vector<HMMInputData>& event_sequences,
                                               const size_t required_match_length,
                                               const size_t max_extension)
{
    std::vector<std::string> candidate_sequences;
    BranchSequence root_elem = { root, -INFINITY };
    std::vector<BranchSequence> branches(1, root_elem);

    size_t round = 0;

    while(round++ < max_extension) {
        
        if(opt::verbose > 2) {
            fprintf(stderr, "\n==== Round %zu ====\n", round);
            fprintf(stderr, "ref: %s\n", ref_subseq.c_str());
        }

        std::vector<BranchSequence> incoming;
        for(size_t branch_idx = 0; branch_idx < branches.size(); ++branch_idx) {

            std::vector<std::string> extension_set = generate_mers(1);

            for(size_t ext_idx = 0; ext_idx < extension_set.size(); ++ext_idx) {
                std::string extended = branches[branch_idx].sequence + extension_set[ext_idx];

                double score = parallel_score(extended, event_sequences);
                BranchSequence new_branch = { extended, score };
                incoming.push_back(new_branch);
            }
        }
        
        // sort by score
        std::sort(incoming.begin(), incoming.end(), BranchSequence::sortByScore);
        branches.clear();
        
        // cull bad branches
        std::vector<double> score_with_best;
        for(size_t branch_idx = 0; branch_idx < incoming.size(); ++branch_idx) {
            BranchSequence& branch = incoming[branch_idx];
            double relative_score = branch.score - incoming[0].score;
            bool is_ref = ref_subseq.find(branch.sequence) != std::string::npos;
            bool selected = false;

            if( branch_idx < 16 || is_ref) {
                selected = true;
            }
         
            bool joined = false;
            // Check if this branch has rejoined the reference
            if(selected) {
                if(is_ref) {
                    // never check the reference for joining
                    branches.push_back(branch);
                } else {
                    std::string branch_suffix = branch.sequence.substr(branch.sequence.size() - required_match_length);
                    if(ref_subseq.rfind(branch_suffix) != std::string::npos) {
                        // branch has joined reference
                        joined = true;
                        candidate_sequences.push_back(branch.sequence);
                    } else {
                        // continue extending
                        branches.push_back(branch);
                    }
                }
            }

            if(opt::verbose > 2) {
                // Debug print
                std::string status = is_ref ? "*" : " ";
                status += selected ? 's' : ' ';
                status += joined ? 'j' : ' ';
                
                fprintf(stderr, "seq: %s %5.2lf %s\n", branch.sequence.c_str(), branch_idx == 0 ? branch.score : relative_score, status.c_str());
                
            }
        }
    }
    return candidate_sequences;
}

// Branch and bound algorithm for calling variants
Haplotype call_variants_for_region_bb(const std::string& contig, int region_start, int region_end)
{
    const int BUFFER = 10;
    const int MAX_EXTEND = 20;

    if(region_start < BUFFER)
        region_start = BUFFER;
    region_end = std::min(region_end, get_contig_length(contig) - BUFFER - 1);
    
    size_t buffered_region_start = region_start - BUFFER;
    size_t buffered_region_end = region_end + BUFFER;

    // load the region, accounting for the buffering
    AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file);
    alignments.load_region(contig, buffered_region_start, buffered_region_end);
    Haplotype derived_haplotype(contig,
                                alignments.get_region_start(),
                                alignments.get_reference());
    
    // What base are we currently calling from
    size_t curr_ref_start = buffered_region_start;

    // Get the reference string over the whole region
    std::string ref_string = alignments.get_reference_substring(contig, buffered_region_start, buffered_region_end);

    // Extend the consensus holding the invariant that the last N bp exactly match the reference sequence
    while(curr_ref_start < region_end) {

        // Call variants over a subregion of size 60
        size_t curr_ref_end = curr_ref_start + 60;
        curr_ref_end = std::min(curr_ref_end, buffered_region_end);

        // Get the event sequences within this window
        std::string ref_subseq = alignments.get_reference_substring(contig, curr_ref_start, curr_ref_end);
        std::vector<HMMInputData> event_sequences = alignments.get_event_subsequences(contig, curr_ref_start, curr_ref_end);
        
        //filter_outlier_data_by_kmers(event_sequences, ref_subseq);

        // Start the extension from the last BUFFER bases of the called sequence
        std::string root = ref_string.substr(curr_ref_start - buffered_region_start, BUFFER);
        if(opt::verbose > 2) {
            fprintf(stderr, "calling region start: %zu brs: %zu root: %s\n", 
                curr_ref_start, 
                buffered_region_start,
                root.c_str());
        }

        // screen out trivial reference extensions quickly
        std::vector<std::string> extension_set = generate_mers(3);
        char ref_base = ref_string[curr_ref_start - buffered_region_start + BUFFER];
        double best_score = -INFINITY;
        std::string best_ext;

        for(size_t ei = 0; ei < extension_set.size(); ei++) {
            std::string extended = root + extension_set[ei];
            double score = parallel_score(extended, event_sequences);

            if(score > best_score) {
                best_score = score;
                best_ext = extension_set[ei];
            }
        }

        /*
        printf("Ref:  %s\n", ref_string.substr(curr_ref_start - buffered_region_start).c_str());
        printf("Root: %s\n", root.c_str());
        printf("[%d] Best: %.2lf base: %s ref: %c delta: %.2lf\n", curr_ref_start + BUFFER, best_score, best_ext.c_str(), ref_base, best_score - 0);
        */
        if(best_ext.front() == ref_base) {
            curr_ref_start += 1;
            continue;
        }

        std::vector<std::string> candidate_sequences = search_for_extensions(root, ref_subseq, event_sequences, BUFFER, MAX_EXTEND);

        // parse variants from the consensus sequences
        std::vector<Variant> candidate_variants;

        for(size_t cs_idx = 0; cs_idx < candidate_sequences.size(); ++cs_idx) {
            std::vector<Variant> variants = 
                generate_variants_from_kLCS(ref_subseq, candidate_sequences[cs_idx], BUFFER);
            candidate_variants.insert(candidate_variants.end(), variants.begin(), variants.end());
        }

        // translate variants into reference coordinates
        for(size_t vi = 0; vi < candidate_variants.size(); ++vi) {
            candidate_variants[vi].ref_name = contig;
            candidate_variants[vi].ref_position += curr_ref_start;
            
            if(opt::verbose > 1) {
                candidate_variants[vi].write_vcf(stderr);
            }
        }

        if(opt::snps_only) {
            filter_out_non_snp_variants(candidate_variants);
        }

        int max_variant_end = -1;

        // select final variant list
        Haplotype subregion_haplotype = derived_haplotype.substr_by_reference(curr_ref_start, curr_ref_end);
        std::vector<Variant> selected_variants = select_variants(candidate_variants, subregion_haplotype, event_sequences);
        for(size_t i = 0; i < selected_variants.size(); ++i) {
            int variant_end = selected_variants[i].ref_position + selected_variants[i].ref_seq.length() - 1;
            max_variant_end = std::max(variant_end, max_variant_end);

            derived_haplotype.apply_variant(selected_variants[i]);
            if(opt::verbose > 0) {
                selected_variants[i].write_vcf(stderr);
            }
        }
        
        if(max_variant_end == -1) {
            // no variants in region
            curr_ref_start += 1;
        } else {
            curr_ref_start = max_variant_end;
        }
    }
    return derived_haplotype;
}


Haplotype call_variants_for_region(const std::string& contig, int region_start, int region_end)
{
    const int BUFFER = 20;
    int STRIDE = 100;

    if(region_start < BUFFER)
        region_start = BUFFER + 20;

    // load the region, accounting for the buffering
    AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file);
    alignments.load_region(contig, region_start - BUFFER, region_end + BUFFER);
    Haplotype derived_haplotype(contig,
                                alignments.get_region_start(),
                                alignments.get_reference());

    for(int subregion_start = region_start;
            subregion_start < region_end; 
            subregion_start += STRIDE)
    {
        int subregion_end = subregion_start + STRIDE;

        int buffer_start = subregion_start - BUFFER;
        int buffer_end = subregion_end + BUFFER;
        buffer_end = std::min(region_end, buffer_end);

        // extract data from alignment database
        std::string ref_string = alignments.get_reference_substring(contig, buffer_start, buffer_end);
        std::vector<std::string> read_strings = alignments.get_read_substrings(contig, buffer_start, buffer_end);
        std::vector<HMMInputData> event_sequences = alignments.get_event_subsequences(contig, buffer_start, buffer_end);
        
        if(opt::verbose > 1) {
            fprintf(stderr, "Calling:\n");
            fprintf(stderr, "%s:%d-%d using buffer range [%d %d]\n", contig.c_str(), subregion_start, subregion_end, buffer_start, buffer_end);
            fprintf(stderr, "%s\n", ref_string.c_str());
        }

        std::vector<Variant> candidate_variants = alignments.get_variants_in_region(contig, buffer_start, buffer_end, 0.1);

        // extract potential variants from read strings
        //std::vector<Variant> candidate_variants = generate_all_snps(ref_string);

        //std::vector<Variant> candidate_variants = generate_variants_from_reads(ref_string, read_strings);
        //filter_variants_by_count(candidate_variants, opt::min_read_evidence);
        if(opt::snps_only) {
            filter_out_non_snp_variants(candidate_variants);
        }

        // remove variants that are inside of the buffer regions
        std::vector<Variant> tmp;
        for(size_t i = 0; i < candidate_variants.size(); ++i) {
            Variant& v = candidate_variants[i];

            int p = v.ref_position;
            if(p - buffer_start >= BUFFER && buffer_end - p >= BUFFER) {
                tmp.push_back(v);
            }
        }
        candidate_variants.swap(tmp);

        // Add variants into the haplotype
        Haplotype subregion_haplotype = derived_haplotype.substr_by_reference(buffer_start, buffer_end);
        std::vector<Variant> selected_variants = select_variants(candidate_variants, subregion_haplotype, event_sequences);
        for(size_t i = 0; i < selected_variants.size(); ++i) {
            derived_haplotype.apply_variant(selected_variants[i]);
            if(opt::verbose > 0) {
                selected_variants[i].write_vcf(stdout);
            }
        }
    }

    return derived_haplotype;
}

void parse_consensus_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'e': arg >> opt::event_bam_file; break;
            case 'w': arg >> opt::window; break;
            case 'o': arg >> opt::output_file; break;
            case 'm': arg >> opt::min_read_evidence; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_SNPS_ONLY: opt::snps_only = 1; break;
            case OPT_PROGRESS: opt::show_progress = 1; break;
            case OPT_VCF: arg >> opt::output_vcf; break;
            case OPT_HELP:
                std::cout << CONSENSUS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONSENSUS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 0) {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } else if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::num_threads <= 0) {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
        die = true;
    }
    
    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << CONSENSUS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int consensus_main(int argc, char** argv)
{
    parse_consensus_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    // Parse the window string
    // Replace ":" and "-" with spaces to make it parseable with stringstream
    std::replace(opt::window.begin(), opt::window.end(), ':', ' ');
    std::replace(opt::window.begin(), opt::window.end(), '-', ' ');

    const int WINDOW_LENGTH = 10000;
    const int WINDOW_OVERLAP = 200;

    std::stringstream parser(opt::window);
    std::string contig;
    int start_base;
    int end_base;
    
    parser >> contig >> start_base >> end_base;
    end_base = std::min(end_base, get_contig_length(contig) - 1);

    FILE* out_fp = NULL;

    if(!opt::output_file.empty()) {
        out_fp = fopen(opt::output_file.c_str(), "w");
    } else {
        out_fp = stdout;
    }

    FILE* out_vcf = NULL;
    if(!opt::output_vcf.empty()) {
        out_vcf = fopen(opt::output_vcf.c_str(), "w");
        Variant::write_vcf_header(out_vcf);
    }

    fprintf(stderr, "TODO: train model\n");
    fprintf(stderr, "TODO: filter data\n");

    Haplotype haplotype = call_variants_for_region_bb(contig, start_base, end_base);

    fprintf(out_fp, ">%s:%d-%d\n%s\n", contig.c_str(), start_base, end_base, haplotype.get_sequence().c_str());

    if(!opt::output_vcf.empty()) {
        std::vector<Variant> variants = haplotype.get_variants();
        for(size_t vi = 0; vi < variants.size(); vi++) {
            variants[vi].write_vcf(out_vcf);
        }
    }

    if(out_fp != stdout) {
        fclose(out_fp);
    }

    if(out_vcf != NULL) {
        fclose(out_vcf);
    }
}
