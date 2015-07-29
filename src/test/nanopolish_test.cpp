//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_test -- test driver program
//
#define CATCH_CONFIG_MAIN
#include <stdio.h>
#include <string>
#include "logsum.h"
#include "catch.hpp"
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"
#include "nanopolish_profile_hmm.h"

// This code needs to be run before any of the program logic
// It sets up pre-computed values and caches
void initialize()
{
    p7_FLogsumInit();
}

TEST_CASE( "string functions", "[string_functions]" ) {

    // kmer rank
    REQUIRE( kmer_rank("AAAAA", 5) == 0 );
    REQUIRE( kmer_rank("GATGA", 5) == 568 );
    REQUIRE( kmer_rank("TTTTT", 5) == 1023 );
    REQUIRE( kmer_rank("GATGA", 5) == rc_kmer_rank("TCATC", 5 ) );

    // lexicographic increment
    std::string str = "AAAAA";
    lexicographic_next(str);
    REQUIRE( str == "AAAAC" );

    str = "AAAAT";
    lexicographic_next(str);
    REQUIRE( str == "AAACA" );

    // complement, reverse complement
    REQUIRE( complement('A') == 'T' );
    REQUIRE( reverse_complement("GATGA") == "TCATC" );

}

TEST_CASE( "math", "[math]") {
    GaussianParameters params;
    params.mean = 4;
    params.stdv = 2;
    params.log_stdv = log(params.stdv);

    REQUIRE( normal_pdf(2.25, params) == Approx(0.1360275) );
    REQUIRE( log_normal_pdf(2.25, params) == Approx(log(normal_pdf(2.25, params))) );
}

std::string event_alignment_to_string(const std::vector<AlignmentState>& alignment)
{
    std::string out;
    for(size_t i = 0; i < alignment.size(); ++i) {
        out.append(1, alignment[i].state);
    }
    return out;
}

TEST_CASE( "hmm", "[hmm]") {

    // read the FAST5
    SquiggleRead sr("test_read", "test/data/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch8_file30_strand.fast5");
    sr.transform();

    // The reference sequence to align to:
    std::string ref_subseq = "ATCAGTAAAATAACGTAGAGCGGTAACCTTGCCATAAAGGTCGAGTTTA"
                             "TTACCATCCTTGTTATAGACTTCGGCAGCGTGTGCTACGTTCGCAGCT";

    // Generate a HMMInputData structure to tell the HMM
    // which part of the read to align
    HMMInputData input[2];
    
    // template strand
    input[0].read = &sr;
    input[0].event_start_idx = 3;
    input[0].event_stop_idx = 88;
    input[0].event_stride = 1;
    input[0].rc = false;
    input[0].strand = 0;
    
    // complement strand
    input[1].read = &sr;
    input[1].event_start_idx = 6788;
    input[1].event_stop_idx = 6697;
    input[1].event_stride = -1;
    input[1].rc = true;
    input[1].strand = 1;

    // Local alignment test
    std::string ref_local = ref_subseq.substr(0, 17);

    HMMInputData local_input[2] = { input[0], input[1] };
    local_input[1].event_start_idx = 6788;
    local_input[1].event_stop_idx = 6770;

    printf("seq: %s\n", ref_subseq.c_str());
    printf("seq: %s %.2lf\n", ref_local.c_str(), profile_hmm_score(ref_local, local_input[1]));

    float max0 = -INFINITY;
    float max1 = -INFINITY;
    std::string ref_extend = ref_local + "AAAA";
    for(size_t n = 0; n < 256; n++) {

        std::vector<AlignmentState> a0 = profile_hmm_align(ref_extend, local_input[0]);
        std::vector<AlignmentState> a1 = profile_hmm_align(ref_extend, local_input[1]);

        double score0 = profile_hmm_score(ref_extend, local_input[0]);
        double score1 = profile_hmm_score(ref_extend, local_input[1]);

        if(score0 > max0) {
            max0 = score0;
        }
        if(score1 > max1) {
            max1 = score1;
        }

        char match = ref_subseq.find(ref_extend) != std::string::npos ? '*' : ' ';

        printf("seq: %s %.2lf %.2lf [%d %d] [%d %d] [%.2lf %.2lf] %c\n", 
                ref_extend.c_str(), score0, score1, 
                a1.front().event_idx, a1.front().kmer_idx,
                a1.back().event_idx, a1.back().kmer_idx,
                max0, max1, match);
        lexicographic_next(ref_extend);
    }

    {
        std::string debug_string = ref_local + "GAGC";
        std::vector<AlignmentState> local_alignment = profile_hmm_align(ref_local + "GAGC", local_input[1]);

        for(size_t i = 0; i < local_alignment.size(); ++i) {
            double lp_e = log_probability_match(sr, 
                                                kmer_rank(debug_string.c_str() + local_alignment[i].kmer_idx, K), 
                                                local_alignment[i].event_idx,
                                                1);

            printf("\t[%d %d]: %.2lf\n", local_alignment[i].event_idx, local_alignment[i].kmer_idx, lp_e);
        }
        std::string ea_str = event_alignment_to_string(local_alignment);
        printf("Local alignment string: %s\n", ea_str.c_str());
    }


    // expected output
    std::string expected_alignment[2];
    expected_alignment[0] = 
        "MMMMMEMKMKMMMMMMMKMMMKMMMKMMMMMMMMMKKMMEEEMMMMMMKMMMM" 
        "MMMKMMMMMKMKMKMEMKKMKMKKMMMMMMEMMMMKMKMEEMMMMKMEEEEEM";

    expected_alignment[1] = 
        "MMKMMMKMEEMMKMKMKMEMMMKMMMKMEMMMKMMMKMMMMMMMMMKKMEMMMM"
        "EMMMMMMMMKMKKMMMMMMMEMMMMMKMMMMMKMEMMMMMKMMMMMEEEEEEEEM";

    double expected_viterbi_last_state[2] = { -237.24, -265.66 };
    double expected_forward[2] = { -237.238434, -265.6596984863 };

    for(int si = 0; si <= 1; ++si) {

        // viterbi align
        std::vector<AlignmentState> event_alignment = profile_hmm_align(ref_subseq, input[si]);
        std::string ea_str = event_alignment_to_string(event_alignment);
    
        // check
        REQUIRE( ea_str == expected_alignment[si]);
        REQUIRE( event_alignment.back().l_fm == Approx(expected_viterbi_last_state[si]));

        // forward algorithm
        double lp = profile_hmm_score(ref_subseq, input[si]);
        REQUIRE(lp == Approx(expected_forward[si]));
    }
}
