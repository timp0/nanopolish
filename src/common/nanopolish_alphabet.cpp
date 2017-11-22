//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alphabet -- support for multiple alphabets
//
#include <cassert>
#include <vector>
#include "nanopolish_alphabet.h"

//
// DNAAlphabet
// 
const uint8_t DNAAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const char* DNAAlphabet::_name = "nucleotide";
const char* DNAAlphabet::_base = "ACGT";
const char* DNAAlphabet::_complement = "TGCA";
const uint32_t DNAAlphabet::_size = 4;

//
// UtoTRNAAlphabet
// 
const uint8_t UtoTRNAAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const char* UtoTRNAAlphabet::_name = "u_to_t_rna";
const char* UtoTRNAAlphabet::_base = "ACGT";
const char* UtoTRNAAlphabet::_complement = "TGCA";
const uint32_t UtoTRNAAlphabet::_size = 4;

//
// methyl-cytosine in CG context
//
const uint8_t MethylCpGAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylCpGAlphabet::_name = "cpg";
const char* MethylCpGAlphabet::_base = "ACGMT";
const char* MethylCpGAlphabet::_complement = "TGCGA";
const uint32_t MethylCpGAlphabet::_size = 5;

const uint32_t MethylCpGAlphabet::_num_recognition_sites = 1;
const uint32_t MethylCpGAlphabet::_recognition_length = 2;
const char* MethylCpGAlphabet::_recognition_sites[] = { "CG" };
const char* MethylCpGAlphabet::_recognition_sites_methylated[] = { "MG" };
const char* MethylCpGAlphabet::_recognition_sites_methylated_complement[] = { "GM" };

//
// Dam methylation: methyl-adenine in GATC context
//
const uint8_t MethylDamAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylDamAlphabet::_name = "dam";
const char* MethylDamAlphabet::_base = "ACGMT";
const char* MethylDamAlphabet::_complement = "TGCTA";
const uint32_t MethylDamAlphabet::_size = 5;

const uint32_t MethylDamAlphabet::_num_recognition_sites = 1;
const uint32_t MethylDamAlphabet::_recognition_length = 4;
const char* MethylDamAlphabet::_recognition_sites[] = { "GATC" };
const char* MethylDamAlphabet::_recognition_sites_methylated[] = { "GMTC" };
const char* MethylDamAlphabet::_recognition_sites_methylated_complement[] = { "CTMG" };

//
// Dcm methylation: methyl-cytosine in CCAGG and CCTGG context
//
const uint8_t MethylDcmAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylDcmAlphabet::_name = "dcm";
const char* MethylDcmAlphabet::_base = "ACGMT";
const char* MethylDcmAlphabet::_complement = "TGCGA";
const uint32_t MethylDcmAlphabet::_size = 5;

const uint32_t MethylDcmAlphabet::_num_recognition_sites = 2;
const uint32_t MethylDcmAlphabet::_recognition_length = 5;
const char* MethylDcmAlphabet::_recognition_sites[] = { "CCAGG", "CCTGG" };
const char* MethylDcmAlphabet::_recognition_sites_methylated[] = { "CMAGG", "CMTGG" };
const char* MethylDcmAlphabet::_recognition_sites_methylated_complement[] = { "GGTMC", "GGAMC" };


//
// Added by Timp 17/10/20
// sin395 methylation: methyl-cytosine in GATC context
//
const uint8_t MethylSin395Alphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylSin395Alphabet::_name = "sin395";
const char* MethylSin395Alphabet::_base = "ACGMT";
const char* MethylSin395Alphabet::_complement = "TGCGA";
const uint32_t MethylSin395Alphabet::_size = 5;

const uint32_t MethylSin395Alphabet::_num_recognition_sites = 2;
const uint32_t MethylSin395Alphabet::_recognition_length = 4;
const char* MethylSin395Alphabet::_recognition_sites[] = { "GATC" };
const char* MethylSin395Alphabet::_recognition_sites_methylated[] = { "GATM" };
const char* MethylSin395Alphabet::_recognition_sites_methylated_complement[] = { "CTAG" };



//
// Added by yfan 17/11/22
// fnu4h methylation: methyl-cytosine in GCNGC context
//
const uint8_t Methylfnu4hAlphabet::_rank[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
  0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* Methylfnu4hAlphabet::_name = "fnu4h";
const char* Methylfnu4hAlphabet::_base = "ACGMT";
const char* Methylfnu4hAlphabet::_complement = "TGCGA";
const uint32_t Methylfnu4hAlphabet::_size = 5;

const uint32_t Methylfnu4hAlphabet::_num_recognition_sites = 4;
const uint32_t Methylfnu4hAlphabet::_recognition_length = 5;
const char* Methylfnu4hAlphabet::_recognition_sites[] = { "GCAGC" , "GCCGC" , "GCGGC" , "GCTGC" };
const char* Methylfnu4hAlphabet::_recognition_sites_methylated[] = { "GMAGC", "GMCGC" , "GMGGC" , "GMTGC" };
const char* Methylfnu4hAlphabet::_recognition_sites_methylated_complement[] = { "CGTCG" , "CGGCG" , "CGCCG" , "CGACG" };


//
// Added by yfan 17/11/22
// sdeaII methylation: methyl-cytosine in CCNGGC context
//
const uint8_t MethylsdeaIIAlphabet::_rank[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
  0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylsdeaIIAlphabet::_name = "sdeaII";
const char* MethylsdeaIIAlphabet::_base = "ACGMT";
const char* MethylsdeaIIAlphabet::_complement = "TGCGA";
const uint32_t MethylsdeaIIAlphabet::_size = 5;

const uint32_t MethylsdeaIIAlphabet::_num_recognition_sites = 4;
const uint32_t MethylsdeaIIAlphabet::_recognition_length = 6;
const char* MethylsdeaIIAlphabet::_recognition_sites[] = { "CCAGGC", "CCCGGC" , "CCGGGC" , "CCTGGC" };
const char* MethylsdeaIIAlphabet::_recognition_sites_methylated[] = { "CCAGGM" , "CCCGGM" , "CCGGGM" , "CCTGGM" };
const char* MethylsdeaIIAlphabet::_recognition_sites_methylated_complement[] = { "GGTCCG" , "GGGCCG", "GGCCCG" , "GGACCG" };


//
// Added by yfan 17/11/22
// hinfI methylation: methyl-adenine in GANTC context
//
const uint8_t MethylhinfIAlphabet::_rank[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
  0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylhinfIAlphabet::_name = "hinfI";
const char* MethylhinfIAlphabet::_base = "ACGMT";
const char* MethylhinfIAlphabet::_complement = "TGCTA";
const uint32_t MethylhinfIAlphabet::_size = 5;

const uint32_t MethylhinfIAlphabet::_num_recognition_sites = 4;
const uint32_t MethylhinfIAlphabet::_recognition_length = 5;
const char* MethylhinfIAlphabet::_recognition_sites[] = { "GAATC" , "GACTC" , "GAGTC" , "GATTC" };
const char* MethylhinfIAlphabet::_recognition_sites_methylated[] = { "GMATC" , "GMCTC" , "GMGTC" , "GMTTC" };
const char* MethylhinfIAlphabet::_recognition_sites_methylated_complement[] = { "CTTAG" , "CTGAG" , "CTCAG" , "CTAAG" };


//
// Added by yfan 17/11/22
// pspjdri methylation: meth cytosine in CCGG context
//
const uint8_t MethylpspjdriAlphabet::_rank[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
  0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylpspjdriAlphabet::_name = "pspjdri";
const char* MethylpspjdriAlphabet::_base = "ACGMT";
const char* MethylpspjdriAlphabet::_complement = "TGCGA";
const uint32_t MethylpspjdriAlphabet::_size = 5;

const uint32_t MethylpspjdriAlphabet::_num_recognition_sites = 1;
const uint32_t MethylpspjdriAlphabet::_recognition_length = 4;
const char* MethylpspjdriAlphabet::_recognition_sites[] = { "CCGG" };
const char* MethylpspjdriAlphabet::_recognition_sites_methylated[] = { "MCGG" };
const char* MethylpspjdriAlphabet::_recognition_sites_methylated_complement[] = { "GGCC" };



//
// T-mods - uracil or mod T substitution
// Added by Timp 17/08/16
//
const uint8_t ModTAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* ModTAlphabet::_name = "tmod";
const char* ModTAlphabet::_base = "ACGMT";
const char* ModTAlphabet::_complement = "TGCAA";
const uint32_t ModTAlphabet::_size = 5;

const uint32_t ModTAlphabet::_num_recognition_sites = 1;
const uint32_t ModTAlphabet::_recognition_length = 1;
const char* ModTAlphabet::_recognition_sites[] = { "T" };
const char* ModTAlphabet::_recognition_sites_methylated[] = { "M" };
const char* ModTAlphabet::_recognition_sites_methylated_complement[] = { "A" };



// Global objects
DNAAlphabet gDNAAlphabet;
MethylCpGAlphabet gMCpGAlphabet;
MethylDamAlphabet gMethylDamAlphabet;
MethylDcmAlphabet gMethylDcmAlphabet;
MethylSin395Alphabet gMethylSin395Alphabet;
Methylfnu4hAlphabet gMethylfnu4hAlphabet;
MethylsdeaIIAlphabet gMethylsdeaIIAlphabet;
MethylhinfIAlphabet gMethylhinfIAlphabet;
MethylpspjdriAlphabet gMethylpspjdriAlphabet;
ModTAlphabet gModTAlphabet;
UtoTRNAAlphabet gUtoTRNAAlphabet;

std::vector<const Alphabet*> get_alphabet_list()
{
    std::vector<const Alphabet*> list = { &gDNAAlphabet, 
                                          &gMCpGAlphabet, 
                                          &gMethylDamAlphabet,
                                          &gMethylDcmAlphabet,
					  &gMethylSin395Alphabet,
					  &gMethylfnu4hAlphabet,
					  &gMethylsdeaIIAlphabet,
					  &gMethylhinfIAlphabet,
					  &gMethylpspjdriAlphabet,
					  &gModTAlphabet,
                                          &gUtoTRNAAlphabet };
    return list;
}

// Select the alphabet that best matches bases
const Alphabet* best_alphabet(const char *bases)
{
    std::vector<const Alphabet*> list = get_alphabet_list();

    for (auto alphabet: list)
        if (alphabet->contains_all(bases))
            return alphabet;

    return nullptr;                
}

// Select the alphabet by name
const Alphabet* get_alphabet_by_name(const std::string& name)
{
    std::vector<const Alphabet*> list = get_alphabet_list();

    for (auto alphabet: list)
        if (alphabet->get_name() == name)
            return alphabet;
    
    fprintf(stderr, "Error, unknown alphabet name: %s\n", name.c_str());
    exit(EXIT_FAILURE);
    return nullptr; 
}

