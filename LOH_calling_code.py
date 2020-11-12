import re,sys,os, glob
from string import *
import math, numpy, scipy, math
from numpy import array
from scipy import stats
#import pvalue_combine

##############
# This code is designed to detect LOH events of interested gene per one cancer patient using two different exome sequencing samples from normal and tumor.


#To run this code, below files need to prepare:

#1. ExAC variant file from http://exac.broadinstitute.org/downloads: ExAC.r0.3.sites.vep.vcf
#2. gene2locus information file: UCSC-2014-10-30_hg19_refGene.txt
#3. cancer sequencing samples (please, see the example files from input directory), (1) normal sample from cancer patient, (2) tumor sample from cancer patient

#############

# Constants
exac_input = "input_files/ExAC.r0.3.sites.vep.vcf" ## download from : http://exac.broadinstitute.org/downloads
neighbor_genes = "input_files/UCSC-2014-10-30_hg19_refGene.txt"
query_KB = 50*1000 # to define neighboring variants. The definition of neighboring variants can be changed.

def extractPASS(chr_query) :
    """Get variants in gNOMAD/ExAC (or similar) that passed all the filters

    Parameters
    ----------
        chr_query : str
            Chromosome number where the user needs to extract the PASS variants

    Returns
    -------
        dict
            Dictionary with the variants that have passed the filters in the population database (column FILTER = PASS)
            Dictionary format: {"chr_number" : ["chr-position-start-end", ...], ...}
    """
    ExAC_PASS = {}
    ExAC_PASS[chr_query] = []
    print("INFO: Reading ExAC input data")
    with open(exac_input, "r") as fi :
        for line in fi:
            line = line.strip()
            field = line.split('\t')
            # Don't read the vcf header
            if '##' in line:
                continue
            elif '#CHROM' in line:
                for i in range(len(field)): # Find the column number that corresponds to the FILTER column
                    if field[i] == 'FILTER':
                        PASS_index = i
            else:
                chr_input = field[0]
                if chr_input != chr_query:
                    continue
                PASS_query = field[PASS_index]
                if PASS_query == 'PASS':
                    variant_query = "{chr}-{pos}-{ref}-{alt}".format(chr = field[0], pos = field[1], ref = field[3], alt = field[4]) # Format (chromosome-position-reference-alterated)
                    ExAC_PASS[chr_input].append(variant_query)

    if len(ExAC_PASS[chr_query]) == 0 :
        print("WARNING: No population variants found. Please check the input")

    return ExAC_PASS

def extractGenes() :
    """Get the length, start, end and chromosome from all the genes

    Returns
    -------
        dict
            Dictionary with all mentioned information in format
            {gene_name : [[length1, length2, ...], [start1, start2, ...], [end1, end2, ...], [chromosome1, chromosome2, ...]]}
    """
    print("INFO: Extracting the gene information")
    gene2locus = {}
    with open(neighbor_genes,'r') as fi :
        for line in fi :
            line = line.strip()
            field = line.split('\t')
            gene_name = field[12]
            start_pos = int(field[4])
            end_pos = int(field[5])
            length = end_pos - start_pos
            chr_info = field[2].replace("chr", "")

            if gene_name not in gene2locus.keys():
                gene2locus[gene_name] = [[],[],[],[]]
                gene2locus[gene_name][0].append(length)
                gene2locus[gene_name][1].append(start_pos)
                gene2locus[gene_name][2].append(end_pos)
                gene2locus[gene_name][3].append(chr_info)
            else:
                gene2locus[gene_name][0].append(length)
                gene2locus[gene_name][1].append(start_pos)
                gene2locus[gene_name][2].append(end_pos)

    return gene2locus

def getNeighborGenes(gene_query, gene2locus) :
    """Finds the gens that are closer than query_KB from the genes of interest

    Gene the closer genes from a list of genes passed as parameter

    Parameters
    ----------
        gene_query : list
            Gene names that the user is interested in finding the neighbors
        gene2locus : dict
            Gene list extracted in extractGenes function

    Returns
    -------
        gene2degree : dict
            Genes of interest and the corresponding neighbors in format: { gene_name : [neighbor1, neighbor2, ...]}
        degree2gene : dict
            Genes of interest' neighbors in format: {neighbor1 : [gene_of_interest], neighbor2 : [gene_of_interest], ...}
        gene_query_total : set
            Unique gene names find that are neighbors of the list of gens passed
    """
    gene2degree = {}
    degree2gene = {}
    gene_query_total = set([]) #A set is an unordered collection of items. Every set element is unique (no duplicates) and must be cannot be changed

    # gene_query has the gene ids the user is interested for
    for gids in gene_query:
        gene2degree[gids] = []
        gene_query_total.add(gids)

        length_query = numpy.median(gene2locus[gids][0])
        start_query = numpy.min(gene2locus[gids][1])
        end_query = numpy.max(gene2locus[gids][2])
        center_query = (start_query + end_query)/2.0
        chr_input = gene2locus[gids][3][0]
        possible_start = center_query - query_KB
        possible_end = center_query + query_KB

        for gene_cand in gene2locus.keys() :
            chr_cand = gene2locus[gene_cand][3][0]
            if chr_input!= chr_cand: # TODO: Canviar este if
                continue
            elif gene_cand == gids:
                continue

            length_cand = numpy.median(gene2locus[gene_cand][0])
            start_cand = numpy.min(gene2locus[gene_cand][1])
            end_cand = numpy.max(gene2locus[gene_cand][2])

            if end_cand > possible_start and end_cand < possible_end :
                gene2degree[gids].append(gene_cand)
                gene_query_total.add(gene_cand)
                if gene_cand in degree2gene.keys():
                    degree2gene[gene_cand].append(gids)
                else:
                    degree2gene[gene_cand] = [gids]
            elif start_cand < possible_end and start_cand > possible_start :
                gene2degree[gids].append(gene_cand)
                gene_query_total.add(gene_cand)
                if gene_cand in degree2gene.keys():
                    degree2gene[gene_cand].append(gids)
                else:
                    degree2gene[gene_cand] = [gids]

    return gene2degree, degree2gene, gene_query_total

def convert2dict(line, header) :
    field = line.strip().split("\t")
    dc = {}
    it = 0
    for h in header :
        if h == "Otherinfo" :
            dc["vcf_chrom"] = field[it]
            dc["vcf_pos"] = field[it+1]
            dc["vcf_id"] = field[it+2]
            dc["vcf_ref"] = field[it+3]
            dc["vcf_alt"] = field[it+4]
            dc["vcf_qual"] = field[it+5]
            dc["vcf_filter"] = field[it+6]
            dc["vcf_info"] = field[it+7]
            header_format = field[it+8].split(":")
            it2 = 0
            values_format = field[it+9].split(":")
            for h2 in header_format :
                dc["vcf_{}".format(h2)] = values_format[it2]
                it2 += 1
        else :
            dc[h] = field[it]
            it += 1

    return dc

def getHghestMAF(vcf) :
    keys = ["AF_fin", "AF_eas", "non_neuro_AF_popmax", "non_topmed_AF_popmax", "controls_AF_popmax", "AF_female", "AF_popmax", "AF_asj", "AF_ami", "AF_male", "AF_nfe", "AF_oth", "AF_amr",
     "non_cancer_AF_popmax", "AF", "AF_afr", "AF_sas", "AF_raw"]
    max = -1
    for f in keys :
        try :
            aux = float(vcf[f])
            if aux > max :
                max = aux
        except ValueError as e : # Control de errors obtained if there is no number in the allele frequency columns
            pass

    if max == -1 :
        max = "NA"
    else :
        print(vcf)
        sys.exit()

    return max

def readGermline(path, chr_query, gene_query_total) :
    print("INFO: Extracting information from germline file")
    germline_variant = {}

    ### collecting germline variants (all possible PASS variants)
    """
    # Input files' format #
        VCF file and ANNOVAR gene_anno information stored in INFO column
        The variannts from the genes of interest and the surrounding genes are stored in germline_variant dict, and in gene2variant dict  after this loop
    """
    header = []
    with open(path, "r") as fi :
        for line in fi :
            if len(header) == 0: # Read the header. Extract the column number for each column
                header = line.strip().split('\t')
            else :
                dcAux = convert2dict(line, header)
                try :
                    chr_input = dcAux["vcf_chrom"]
                    if chr_input == chr_query:
                        pos_query = dcAux["vcf_pos"]
                        # variant format: chr-pos-ref-alt
                        variant = '{chr}-{pos}-{ref}-{alt}'.format(chr = dcAux["vcf_chrom"], pos = dcAux["vcf_pos"], ref = dcAux["vcf_ref"], alt = dcAux["vcf_alt"])
                        PASS_index = dcAux["vcf_filter"]
                        gene_check = 0
                        if PASS_index == 'PASS':# to collect high-quality variant
                            if dcAux["Gene.refGene"] in gene_query_total : #collect variants of query genes and neighboring genes
                                gene_name = dcAux["Gene.refGene"]
                                exon = dcAux["Func.refGene"]
                                mut_type = dcAux["ExonicFunc.refGene"]
                                MAF = getHghestMAF(dcAux)
                                # Store the information regarding the genotype (GT) and the ref/alt depth (AD)
                                GQX_info = "{}:{}".format(dcAux["vcf_GT"], dcAux["vcf_AD"])
                                genotype = dcAux["vcf_GT"]
                                if genotype == '0/1' :
                                    germline_variant[variant] = ["{gene}\t{mutation}\t{maf}\t{exon}\t{GQX}".format(gene = gene_name, mutation = mut_type, maf = MAF, exon = exon, GQX = GQX_info)]
                except KeyError as e :
                    print("ERROR: Key {} not found in dict. The keys found were: {}".format(e, dcAux.keys()))

    return germline_variant

def germline2somatic_variant_mapping_LOHcalling (germline_sample, somatic_sample, chr_query, gene_query):

    #### step 1: collecting ExAC PASS variants in format "chr_number-position-reference-alterated"
    # ExAC_PASS = extractPASS(chr_query)

    #### step 2: neighboring gene
    gene2locus = extractGenes()

    ### Finding neighboring genes
    gene2degree, degree2gene, gene_query_total = getNeighborGenes(gene_query, gene2locus)



    """
    # Variables filled after this loop #
        * germline_variant dict format: { "chr-pos-ref-alt" : geneName\tmutation_type, MAF, exonic_mutation_type, 0/1}
        * gene2variant dict format : { "gene_name" :  "chr-pos-ref-alt"} NOT USED
    """
    fsomatic = open('%s'%(somatic_sample), 'r')
    somatic_variant = {}
    ### mapping germline variants from normal sample to matched somatic variants from tumor sample
    for line in fsomatic.xreadlines():
        line = line.strip()
        field = line.split('\t')
        if '#' in line: # Get the indexes of the interesting columns
            if '#CHROM' in line:
                for i in range(len(field)):
                    if field[i] == '#CHROM':
                        chr_index = i
                    elif field[i] == 'Sample_index':
                        sample_index = i
                    elif field[i] == 'FILTER':
                        filter_index = i
        else:  # Extract the variants information for the vcf body
            chr_input = field[chr_index]
            if chr_input!= chr_query:
                continue
            # profile format : chr-pos-ref-alt
            profile = '%s-%s-%s-%s'%(field[chr_index], field[chr_index+1], field[chr_index+3], field[chr_index+4])
            # Remove the variants in the gens of interest that are in the germline file
            if profile in germline_variant.keys():
                somatic_GQX_info = field[sample_index]
                PASS_index = field[filter_index]
                somatic_variant[profile] = ['%s\t%s'%(PASS_index, somatic_GQX_info)]

    ############################################################
    # Up to here we have:
    # * ExAC_PASS -> List with all the SNPs from ExAC, in the chromosome passed as parameter, that passed all the filters
    # * gene2locus -> Dict with the name of the gene as key, followed by lists with the lengths, start positions, end positions, and chromosome
    # * gene_query_total -> List of genes of interested passed as parameter plus the genes that are closer to them
    # * germline_variant -> Dict with the variants of interest from the germline vcf. Dict format specified in line 208
    # * somatic_variant -> Dict with the variants of interest from the Somatic vcf. Dict format in line 208
    #############################################################
    gene_info = {}
    for ids in gene_query: # gene_query is the list of genes passed as parameter to the function
        gene_info[ids] = [[]]

    for ids in germline_variant.keys():
        if ids in somatic_variant.keys():
            variant = ids
            ##############
            germline_mapping_info = germline_variant[ids][0].split('\t')

            gene = germline_mapping_info[0]
            mut_type = germline_mapping_info[1]
            MAF = germline_mapping_info[2]
            exon = germline_mapping_info[3]
            chr_input = ids.split('-')[0]

            ##############
            germline_GQX = germline_mapping_info[-1]
            germline_read_count = germline_GQX.split(':')[-1].split(',')
            germline_score1 = int(germline_read_count[0])
            germline_score2 = int(germline_read_count[1])

            somatic_mapping_info = somatic_variant[ids][0].split('\t')
            somatic_GQX = somatic_mapping_info[-1]
            somatic_read_count = somatic_GQX.split(':')[-1].split(',')
            somatic_score1 = int(somatic_read_count[0])
            somatic_score2 = int(somatic_read_count[1])

            if  float(somatic_score1 + somatic_score2) == 0.0:
                continue

            fisher_pvalue = stats.fisher_exact([[germline_score1, germline_score2], [somatic_score1, somatic_score2]])[1]
            variant_info = '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(variant,fisher_pvalue, mut_type, MAF, exon, germline_GQX, somatic_GQX)

            if gene in gene_info.keys():# collecting all possible variants of query gene
                gene_info[gene][0].append(variant_info)

            elif gene in degree2gene.keys():# collecting all possible variants of neighboring genes
                neighboring_list = degree2gene[gene]

                for friend_query in neighboring_list:
                    ###########################
                    start_pos_query = numpy.min(gene2locus[friend_query][1])
                    end_pos_query = numpy.max(gene2locus[friend_query][2])

                    center_pos_query = (start_pos_query + end_pos_query)/2.0
                    possible_start = center_pos_query - query_KB
                    possible_end = center_pos_query + query_KB
                    ###########################
                    if position >= possible_start and position <= possible_end:
                        gene_info[friend_query][0].append(variant_info)

                    elif position <= possible_end and position >= possible_start:
                        gene_info[friend_query][0].append(variant_info)

    #################################
    effect_size_upper = 0.7 ## cut-off of effect size can be changed.
    effect_size_lower = 0.3
    gene_LOH = {}
    for gids in gene_info.keys():
        feature_list = gene_info[gids][0]
        gene_LOH[gids] = []

        #0: variant,1: fisher_pvalue, 2: mut_type, 3: MAF, 4: exon type, 5: germline_VAF, 6: somatic_VAF

        for sids in feature_list:
            mut = sids.split('\t')[2]
            exon = sids.split('\t')[4]

            germline_read = sids.split('\t')[-2]
            somatic_read = sids.split('\t')[-1]

            germline_info = germline_read.split(':')[-1].split(',')
            somatic_info = somatic_read.split(':')[-1].split(',')
            NormalVAF = int(germline_info[1]) / float(int(germline_info[0]) + int(germline_info[1]))
            TumorVAF = int(somatic_info[1]) / float(int(somatic_info[0]) + int(somatic_info[1]))
            fisher_pvalue = float(sids.split('\t')[1])
            if TumorVAF >= effect_size_upper or TumorVAF <= effect_size_lower:#effect size threshold to define LOH
                gene_LOH[gids].append(fisher_pvalue)

    fout = open('./output/LOH_mapping_output.txt', 'w')
    fout.write('Gene\tLOH_type\tNum_mapping_variant\n')

    for gene_name in gene_info.keys():
        pvalue_combination = pvalue_combine.combine_pvalues(gene_LOH[gene_name])[1]
        if pvalue_combination <= 0.05:
            fout.write('%s\t%s\t%s\n'%(gene_name, 'LOH', len(gene_LOH[gene_name])))
        else:
            fout.write('%s\t%s\t%s\n'%(gene_name, 'noLOH', len(gene_LOH[gene_name])))
    fout.close()
    print(fout)
    return 'LOH detection'

# Unit tests for all the functions
chr_query = 'chr17'
gene_query = ['BRCA1']
# test = extractPASS('17')
gene2locus = extractGenes()

gene2degree, degree2gene, gene_query_total = getNeighborGenes(gene_query, gene2locus)
germline_path = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6-6a6e-4e2c-a704-90952afeef25/strelkaGerm/results/variants/strelka.hg38_multianno.txt"
test2 = readGermline(germline_path, chr_query, gene_query_total)
print(test2)


#################################################

# TODO uncomment this to do the tests
# First dummy example

# germline_sample = 'input_files/normal_ex2.vcf'
# somatic_sample = 'input_files/tumor_ex2.vcf'
# gene_query = ['BRCA1']
# chr_query = '17'
#
# print germline2somatic_variant_mapping_LOHcalling(germline_sample, somatic_sample, chr_query, gene_query)
