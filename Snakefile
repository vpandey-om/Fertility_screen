
import sys
# #input files
# # these are the files where all vectors and mutant genes can be found
# # based on this information we will find out how many are dropouts or
# #there is contaminations
# input_vectors=['data/input_vector.txt','data/input_pool2.txt',
# 'data/input_pool3.txt','data/input_pool4.txt','data/input_pool6.txt',
# 'data/input_pool5.txt','data/input_pool7.txt','data/input_pool5_2.txt',
# 'data/input_pool7_2.txt']
# ## these are the total counts we found it after barseq experiment
# # these files are genrated based on raw fastq files form  sequencing
# barcode_count_files=['data/result_240120_barcode_counts_table.txt',
# 'data/barcode_counts_table_170620_pool2.txt',
# 'data/barcode_counts_table_200720_pool3.txt',
# 'data/barcode_counts_table_170620_pool4.txt',
# 'data/barcode_counts_table_280720_pool6.txt',
# 'data/barcode_counts_table_20920_pool5.txt',
# 'data/barcode_counts_table_20920_pool7.txt',
# 'data/barcode_counts_table_20920_pool5.txt',
# 'data/barcode_counts_table_20920_pool7.txt',
# ]
#
# ## these are experimental designs. This will help us to do variability analysis
# manifest_files=['data/manifests_pool1.txt','data/manifest_pool2.txt',
# 'data/manifest_pool3.txt','data/manifest_pool4.txt',
# 'data/manifest_pool6.txt','data/manifest_pool5_1.txt',
# 'data/manifest_pool7_1.txt','data/manifest_pool5_2.txt',
# 'data/manifest_pool7_2.txt'
# ]
#
# ### params for motility screen
# motility_vector=['data/input_male_pool_dis.txt','data/input_male_pool_slow.txt']
# motility_barcode=['data/barcode_counts_table_050221_male_pool.txt']
# motility_manifest=['data/manifest_male_pool_dis.txt','data/manifest_male_pool_slow.txt']
# ##
#
# config={'vectors': input_vectors,'barcode':barcode_count_files,
# 'manifest':manifest_files,'motility_vectors':motility_vector,
# 'motility_barcode':motility_barcode, 'motility_manifest':motility_manifest,
# 'input_file':"barseq-raw/testdata/sequence",
# "barcode_file":"barseq-raw/barcode_to_gene_210920.csv",
#  "output_dir":"barseq-raw/rawCountResult"}
configfile: "config.yaml"
# input_file=config["input_file"]
# barcode_file=config["barcode_file"]
# output_dir=config["output_dir"]
print(config['input_file'])
# sys.exit("Exiting Snakemake due to a specific condition")

# Rule to run script1 on input files
# this script is for combining all pools



rule raw_barseq:
    params:
        input_file=config['input_file'],
        barcode_file=config['barcode_file'],
        output_dir=config['output_dir']
    conda:
        "rawbarseq_environment.yaml"
    shell:
        # "python barseq-raw/barseq.py -i {input[0]} -b {input[1]} -r {output}"
        "python barseq-raw/barseq.py -i {params.input_file} -b {params.barcode_file} -r {params.output_dir}"
    # script:
    #     "barseq-raw/barseq.py"

rule remove_zeros:
    params:
        output_dir=config['output_dir']
    conda:
        "rawbarseq_environment.yaml"
    shell:
        # "python barseq-raw/barseq.py -i {input[0]} -b {input[1]} -r {output}"
        "python barseq-raw/remove_zero_barseq_count.py {params.output_dir}/barcode_counts_table.csv {params.output_dir}/removed_zeros_barcode_counts_table.csv"

# Define a default rule to execute run_script1 and run_script2 in sequence
# rule all:
#     params:
#         output_dir=config['output_dir']
#     input:
#         "{params.output_dir}/barcode_counts_table.csv"

rule combine_pools:
    conda:
        "environment.yaml"
    priority: 10
    params:
        vectors=config['input_vectors'],
        barcode=config['barcode_count_files'],
        manifest=config['manifest_files']
    output:
        "output/final_phenocall.txt"
    script:
        "codes/combine_all_pool.py"

## plot male and female fertility 1) male and female fertilty rate with ranking order
# 2) scatter plot of male and female using half circle filled with diffrent color
rule plot_fertility:
    input:
        "data/phenotype/Phenotype_call_final_100621.xlsx"
    priority: 10
    conda:
        "environment.yaml"
    output:
        ["output/s_curve_male.pdf","output/s_curve_female.pdf","output/half_circle_male_female_scatter.pdf"]
    script:
        "codes/s_curve_plot_Fig1.py"


# This rule is for mpmp analysis
# enrichment analysis for male_only, female_only, reduced motality
rule mpmp_enrichment:
    input:
        ["data/mpmp/mpmp_pbe.pickle","data/mpmp/MPMP.xlsx","data/mpmp/male_only_genes.txt","data/mpmp/female_only_genes.txt","data/mpmp/notreduced_motality.txt"]
    priority: 10
    conda:
        "environment.yaml"
    output:
        ["output/male_mpmp_enrich.txt","output/female_mpmp_enrich.txt","output/notreduced_motality_enrich.txt"]
    script:
        "codes/mpmp_enrichment_table.py"

## run mpmp results in violin most significant pathways
rule mpmp_violin:
    input:
        ["data/mpmp/mpmp_pbe.pickle","data/mpmp/MPMP.xlsx","data/mpmp/male_mpmp_enrich.txt",
        "data/mpmp/female_mpmp_enrich.txt","data/phenotype/Phenotype_call_final_100621.xlsx"]
    conda:
        "plotnine_environment.yaml"
    output:
        ["output/male_female_mpmp_violin.pdf"]
    script:
        "codes/mpmp_violin_plot.py"


## run motility screen
rule motility_screen:
    input:
        ['data/barcode_counts_table_050221_male_pool.txt']
    output:
        ["output/motility_disp.xlsx","output/motility_slow.xlsx"]
    conda:
        "environment.yaml"
    priority: 1
    params:
        vectors=config['motility_vectors'],
        barcode=config['motility_barcode'],
        manifest=config['motility_manifest']
    script:
        "codes/motility_screen.py"

## combine error analysis


rule error_plot:
    input:
        ["data/errorP/pilot.p","data/errorP/pool1.p",
        "data/errorP/pool2.p","data/errorP/pool3.p",
        "data/errorP/pool4.p","data/errorP/pool6.p",
        "data/errorP/pool5_1.p","data/errorP/pool7_1.p",
        "data/errorP/pool5_2.p","data/errorP/pool7_2.p"]
    priority: 10
    conda:
        "environment.yaml"
    output:
        ["output/error_analysis.pdf"]
    script:
        "codes/error_analysis.py"


rule mca_gene_plot:
    input:
        ["data/ginny_gene_mca.xlsx",
        "data/phenotype/Phenotype_call_final_100621.xlsx",
        "data/mpmp/male_only_genes.txt",
        "data/mpmp/female_only_genes.txt",
        ]
    priority: 10
    conda:
        "environment.yaml"
    output:
        ["output/male_mca.pdf","output/female_mca.pdf"]
    script:
        "codes/mca_gene_cluster.py"




## run_script2 is for generating ordered fertility rate(s-curve)
