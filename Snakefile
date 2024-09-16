# Snakefile

# Define configuration file that lists samples and metadata
# That is the file you should edit prior to analysis
configfile: "config.yaml"

# Get the software container from config file
singularity: config["singularity"]

# Get list of sample names from config file
SAMPLES = config["samples"]

# Get input file names from config file
ENVFACTORS = config["input_files"]["envfactors"]
ENVFACTORSNAMES = [
    "hillshade",
    "aspect.eastness",
    "aspect.northness",
    "topographic.exposure",
    "horizontal.curvature",
    "vertical.curvature",
    "vrm",
    "wetness.index.dinf",
    "diffus.solar.radiation",
    "direct.solar.radiation",
    "total.solar.radiation",
    "bio1",
    "bio2",
    "bio3",
    "bio4",
    "bio5",
    "bio6",
    "bio7",
    "bio8",
    "bio9",
    "bio10",
    "bio11",
    "bio12",
    "bio13",
    "bio14",
    "bio15",
    "bio16",
    "bio17",
    "bio18",
    "bio19",
    "slope.degrees"
]
POOLSIZES = config["input_files"]["poolsizes"]

# Get parameters from config file
N_SUBS = config["parameters"]["n_subs"]
N_CORE_REPLICATES = config["parameters"]["n_core_replicates"]
MIN_HAPLOID_POOL_SIZE = config["parameters"]["min_haploid_pool_size"]
N_PILOT = config["parameters"]["n_pilot"]

# Get resources from config file
RESOURCES = config["resources"]

# function to calculate the number of populations
import os

def count_words_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    return len(content.split())

rule all:
    input:
        expand("data/{sample}_efile_{envfactorname}", sample=SAMPLES, envfactorname=ENVFACTORSNAMES),
        expand("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out", sample=SAMPLES, i=range(1, N_CORE_REPLICATES+1)),    
        expand("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_summary_betai_reg.out", 
            sample=SAMPLES, envfactorname=ENVFACTORSNAMES, i=range(1, N_SUBS+1)),
    resources:
        runtime=RESOURCES["all"]["runtime"],
        mem_mb=RESOURCES["all"]["mem_mb"],
        slurm_partition=RESOURCES["all"]["slurm_partition"]

rule generate_complementary_inputs:
    input:
        envfactors=ENVFACTORS,
        poolsizes=POOLSIZES,
    resources:
        runtime=RESOURCES["generate_complementary_inputs"]["runtime"],
        mem_mb=RESOURCES["generate_complementary_inputs"]["mem_mb"],
        slurm_partition=RESOURCES["generate_complementary_inputs"]["slurm_partition"]
    output:
        poolsizefile="data/{sample}_poolsizes",
    script: "scripts/generate_complementary_inputs.R"

rule vcf2genobaypass:
    input:
        vcf="data/{sample}.vcf",
        poolsizes="data/{sample}_poolsizes",
        poolnames="data/{sample}_poolnames",
    params:
        prefix="{sample}",
        subs=N_SUBS,
    resources:
        runtime=RESOURCES["vcf2genobaypass"]["runtime"],
        mem_mb=RESOURCES["vcf2genobaypass"]["mem_mb"],
        slurm_partition=RESOURCES["vcf2genobaypass"]["slurm_partition"]
    output:
        genobaypass=["results/subsets/{{sample}}.genobaypass.sub{}".format(i) for i in range(1, N_SUBS+1)],
        snpdet=["results/subsets/{{sample}}.snpdet.sub{}".format(i) for i in range(1, N_SUBS+1)]
    script: "scripts/poolfstat_subsample.R"

#########################################
## Detecting outlier loci with baypass ##
#########################################

## Notes: 
## (1) This part of the code is run directly from the command line in the Terminal.
## (2)-d0yij = 1/5 of min haploid pool size; use npilot = 100 for final analysis (see the Manual for detailed parameter descriptions).
## (3) baypass is not scaling linearly with the no of threads; analyze subsets of data on single threads in parallel; 
## pooldata.subset() function in poolfstat can be used.
## (4) Option 3: Contrast analysis is applied on the subsetted ACORN data!

## Option 1. Scanning the genome for differentiation (without covariate) 
# running core model for scanning the genome for differentiation using the XtX statistics
rule run_baypass_core:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}"
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: count_words_in_file(f"data/{wildcards.sample}_poolnames"),
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    resources:
        runtime=RESOURCES["run_baypass_core"]["runtime"],
        mem_mb=RESOURCES["run_baypass_core"]["mem_mb"],
        slurm_partition=RESOURCES["run_baypass_core"]["slurm_partition"]
    output:
        mat_omega = protected("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out"),
    shell:
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.sub} -poolsizefile {params.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -outprefix results/{wildcards.sample}_baypassSplitOut_core/core_{wildcards.i}
        """

rule generate_complementary_inputs2:
    input:
        envfactors=ENVFACTORS,
        poolsizes=POOLSIZES,
    resources:
        runtime=RESOURCES["generate_complementary_inputs"]["runtime"],
        mem_mb=RESOURCES["generate_complementary_inputs"]["mem_mb"],
        slurm_partition=RESOURCES["generate_complementary_inputs"]["slurm_partition"]
    output:
        efile="data/{sample}_efile_{envfactorname}",
    script: "scripts/generate_complementary_inputs.R"

## Option 2. Identifying SNPs associated with population covariate data
rule run_baypass_covariate:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
        efile="data/{sample}_efile_{envfactorname}",
        poolsizefile="data/{sample}_poolsizes",
    params:
        threads=1,
        npop=lambda wildcards: count_words_in_file(f"data/{wildcards.sample}_poolnames"),
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    resources:
        runtime=RESOURCES["run_baypass_covariate"]["runtime"],
        mem_mb=RESOURCES["run_baypass_covariate"]["mem_mb"],
        slurm_partition=RESOURCES["run_baypass_covariate"]["slurm_partition"]
    output:
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_covariate.std"),
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_DIC.out"),
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_summary_beta_params.out"),
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_summary_betai_reg.out"),
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_summary_pi_xtx.out"),
        protected("results/{sample}_baypassSplitOut_covariate/{envfactorname}/covariate_{i}_summary_yij_pij.out")
    shell:    
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.sub} -poolsizefile {input.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -omegafile {input.omegafile} \
        -efile {input.efile} \
        -outprefix results/{wildcards.sample}_baypassSplitOut_covariate/{wildcards.envfactorname}/covariate_{wildcards.i}
        """

# rule baypass_covariate_diagnostics:
#     input:
#         summary_betai_reg = "results/{sample}_baypassSplitOut_covariate/covariate_1_summary_betai_reg.out",
#         envfactor_names = "data/{sample}_efile_envfactor_names"
#     resources:
#         runtime=RESOURCES["baypass_covariate_diagnostics"]["runtime"],
#         mem_mb=RESOURCES["baypass_covariate_diagnostics"]["mem_mb"],
#         slurm_partition=RESOURCES["baypass_covariate_diagnostics"]["slurm_partition"]
#     output:
#         BFis = "results/{sample}_std_IS_model_BFis.png",
#         eBPis = "results/{sample}_std_IS_model_eBPis.png",
#         Beta_is = "results/{sample}_std_IS_model_Beta_is.png",
#     script: "scripts/model_diagnostics_covariate.R"

# rule concatenate_results_covariate:
#     input:
#         envfactor_names = "data/{sample}_efile_envfactor_names",
#         covariate = expand("results/{{sample}}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out", 
#         i=range(1, N_SUBS+1)),
#     params:
#         prefixcovariate = "results/{sample}_baypassSplitOut_covariate/covariate",
#         subs=N_SUBS,
#         snpdetprefix = "results/subsets/{sample}.snpdet.sub",
#         retrieve_c2= False
#     resources:
#         runtime=RESOURCES["concatenate_results_covariate"]["runtime"],
#         mem_mb=RESOURCES["concatenate_results_covariate"]["mem_mb"],
#         slurm_partition=RESOURCES["concatenate_results_covariate"]["slurm_partition"]       
#     output:
#         covariateresults = "results/{sample}_concatenated_res_covariate.csv",
#         manhattanplot = "results/{sample}_manhattanplot_covariate.png"
#     script: "scripts/concatenate_res.R"