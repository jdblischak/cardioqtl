# Snakefile
#
# This Snakefile runs the eQTL analysis starting from the fastq files.
#
# To configure the paths to data files and other settings, edit
# config.yaml.
#
# To configure job submission settings for your cluster, edit
# cluster.json and submit-snakemake.sh.
#
# To run on RCC Midway2 use `bash submit-snakemake.sh`

import glob
import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

configfile: "config.yaml"

# Specify Ensembl release for genome sequence and annotation
ensembl_archive = config["ensembl_archive"]
ensembl_rel = config["ensembl_rel"]
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-" + \
              str(ensembl_rel) + \
              "/fasta/homo_sapiens/dna/"
ensembl_exons = "exons-ensembl-release-" + str(ensembl_rel) + ".saf"
ensembl_genome = config["ensembl_genome"]

# Paths to data (must end with forward slash)
dir_proj = config["dir_proj"]
dir_data = dir_proj + "data/"
dir_fq = dir_data + "fastq/"
scratch = config["scratch"]
dir_genome = scratch + "genome-ensembl-release-" + str(ensembl_rel) + "/"
dir_fq_tmp = scratch + "cardioqtl-fastq/"
dir_bam = scratch + "cardioqtl-bam/"
dir_counts = scratch + "cardioqtl-counts/"
dir_pca = dir_data + "pca/"
dir_vcf = dir_data + "vcf/"
dir_id = dir_data + "id/"
dir_pheno = dir_data + "pheno/"
dir_plink = dir_data + "plink/"
dir_gemma = dir_data + "gemma/"
dir_gwas = dir_data + "gwas/"

assert os.path.exists(dir_proj), "Project directory exists"
assert os.path.exists(scratch), "Scratch directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)

# Specify chromosomes
#
# For quantifying gene expression
chr_genes = config["chr_genes"]
# For mapping eQTLs
chr_snps = config["chr_snps"]

# Specify cis-window around TSS for testing eQTLs (TSS +/- window)
window = config["window"]

# Number of PCs to regress in linear model for testing eQTLs
n_pcs = config["n_pcs"]

# GWAS phenotypes
phenos = config["phenos"]

# Input samples ----------------------------------------------------------------

samples = set(glob_wildcards(dir_fq + "{samples}-H{flow_cell}-l{lane,[1-8]}.fastq.gz").samples)
# For testing:
#files_fq = glob.glob(dir_fq + "*fastq.gz")
#samples = [os.path.basename(x).split("-")[0] for x in files_fq]

# Constrain sample wildcard to always be all numbers. Necessary to
# resolve some of the complex rules.
wildcard_constraints: sample = "\d+"

# Functions --------------------------------------------------------------------

# Find all fastq.gz files for a given sample.
# Inspired by this post on the Snakemake Google Group:
# https://groups.google.com/forum/#!searchin/snakemake/multiple$20input$20files%7Csort:relevance/snakemake/bpTnr7FgDuQ/ybacyom6BQAJ
def merge_fastq(wc):
    unknowns = glob_wildcards(dir_fq +
               "{s}-H{{flow_cell}}-l{{lane,[1-8]}}.fastq.gz".format(s = wc.sample))
    files = expand(dir_fq + "{s}-H{{flow_cell}}-l{{lane}}.fastq.gz".format(s = wc.sample),
                   zip, flow_cell = unknowns.flow_cell, lane = unknowns.lane)
    return files

# Targets ----------------------------------------------------------------------

rule run_gemma:
    input: expand(dir_gemma + "top-pca-{pc}.txt", pc = [x for x in range(n_pcs + 1)] + [15, 20, 30, 40])

rule prepare_gemma:
    input: dir_data + "tss.txt",
           dir_plink + "cardioqtl.bed",
           dir_data + "counts-normalized.txt",
           expand(dir_pca + "pca-{pc}.txt", pc = [x for x in range(n_pcs + 1)])

rule run_pca:
    input: expand(dir_pca + "pca-{pc}.txt", pc = [x for x in range(n_pcs + 1)])

rule counts_for_gemma:
    input: dir_data + "counts-normalized.txt"

rule run_featurecounts:
    input: expand(dir_counts + "{sample}.genecounts.txt", sample = samples)

rule prepare_featurecounts:
    input: dir_genome + ensembl_exons

rule run_verifyBamID:
    input: dir_data + "verify.txt"

rule run_subjunc:
    input: expand(dir_bam + "{sample}-sort.bam.bai", sample = samples)

rule prepare_subjunc:
    input: dir_genome + ensembl_genome + ".reads"

# Quanitify expression with Subjunc/featureCounts ------------------------------

rule download_genome:
    output: dir_genome + "Homo_sapiens." + ensembl_genome + \
            ".dna_sm.chromosome.{chr}.fa.gz"
    params: chr = "{chr}", build = ensembl_genome
    shell: "wget -O {output} {ensembl_ftp}Homo_sapiens.{params.build}.dna_sm.chromosome.{params.chr}.fa.gz"

rule unzip_chromosome_fasta:
    input: dir_genome + "Homo_sapiens." + ensembl_genome + \
           ".dna_sm.chromosome.{chr}.fa.gz"
    output: temp(dir_genome + "Homo_sapiens." + ensembl_genome + \
                 ".dna_sm.chromosome.{chr}.fa")
    shell: "zcat {input} > {output}"

rule subread_index:
    input: expand(dir_genome + "Homo_sapiens." + ensembl_genome + \
                  ".dna_sm.chromosome.{chr}.fa", \
                  chr = chr_genes)
    output: dir_genome + ensembl_genome + ".reads"
    params: prefix = dir_genome + ensembl_genome
    shell: "subread-buildindex -o {params.prefix} {input}"

rule combine_fastq:
    input: merge_fastq
    output: temp(dir_fq_tmp + "{sample}.fastq")
    shell: "zcat {input} > {output}"

rule subjunc:
    input: read = dir_fq_tmp + "{sample}.fastq",
           index = dir_genome + ensembl_genome + ".reads"
    output: temp(dir_bam + "{sample}.bam")
    params: prefix = dir_genome + ensembl_genome
    threads: 8
    shell: "subjunc -i {params.prefix} -r {input.read} -u -T {threads} > {output}"

rule sort_bam:
    input: dir_bam + "{sample}.bam"
    output: dir_bam + "{sample}-sort.bam"
    shell: "samtools sort -o {output} {input}"

rule index_bam:
    input: dir_bam + "{sample}-sort.bam"
    output: dir_bam + "{sample}-sort.bam.bai"
    shell: "samtools index {input}"

rule create_exons_saf:
    output: dir_genome + ensembl_exons
    shell: "Rscript scripts/create-exons.R {ensembl_archive} > {output}"

rule feauturecounts:
    input: bam = dir_bam + "{sample}-sort.bam",
           bai = dir_bam + "{sample}-sort.bam.bai",
           exons = dir_genome + ensembl_exons
    output: dir_counts + "{sample}.genecounts.txt"
    threads: 8
    shell: "featureCounts -a {input.exons} -F SAF -T {threads} -o {output} {input.bam}"

rule combine_featurecounts:
    input: expand(dir_counts + "{sample}.genecounts.txt", sample = samples)
    output: dir_data + "counts-subread.txt"
    shell: "Rscript scripts/combine-featurecounts.R {input} > {output}"

rule clean_counts:
    input: dir_data + "counts-subread.txt"
    output: dir_data + "counts-clean.txt"
    shell: "Rscript scripts/clean-counts.R {input} > {output}"

rule normalize_counts:
    input: dir_data + "counts-clean.txt"
    output: dir_data + "counts-normalized.txt"
    shell: "Rscript scripts/normalize-counts.R {input} > {output}"

# Map eQTLs with GEMMA ---------------------------------------------------------

localrules: create_filter_file_for_plink

rule pca:
    input: dir_data + "counts-normalized.txt"
    output: dir_pca + "pca-{pc}.txt"
    params: pc = "{pc}"
    shell: "Rscript scripts/compute-pca.R {params.pc} {input} > {output}"

rule get_tss:
    input: dir_data + "counts-clean.txt"
    output: dir_data + "tss.txt"
    shell: "Rscript scripts/get-tss.R {ensembl_archive} {input} > {output}"

rule create_filter_file_for_plink:
    input: exp = dir_data + "counts-clean.txt"
    output: filter = dir_plink + "filter-individuals.txt"
    run:
        in_handle = open(input.exp, "r")
        header = in_handle.readline()
        individuals = header.strip().split("\t")
        out_handle = open(output.filter, "w")
        for i in individuals:
            line = "HUTTERITES\t%s\n"%(i)
            out_handle.write(line)
        out_handle.close()
        in_handle.close()

# https://www.cog-genomics.org/plink/1.9/input#bed
rule format_plink_eqtl:
    input: bed = dir_plink + "hutt.imputed.rename.bed",
           bim = dir_plink + "hutt.imputed.rename.bim",
           fam = dir_plink + "hutt.imputed.rename.fam",
           filter = dir_plink + "filter-individuals.txt"
    output: bed = dir_plink + "cardioqtl.bed",
            bim = dir_plink + "cardioqtl.bim",
            fam = dir_plink + "cardioqtl.fam"
    params: prefix_in = dir_plink + "hutt.imputed.rename",
            prefix_out = dir_plink + "cardioqtl"
    shell: "plink2 -bfile {params.prefix_in} \
           --make-bed \
           --keep {input.filter} \
           --indiv-sort natural \
           --out {params.prefix_out}"

rule gemma:
    input: counts = dir_data + "counts-normalized.txt",
           tss = dir_data + "tss.txt",
           pca = dir_pca + "pca-{pc}.txt",
           relat = dir_data + "relatedness-matrix-all.txt",
           bed = dir_plink + "cardioqtl.bed",
           bim = dir_plink + "cardioqtl.bim",
           fam = dir_plink + "cardioqtl.fam"
    output: top = dir_gemma + "top-pca-{pc}.txt"
    params: prefix_plink = dir_plink + "cardioqtl"
    shell: "Rscript scripts/run-gemma.R {input.counts} {input.tss} {input.pca} {input.relat} {params.prefix_plink} {window} {dir_pheno} {dir_plink} {dir_gemma}"

# Verify identity with verifyBamID ---------------------------------------------

# Convert exons in SAF format to BED format. Duplicate exons are maintained.
rule convert_to_bed:
    input: saf = dir_genome + ensembl_exons
    output: bed = dir_genome + "exons.bed"
    run:
        saf = open(input.saf, "r")
        bed = open(output.bed, "w")
        # Discard header
        saf.readline()

        for line in saf:
            cols = line.strip().split("\t")
            id = cols[0]
            chr = cols[1]
            start = str(int(cols[2]) - 1)
            end = cols[3]
            strand = cols[4]
            entry = "%s\t%s\t%s\t%s\t%s\t%s\n"%(chr, start, end, id, 0, strand)
            bed.write(entry)

        saf.close()
        bed.close()


# Add "chr" to chromosome name
rule add_chr_to_vcf:
    input: dir_vcf + "dox-hg38-chr{CHR}.vcf.gz"
    output: temp(dir_vcf + "dox-hg38-chr{CHR}-chr.vcf")
    run:
        import gzip

        vcf_in = gzip.open(input[0], "rt")
        vcf_out = open(output[0], "w")
        for line in vcf_in:
            if line[0] == "#":
                vcf_out.write(line)
            else:
                vcf_out.write("chr" + line)

        vcf_in.close()
        vcf_out.close()

# Select only those SNPs in annotated exons of protein-coding genes
rule select_exonic_snps:
    input: vcf = dir_vcf + "dox-hg38-chr{CHR}-chr.vcf",
           exons = dir_genome + "exons.bed"
    output: vcf = temp(dir_vcf + "dox-hg38-chr{CHR}-exons.vcf")
    shell: "bedtools intersect -a {input.vcf} -b {input.exons} -u -header > {output.vcf}"

# Combine exonic SNPs into one file
rule combine_snps:
    input: vcf = expand(dir_vcf  + "dox-hg38-chr{CHR}-exons.vcf", CHR = chr_snps[:-3])
    output: vcf = dir_vcf + "dox-hg38-exons.vcf"
    shell: "cat <(grep CHROM {input.vcf[0]}) <(cat {input.vcf} | grep -v '#') > {output.vcf}"

# Run verifyBamID to obtain the best individual match for the BAM file
rule verify_bam:
    input: vcf =  dir_vcf + "dox-hg38-exons.vcf",
           bam = dir_bam + "{sample}-sort.bam",
           index = dir_bam + "{sample}-sort.bam.bai"
    output: selfSM = dir_id + "{sample}.selfSM",
            bestSM = dir_id + "{sample}.bestSM",
            depthSM = dir_id + "{sample}.depthSM"
    params: prefix = dir_id + "{sample}",
            individual = "{sample}"
    shell: "verifyBamID --vcf {input.vcf} --bam {input.bam} --best --ignoreRG --smID {params.individual} --out {params.prefix}"

rule parse_verify:
    input: selfSM = dir_id + "{sample}.selfSM",
           bestSM = dir_id + "{sample}.bestSM",
           depthSM = dir_id + "{sample}.depthSM"
    output: dir_id + "{sample}-results.txt"
    run:
        selfSM = open(input.selfSM, "rt")
        bestSM = open(input.bestSM, "rt")
        depthSM = open(input.depthSM, "rt")
        results = open(output[0], "w")

        for line in selfSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[0] == "#SEQ_ID" and cols[2] == "CHIP_ID" and \
                       cols[3] == "#SNPS" and cols[4] == "#READS" and \
                       cols[5] == "AVG_DP" and cols[6] == "FREEMIX" and \
                       cols[11] == "CHIPMIX", "selfSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                seq_id = cols[0]
                chip_id = cols[2]
                assert chip_id == seq_id, "selfSM compares to self"
                snps = cols[3]
                reads = cols[4]
                avg_dp = cols[5]
                selfSM_freemix = cols[6]
                selfSM_chipmix = cols[11]

        for line in bestSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[0] == "#SEQ_ID" and cols[2] == "CHIP_ID" and \
                       cols[3] == "#SNPS" and cols[4] == "#READS" and \
                       cols[5] == "AVG_DP" and cols[6] == "FREEMIX" and \
                       cols[11] == "CHIPMIX", "bestSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                seq_id = cols[0]
                bestSM_chip_id = cols[2]
                match = str(bestSM_chip_id == seq_id).upper()
                bestSM_freemix = cols[6]
                bestSM_chipmix = cols[11]

        # Report the number of SNPs that had more than a minimum read depth
        depth_min = 10
        snps_w_min = 0
        for line in depthSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[1] == "DEPTH" and cols[2] == "#SNPs", \
                       "depthSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                depth = int(cols[1])
                n_snps = int(cols[2])
                if depth >= depth_min:
                    snps_w_min = snps_w_min + n_snps

        out_header = ["id", "best", "match",
                      "self_freemix", "self_chipmix",
                      "best_freemix", "best_chipmix",
                      "snps", "reads", "avg_dp",
                      "min_dp", "snps_w_min"]
        out_cols = [seq_id, bestSM_chip_id, match,
                    selfSM_freemix, selfSM_chipmix,
                    bestSM_freemix, bestSM_chipmix,
                    snps, reads, avg_dp,
                    str(depth_min), str(snps_w_min)]
        results.write("\t".join(out_header) + "\n")
        results.write("\t".join(out_cols) + "\n")

        selfSM.close()
        bestSM.close()
        depthSM.close()
        results.close()

rule combine_verify:
    input: expand(dir_id + "{sample}-results.txt", sample = samples)
    output: dir_data + "verify.txt"
    shell:
        "head -n 1 {input[0]} > {output};"
        "cat {input} | grep -v \"id\" | sort -k1n >> {output}"

# Map phenotype QTLs with GEMMA ------------------------------------------------

rule run_gwas:
    input: expand(dir_gwas + "{pheno}/gwas-{pheno}.html", pheno = phenos)

# Remove any eQTL individuals from GWAS and sort by IID
rule filter_gwas:
    input: raw = dir_gwas + "{pheno}/raw-{pheno}.txt",
           filter = dir_plink + "filter-individuals.txt"
    output: clean = dir_gwas + "{pheno}/clean-{pheno}.txt"
    shell: "Rscript scripts/filter-gwas.R {input.raw} {input.filter} > {output.clean}"

# Split GWAS data into 3 files (none have header rows):
#
# 1. keep: FID/IID to filter individuals with PLINK
# 2. pheno: FID/IID/phenotype to add phenotype with PLINK
# 3. cov: covariates to model covariates with GEMMA
rule split_gwas:
    input: clean = dir_gwas + "{pheno}/clean-{pheno}.txt"
    output: keep = dir_gwas + "{pheno}/keep-{pheno}.txt",
            pheno = dir_gwas + "{pheno}/pheno-{pheno}.txt",
            cov = dir_gwas + "{pheno}/cov-{pheno}.txt"
    run:
        clean = open(input.clean, "r")
        header = clean.readline()
        assert header[:6] == "findiv", \
            "IID is in column 1"
        keep = open(output.keep, "w")
        pheno = open(output.pheno, "w")
        cov = open(output.cov, "w")
        for line in clean:
            cols = line.strip("\n").split("\t")
            iid = cols[0]
            trait = cols[1]
            covariates = cols[2:]
            keep.write("HUTTERITES\t%s\n"%(iid))
            pheno.write("HUTTERITES\t%s\t%s\n"%(iid, trait))
            cov.write("\t".join(covariates) + "\n")
        clean.close()
        keep.close()
        pheno.close()
        cov.close()


# Create PLINK file for GWAS
rule create_gwas_plink:
    input: bed = dir_plink + "hutt.imputed.rename.bed",
           bim = dir_plink + "hutt.imputed.rename.bim",
           fam = dir_plink + "hutt.imputed.rename.fam",
           keep = dir_gwas + "{pheno}/keep-{pheno}.txt",
           pheno = dir_gwas + "{pheno}/pheno-{pheno}.txt"
    output: bed = dir_gwas + "{pheno}/plink-{pheno}.bed",
            bim = dir_gwas + "{pheno}/plink-{pheno}.bim",
            fam = dir_gwas + "{pheno}/plink-{pheno}.fam"
    params: prefix_in = dir_plink + "hutt.imputed.rename",
            prefix_out = dir_gwas + "{pheno}/plink-{pheno}"
    shell: "plink2 --bfile {params.prefix_in} --make-bed \
           --keep {input.keep} --indiv-sort natural \
           --pheno {input.pheno}  --out {params.prefix_out}"

# Run GEMMA on each phenotype. Minor allele frequency cutoff relaxed
# from default of 0.01 to 0.05 for consistency with Cusanovich et al.,
# 2016.
rule gwas_gemma:
    input: bed = dir_gwas + "{pheno}/plink-{pheno}.bed",
           bim = dir_gwas + "{pheno}/plink-{pheno}.bim",
           fam = dir_gwas + "{pheno}/plink-{pheno}.fam",
           cov = dir_gwas + "{pheno}/cov-{pheno}.txt",
           relat = dir_data + "relatedness-matrix-all.txt"
    output: assoc = dir_gwas + "{pheno}/gemma-{pheno}.assoc.txt",
            log = dir_gwas + "{pheno}/gemma-{pheno}.log.txt"
    params: prefix_in = dir_gwas + "{pheno}/plink-{pheno}",
            prefix_out = "gemma-{pheno}",
            outdir = dir_gwas + "{pheno}"
    shell: "gemma -bfile {params.prefix_in} \
           -k {input.relat} -km 2 \
           -c {input.cov} -lmm 4 \
           -maf 0.05 \
           -o {params.prefix_out} \
           ; mv output/{params.prefix_out}* {params.outdir}"

rule gwas_report:
    input: assoc = dir_gwas + "{pheno}/gemma-{pheno}.assoc.txt",
           template = "scratch/gwas.Rmd"
    output: dir_gwas + "{pheno}/gwas-{pheno}.html"
    shell: "Rscript scripts/report-gwas.R {input.template} {wildcards.pheno} {output}"
