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
dir_vcf = dir_data + "vcf/"
dir_id = dir_data + "id/"
dir_tss = dir_data + "tss/"
dir_pheno = dir_data + "pheno/"
dir_plink = dir_data + "plink/"
dir_gemma = dir_data + "gemma/"

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

rule all:
    input: dir_gemma + "gemma-top.txt"

rule gemma:
    input: dynamic(dir_gemma + "{gene}.assoc.txt")

rule prepare_tss:
    input: dynamic(dir_tss + "tss-{gene}.txt")

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

rule get_tss:
    input: dir_data + "counts-clean.txt"
    output: dir_tss + "tss.txt"
    shell: "Rscript scripts/get-tss.R {ensembl_archive} {input} > {output}"

rule get_tss_gene:
    input: dir_tss + "tss.txt"
    output: dynamic(dir_tss + "tss-{gene}.txt")
    params: outdir = dir_tss
    run:
        i = 0
        f = open(input[0], "r")
        for line in f:
            cols = line.strip().split("\t")
            gene = cols[0]
            assert gene[:4] == "ENSG", "Proper gene name"
            chrom = cols[1]
            if chrom not in chr_snps:
                continue
            tss = int(cols[2])
            start = tss - window
            if start < 0:
                start = str(0)
            else:
                start = str(start)
            end = tss + window
            end = str(end)
            fname = params.outdir + "tss-" + gene + ".txt"
            handle = open(fname, "w")
            cols_new = [gene, chrom, start, end] + cols[3:]
            handle.write("\t".join(cols_new) + "\n")
            handle.close()
            i += 1
            if i == 15145:
                break
        f.close()

rule subset_relatedness:
    input: relat = dir_data + "relatedness-matrix-all.txt",
           exp = dir_data + "counts-clean.txt"
    output: dir_data + "relatedness-matrix-sub.txt"
    shell: "Rscript scripts/subset-relatedness.R {input.exp} {input.relat} > {output}"

# Merge VCF files (Plink only accepts one as input)
rule merge_vcf:
    input: vcf = expand(dir_vcf  + "dox-hg38-chr{CHR}.vcf.gz", CHR = chr_snps)
    output: vcf = dir_vcf + "dox-hg38.vcf.gz"
    shell: "cat <(zcat {input.vcf[0]} | grep CHROM) \
                <(zcat {input.vcf} | grep -v '#') | \
            gzip -c > {output.vcf}"

rule create_filter_file_for_plink:
    input: exp = dir_data + "counts-clean.txt"
    output: filter = dir_plink + "filter-individuals.txt"
    shell: """
           head -n 1 {input.exp} | tr '\\t' '\\n' > /tmp/tmp-filter.txt
           paste /tmp/tmp-filter.txt /tmp/tmp-filter.txt > {output.filter}
           rm /tmp/tmp-filter.txt
           """

# https://www.cog-genomics.org/plink/1.9/input#vcf
# https://www.cog-genomics.org/plink/1.9/input#bed
rule convert_to_binary_plink:
    input: vcf = dir_vcf + "dox-hg38.vcf.gz",
           filter = dir_plink + "filter-individuals.txt"
    output: bed = dir_plink + "dox-hg38.bed",
            bim = dir_plink + "dox-hg38.bim",
            fam = dir_plink + "dox-hg38.fam"
    params: out = dir_plink + "dox-hg38"
    shell: "plink2 --vcf {input.vcf} \
           --double-id \
           --vcf-half-call missing \
           --biallelic-only \
           --make-bed \
           --keep {input.filter} \
           --indiv-sort natural \
           --out {params.out}"

# Format is:
# Column 1/2 = IID/FID
# Column 3 = Phenotype (i.e. gene expression levels)
rule pheno_file:
    input: tss = dir_tss + "tss-{gene}.txt",
           counts = dir_data + "counts-normalized.txt"
    output: dir_pheno + "{gene}.pheno"
    params: gene = "{gene}"
    shell: """
           paste \
           <(head -n 1 {input.counts} | tr '\\t' '\\n') \
           <(head -n 1 {input.counts} | tr '\\t' '\\n') \
           <(grep {params.gene} {input.counts} | tr '\\t' '\\n' | sed -e '1d') \
           > {output}
           """

rule plink_per_gene:
    input: bed = dir_plink + "dox-hg38.bed",
           bim = dir_plink + "dox-hg38.bim",
           fam = dir_plink + "dox-hg38.fam",
           tss = dir_tss + "tss-{gene}.txt",
           pheno = dir_pheno + "{gene}.pheno"
    output: bed = dir_plink + "{gene}.bed",
            bim = dir_plink + "{gene}.bim",
            fam = dir_plink + "{gene}.fam"
    params: prefix = dir_plink + "dox-hg38",
            out = dir_plink + "{gene}"
    shell: "plink2 --bfile {params.prefix} --make-bed \
           --pheno {input.pheno} \
           --chr `cut -f2 {input.tss}` \
           --from-bp `cut -f3 {input.tss}` \
           --to-bp  `cut -f4 {input.tss}` \
           --out {params.out}"

rule run_gemma:
    input: bed = dir_plink + "{gene}.bed",
           bim = dir_plink + "{gene}.bim",
           fam = dir_plink + "{gene}.fam",
           relat = dir_data + "relatedness-matrix-sub.txt"
    output: assoc = dir_gemma + "{gene}.assoc.txt",
            log = dir_gemma + "{gene}.log.txt"
    params: pre_plink = dir_plink + "{gene}",
            pre_gemma = "{gene}",
            outdir = dir_gemma
    shell: """
           gemma -bfile {params.pre_plink} \
           -k {input.relat} -km 2 \
           -lmm 2 \
           -o {params.pre_gemma}

           mv output/{params.pre_gemma}* {params.outdir}
           """

rule parse_gemma:
    input: assoc = dir_gemma + "{gene}.assoc.txt"
    output: dir_gemma + "{gene}.top.txt"
    params: gene = "{gene}"
    shell: "Rscript scripts/parse-gemma.R {params.gene} {input} > {output}"

rule combine_gemma:
    input: top = dynamic(dir_gemma + "{gene}.top.txt")
    output: dir_gemma + "gemma-top.txt"
    shell: "cat <(head -n 1 {input.top[0]}) \
                <(cat {input} | grep -v mle) > {output}"

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
