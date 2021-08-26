import os.path
import subprocess
import json
import csv

configfile: "SampleList.test"
current_dir = os.getcwd() + "/"

sequence_fasta = current_dir+ "path/to/reference.fasta"

# Enable or disable individual rules.
run_rmdup_bams = True
run_fastqc_original_samples = True
run_fastqc_rmdup_bams = True
run_mapdamage = True
run_qualimap = True
run_haplogrep = False
run_haplocheck = False
run_ind_mapping_report = True
run_agr_mapping_report = True
run_multiqc_report = False
run_bam_to_fasta = False

# folder paths
samples = current_dir + "Data/samples/"
merged_trimmed = current_dir + "Data/merged_trimmed/"
aligned_mito = current_dir + "Data/aligned_mito/"
bam_folder = current_dir + "Data/bams/"
filtered = current_dir + "Data/filtered_bams/"
q30 = current_dir + "Data/q30/"
sort_to_left_dir = current_dir + "Data/sort/"
rem_dups = current_dir + "Data/rmdup/"
min35_folder = current_dir + "Data/min_35/"
clipped_folder = current_dir + "Data/clipped/"
fqc_folder = current_dir + "Data/fqc_samples/"  
fqctrim_folder = current_dir + "Data/fqc_trimmed/" 
md_folder = current_dir + "Data/mapdamage_reports/"
hg_folder = current_dir + "Data/haplogrep/"
hc_folder = current_dir + "Data/haplocheck/"
rmdupf_folder = current_dir + "Data/rmdup_fasta/"
qualimap_folder = current_dir + "Data/qualimap/"
mapping_folder = current_dir + "Data/mapping_reports/"

# app paths
# leeHom = current_dir + "apps/leeHom/src/leeHomMulti"
# schmutzi_endocaller = current_dir + "apps/schmutzi/src/endoCaller"
leeHom = "leeHom"
schmutzi_endocaller = current_dir + "apps/schmutzi/src/endoCaller"
realign_sam_prg = current_dir + "apps/realign.jar" 

g_threads = 12

# gets all the files for our pipeline
def get_files():
    files = config["samples"]
    ids = []
    for file in files:
        ids.append(file["id"])
    return ids

def target_files():
  target_list = []

  if run_rmdup_bams:
    target_list += expand(rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam", sample=get_files())

  if run_fastqc_original_samples:
    target_list += expand(fqc_folder + "{sample}_{part}_001_fastqc.html", sample=get_files(), part=["R1", "R2"])

  if run_fastqc_rmdup_bams:
    target_list += expand(fqctrim_folder + "{sample}_trimmed_merged_fastqc.html", sample=get_files())

  if run_mapdamage:
    target_list += expand(md_folder + "{sample}/Fragmisincorporation_plot.pdf", sample=get_files())

  if run_qualimap:
    target_list += expand(qualimap_folder + "{sample}/", sample=get_files())

  if run_haplogrep:
    target_list += expand(hg_folder + "{sample}_trimmed_mapped_realigned_f4_q30_sort_rmdup_fasta_haplogrep.out", sample=get_files())

  if run_haplocheck:
    target_list += expand(hc_folder + "{sample}_haplochecked", sample=get_files())

  if run_ind_mapping_report:
    target_list += expand(mapping_folder + "{sample}_mapping_report.json", sample=get_files())

  if run_agr_mapping_report:
    target_list.append(mapping_folder + "aggregated_mapping_report_mqc.tsv")

  if run_multiqc_report:
    target_list.append(current_dir + "Data/multiqc/")

  if run_bam_to_fasta:
    target_list += expand(rmdupf_folder + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.fasta", sample=get_files())

  return target_list

def mk_dirs():
  #mk_dirs if they don't exist
  dir_list = [
  merged_trimmed,
  aligned_mito,
  bam_folder,
  filtered,
  q30,
  sort_to_left_dir,
  rem_dups,
  min35_folder,
  clipped_folder,
  fqc_folder,
  fqctrim_folder,
  md_folder,
  hg_folder,
  hc_folder,
  rmdupf_folder,
  qualimap_folder,
  mapping_folder,
  ]

  for dr in dir_list:
    subprocess.check_output("mkdir -p " + dr, shell=True, text=True)

def parse_leeHom(text):
    text = text.strip().split(";")
    parsed = {
        "total": text[0].split(" ")[-1],
        "merged_trim": text[1].split(" ")[-1],
        "merged_overlap": text[2].split(" ")[-1],
        "kept": text[3].split(" ")[-1],
        "trimmed_sr": text[4].split(" ")[-1],
        "adapter_dimers_chimeras": text[5].split(" ")[-1],
        "failed_key": text[6].split(" ")[-1],
        "umi_problems": text[7].split(" ")[-1],
    }
    return parsed


def parse_markdup(text):
    text = text.strip().split("\n")
    parsed = {}
    for line in text:
        splitter = line.split(":")
        parsed[splitter[0].lower().replace(" ", "_")] = splitter[1].strip()
    return parsed

def parse_qualimap(text):
  text = text.strip().splitlines()
  parsed = {
    "mean": -1,
    "std": -1,
    "reads": -1,
  }

  for line in text:
    if 'mean coverageData' in line:
        parsed["mean"] = line.split(" ")[-1]
    if 'std coverageData' in line:
        parsed["std"] = line.split(" ")[-1]
    if 'number of reads' in line:
        parsed["reads"] = line.split(" ")[-1]

  return parsed


def parse_it(sample, trim_report, q_bam, qmapdir):
    trim_report = open(trim_report, "r")
    # mdup_report = open(mdup_report, "r")

    quality_bam = subprocess.check_output("samtools view -c " + q_bam, shell=True, text=True)

    if os.path.exists(qmapdir + "/genome_results.txt"):
      qmap_report = open(qmapdir + "/genome_results.txt")
      qmap_dict = parse_qualimap(qmap_report.read())
    else:
      qmap_dict = parse_qualimap("")

    leehom_dict = parse_leeHom(trim_report.read())
    #markdup_dict = parse_markdup(mdup_report.read())

    trim_report.close()
    #mdup_report.close()

    parsed = {
        "sample": sample,
        "leeHom": leehom_dict,
        # "markdup": markdup_dict,
        "q30_bam": quality_bam.strip(),
        "qualimap": qmap_dict,
    }
    return parsed

rule target:
    input:
      target_files()

rule trim_merge:
    input:
         fastq1 = samples + "{sample}_R1_001.fastq.gz",
         fastq2 = samples + "{sample}_R2_001.fastq.gz"
    output:
          fastq = merged_trimmed + "{sample}_trimmed_merged.fq.gz",
          report = merged_trimmed + "{sample}_report.txt"

    params:
        threads = g_threads,
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        second_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
        sample = merged_trimmed + "{sample}_trimmed_merged",
        leeHom = leeHom

    shell:
        "{params.leeHom} -t {params.threads} -f {params.forward_adapter} -s {params.second_adapter} "
              "--ancientdna -fq1 {input.fastq1} -fq2 {input.fastq2} -fqo {params.sample} 2> {output.report}"


rule index_ref:
    input:
        ref = sequence_fasta
    output:
        "sequence.fasta.amb",
         "sequence.fasta.ann",
         "sequence.fasta.bwt",
         "sequence.fasta.pac",
         "sequence.fasta.sa"
    run:
        shell("samtools faidx {input.ref}")

rule alignment:
    input:
         fastq = merged_trimmed + "{sample}_trimmed_merged.fq.gz",
         ref = sequence_fasta

    output:
          sam = aligned_mito + "{sample}_trimmed_merged_mapped.sam"
    run:
        shell("bwa mem -t {g_threads} {input.ref} {input.fastq} > {output.sam}")


rule awk_rm_softclip:
  input:
    sam = aligned_mito + "{sample}_trimmed_merged_mapped.sam",
    ref = sequence_fasta

  output:
    clipped = clipped_folder + "{sample}_trimmed_merged_mapped_clipped.sam"

  run:
    shell("samclip --ref {input.ref} {input.sam} > {output.clipped}")

  # shell:
  #   """
  #     awk 'BEGIN {{OFS="t"}} {{split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {{$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}}; if (C[length(C)]=="S") {{L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }}; gsub(/[0-9]*S/,"",$6); print}}' {input.sam} > {output.clipped}
  #   """

# rule realign_sam_to_bam:
#   input:
#     sam = clipped_folder + "{sample}_trimmed_merged_mapped_clipped.sam",
#     ref = sequence_fasta

#   output:
#     bam = bam_folder + "{sample}_trimmed_merged_mapped_realigned.bam"

#   params:
#     elongate = 500

  # run:
  #   shell("java -jar {realign_sam_prg} -e {params.elongate} -i {input.sam} -r {input.ref};mv {aligned_mito}{wildcards.sample}_trimmed_merged_mapped_realigned.bam {bam_folder}")

rule sam_to_bam:
    input:
          sam = clipped_folder + "{sample}_trimmed_merged_mapped_clipped.sam",
    output:
          bam = bam_folder + "{sample}_trimmed_merged_mapped_realigned.bam"
    run:
        shell("samtools view -@ {g_threads} -bSh {input.sam} > {output.bam}")

rule read_filter_f4:
    input:
         bam = bam_folder + "{sample}_trimmed_merged_mapped_realigned.bam"
    output:
          filtered = filtered + "{sample}_trimmed_merged_mapped_realigned_f4.bam"
    run:
        shell("samtools view -@ {g_threads} -bh -F4 {input.bam} > {output.filtered}")

rule min_35:
    input:
        filtered = filtered + "{sample}_trimmed_merged_mapped_realigned_f4.bam"
    output:
        filtered = min35_folder + "{sample}_trimmed_merged_mapped_realigned_f4_min35.bam"
    run:
        shell("samtools view -@ {g_threads} -h {input.filtered} | awk 'length($10) > 30 || $1 ~ /^@/' | samtools view -bS - > {output.filtered}")

rule read_filter_q30:
    input:
         filtered = min35_folder + "{sample}_trimmed_merged_mapped_realigned_f4_min35.bam"
    output:
          quality = q30 + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30.bam"
    params:
        filter_amount = 30
    run:
        shell('samtools view -@ {g_threads} -bh -q {params.filter_amount} {input.filtered} > {output.quality}')

rule sort_to_left:
    input:
         quality = q30 + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30.bam"
    output:
          sort_to_left =  sort_to_left_dir + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort.bam"
    run:
        shell("samtools sort -@ {g_threads} {input.quality} -o {output.sort_to_left}")

rule remove_duplicates:
    input:
         sort_to_left =  sort_to_left_dir + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort.bam"
    output:
          rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam",
    run:
        shell("bam-rmdup -o {output.rmdupss} {input.sort_to_left} --verbose")

rule bam_to_fasta:
  input:
    rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam",
    ref = sequence_fasta,

  output:
    fasta = rmdupf_folder + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.fasta"

  shell:
    "{schmutzi_endocaller} {input.rmdupss} {input.ref} {output.fasta}"

rule index_rmdups:
    input:
         rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam"
    output:
        indexed_rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam.bai"
    run:
        shell("samtools index -@ {g_threads} {input.rmdupss}")

# reporting shizzle
rule fqc_sample:
  input:
    fastq = samples + "{sample}_{part}_001.fastq.gz"

  output:
    folder = fqc_folder + "{sample}_{part}_001_fastqc.html"

  wildcard_constraints:
    part="R1|R2"

  shell:
    "fastqc -o {fqc_folder} {input.fastq}"

rule fqc_trimmed:
  input:
    fastq = merged_trimmed + "{sample}_trimmed_merged.fq.gz"

  output:
    fastqc = fqctrim_folder + "{sample}_trimmed_merged_fastqc.html"

  shell:
    "fastqc -o {fqctrim_folder} {input.fastq}"

rule mapDamage:
    input:
         rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam",
         indexed_rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam.bai"

    output:
         sample = md_folder + "{sample}/Fragmisincorporation_plot.pdf"

    params:
          settings = "--rescale",
          outfolder = md_folder + "{sample}"

    shell:
        "mapDamage -i {input.rmdupss} -r {sequence_fasta} -d {params.outfolder} {params.settings}"

rule qualiMap:
  input:
      rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam"

  output:
    results = directory(qualimap_folder + "{sample}/"),

  shell:
    "qualimap bamqc -bam {input} -outdir {output.results} -outformat HTML || mkdir -p {output.results}"

rule mapping_report:
  input:
    quality = q30 + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30.bam",
    rmdup = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam",
    trim_report = merged_trimmed + "{sample}_report.txt",
    qualimap = qualimap_folder + "{sample}/"

  output:
    mapping_report = mapping_folder + "{sample}_mapping_report.json"

  run:
    f = open(output.mapping_report, "w")
    f.write(json.dumps(parse_it(wildcards.sample, input.trim_report, input.quality, input.qualimap), indent=4, sort_keys=True))
    f.close()

rule agr_mapping_report:
  input:
    expand(mapping_folder + "{sample}_mapping_report.json", sample=get_files())

  output:
    tsv = mapping_folder + "aggregated_mapping_report_mqc.tsv"

  run:
    aggr_report = open(output.tsv, "w")
    tsv = csv.writer(aggr_report,  delimiter='\t')
    tsv.writerow(['sample', 'lh_total', 'lh_mergTrim', 'lh_mergOver', 'lh_kept', 'lh_chimer', 'q30_bam',
              'qm_cov_mean', 'qm_cov_std', 'qm_reads'])

    for f in input:
      report = open(f, "r")
      data = json.load(report)
      tsv.writerow([
    data['sample'],data['leeHom']['total'], data['leeHom']['merged_trim'], data['leeHom']['merged_overlap'],
    data['leeHom']['kept'], data['leeHom']['adapter_dimers_chimeras'], data['q30_bam'],
    data['qualimap']['mean'], data['qualimap']['std'], data['qualimap']['reads']])
      report.close()
    aggr_report.close()

rule multiqc:
  output:
    directory(current_dir + "Data/multiqc/")
  shell:
    "multiqc -o {output} Data/"


rule haplogrep:
    input:
         fasta = rmdupf_folder + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.fasta"
    output:
          haplogrep = hg_folder + "{sample}_trimmed_mapped_realigned_f4_q30_sort_rmdup_fasta_haplogrep.out"
    shell:
        "haplogrep --in {input.fasta} --format fasta --out {output.haplogrep}"

rule haplocheck:
  input:
    rmdupss = rem_dups + "{sample}_trimmed_merged_mapped_realigned_f4_min35_q30_sort_rmdup.bam"

  output:
    haplocheck = directory(hc_folder + "{sample}/")

  shell:
    "cloudgene run haplocheck@1.3.2 --files {input.rmdupss} --format bam --output {output.haplocheck}"


mk_dirs()

























