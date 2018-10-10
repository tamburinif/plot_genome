# Fiona Tamburini
# Sept 2018

# Finds contigs aligning to a reference genome, uses nucmer to align/scaffold
# contigs against the reference and rotate genome to same starting coordinates

localrules: nucmer, filter, show_coords, scaffold_1, scaffold_2, stitch, done, combine

import subprocess
from Bio.Seq import Seq
from Bio import SeqIO

# read in list of samples
samples = []
try:
	f = open(config['samples'])
	for line in f:
		sample = line.rstrip('\n')
		samples += [sample]
except (IOError, OSError) as e:
	print("File {f} does not exist".format(f=config['samples']))
else:
	f.close()

ref = config['ref']
refname = config['refname']

# index blast database
rule makeblastdb:
    input: config['ref']
    output:
        "{ref}.nhr".format(ref=config['ref']),
        "{ref}.nin".format(ref=config['ref']),
        "{ref}.nsq".format(ref=config['ref'])
    resources:
        mem=2,
        time=1
    shell:
        "makeblastdb -in {input} -dbtype nucl"

# blast align contigs to reference
rule blast:
    input:
        ref=config['ref'],
        idx=rules.makeblastdb.output,
        assembly=config['sample_dir'] + "{sample}.fasta"
    output: "blast/{sample}.blast.out"
    resources:
        mem=8,
        time=2
    shell:
        "blastn -db {input.ref} -query {input.assembly} -outfmt 6 -out {output}"

rule blast_all:
    input: expand("blast/{sample}.blast.out", sample=samples)
    output: "blast_done"
    resources:
        mem=1,
        time=1
    shell:
        "touch {output}"

# get a list of contigs over 1kb that align to ref genome
# changed - aligned length must be over 1kb
rule contig_list:
    input: rules.blast.output
    output: "contig_list/{sample}.crass.list"
    resources:
        mem=1,
        time=1
    shell:
        "awk '$4 > 1000 {{print $1}}' {input} | sort -u > {output}"

# extract sequences that align to ref genome
rule get_seqs:
    input:
        contigs=rules.contig_list.output,
        assembly=config['sample_dir'] + "{sample}.fasta"
    output: "crass_genomes/{sample}.crass.fasta"
    resources:
        mem=1,
        time=1
    shell:
        "seqtk subseq {input.assembly} {input.contigs} > {output}"

# concatenate all sequences and rename fasta header
# rule format_fasta:
#     input: rules.get_seqs.output
#     output: "format_crass_genomes/{sample}.fasta"
#     resources:
#         mem=1,
#         time=1
#     shell:
#         "echo '>{wildcards.sample}' > {output}; "\
#         "grep -v '>' {input} >> {output}"

# align contigs to reference genome with nucmer
# rule nucmer:
#     input:
#         ref,
#         rules.get_seqs.output
#     output: "nucmer/{sample}.delta"
#     resources:
#         mem=4,
#         time=1
#     shell:
#         "nucmer {input[0]} {input[1]} -p nucmer/{wildcards.sample}"

# filter for primary alignments
# rule filter:
#     input: rules.nucmer.output
#     output: "nucmer/{sample}.filtered.delta"
#     resources:
#         mem=1,
#         time=1
#     shell:
#         "delta-filter -1 {input} > {output}"

# get coordinates from nucmer delta file
# rule show_coords:
# 	input: rules.filter.output
# 	output: "nucmer/{sample}.coords"
# 	resources:
# 		mem=4,
# 		time=1
# 	shell:
# 		"show-coords {input} > {output}"

# scaffold with medusa
# rule medusa:
# 	input: rules.get_seqs.output
# 	output: "medusa_scaffolds/{sample}.crass.fastaScaffold.fasta"
# 	run:
# 		nseqs = int(subprocess.check_output("grep -c '>' {f}".format(f=input), shell=True).decode("utf-8").rstrip())
# 		if nseqs == 1:
# 			shell("ln -s $PWD/{i} {o}".format(i=input, o=output))
# 		else:
# 			shell("java -jar ~/fiona/tools/medusa/medusa.jar -i {input} -f medusa_drafts "\
# 			"-scriptPath ~/fiona/tools/medusa/medusa_scripts -o {output}")
# 		shell("sed -E -i 's/^>.+/>{s}/g' {o}".format(s=wildcards.sample, o=output))

# detect strand with mummer
rule nucmer:
	input:
		ref,
		rules.get_seqs.output
	output: "nucmer/{sample}.delta"
	resources:
		mem=4,
		time=1
	shell:
		"nucmer {input[0]} {input[1]} -p nucmer/{wildcards.sample}"

# filter for primary alignments
rule filter:
    input: rules.nucmer.output
    output: "nucmer/{sample}.filtered.delta"
    resources:
        mem=1,
        time=1
    shell:
        "delta-filter -q {input} > {output}"

rule show_coords:
    input: rules.filter.output
    output: "nucmer/{sample}.filtered.coords"
    resources:
        mem=1,
        time=1
    shell:
        "show-coords -THrd {input} > {output}"

rule scaffold_1:
	input:
		rules.show_coords.output
	output:
		"nucmer/{sample}.coords.small"
	shell:
		"cut -f9,11 {input} | awk '!seen[$0]++' > {output}"

rule scaffold_2:
	input:
		rules.scaffold_1.output,
		rules.get_seqs.output
	output:
		"scaffold/{sample}.fasta"
	run:
		with open(input[0], "r") as contigs:
			for line in contigs:
				strand, contig = line.rstrip().split('\t')
				if strand == "1":
					shell("samtools faidx {i} {c} >> {o}".format(i=input[1], c=contig, o=output))
				else:
					shell("samtools faidx --reverse-complement " \
					"--mark-strand sign {i} {c} >> {o}".format(i=input[1], c=contig, o=output))

rule stitch:
	input: rules.scaffold_2.output
	output: "scaffold/{sample}.fastaScaffold"
	shell:
		"echo '>{wildcards.sample}' > {output}; " \
		"sed -E 's/^>.+/NNNNNNNNNN/g' {input} | sed '1d' >> {output}"

rule done:
    input: expand("scaffold/{sample}.fastaScaffold", sample=samples)
    output: "done"
    resources:
        mem=1,
        time=1
    shell:
        "touch {output}"

# find longest match between crassphage ref and genome
# rule rev_comp:
# 	input:
# 		rules.nucmer.output,
# 		rules.medusa.output
# 	output:	"flipped/{sample}.flipped.fasta"
# 	run:
# 		strand = subprocess.check_output("show-coords -B {d} | "\
# 		"tr '[:blank:]' '\t' | sort -k11,11n | head -1 | "\
# 		"cut -f20".format(d=input[0]), shell=True).decode("utf-8").rstrip()
#
# 		# print(wildcards.sample + ': ' + str(strand))
# 		# print(str(strand) == "Plus")
#
# 		# if on the plus strand, do nothing
# 		if strand == 'Plus':
# 			shell("ln -s $PWD/{i} {o}".format(i=input[1], o=output))
#
# 		# otherwise rev comp
# 		else:
# 			record = SeqIO.read(input[1], "fasta")
# 			with open(output[0], "w") as output_handle:
# 				output_handle.write('>' + wildcards.sample + '\n')
# 				output_handle.write(str(record.seq.reverse_complement()) + '\n')

# get fasta file from alignments
# rule fasta:
#     input: rules.filter.output
#     output: "final/{sample}.rotated.fasta"
#     resources:
#         mem=4,
#         time=1
#     shell:
#         "show-aligns {input} " + refname + " {wildcards.sample} | "\
#         "./nucmerAlignsToFasta.py {wildcards.sample} | sed 's/\.//g'> {output}"

rule combine:
    input: expand("scaffold/{sample}.fastaScaffold", sample=samples)
    output: "all_samples.fasta"
    resources:
        mem=1,
        time=1
    shell:
        "cat " + ref + " {input} > {output}"

# rule rotate:
# 	input: rules.combine.output
# 	output: "all_samples-Rotated.fasta"
# 	shell:
# 		"~/fiona/tools/csa_r/CSA R {input}"

rule prodigal:
    input: rules.combine.output
    output:
        "prodigal/all_samples.faa",
        "prodigal/all_samples.gff"
    resources:
        mem=8,
        time=1
    shell:
        "prodigal -a {output[0]} -f gff -o {output[1]} -i {input} -t /home/tamburin/fiona/tools/annotate_crassphage/podoviridae.trn; "\
        "sed -i 's/ /_/g' {output[0]}"

# rule glimmer:
#     input: rules.combine.output
#     output:
#         "prodigal/all_samples.faa",
#         "prodigal/all_samples.gff"
#     resources:
#         mem=8,
#         time=1
#     shell:
#         "prodigal -a {output[0]} -f gff -o {output[1]} -i {input}; "\
#         "sed -i 's/ /_/g' {output[0]}"

rule cdhit:
    input: "prodigal/all_samples.faa"
    output: "cdhit/all_samples.faa.clstr"
    resources:
        mem=8,
        time=1
    shell:
        "cd-hit -c 0.75 -s .9 -d 0 -i {input} -o cdhit/all_samples.faa"

rule parse_clusters:
    input:
        "cdhit/all_samples.faa.clstr",
        "prodigal/all_samples.gff"
    output: "all_samples.genes"
    resources:
        mem=1,
        time=1
    script:
        "parseClusters.py"
