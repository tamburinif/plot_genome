# Fiona Tamburini
# Sept 2018

# Finds contigs aligning to a reference genome, uses nucmer to align/scaffold
# contigs against the reference and rotate genome to same starting coordinates

localrules: makeblastdb, contig_list, get_seqs, nucmer, filter, show_coords, scaffold_1, scaffold_2, stitch, done, combine

import subprocess

##### prefix for output files #####
if config['prefix'] is None:
	prefix = "all_samples"
else:
	prefix = config['prefix']

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

# reference configuration
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
		mem=2,
		time=2,
		threads=8
	shell:
		"blastn -db {input.ref} -query {input.assembly} -num_threads {resources.threads} \
		-outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart \
		send evalue bitscore\" -out {output}"

# aggregate rule to execute all blast searches
rule blast_all:
    input: expand("blast/{sample}.blast.out", sample=samples)
    output: "blast_done"
    resources:
        mem=1,
        time=1
    shell:
        "touch {output}"

# rule all:
#     input: expand("crass_genomes/{sample}.crass.fasta", sample=samples)
#     output: "done"
#     resources:
#         mem=1,
#         time=1
#     shell:
#         "touch {output}"

# get a list of contigs over 200 bb that align to ref genome
# changed - aligned length must be over 1kb
rule contig_list:
	input: rules.blast.output
	output: "contig_list/{sample}.crass.list"
	resources:
		mem=1,
		time=1
	params:
		min_length_bp=200,
		min_p=0.1
	shell:
		"awk '$5 > {params.min_length_bp} && $5/$4 > {params.min_p} {{print $1}}' \
		{input} | sort -u > {output}"

# extract contigs that align to reference genome
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

# filter nucmer output for primary alignments
rule filter:
    input: rules.nucmer.output
    output: "nucmer/{sample}.filtered.delta"
    resources:
        mem=1,
        time=1
    shell:
        "delta-filter -q {input} > {output}"

# get nucmer coordinates to scaffold contigs
rule show_coords:
    input: rules.filter.output
    output: "nucmer/{sample}.filtered.coords"
    resources:
        mem=1,
        time=1
    shell:
        "show-coords -THrd {input} > {output}"

# scaffold contigs
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

# concatenate contigs
rule stitch:
	input: rules.scaffold_2.output
	output: "scaffold/{sample}.fastaScaffold"
	shell:
		"(grep -v '>' {input} | awk 'BEGIN {{ ORS=\"\"; print \">{wildcards.sample}\\n\" }} {{ print }}'; printf '\\n') > {output}"

# rule done:
#     input: expand("scaffold/{sample}.fastaScaffold", sample=samples)
#     output: "done"
#     resources:
#         mem=1,
#         time=1
#     shell:
#         "touch {output}"

# combine all genomes into multifasta
rule combine:
    input: expand("scaffold/{sample}.fastaScaffold", sample=samples)
    output: "{name}.fasta".format(name = prefix)
    resources:
        mem=1,
        time=1
    shell:
        "cat " + ref + " {input} > {output}"

# find genes in multifasta genomes
rule prodigal:
	input:
		fasta=rules.combine.output,
		training=config['training_data']
	output:
		"prodigal/{name}.faa".format(name = prefix),
		"prodigal/{name}.gff".format(name = prefix)
	resources:
		mem=8,
		time=1
	shell:
		"prodigal -a {output[0]} -f gff -o {output[1]} -i {input.fasta} -t {input.training}; "\
		"sed -i 's/ /_/g' {output[0]}"

# cluster protein sequences
rule cdhit:
	input: "prodigal/{name}.faa".format(name = prefix)
	output:
		"cdhit/{name}.faa.clstr".format(name = prefix),
		"cdhit/{name}.faa".format(name = prefix)
	resources:
		mem=8,
		time=1
	params:
		pid=config['pid']
		length=config['length']
	shell:
		"cd-hit -c {params.pid} -s {params.length} -d 0 -i {input} -o {output[1]}"

# parse clusters into gff-like gene format
rule parse_clusters:
    input:
        "cdhit/{name}.faa.clstr".format(name = prefix),
        "prodigal/{name}.gff".format(name = prefix),
    output: "{name}.genes".format(name = prefix)
    resources:
        mem=1,
        time=1
    script:
        "parseClusters.py"

# plot genomes, color by homolog
# rule plot:
	# input: rules.parse_clusters.output
	# output: "{name}.pdf".format(name = prefix)
	# script:
	# 	"plotGenes.R"
