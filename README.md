# Plot genomes workflow
Identify open reading frames and plot orthologs between small genomes

1. Scaffold assemblies against reference
2. Identify open reading frames with Prodigal
3. Identify orthologous genes and plot genome

## Config file
```
# reference genome for blast alignment
ref: /path/to/reference.fasta

# prefix for output files:
prefix: something

# path to assemblies (.fasta or .fa)
sample_dir: /path/to/assemblies

# min length of contigs to include in scaffold
min_len: 200

# prodigal training data
# podoviridae: /labs/asbhatt/fiona/tools/annotate_crassphage/podoviridae.trn
training_data: /path/to/training_data

# parameters for CDHIT clustering
pid: 0.75
length: 0.9
```

## Usage

source activate crassviz

snakemake --configfile config.yaml --snakefile /labs/asbhatt/fiona/tools/annotate_crassphage/Snakefile --profile scg -j 999 --rerun-incomplete --restart-times 0 --latency-wait 20
