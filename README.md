# Plot genomes workflow
Identify open reading frames and plot orthologs between small genomes

1. Scaffold assemblies against reference
2. Identify open reading frames with Prodigal
3. Identify orthologous genes and plot genome

## Config file

  ref: reference.fasta
  
  refname: reference
  
  sample_dir: /path/to/assemblies
  
  samples: samples.list
  

## Usage

source activate annotate

snakemake foo.genes -s /path/to/Snakefile --configfile config.yaml
