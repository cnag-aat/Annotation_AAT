from datetime import datetime
import re
import os

rule get_training_candidates:
  input:
    transdecoder = "step03_annotation_pipeline.V01/protein_and_transcript_mappings/pasa/pasadb.assemblies.fasta.transdecoder.genome.gff3",
    Repeats = "step01_Repeat_Annotation.V01/Repeats.4jb.gff3",
    genome = "genome.fa"
  output:
    out =  "good_candidates.1trans.gff3"
  params:
    outdir = 'step03_annotation_pipeline.V01/train_gene_predictors/get_training_candidates',
    scripts_dir = "../scripts/",
    project = "test",
    blastdb = "/scratch/project/devel/aateam/blastdbs/swissprot",
    eval = 0.000001
  conda:
    "../envs/blast_bedtools.yaml"
  threads: 4
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir};"
    "ln -s {input.transdecoder};"
    "ln -s {input.Repeats};"
    "export PATH={params.scripts_dir}:$PATH;"
    "gawk \'$3==\"mRNA\"\' {input.transdecoder} | grep complete | perl -ane \'chomp; if ($_ =~m/ID=([^;]+)/){{print \"$1\\n\";}}\' > transdecoder.genome.complete.ids;"
    "{params.scripts_dir}filtering_transcriptlist.pl  -gff {input.transdecoder} -l transdecoder.genome.complete.ids -a > transdecoder.genome.complete.gff3;"
    "{params.scripts_dir}assignIDs_generic_fast.pl -f transdecoder.genome.complete.gff3  -project {params.project} -v 1;"
    "cut -f 1,4,5 {input.Repeats} | sort -k1,1 -k2,2n | perl -ane \'chomp; $len = $F[1]-1; print \"$F[0]\\t$len\\t$F[2]\\n\";\' > Repeats.sorted.bed;"
    "bedtools merge -i Repeats.sorted.bed > Repeats.sorted.merged.bed;"
    "mv *.{params.project}.1.gff3 {params.project}.1.gff3;"
    "bedtools intersect -wa -u -a {params.project}.1.gff3 -b Repeats.sorted.merged.bed | gawk \'$3==\"CDS\"\' | perl -ane \'chomp; if ($_ =~m/Parent=([^;]+)/){{print \"$1\\n\";}}\' | uniq > {params.project}.1.Repeats.ids;"
    "{params.scripts_dir}filtering_transcriptlist.pl -gff {params.project}.1.gff3 -l {params.project}.1.Repeats.ids -b > {params.project}.1.no_repeats.gff3;"
    "{params.scripts_dir}CDS2seq.v2.pl -gff {params.project}.1.no_repeats.gff3 -seq {input.genome};"
    "get_longest_peptide.pl {params.project}.1.no_repeats.pep.fa > {params.project}.1.no_repeats.longest_peptide.fa;"
    "blastp -outfmt \"6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore\" -query {params.project}.1.no_repeats.longest_peptide.fa -db {params.blastdb} -num_threads {threads} > blast.out;"
    "gawk \'$3>50 && $13<{params.eval}\' blast.out | cut -f 1 | sort | uniq > good_candidates.ids;"
    "filtering_protein_list.pl -gff {params.project}.1.gff3 -l good_candidates.ids -a | grep -v \'#\' > good_candidates.gff3;"
    "gawk \'$3==\"transcript\"\' good_candidates.gff3  |  perl -ane \'chomp; if ($_ =~ m/ID=([^;]+)/){{$t = $1;}} @g = split /T/, $t; print \"$t\\n\" if (!exists $ids{{$g[0]}}); $ids{{$g[0]}}++;\'" 
    "| shuf | head -n 1000 > good_candidates.trans.ids;"
    "filtering_transcriptlist.pl -gff good_candidates.gff3 -l good_candidates.trans.ids -a > {output.out};"

rule train_augustus:
  input:
    candidates = "step03_annotation_pipeline.V01/train_gene_predictors/get_training_candidates/good_candidates.1trans.gff3",
    genome = "genome.masked.fa"
  output:
    trained = "human.trained"
  params:
    outdir = 'step03_annotation_pipeline.V01/train_gene_predictors/train_augustus/',
    species = "human",
    config_path="/software/assembly/conda/augustus3.5.0/config/"
  conda:
    "../envs/augustus3.5.0.yaml"
  threads: 24
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir};"
    "ln -s {input.candidates};"
    "export AUGUSTUS_CONFIG_PATH={params.config_path};"
    "$CONDA_PREFIX/bin/gff2gbSmallDNA.pl {input.candidates} {input.genome} 10000 candidates.gb;"
    "$CONDA_PREFIX/bin/randomSplit.pl candidates.gb 600;"
    "$CONDA_PREFIX/bin/new_species.pl --species={params.species};"
    "augustus --species={params.species} candidates.gb.test | tee firsttest.out;"
    "$CONDA_PREFIX/bin/optimize_augustus.pl --species={params.species} --cpus={threads} candidates.gb.train --kfold={threads};"
    "augustus --species={params.species} candidates.gb.test | tee secondtest.out;"
    "touch {output.trained};"
