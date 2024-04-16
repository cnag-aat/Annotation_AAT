from datetime import datetime
import re
import os

rule repeat_masker:
  input: 
    genome = "genome.fa"
  output:
    masked = "genome.fa.repeatmasker",
    out= "genome.fa.out"
  params:
    db = " -species human",
    repdir = "step01_RepeatAnnotation.V01/",
    inbasename = "genome"
  conda:
    "../envs/repeatmasker-4.1.5-0.yaml"
  shell:
    "cd $TMPDIR;"
    "RepeatMasker -nolow -gff -dir {params.repdir} -pa {threads} {params.db} {input.genome};"
    "mv {params.repdir}{params.inbasename}.masked {output.masked};"
    "mv {params.repdir}{params.inbasename}.out {output.out};"

rule redmask:
  input:
    genome = "genome.fa"
  output:
    fasta = "step01_RepeatAnnotationPipeline.V01/redmask.msk.fa",
    bed = "step01_RepeatAnnotationPipeline.V01/redmask.rep.bed",
  params:
    word_length = 15,
    min_kmer = 3,
    additional_redmask_opts = "",
    run_dir = "step01_RepeatAnnotationPipeline.V01/redmask_run/"
  conda:
    "../envs/redmask.yaml"
  shell:
    "mkdir -p {params.run_dir}/redmask_in;"
    "cd {params.run_dir}/redmask_in;"
    "ln -s {input.genome} input.fa;"
    "Red -gnm {params.run_dir}/redmask_in -msk {params.run_dir} -rpt {params.run_dir} "
    " -len {params.word_length} -min {params.min_kmer} {params.additional_redmask_opts} -frm 2;"
    "ln -s {params.run_dir}/input.msk {output.fasta};"
    "cat  {params.run_dir}/input.rpt | sed 's/>//g' > {output.bed};"

rule filter_prot:
  input:
    genome = "genome.fa",
    repeat_bed = "repeats.bed",
  output:
    repeat_fasta = "repeats.fa",
    blast = "blast.out",
    masked_genome = "genome.masked.fa",
    redmask_bed = "repeats.noblast.bed"
  params:
    blastdb = "swissprot",
    evalue = 0.000001, 
    repdir = "step01_Repeat_Annotation.V01/"
  conda:
    "../envs/blast_bedtools.yaml"
  shell:
    "bedtools getfasta -fi {input.genome} -bed {input.repeat_bed} -fo {output.repeat_fasta};"
    "blastx -query {output.repeat_fasta} -db {params.blastdb} -num_threads {threads} -outfmt 6 -evalue {params.evalue} -out {output.blast};"
    "cut -f 1 {output.blast} | sort | uniq | perl -ane 'if ($_ =~ m/([^:]+):([^-]+)-([^\n]+)/){{print \"$1\t$2\t$3\n\";}}' > {params.repdir}repeats.blast.bed ;"
    "bedtools intersect -f 1 -r -v -a {input.repeat_bed} -b {output.redmask_bed} | bedtools merge -i - > {params.repdir}repeats.noblast.bed;"
    "bedtools maskfasta -fi {input.genome} -bed {params.repdir}repeats.noblast.bed -fo {output.masked_genome};"

rule get_repeats_gff:
  input:
    redmask_bed = "repeats.noblast.bed",
    repeatmasker_out = "genome.fa.out"
  output:
    intermediate_gffs = ["repeats.noblast.gff3", "repeatmasker.gff3"],
    repeats_gff = "Repeats.4jb.gff3"
  params:
    repdir = "step01_Repeat_Annotation.V01/"
  threads: 1
  shell:
    "if [ -f {input.redmask_bed} ]; then cat {input.redmask_bed} | perl -ane \'chomp; $start = $F[1]+1; " +\
    "print \"$F[0]\\tredmask\\tInterspersed_Repeat\\t$start\\t$F[2]\\t.\\t.\\t.\\tID=$F[0]:$start-$F[2]\\n\";\'"
    " > {params.repdir}/repeats.noblast.gff3; fi;"
    "if [ -f {input.repeatmasker_out} ]; then cat {input.repeatmasker_out} | perl -ne \'chomp; next if /^$/o; $_ =~ s/^\\s+|\\s+$//g; next if /SW/; " +\
    "next if /score/;@line = split /\\s+/, $_; @class = split /\//, $line[10]; if ($line[8] eq \"+\")" +\
    " {{print \"$line[4]\\tRepeatMasker\\tInterspersed_Repeat\\t$line[5]\\t$line[6]\\t$line[0]\\t$line[8]\\t.\\t" +\
    "Name=$line[9];Class=$class[0];Family=$class[1];Divergence=$line[1]%;Deletions=$line[2]%;" +\
    "Insertions=$line[3]%;Target=$line[9]:$line[12]-$line[11]\\n\";}} " +\
    "else {{print \"$line[4]\\tRepeatMasker\\tInterspersed_Repeat\\t$line[5]\\t$line[6]\\t$line[0]\\t-\\t.\\t"+\
    "Name=$line[9];Class=$class[0];Family=$class[1];Divergence=$line[1]%;Deletions=$line[2]%;Insertions=$line[3]%;Target=$line[9]:$line[13]-$line[12]\\n\";}}\'" +\
    "> {params.repdir}/repeatmasker.gff3; fi;"
    "cat {output.intermediate_gffs} | sort -k1,1 -k4,4n > {output.repeats_gff};"

