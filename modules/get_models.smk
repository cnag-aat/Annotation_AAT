rule stringtie:
  input:
    bam = "star/readsAligned.sortedByCoord.out.bam"
  output: 
    models = "stringtie/reads.stringtie.gtf"
  params:
    stringtie_opts = "",
    basename = "reads",
  threads: 4
  conda:
    "../envs/stringtie2.2.1.yaml"
  shell:
    "stringtie {input.bam} {params.stringtie_opts} -l {params.basename} -p {threads} -o {output.models};"

rule join_models:
  input:
    genome = "assemvbly.fa",
    models = "stringtie/reads.stringtie.gtf" 
  output:
    joined_models = "TACO_assembled.gtf",
  params:
    TACO_opts = "",
    outdir = "TACO_output/"
  threads: 4
  conda:
    "../envs/taco0.7.3.yaml"
  shell:
    "cd {params.outdir};"
    "echo {input.models}| sed 's/\s/\\n/g' > gtf_files.txt;"
    "taco_run -o {params.outdir}TACO_dir -p {threads} --filter-splice-juncs {params.TACO_opts} --ref-genome-fasta {input.genome} gtf_files.txt;"
    "cat {params.outdir}TACO_dir/assembly.gtf | sed 's/expr/FPKM/g' > {output.joined_models};"

rule portcullis:
  input:
    genome = "assembly.fa",
    bams =  "star/readsAligned.sortedByCoord.out.bam" ,
  output:
    junctions = "junctions.gff"
  params:
  threads: 4
  conda:
    "../envs/portcullis1.2.4.yaml"
  shell:
    "portcullis full --intron_gff -t {threads} {input.genome} {input.bams};"
    "gawk \'$7!=\"?\"\' portcullis_out/3-filt/portcullis_filtered.pass.junctions.intron.gff3 > {output.junctions};"

rule ESPRESSO:
  input:
    genome = "assembly.fa",
    sams =  "star/readsAligned.sortedByCoord.out.sam",
    tsv = "samples.tsv"
  output:
    junctions = "junctions.gff",
  params:
    ESPRESSO_path = "/software/assembly/src/ESPRESSO/espresso_v_1_3_0_beta/src/",
    out_dir = "ESPRESSO_out"
  threads: 4
  conda:
    "../envs/ESPRESSO1.3.0.yaml"
  shell:
    "perl {params.ESPRESSO_path}ESPRESSO_S.pl -L {input.tsv} -F {input.genome} -O {params.out_dir} -T {threads};"
    "grep -v SJ_cluster {params.out_dir}/other_SJ_simplified.list | gawk \'$5==0 || $5==1\' | " +\
    " perl -ane \'chomp; @scaff = split /\:/, $F[1]; $start = $F[2] + 1; $freq=$F[5]+1; " +\
    " if ($F[4] == 0){{$strand = \"+\"}}else {{$strand=\"-\"}};print " +\
    "\"$scaff[0]\tESPRESSO\tintron\t$start\t$F[3]\t$freq\t$strand\t.\tgrp=$F[0];src=E\n\";\' > {output.junctions};"
    "rm {input.sams};"
