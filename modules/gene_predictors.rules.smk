from datetime import datetime
import re
import os

rule augustus_jobarray:
  input:
    fasta = "assembly.masked.1.fa"
  output:
    touch = "augustus_run.ok" 
  params:
    prefix = "assembly.masked",
    out_prediction_augustus = "assembly.masked.augustus_gene_predictions.gff3",
    species = "human",
    alternatives_from_sampling =  "true",
    alternatives_from_evidence =   "true",
    uniqueGeneId =   "true",
    gff3 =  "on",
    sample =  60,
    noInFrameStop =   "true",
    maxtracks =  2,
    singlestrand =   "false",
    strand = "both",
    min_intron_len = 30,
    additional_aug_opts = ""
  conda: 
    "../envs/augustus3.5.0.yaml"
  threads: 1
  shell:
      "export AUGUSTUS_CONFIG_PATH=/software/assembly/conda/augustus3.5.0/config/;"
      "augustus --species={params.species} --alternatives-from-sampling={params.alternatives_from_sampling} " +\
      "--alternatives-from-evidence={params.alternatives_from_evidence} --sample={params.sample} --gff3={params.gff3}" +\
      " --noInFrameStop={params.noInFrameStop} --uniqueGeneId={params.uniqueGeneId} --maxtracks={params.maxtracks} --strand={params.strand}" +\
      " --singlestrand={params.singlestrand} --min_intron_len={params.min_intron_len} {params.additional_aug_opts} " +\
      "{params.prefix}.$SLURM_ARRAY_TASK_ID.fa  > {params.out_prediction_augustus}.$SLURM_ARRAY_TASK_ID; "
      "[ ! -f {output.touch} ]; touch {output.touch};"

rule augustus:
  input:
    fasta = "assembly.masked.fa"
  output:
    predictions = "assembly.masked.augustus_gene_predictions.gff3" 
  params:
    species = "human",
    alternatives_from_sampling =  "true",
    alternatives_from_evidence =   "true",
    uniqueGeneId =   "true",
    gff3 =  "on",
    sample =  60,
    noInFrameStop =   "true",
    maxtracks =  2,
    singlestrand =   "false",
    strand = "both",
    min_intron_len = 30,
    additional_aug_opts = ""
  conda: 
    "../envs/augustus3.5.0.yaml"
  threads: 1
  shell:
      "export AUGUSTUS_CONFIG_PATH=/software/assembly/conda/augustus3.5.0/config/;"
      "augustus --species={params.species} --alternatives-from-sampling={params.alternatives_from_sampling} " +\
      "--alternatives-from-evidence={params.alternatives_from_evidence} --sample={params.sample} --gff3={params.gff3}" +\
      " --noInFrameStop={params.noInFrameStop} --uniqueGeneId={params.uniqueGeneId} --maxtracks={params.maxtracks} --strand={params.strand}" +\
      " --singlestrand={params.singlestrand} --min_intron_len={params.min_intron_len} {params.additional_aug_opts} " +\
      "{input.fasta}  > {output.predictions}; "

rule augustus_hints_jobarray:
  input:
    fasta = "assembly.masked.1.fa",
    hints = "junctions.gff3",
    extrinsic_file = "extrinsic.E.cfg"
  output:
    touch ="augustus_hints_run.ok"
  params:
    scripts_dir = "scripts/",
    prefix = "assembly.masked",
    out_prediction_augustus = "assembly.masked.augustus_hints_gene_predictions.gff3",
    species = "human",
    alternatives_from_sampling =  "true",
    alternatives_from_evidence =   "true",
    uniqueGeneId =   "true",
    gff3 =  "on",
    sample =  60,
    noInFrameStop =   "true",
    maxtracks =  2,
    singlestrand =   "false",
    strand = "both",
    min_intron_len = 30,
    additional_aug_opts = "",
    rmcmd = "cd ../..; rm -r $TMPDIR/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID;"
  conda: 
    "../envs/augustus3.5.0.yaml"  
  threads: 1
  shell:
    "export AUGUSTUS_CONFIG_PATH=/software/assembly/conda/augustus3.5.0/config/;"
    "dir=$TMPDIR/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID;"
    "echo $dir;"
    "mkdir -p $dir;"
    "cd $dir;"
    "cp {input.hints} hints.gff;"
    "sort -k1,1 -k4,5n -k7,7 hints.gff > hints.sorted.gff;"
    "ln -s {params.prefix}.$SLURM_ARRAY_TASK_ID.fa masked_genome_chunk.fa;"
    "cat hints.sorted.gff | cut -f 1 | uniq > seqs_with_hints.ids;"
    "{params.scripts_dir}/filter_fasta_from_ids.pl -f masked_genome_chunk.fa -l seqs_with_hints.ids -a > masked_genome_chunk.with_hints.fa;"
    "cat hints.sorted.gff | {params.scripts_dir}/run_augustus_with_hints.v3.pl -f masked_genome_chunk.with_hints.fa " +\
    " -s {params.species} -alt {params.alternatives_from_sampling} -sample {params.sample} -gff {params.gff3} -iF {params.noInFrameStop} " +\
    " -uniq {params.uniqueGeneId} -max {params.maxtracks} -str {params.strand} -ss {params.singlestrand} -int {params.min_intron_len} " +\ 
    " -ef {input.extrinsic_file} -ad \"{params.additional_aug_opts}\";" 
    "cat *augustus_introns.gff3 > {params.out_prediction_augustus}.$SLURM_ARRAY_TASK_ID;"
    "{params.rmcmd};"
    "[ ! -f {output.touch} ]; touch {output.touch};"

rule augustus_hints:
  input:
    fasta = "assembly.masked.fa",
    hints = "junctions.gff3",
    extrinsic_file = "extrinsic.E.cfg"
  output:
    predictions = "assembly.masked.augustus_hints_gene_predictions.gff3" 
  params:
    scripts_dir = "scripts/",
    species = "human",
    alternatives_from_sampling =  "true",
    alternatives_from_evidence =   "true",
    uniqueGeneId =   "true",
    gff3 =  "on",
    sample =  60,
    noInFrameStop =   "true",
    maxtracks =  2,
    singlestrand =   "false",
    strand = "both",
    min_intron_len = 30,
    additional_aug_opts = "",
    rmcmd = "cd ../..; rm -r $TMPDIR/augustus_hints;"
  conda: 
    "../envs/augustus3.5.0.yaml"  
  threads: 1
  shell:
    "export AUGUSTUS_CONFIG_PATH=/software/assembly/conda/augustus3.5.0/config/;"  
    "mkdir -p $TMPDIR/augustus_hints;"
    "cd $TMPDIR/augustus_hints;"
    "cp {input.hints} hints.gff;"
    "sort -k1,1 -k4,5n -k7,7 hints.gff > hints.sorted.gff;"
    "ln -s {input.fasta} masked_genome_chunk.fa;"
    "cat hints.sorted.gff | cut -f 1 | uniq > seqs_with_hints.ids;"
    "{params.scripts_dir}/filter_fasta_from_ids.pl -f masked_genome_chunk.fa -l seqs_with_hints.ids -a > masked_genome_chunk.with_hints.fa;"
    "cat hints.sorted.gff | {params.scripts_dir}/run_augustus_with_hints.v3.pl -f masked_genome_chunk.with_hints.fa " +\
    " -s {params.species} -alt {params.alternatives_from_sampling} -sample {params.sample} -gff {params.gff3} -iF {params.noInFrameStop} " +\
    " -uniq {params.uniqueGeneId} -max {params.maxtracks} -str {params.strand} -ss {params.singlestrand} -int {params.min_intron_len} " +\ 
    " -ef {input.extrinsic_file} -ad \"{params.additional_aug_opts}\";"  
    "cat *augustus_introns.gff3  > {output.predictions};"
    "{params.rmcmd};"

rule geneid_jobarray:
  input:      
    fasta = "assembly.masked.1.fa",
    geneid_parameters = "human.params.txt"
  output:
    touch = "geneid_run.ok"  
  params:
    geneid_options = " -3U ",
    path = "/software/assembly/src/geneid/",
    prefix = "assembly.masked",
    out_prediction_geneid = "assembly.masked.geneid_gene_predictions.gff3",
  threads: 2
  shell:
    "{params.path}bin/geneid -P {input.geneid_parameters} {params.geneid_options} "+\
    "{params.prefix}.$SLURM_ARRAY_TASK_ID.fa  > {params.out_prediction_geneid}.$SLURM_ARRAY_TASK_ID; "
    "[ ! -f {output.touch} ]; touch {output.touch};"

rule geneid:
  input:      
    fasta = "assembly.masked.fa",
    geneid_parameters = "human.params.txt"
  output: 
    out_prediction_geneid = "assembly.masked.geneid_gene_predictions.gff3",
  params:
    geneid_options = " -3U ",
    path = "/software/assembly/src/geneid/"
  threads: 2
  shell:
    "{params.path}bin/geneid -P {input.geneid_parameters} {params.geneid_options} "+\
    "{input.fasta}  > {output.out_prediction_geneid} "

rule geneid_introns_jobarray:
  input:  
    fasta = "assembly.masked.1.fa",
    geneid_parameters = "human.params.txt",
    junctions = "junctions.gff3"
  output:
    touch = "geneid_introns_run.ok" 
  params:
    scripts_dir = "scripts/",
    geneid_options = " -3nU ",
    path = "/software/assembly/src/geneid/",
    prefix = "assembly.masked",
    out_prediction_geneid_introns = "assembly.masked.geneid_introns_gene_predictions.gff3",
    rmcmd = "cd ../..; rm -r $TMPDIR/geneid_introns;"
  threads: 2
  shell:
    "[ ! -f {output.touch} ]; touch {output.touch};"
    "dir=$TMPDIR/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID;"
    "echo $dir;"
    "mkdir -p $dir;"
    "cd $dir;"
    "cat {input.junctions} | sed 's/intron/Intron/g' | sort -k1,1 -k4,5n -k7,7 > hints.sorted.gff ;"  
    "ln -s {params.prefix}.$SLURM_ARRAY_TASK_ID.fa masked_genome_chunk.fa;"
    "cat hints.sorted.gff | cut -f 1 | uniq > seqs_with_hints.ids;"
    "{params.scripts_dir}/filter_fasta_from_ids.pl -f masked_genome_chunk.fa -l seqs_with_hints.ids -a > masked_genome_chunk.with_hints.fa;"
    "{params.scripts_dir}/rungeneidwithhints.pl masked_genome_chunk.with_hints.fa hints.sorted.gff " +\
    " {input.geneid_parameters} {params.geneid_options} {params.path};"
    "cat *geneid_introns.gff3 > {params.out_prediction_geneid_introns}.$SLURM_ARRAY_TASK_ID;"
    "{params.rmcmd};"    

rule geneid_introns:
  input:  
    fasta = "assembly.masked.fa",
    geneid_parameters = "human.params.txt",
    junctions = "junctions.gff3"
  output:
    predictions = "assembly.masked.geneid_introns_gene_predictions.gff3",
  params:
    scripts_dir = "scripts/",
    geneid_options = " -3nU ",
    path = "/software/assembly/src/geneid/",
    prefix = "assembly.masked",
    rmcmd = "cd ../..; rm -r $TMPDIR/geneid_introns;"
  threads: 2
  shell:
    "mkdir -p $TMPDIR/geneid_introns;"
    "cd $TMPDIR/geneid_introns;"
    "cp {input.junctions} hints.gff;"
    "sort -k1,1 -k4,5n -k7,7 hints.gff > hints.sorted.gff;"
    "ln -s {input.fasta} masked_genome_chunk.fa;"
    "cat hints.sorted.gff | cut -f 1 | uniq > seqs_with_hints.ids;"
    "{params.scripts_dir}/filter_fasta_from_ids.pl -f masked_genome_chunk.fa -l seqs_with_hints.ids -a > masked_genome_chunk.with_hints.fa;"
    "{params.scripts_dir}/rungeneidwithhints.pl masked_genome_chunk.with_hints.fa hints.sorted.gff " +\
    " {input.geneid_parameters} {params.geneid_options} {params.path};"
    "cat *geneid_introns.gff3  > {output.predictions};"
    "{params.rmcmd};"
    # "{params.path}bin/geneid -R {input.junctions} -P {input.geneid_parameters} {params.geneid_options} "+\
    # "{input.fasta}  > {output.predictions}; "

rule genemark:
  input:
    fasta = "assembly.masked.fa"
  output:
    out = "step03_annotation_pipeline.V01/gene_predictions/genemark.gtf",
    EVM_out = "step03_annotation_pipeline.V01/gene_predictions/genemark_preEVM.gff3"
  params:
    scripts_dir = "scripts/",
    maxgap = 5000, 
    mincontig = 50000, 
    maxcontig = 500000,
    add_opts = "",
    EVM_dir = "step04_EVM.V01",
    link_out = "step04_EVM.V01/genemark_predictions.gff3",
    create_weights_gmk = "echo \"PREDICTION\tGeneMark.hmm3\t3\">> weights_1.txt;",
    rmcmd = "cd ../..; rm -r $TMPDIR;"
  envmodules:
    "GeneMark-ET"
  threads: 8,
  shell:  
    "TMPDIR=$TMPDIR/genemark; mkdir -p $TMPDIR; cd $TMPDIR;"
    "cp {input.fasta} genome_masked.fa;"
    "gmes_petap.pl --sequence genome_masked.fa --ES --max_gap {params.maxgap} "+\
    " --min_contig {params.mincontig} --cores {threads} --max_contig {params.maxcontig} {params.add_opts};"
    "cp genemark.gtf {output.out};"
    "cat {output.out} | {params.scripts_dir}/gtf_to_gff2or3.pl -noaddstop -v 2.5 -ov 3 -mrna "+\
    " | grep -v codon | grep -v \'#\' > {output.EVM_out};"
    "cd {params.EVM_dir};"
    "ln -s {output.EVM_out} genemark_predictions.gff3;"
    "{params.create_weights_gmk}"
    "{params.rmcmd}"

rule genemark_ET:
  input:    
    fasta = "assembly.masked.fa",
    hints = "junctions.gff3",
  output:
    out = "step03_annotation_pipeline.V01/gene_predictions/genemark-ET.gtf",
    EVM_out = "step03_annotation_pipeline.V01/gene_predictions/genemark-ET_preEVM.gff3"
  params:
    scripts_dir = "scripts/",
    maxgap = 5000, 
    mincontig = 50000, 
    maxcontig = 500000,
    add_opts = "",
    EVM_dir = "step04_EVM.V01",
    link_out = "step04_EVM.V01/genemark-ET_predictions.gff3",
    create_weights_gmk = "echo \"PREDICTION\tGeneMark-ET\t3\">> weights_1.txt;",
    rmcmd = "cd ../..; rm -r $TMPDIR;",
  envmodules:
    "GeneMark-ET"
  threads: 8,
  shell:  
    "TMPDIR=$TMPDIR/genemark-ET; mkdir -p $TMPDIR; cd $TMPDIR;"
    "cp {input.fasta} genome_masked.fa;"
    "cp {input.hints} junctions.gff;"
    "gmes_petap.pl --sequence genome_masked.fa --ET junctions.gff --max_gap {params.maxgap} "+\
    " --min_contig {params.mincontig} --cores {threads} --max_contig {params.maxcontig} " +\
    " {params.add_opts};"
    "cp genemark.gtf {output.out};"
    "cat {output.out} | {params.scripts_dir}/gtf_to_gff2or3.pl -noaddstop -v 2.5 -ov 3 -mrna | "+\
    " grep -v codon | grep -v '#' | sed 's/GeneMark.hmm3/GeneMark-ET/g' > {output.EVM_out};"
    "cd {params.EVM_dir};"
    "ln -s {output.EVM_out} genemark-ET_predictions.gff3;"
    "{params.create_weights_gmk}"
    "{params.rmcmd}"
