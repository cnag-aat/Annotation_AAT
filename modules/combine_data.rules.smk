from datetime import datetime
import re
import os

rule get_coding_junctions:
  input:
    junctions = "step02_RNApipeline.V01/allhints.gff3", 
    coding = "proteins.gff3"
  output:
    coding_junctions = "allhints.coding.gff3"
  params:
  threads: 1
  conda: 
    "../envs/bedtools2.30.0.yaml"
  shell:
    "dir=$TMPDIR/get_coding_junctions; mkdir -p $dir;"
    "cd $dir;"
    "gawk '$3==\"mRNA\"' {input.coding} | sort -k1,1 -k4,4n > coding.trans.sorted.gff3;"
    "sort -k1,1 -k4,4n {input.junctions} > allhints.sorted.gff3;"
    "bedtools intersect -u -a allhints.sorted.gff3 -b coding.trans.sorted.gff3 > {output.coding_junctions};"
    "cd ..; rm -r $dir;"

rule prepare_evm:
  input:
    genome = 'assembly.fa',
    predictions = ["geneid_predictions.gff3", "augustus_predictions", "genemark_predictions"],
    other_evm_inputs = ["transcripts.gff3", "proteins.gff3"],
  output:
    checkpoint = "step04_EVM.V01/evm_ready",
    predictions_out = "step04_EVM.V01/predictions.gff3",
    genome_out = "step04_EVM_V01/genome.fa",
    weights = expand( "step04_EVM.V01/weights_{w}.clean.txt", w=range(1, 3)),
  params:
    EVM_dir = "step04_EVM.V01",
    sortcommand = "sort weights_1.txt | uniq > weights_1.clean.txt;"
  threads: 1
  shell:
    "cd {params.EVM_dir};"
    "ln -s {input.genome} {output.genome_out};"
    "cat {input.predictions} > {output.predictions_out};"
    "{params.sortcommand}"
    "touch {output.checkpoint};"

rule EVM:
  input:
    weight_file = "step04_EVM.V01/weights_1.clean.txt",
    lgenome = 'assembly.fa',
    predictions = "step04_EVM.V01/predictions.gff3",
    lproteins = "step04_EVM.V01/proteins.gff3",
    ltranscripts = "step04_EVM.V01/transcripts.gff3",
  output:
    evm_models = "step04_EVM.V01/evm_weights_1.gff3",
    out = "step04_EVM.V01/evm_weights_1.out",
  params:
    set = 1,
    evidence_opts = "",
    additional_evm_opts = "",
    scripts_dir = "../scripts/"
  threads: 16
  conda:
    "../envs/evm.yaml"
  shell:
    "dir=$TMPDIR/evm.{params.set}; mkdir -p $dir;"
    "cd $dir;"
    "echo \"Running EVM step. Partitioning\";"
    "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome {input.lgenome} --gene_predictions {input.predictions} "+\
    " {params.evidence_opts} --segmentSize 2000000 --overlapSize 1000000 --partition_listing partitions_list.out;"
    "echo writing commands;"
    "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome {input.lgenome} --weights {input.weight_file} " +\
    " --gene_predictions {input.predictions} {params.evidence_opts} --output_file_name evm_{params.set}.out --partitions partitions_list.out {params.additional_evm_opts} > evm_{params.set}.cmd;"
    "echo running evm;"
    "{params.scripts_dir}/paralellize_evm.sh {wildcards.w} {params.scripts_dir}split_commands_file_evm.pl {threads} $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/execute_EVM_commands.pl;"
    "echo collecting outputs;"
    "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm_{wildcards.w}.out;"
    "echo converting to gff3;"
    "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm_{wildcards.w}.out --genome {input.lgenome};" 
    "cat */evm_{wildcards.w}.out > {output.out};"
    "cat */evm_{wildcards.w}.out.gff3 > {output.evm_models};"
   # "find . -name {wildcards.w}.out | xargs cat > {output.out};"
   # "find . -name {wildcards.w}.out.gff3 | xargs cat > {output.evm_models};"
    "cd ..; rm -r $dir;"
    "echo 'EVM done!';"

rule EVM2:
  input:
    weight_file = "step04_EVM.V01/weights_1.clean.txt",
    lgenome = "step04_EVM.V01/genome.fa",
    predictions = "step04_EVM.V01/predictions.gff3",
    lproteins = "step04_EVM.V01/proteins.gff3",
    ltranscripts = "step04_EVM.V01/transcripts.gff3",
    repeats = "step01_RepeatAnnotationPipeline/Repeats.4jb.gff3"
  output:
    evm_cds = "step04_EVM2.V01/{params.sample_id}.EVM.cds",
    evm_pep= "step04_EVM2.V01/{params.sample_id}.EVM.pep",
    evm_bed= "step04_EVM2.V01/{params.sample_id}.EVM.bed",
    evm_models = "step04_EVM2.V01/{params.sample_id}.EVM.gff3"
  params:
    sample_id = "",
    software_path = "/software/assembly/conda/EVM2.1.0/EVidenceModeler-v2.1.0/",
    additional_evm2_opts = "",
    EVM_dir = "step04_EVM2.V01",
    evidence_opts = "--protein_alignments proteins.gff3 --transcript_alignments transcripts.gff3",
    segmentSize = "100000",
    overlapSize = "10000",
    rmcmd = "rm -r $TMPDIR/evm."
  threads: 16
  conda:
    "../envs/evm2.1.yaml"
  shell:
    "dir=$TMPDIR/evm.{params.sample_id}; mkdir -p $dir;"
    "cd $dir;"
    "export PATH=$PATH:{params.software_path};"
    "EVidenceModeler --sample_id {params.sample_id} --genome {input.lgenome} --weights {input.weight_file} --gene_predictions {input.predictions} " +\
    " {params.evidence_opts} --repeats {input.repeats}" +\
    " {params.additional_evm2_opts} --segmentSize {params.segmentSize} --overlapSize {params.overlapSize} --CPU {threads};"
    "cp {params.sample_id}.EVM.cds {output.evm_cds};"
    "cp {params.sample_id}.EVM.pep {output.evm_pep};"
    "cp {params.sample_id}.EVM.bed {output.evm_bed};"
    "cp {params.sample_id}.EVM.gff3 {output.evm_models};"
    "{params.rmcmd};"

rule select_EVM:
  input:
    models =  "step04_EVM.V01/evm_weights_1.gff3",
    base_trans = "step02_RNApipeline.V01/transcripts.gtf"
  output:
    EVM_out = "step04_EVM_V01/evm.best.gff3"
  params:
    EVM_dir = "step04_EVM.V01",
    scripts_dir = "../scripts/",
    total_weight_files = 1,
    reference_field = "exon"
  threads: 1
  conda:
    "../envs/bedtools2.30.0.yaml"
  shell: 
    "cd {params.EVM_dir};"
    "gawk \'$3==\"{params.reference_field}\"\' {input.base_trans} > reference_exons.gtf;"
    "{params.scripts_dir}select_best_evm.V04.pl reference_exons.gtf {input.models};"