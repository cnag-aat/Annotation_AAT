from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:
module jbrowse_workflow:
  snakefile: "../modules/get_jbrowse.rules.smk"

##Prepare runs
jbrowse_targets["tracks"] = []
def create_browse_dictionary(program, gff, cmd="", classname = "pred", type = "mRNA"):
    tracks.append(jbrowse_dir + "tracks/" + program)
    browse[program] = {}
    browse[program]["input"] = gff
    browse[program]["type"] = type
    browse[program]["class"] = classname
    browse[program]["new_input"] = gff
    browse[program]["setup_cmd"] = cmd
    jbrowse_targets["tracks"].append(jbrowse_dir + "browse." + program + ".ok")
    if cmd:
      browse[program]["new_input"] = jbrowse_dir + program + ".gff3"

tracks = []
browse = {}
if config["Parameters"]["run_augustus"]:
  create_browse_dictionary("Augustus", config["Outputs"]["augustus_preEVM"])

if config["Parameters"]["run_augustus_hints"]:
  create_browse_dictionary("Augustus_Hints", config["Outputs"]["augustus_hints_preEVM"])

if config["Parameters"]["run_genemark"]:
  create_browse_dictionary("Genemark", config["Outputs"]["genemark_preEVM"])

if config["Parameters"]["run_genemark-ET"]:
  create_browse_dictionary("Genemark-ET", config["Outputs"]["genemark_ET_preEVM"])

if config["Parameters"]["run_transdecoder"]:
  create_browse_dictionary("Transdecoder", EVM_dir + "transdecoder_predictions.gff3")

if config["Parameters"]["run_geneid"]:
  cmd = "cat " + config["Outputs"]["geneid_preEVM"] +\
        "| perl -ane 'chomp; if ($F[2] =~ m/CDS|exon/){{my @info = split /\;/, $_; print \"$info[0];$info[1]\n\";}}else {{print \"$_\n\";}}' > " +\
        jbrowse_dir + "Geneid.gff3;"
  create_browse_dictionary("Geneid", config["Outputs"]["geneid_preEVM"], cmd)

if config["Parameters"]["run_geneid_introns"]:
  cmd = "cat " + config["Outputs"]["geneid_introns_preEVM"] +\
        "| perl -ane 'chomp; if ($F[2] =~ m/CDS|exon/){{my @info = split /\;/, $_; print \"$info[0];$info[1]\n\";}}else {{print \"$_\n\";}}' > " +\
        jbrowse_dir + "Geneid_introns.gff3;"
  create_browse_dictionary("Geneid_introns", config["Outputs"]["geneid_introns_preEVM"], cmd)

if config["Parameters"]["run_miniprot"]:
  cmd = "grep -v '#' " + config["Outputs"]["miniprot_gene"] +\
        "| sed 's/mRNA/Protein_Alignment/g' > " + jbrowse_dir + "Miniprot.gff3;"
  create_browse_dictionary("Miniprot", config["Outputs"]["miniprot_gene"], cmd, "protali", "Protein_Alignment")

if config["Parameters"]["run_pasa"]:
  cmd = "cat " + config["Outputs"]["EVM_dir"] + "transcripts.gff3 |"+\
  scripts_dir + "pasaassemblies_2_jb.pl > PASA_assemblies.gff3;"
  create_browse_dictionary("PASA_assemblies", config["Outputs"]["EVM_dir"] + "transcripts.gff3", cmd, "transali" )

if config["Outputs"]["GTF models"]:
  cmd = "gawk '$3==\"exon\"' " + config["Outputs"]["GTF models"] + "|" + scripts_dir + "gtf_onlyexons_2transcripts.pl | " +\ 
        "perl -ane 'chomp; @F = split /\t/, $_; if ($F[2] eq \"transcript\") {{if ($_ =~m/transcript_id \"([^\";]+)/){{$F[8]= \"ID=$1\";}}}} " +\
        " elsif ($F[2] eq \"exon\"){{if ($_ =~m/transcript_id \"([^\";]+)/){{$t = $1;}}if ($_ =~ m/exon_number \"([^\"]+)/){{$F[8] = \"Parent=$t;ID=$t.$1\";}}}}$p = join \"\t\", @F;print \"$p\n\";' " +\ 
        " > Stringtie-TACO.gff3;"
  create_browse_dictionary("Stringtie-TACO", config["Outputs"]["GTF models"], cmd, "transali", "transcript" )

if config["Parameters"]["run_geneid_introns"] or config["Parameters"]["run_genemark-ET"] or config["Parameters"]["run_augustus_hints"]:
  cmd = "cat " + config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3| sed 's/intron/Intron/g' > junctions.gff3;"
  create_browse_dictionary("junctions", config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3", cmd, "junctions", "Intron" )

if config["Parameters"]["run_EVM"]:
  cmd = ""
  create_browse_dictionary("EVM-Models", config["Outputs"]["evm_out"], cmd, "protali")

if config["Parameters"]["run_update"]:
  gff = config["Outputs"]["update_dir"] + config["Parameters"]["project_name"][0] + config["Parameters"]["project_name"][1] + ".gff3" 
  cmd = "grep -v '#' " + gff + " | sed 's/transcript/mRNA/g' > Annotation" + config["Parameters"]["project_name"][1] + ".gff3;"
  create_browse_dictionary("Annotation"+ config["Parameters"]["project_name"][1], gff, cmd, "transcript")

if config["Parameters"]["run_non_coding"]:
  gff = config["Outputs"]["ncRNA_annotation_dir"] + config["Parameters"]["project_name"][0] + "nc" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3"
  cmd = "grep -v '#' " + gff + "| sed 's/transcript/ncRNA/g' > ncAnnotation" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3;"
  create_browse_dictionary("ncAnnotation"+ config["ncRNA_annotation"]["ncRNA_version"], gff, cmd, "ncRNA", "ncRNA")

if config["Outputs"]["Repeat_gff"]:
  create_browse_dictionary("Repeats", config["Outputs"]["Repeat_gff"], "", "feature3", "Interspersed_Repeat")

jbrowse_targets["seq"] = jbrowse_dir + "seq.ok"

##Execute rules
use rule browse_seq from jbrowse_workflow with:
  input:
    genome= genome
  output:
    checkpoint = jbrowse_dir + "seq.ok"
  conda:
    "../envs/JBrowse1.16.yaml"
  log:
    jlogs_dir + str(date) + ".j%j.browse_seq.out",
    jlogs_dir + str(date) + ".j%j.browse_seq.err",
  benchmark:
    jlogs_dir + str(date) + ".browse_seq.benchmark.txt" 
  threads: 2

use rule get_GCcontent from jbrowse_workflow with:
  input:
    genome = genome,
    lengths = config["Inputs"]["genome_lengths"] if config["Inputs"]["genome_lengths"] else ""
  output:
    wig = jbrowse_dir + "{name}.GC.wig",
    bw = jbrowse_dir + "{name}.GC.bw"
  params:
    scripts_dir = scripts_dir,
    window = 50,
    jbrowse_dir = jbrowse_dir
  log:
    jlogs_dir + str(date) + ".j%j.{name}.getGCcontent.out",
    jlogs_dir + str(date) + ".j%j.{name}.getGCcontent.err",
  benchmark:
    jlogs_dir + str(date) + ".{name}.getGCcontent.benchmark.txt"

use rule browse_tracks from jbrowse_workflow with:
  input:
    gff = lambda wildcards: browse[wildcards.prog]["input"] 
  output:
    checkpoint = jbrowse_dir + "browse.{prog}.ok"
  params:
    jbrowse_dir = jbrowse_dir,
    processed_gff = lambda wildcards: browse[wildcards.prog]["new_input"],
    type = lambda wildcards: browse[wildcards.prog]["type"] ,
    label = "{prog}",
    className = lambda wildcards: browse[wildcards.prog]["class"], 
    setup = lambda wildcards: browse[wildcards.prog]["setup_cmd"]
  conda:
    "../envs/JBrowse1.16.yaml"
  log:
    jlogs_dir + str(date) + ".j%j.browse_tracks.{prog}.out",
    jlogs_dir + str(date) + ".j%j.browse_tracks.{prog}.err",
  benchmark:
    jlogs_dir + str(date) + ".browse_tracks.{prog}.benchmark.txt" 
  threads: 8

use rule get_bw from jbrowse_workflow with:
  input:
    bam = config["Outputs"]["RNA_outdir"] + "{dir}/{sample}.bam",
    glen = config["Inputs"]["genome_lengths"] if config["Inputs"]["genome_lengths"] else ""
  output:
    bw = jbrowse_dir + "RNA/{dir}/{sample}.bw"
  params:
    opt = " -split ",
    dir = jbrowse_dir + "RNA/{dir}",
    bg = jbrowse_dir + "RNA/{dir}/{sample}.bg",
  log:
    jlogs_dir + str(date) + ".j%j.{dir}_{sample}.get_bw.out",
    jlogs_dir + str(date) + ".j%j.{dir}_{sample}.get_bw.err",
  benchmark:
    jlogs_dir + str(date) + ".get_bw.{dir}_{sample}.benchmark.txt" 
  conda:
    "../envs/bedtools2.30.0.yaml"
  threads: 2

use rule get_tar from jbrowse_workflow with:
  input:        
    checkpoints = lambda wildcards: jbrowse_targets[wildcards.indir]
  output:
    tar = jbrowse_dir + "{indir}.tar.gz"
  params:
    jbdir = jbrowse_dir,
    dir = "{indir}",
  log:
    jlogs_dir + str(date) + ".j%j.{indir}.get_tar.out",
    jlogs_dir + str(date) + ".j%j.{indir}.get_tar.err",
  benchmark:
    jlogs_dir + str(date) + ".get_tar.{indir}.benchmark.txt" 
  threads: 1