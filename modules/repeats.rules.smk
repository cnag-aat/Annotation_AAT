from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule redmask:
  input:
    genome = "inputs/"
  output:
    dir_output = "step01_RepeatAnnotationPipeline_redMask/"
  params:
    word_length = 15,
    min_kmer = 3,
    additional_redmask_opts = "",
    scripts_dir = "../scripts/"
  threads: 1
  log:
    "logs/" + str(date) + ".j%j.red_rule.out",
    "logs/" + str(date) + ".j%j.red_rule.err",
  conda:
    "../envs/redmask.yaml"
  shell:
    "mkdir -p {output.dir_output};"
    "Red -gnm {input.genome} -msk {output.dir_output} -rpt {output.dir_output} -cnd {output.dir_output}"
    " -len {params.word_length} -min {params.min_kmer} {params.additional_redmask_opts} -frm 2;"

