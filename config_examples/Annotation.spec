{
  "all": {
    "name": "{rule}_{base}_annotation_pipeline",
    "qos": "test",
    "time": "00:05:00",
    "queue": "genD"
  },
  "get_chunks_fasta": {
    "name": "{rule}_{base}_{input.fasta}",
    "qos": "test",
    "time": "00:10:00",
    "queue": "genD",
    "mem": "10G"
  },
  "get_coding_junctions": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "01:00:00",
    "queue": "genD",
    "mem": "1000"
  },
  "augustus": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "10:00:00",
    "queue": "genD",
    "mem": "50G",
    "array": "1-{masked_chunks}%15"
  },
  "augustus_hints": {
    "name": "{rule}_{base}",
    "qos": "long",
    "time": "24:00:00",
    "queue": "genD",
    "mem": "50G",
    "array": "1-{masked_chunks}%15"
  },
  "geneid": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "1:00:00",
    "queue": "genD",
    "mem": "15G"
  },
  "geneid_introns": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "3:00:00",
    "queue": "genD",
    "mem": "50G"
  },
  "genemark": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "12:00:00",
    "queue": "genD",
    "mem": "20G"
  },
  "genemark_ET": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "12:00:00",
    "queue": "genD",
    "mem": "50G"
  },
  "pasa": {
    "name": "{rule}_{base}",
    "qos": "long",
    "time": "16:05:00",
    "queue": "genD",
    "mem": "20G"
  },
  "Transdecoder": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "10:00:00",
    "queue": "genD",
    "mem": "10G"
  },
  "miniprot": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "10:00:00",
    "queue": "genD",
    "mem": "20G"
  },
  "merge_gffs": {
    "name": "{rule}_{base}_{wildcards.name}",
    "qos": "test",
    "time": "00:05:00",
    "queue": "genD",
    "mem": "100"
  },
  "predictions4EVM": {
    "name": "{rule}_{base}_{wildcards.name}",
    "qos": "test",
    "time": "00:10:00",
    "queue": "genD",
    "mem": "1000"
  },
  "prepare_evm": {
    "name": "{rule}_{base}",
    "qos": "test",
    "time": "00:10:00",
    "queue": "genD",
    "mem": "1000"
  },
  "EVM": {
    "name": "{rule}_{base}_{wildcards.w}",
    "qos": "long",
    "time": "24:00:00",
    "queue": "genD",
    "mem": "30G"
  },
  "select_EVM": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "01:00:00",
    "queue": "genD",
    "mem": "1000"
  },
  "PASA_update": {
    "name": "{rule}_{wildcards.round}_{base}",
    "qos": "long",
    "time": "24:00:00",
    "queue": "genD",
    "mem": "15G"
  },
  "process_update": {
    "name": "{rule}_{base}_{params.project}{params.version}",
    "qos": "normal",
    "time": "3:00:00",
    "queue": "genD",
    "mem": "5G"
  },
  "cmsearch": {
    "name": "{rule}_{base}",
    "qos": "vlong",
    "time": "48:00:00",
    "queue": "genD",
    "mem": "100G"
  },
  "tRNAscan": {
    "name": "{rule}_{base}",
    "qos": "normal",
    "time": "6:00:00",
    "queue": "genD",
    "mem": "5G"
  },
  "lncRNAannotation": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "0:30:00",
    "queue": "genD",
    "mem": "10G"
  },
  "Blast_prot": {
    "name": "{rule}_{base}_{wildcards.dirs}{wildcards.i}",
    "qos": "normal",
    "time": "3:30:00",
    "queue": "genD",
    "mem": "10G"
  },
  "ncAnnotation": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "0:30:00",
    "queue": "genD",
    "mem": "10G"
  },
  "get_GCcontent": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "3:00:00",
    "queue": "genD",
    "mem": "100G"
  },
  "browse_seq": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "3:00:00",
    "queue": "genD",
    "mem": "10G"
  },
  "browse_tracks": {
    "name": "{rule}_{base}",
    "qos": "vshort",
    "time": "1:00:00",
    "queue": "genD",
    "mem": "10G"
  },
  "get_tar": {
    "name": "{rule}_{base}",
    "qos": "short",
    "time": "3:00:00",
    "queue": "genD",
    "mem": "10G"
  }
}
