#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re
from pathlib import Path

##Author: Jessica Gomez-Garrido
##CNAG
##email: jessica.gomez@cnag.eu

def get_wildcards(dir, wildcards, ext, odir, sam):
  wsamples = []
  for r, d, f in os.walk(dir):
    for file in f:
      if file.endswith(ext):
        a = file.replace(ext,'')
        if wildcards == None:
          wildcards = a
        else:
          wildcards += "," + a
        wsamples.append(odir + a + sam + "\t" + a)

  return wildcards, wsamples

class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.configFile = "Annotation.config"	           #Name of the json configuration file with the pipeline parameters to be created
        self.specFile = "Annotation.spec"                  #Name of the spec file to be created
        self.base_name = None                              #Assembly base_name        
        self.pipeline_workdir = os.getcwd() + "/" 
        self.annotation_version = "01"	#Version of the annotation process.
        self.project_name = None                           #Name of the project and version, to give the names to the final annotation output. 
        self.run_redmask = True                           #By default run redmask step
        self.run_pasa = True                               #By default run pasa step
        self.run_transdecoder = True                       #By default run transdecoder step
        self.run_miniprot = True                           #By default run miniprot step
        self.run_augustus = True                           #By default run augustus step
        self.run_augustus_hints = True                     #By default run augustus with hints step
        self.run_geneid = True                             #By default run geneid step
        self.run_geneid_introns = True                     #By default run geneid with introns step
        self.run_genemark = True                           #By default run genemark step
        self.run_genemark_ET = True                        #By default run genemark-ET step
        self.run_EVM = True                                #By default run EVM 
        self.run_update = True                             #By default run the Pasa Update step
        self.run_non_coding = True                         #By default run non_coding step
        self.get_Jbrowse = True                            #By default run get JBrowse step        
        self.keep_intermediate = False

        #ALL SPEC PARAMETERS
        self.all_qos = "test"
        self.all_time = "00:05:00"
        self.all_queue = "genD"

        #INPUT PARAMETERS
        self.genome = None 
        self.scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"                          #Directory with the different scripts for the pipeline
        self.illumina_dir = None                      #Directory with the RNAseq reads
        self.cDNA_dir = None                          #Directory with the cDNA reads
        self.dRNA_dir = None                          #Directory with the dRNA reads
        self.pb_dir = None                            #Directory with the pacbio isoseq reads
        self.rna_samples_tsv = "samples.tsv"          #Samples tsv file to use with ESPRESSO
        self.transcripts = None                       #Path to the fasta with transcript evidence.
        self.pasa_config = None                       #Path to the Pasa configuration file.
        self.update_config = None                     #Path to the Pasa update configuration file.
        self.proteins = None                          #Path to the fasta with protein evidence.
        self.glen = None                              #Path to the assembly.genome file

        #OUTPUT PARAMETERS
        self.genome_masked = None	                       #Path to the masked genome assembly in fasta format.
        self.repeat_annotation_dir = "step01_Repeat_Annotation.V" + self.annotation_version + "/"	#Directory to run the repeat annotation on
        self.Repeat_gff = self.repeat_annotation_dir + "Repeats.4jb.gff3" #Path to the Repeat Annotation gff output.
        self.rna_outdir = "step02_RNAprocessing.V" + self.annotation_version + "/"  #Directory to run the rna step on 
        self.TACO_dir = self.rna_outdir                 #Directory to run TACO
        self.junctions = "alljunctions.final.gff3"      #Final file with the junctions
        self.RNAmodels = "TACO_assembled.gtf"             #Final file with the TACO results
        self.annotation_basedir = "step03_annotation_pipeline.V" + self.annotation_version + "/"  #Directory to run the base annotation steps on
        self.dir_masked_chunks = "chunks/masked_reference/"	#Directory to keep the chunks of the masked genome
        self.pasa_dir =  ""     #Directory to keep all the pasa outputs.   
        self.miniprot_cds =  self.annotation_basedir + "/protein_and_transcript_mappings/miniprot/proteins_miniprot_cds.gff3" #Output file for the miniprot output in a cds gff3 format.
        self.miniprot_gene =  self.annotation_basedir + "/protein_and_transcript_mappings/miniprot/proteins_miniprot_gene.gff3"    #Output file for the miniprot output in a gene gff3 format.
        self.EVM_dir = "step04_EVM.V" + self.annotation_version	#Directory to keep the files for the EVM step
        self.augustus_prediction = self.annotation_basedir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"  #Output file for the augustus predictions.
        self.augustus_preEVM = self.annotation_basedir + "/gene_predictions/augustus/augustus_preEVM.gff3"  #Output file for the augustus predictions converted for EVM.
        self.geneid_prediction = self.annotation_basedir + "/gene_predictions/geneid_gene_prediction.gff3"  #Output file for the geneid predictions.
        self.geneid_preEVM = self.annotation_basedir + "/gene_predictions/geneid_preEVM.gff3"  #Output file for the geneid predictions converted for EVM.
        self.genemark_prediction = self.annotation_basedir + "/gene_predictions/genemark.gtf"  #Output file for the genemark predictions.
        self.genemark_preEVM = self.annotation_basedir + "/gene_predictions/genemark_preEVM.gff3"  #Output file for the genemark predictions converted for EVM.
        self.geneid_introns_prediction = self.annotation_basedir + "/gene_predictions/geneid_introns_gene_prediction.gff3"  #Output file for the geneid with introns predictions.
        self.geneid_introns_preEVM = self.annotation_basedir + "/gene_predictions/geneid_introns_preEVM.gff3"  #Output file for the geneid with introns predictions converted for EVM.
        self.genemark_ET_prediction = self.annotation_basedir + "/gene_predictions/genemark-ET.gtf"  #Output file for the genemark-ET predictions.
        self.genemark_ET_preEVM = self.annotation_basedir + "/gene_predictions/genemark-ET_preEVM.gff3"  #Output file for the genemark-ET predictions converted for EVM.
        self.augustus_hints_prediction = self.annotation_basedir + "/gene_predictions/augustus_hints/augustus_hints_gene_prediction.gff3"  #Output file for the augustus with hints predictions.
        self.augustus_hints_preEVM = self.annotation_basedir + "/gene_predictions/augustus_hints/augustus_hints_preEVM.gff3"  #Output file for the augustus with hints predictions converted for EVM.
        self.evm_out = self.EVM_dir + "/evm.best.gff3"                                  #File with the final EVM models
        self.update_dir = "step05_annotation_update.V" + self.annotation_version  + "/"   #Directory to keep the files for annotation update step.
        self.ncRNA_annotation_dir = "step06_ncRNA_annotation.V" + self.annotation_version  + "/"   #Directory to keep the files of the ncRNA annotation step.
        
        #CHUNKS PARAMETERS
        self.masked_chunks = "23"	#Number of chunks of the masked genome for parallelizing some gene predictors run.
        self.protein_chunks = "20"	#Number of chunks to split the protein files for running blast and classify the lncRNAs. 
        
        #GETCHUNKS SPEC PARAMETERS
        self.chunks_qos = "test"
        self.chunks_time = "00:10:00"
        self.chunks_queue = "genD"
        self.chunks_mem = "10G"

        #REPEAT ANNOTATION PARAMETERS
        self.repeat_library = None                                                 #fasta file containing a pre-existant library of repeats 
        self.species_rdatabase = None                                            #Existant database to run a first time Repeat Masker
        self.rmaskCores = 8                                                         #Number of threads to run RepeatMasker
        
        #REPEAT MASKER SPEC PARAMETERS
        self.rmask_qos = "normal"
        self.rmask_time = "12:00:00"
        self.rmask_queue = "genD"
        self.rmask_mem = "100G"

        #REDMASK PARAMETERS
        self.red_wordlen = 15
        self.red_minkmer = 3
        self.red_add_opts = ""

        #REDMASK SPEC PARAMETERS
        self.redmask_qos = "normal"
        self.redmask_time = "12:00:00"
        self.redmask_queue = "genD"
        self.redmask_mem = "100G"

        #BLAST PARAMETERS
        self.blastdb = "/scratch/project/devel/aateam/blastdbs/swissprot"          #Blast database to check presence of protein families in RedMask library
        self.evalue = 0.000001
        self.blastCores = 24

        #BLAST SPEC PARAMETERS
        self.blast_qos = "long"
        self.blast_time = "24:00:00"
        self.blast_queue = "genD"
        self.blast_mem = "100G"

        #GET_REPEATS_GFF SPEC PARAMETERS
        self.repgff_qos = "short"
        self.repgff_time = "1:00:00"
        self.repgff_queue = "genD"
        self.repgff_mem = "10G"        

        #TRIMGALORE PARAMETERS
        self.trim_galore_opts = "--max_n 0 --gzip -q 20 --paired --retain_unpaired"
        self.Trim_Illumina_cores = 4    

        #ILLUMINA RNA PARAMETERS
        self.star_genome_dir = self.rna_outdir + "star_genome"                    #Directory for the genome index   
        self.illumina_PE = True
        self.star_dir = self.rna_outdir + "star"                        #Directory for mapping illumina
        self.starCores = 4
        self.indexstar_additional_options = ""             #Additional options to run star index with
        self.star_additional_options = ""             #Additional options to run star with
        self.stringtie_illumina_opts = ""             #Options to run stringtie in illumina mappings
        self.TACO_illumina_opts = ""                  #Options to run TACO in illumina models

        #STAR INDEX SPEC PARAMETERS
        self.star_index_qos = "normal"
        self.star_index_time = "12:00:00"
        self.star_index_queue = "genD"
        self.star_index_mem = "100G"

        #TRIMGALORE SPEC PARAMETERS
        self.trimgalore_qos = "normal"
        self.trimgalore_time = "3:00:00"
        self.trimgalore_queue = "genD"
        self.trimgalore_mem = "50G"

        #STAR SPEC PARAMETERS
        self.star_qos = "normal"
        self.star_time = "8:00:00"
        self.star_queue = "genD"
        self.star_mem = "500G"

        #cDNA RNA PARAMETERS
        self.cdna_minimap_dir = self.rna_outdir + "cDNA"         #Directory for the cDNA mappimgs
        self.stringtie_cDNA_opts = "--conservative -R" #Options to run stringtie in cDNA mappings
        self.minimap2_cDNA_opts = ""
        self.TACO_cDNA_opts = "--isoform-frac 0.01"    #Options to run TACO in cDNA models

        #dRNA RNA PARAMETERS
        self.drna_minimap_dir = self.rna_outdir + "dRNA"                                       #Directory for the dRNA mappimgs 
        self.stringtie_dRNA_opts = ""                                        #Options to run stringtie in dRNA mappings
        self.minimap2_dRNA_opts = " -uf -k14"
        self.TACO_dRNA_opts = "--isoform-frac 0.01 --filter-min-expr 0.2"    #Options to run TACO in dRNA models

        #PB RNA PARAMETERS
        self.pb_minimap_dir = "isoseq"                 #Directory for the pacbio isoseq mappimgs
        self.stringtie_pb_opts = "--conservative -R" #Options to run stringtie in pb isoseq mappings
        self.minimap2_pb_opts = ":hq -uf"
        self.TACO_pb_opts = "--isoform-frac 0.01"    #Options to run TACO in pb isoseq models

        #MINIMAP SPEC PARAMETERS
        self.minimap_qos = "long"
        self.minimap_time = "20:00:00"
        self.minimap_queue = "genD"
        self.minimap_mem = "100G"

        #MODEL RNA PARAMETERS
        self.minimapCores = 8
        self.TACO_global_opts = "--isoform-frac 0 --filter-min-expr 0" 
        self.tacoCores = 4
        self.espressoCores = 4
        self.espresso_path = "/software/assembly/src/ESPRESSO/espresso_v_1_3_0_beta/src/"
        self.espresso_outdir = "ESPRESSO_out/"

        #STRINGTIE SPEC PARAMETERS
        self.stringtie_qos = "normal"
        self.stringtie_time = "6:00:00"
        self.stringtie_queue = "genD"
        self.stringtie_mem = "50G"

        #TACO SPEC PARAMETERS
        self.taco_qos = "normal"
        self.taco_time = "12:00:00"
        self.taco_queue = "genD"
        self.taco_mem = "150G"

        #BAM2SAM SPEC PARAMETERS
        self.bam2sam_qos = "short"
        self.bam2sam_time = "3:00:00"
        self.bam2sam_queue = "genD"
        self.bam2sam_mem = "15G"

        #ESPRESSO SPEC PARAMETERS
        self.espresso_qos = "normal"
        self.espresso_time = "10:00:00"
        self.espresso_queue = "genD"
        self.espresso_mem = "500G"

        #GET CODING JUNCTIONS SPEC PARAMETERS
        self.codingHints_qos = "short"
        self.codingHints_time = "01:00:00"
        self.codingHints_queue = "genD"
        self.codingHints_mem = "10G"

        #PASA PARAMETERS
        self.pasadb = ""                                 #Name of the pasa database, it must coincide with the name given in pasa-config.
        self.pasa_weights = [8, 10, 8]                     #Weights given to pasa mappings when running EVM.
        self.create_database = False                       #By default do not create pasa database.
        self.aligners = 'gmap'                             #Program to use to align the transcripts     
        self.add_option = ""                               #Option used to add extra options for PASA
        self.update_rounds = 2                             #Number of rounds of PASA Updates to run

        #PASA SPEC PARAMETERS
        self.PASA_qos = "long"
        self.PASA_time = "16:05:00"
        self.PASA_queue = "genD"
        self.PASA_mem = "100G"       

        #TRANSDECODER PARAMETERS
        self.transdecoder_weights = [3, 2, 3]              #Weights given to the pasa transdecoder output gff3 file.

        #TRANSDECODER SPEC PARAMETERS
        self.transdecoder_qos = "normal"
        self.transdecoder_time = "10:00:00"
        self.transdecoder_queue = "genD"
        self.transdecoder_mem = "50G"  

        #MINIPROT PARAMETERS
        self.miniprot_path = "/software/assembly/src/miniprot/miniprot/"    ##path to miniprot installation
        self.miniprot_cores = 8                               #Number of threads.
        self.miniprot_weights = [10, 8, 10]                   #Weights given to miniprot mappings when running EVM.
        self.additional_miniprot_options = ""               #Additional miniprot options to run it, see miniprot help for more information. 

        #MINIPROT SPEC PARAMETERS
        self.miniprot_qos = "normal"
        self.miniprot_time = "10:00:00"
        self.miniprot_queue = "genD"
        self.miniprot_mem = "20G" 

        #GET_TRAINING_CANDIDATES SPEC PARAMETERS
        self.getcand_qos = "normal"
        self.getcand_time = "12:00:00"
        self.getcand_queue = "genD"
        self.getcand_mem = "100G"  

        #TRAIN AUGUSTUS SPEC PARAMETERS
        self.trainaug_qos = "long"
        self.trainaug_time = "24:00:00"
        self.trainaug_queue = "genD"
        self.trainaug_mem = "100G"  

        #AUGUSTUS PARAMETERS
        self.aug_species = None	#Species name to run augustus with its trained parameters. For augustus and augustus with hints steps.
        self.aug_config_path = "/software/assembly/conda/augustus3.5.0/config/"
        self.aug_optimize_threads = 24   #Number of threads to run the parameter optimization during the training step.
        self.aug_alternatives_from_evidence = "true"	#Report alternative transcripts when they are suggested by hints.
        self.aug_alternatives_from_sampling = "true"	#Report alternative transcripts generated through probabilistic sampling.
        self.aug_uniqueGeneId = "true"	#If true, output gene identifyers like this: seqname.gN. For augustus and augustus with hints.  
        self.aug_gff3 = "on"	#Output in gff3 format. For augustus and augustus with hints.
        self.aug_sample = 60	#For augustus and augustus with hints.
        self.aug_noInFrameStop = "true"	#Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. For augustus and augustus with hints.
        self.aug_maxtracks = 2	#Maximum number of tracks allowed. For augustus and augustus with hints.
        self.aug_singlestrand = "false"	#Predict genes independently on each strand, allow overlapping genes on opposite strands. For augustus and augustus with hints.
        self.aug_strand= "both"	#For augustus and augustus with hints.
        self.aug_min_intron_len= 30	#Minimum predicted intron length. For augustus and augustus with hints.                         
        self.augustus_weights = [2, 2, 1]	#Weights given to augustus predictions when running EVM.
        self.additional_augustus_options = ""	#Additional augustus options to run it, see augustus help for more information.

        #AUGUSTUS SPEC PARAMETERS
        self.augustus_qos = "normal"
        self.augustus_time = "10:00:00"
        self.augustus_queue = "genD"
        self.augustus_mem = "15G"

        #AUGUSTUS hints PARAMETERS
        self.extrinsic_file_augustus_hints = "/software/assembly/conda/augustus3.5.0/config/extrinsic/extrinsic.E.cfg"          #Extrinsic file to use when running augustus with hints. For more information read augustus documentation.
        self.augustus_hints_weights = [3,3,3]            #Weights given to augustus with intron predictions when running EVM.
        self.additional_augustus_hints_options = " --softmasking=0 "    #Additional augustus options to run augustus with hints, see augustus help for more information.

        #AUGUSTUS hints SPEC PARAMETERS
        self.augustus_hints_qos = "long"
        self.augustus_hints_time = "24:00:00"
        self.augustus_hints_queue = "genD"
        self.augustus_hints_mem = "15G"

        #GENEID PARAMETERS
        self.geneid_parameters = None	   #Path to the geneid parameters file. For geneid, geneid with introns and framefixing (part of annotation update) steps.
        self.geneid_path = "/software/assembly/src/geneid/"
        self.geneid_weights = [2,1,2]                      #Weights given to geneid predictions when running EVM.
        self.geneid_options = " -3U "                        #Desired geneid options to run it, see geneid documentation for more information.

        #GENEID SPEC PARAMETERS
        self.geneid_qos = "short"
        self.geneid_time = "1:00:00"
        self.geneid_queue = "genD"
        self.geneid_mem = "15G"

        #GENEID INTRONS PARAMETERS
        self.geneid_introns_weights = [3,3,3]                #Weights given to geneid with intron predictions when running EVM.
        self.geneid_introns_options = " -3nU"                #Desired geneid options to run geneid with introns, see geneid documentation for more information.

        #GENEID INTRONS SPEC PARAMETERS
        self.geneid_introns_qos = "short"
        self.geneid_introns_time = "3:00:00"
        self.geneid_introns_queue = "genD"
        self.geneid_introns_mem = "50G"

       #GENEMARK PARAMETERS
        self.gmk_min_contig = 50000                        #Will ignore contigs shorter then min_contig in training
        self.gmk_max_contig = 5000000                      #will split input genomic sequence into contigs shorter than max_contig.
        self.gmk_max_gap = 5000                            #Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are       interpreted as standing within gaps 
        self.gmk_cores = 24                                #Number of threads for running genemark. 
        self.additional_genemark_options = ""              #Additional genemark options to run it, see genemark documentation for more information.
        self.genemark_weights = [1, 1, 1]                  #Weights given to genemark predictions when running EVM.

        #GENEMARK SPEC PARAMETERS
        self.genemark_qos = "normal"
        self.genemark_time = "12:00:00"
        self.genemark_queue = "genD"
        self.genemark_mem = "20G"   

        #GENEMARK-ET PARAMETERS
        self.genemark_ET_weights = [3,3,3]                 #Weights given to augustus with intron predictions when running EVM.
        self.additional_genemark_ET_options = ""           #Additional genemark-ET options to run it, see genemark documentation for more information.

        #GENEMARK ET SPEC PARAMETERS
        self.genemark_et_qos = "normal"
        self.genemark_et_time = "12:00:00"
        self.genemark_et_queue = "genD"
        self.genemark_et_mem = "50G"  

        #MERGEGFF SPEC PARAMETERS
        self.mergegff_qos = "test"
        self.mergegff_time = "00:10:00"
        self.mergegff_queue = "genD"
        self.mergegff_mem = "1G"

        #PRED4EVM SPEC PARAMETERS
        self.pred4evm_qos = "test"
        self.pred4evm_time = "00:10:00"
        self.pred4evm_queue = "genD"
        self.pred4evm_mem = "10G"

        #EVM PARAMETERS
        self.evm_path = "/software/assembly/src/EVM2.1.0/EVidenceModeler-v2.1.0/"       #Path to the EVM software installation
        self.evm_segmentsize = 2000000                                                  #Size of the genome partitions for EVM
        self.evm_overlapsize = 1000000                                                  #Size of the overlap between the different EVM partitions
        self.evm_cores = 24                                                             #Number of threads to run EVM 
        self.additional_evm_options = ""	#Additional augustus options to run it, see evm help for more information.
       
        #PREPARE EVM SPEC PARAMETERS
        self.prepevm_qos = "test"
        self.prepevm_time = "00:10:00"
        self.prepevm_queue = "genD"
        self.prepevm_mem = "5G"

        #EVM SPEC PARAMETERS
        self.evm_qos = "long"
        self.evm_time = "24:00:00"
        self.evm_queue = "genD"
        self.evm_mem = "50G"

        #SELECT EVM SPEC PARAMETERS
        self.selectevm_qos = "short"
        self.selectevm_time = "01:00:00"
        self.selectevm_queue = "genD"
        self.selectevm_mem = "5G"

        #PASA UPDATE SPEC PARAMETERS
        self.pasaupdate_qos = "long"
        self.pasaupdate_time = "24:00:00"
        self.pasaupdate_queue = "genD"
        self.pasaupdate_mem = "15G"

        #PROCESS UPDATE SPEC PARAMETERS
        self.processupdate_qos = "normal"
        self.processupdate_time = "3:00:00"
        self.processupdate_queue = "genD"
        self.processupdate_mem = "5G"      

        #ncRNA ANNOTATION PARAMETERS
        self.ncRNA_version = "A"                           #Version for the non-coding RNA annotation
        self.cmsearch_CPUs = 32                            #Number of CPUs to run cmsearch
        self.Rfam = "/scratch_isilon/groups/assembly/data/databases/RFAM/Rfam.cm"    #CM file with the Rfam library.

        #CMSEARCH SPEC PARAMETERS
        self.cmsearch_qos = "long"
        self.cmsearch_time = "24:00:00"
        self.cmsearch_queue = "genD"
        self.cmsearch_mem = "100G"

        #tRNAScan SPEC PARAMETERS
        self.tRNAscan_qos = "normal"
        self.tRNAscan_time = "6:00:00"
        self.tRNAscan_queue = "genD"
        self.tRNAscan_mem = "5G"

        #lncRNA SPEC PARAMETERS
        self.lncRNA_qos = "short"
        self.lncRNA_time = "0:30:00"
        self.lncRNA_queue = "genD"
        self.lncRNA_mem = "10G"        

        #BLAST PROT SPEC PARAMETERS
        self.blast_prot_qos = "normal"
        self.blast_prot_time = "3:30:00"
        self.blast_prot_queue = "genD"
        self.blast_prot_mem = "10G"

        #ncRNA SPEC PARAMETERS
        self.ncRNA_qos = "short"
        self.ncRNA_time = "0:30:00"
        self.ncRNA_queue = "genD"
        self.ncRNA_mem = "10G"  

        #get_GC SPEC PARAMETERS
        self.getGC_qos= "short"
        self.getGC_time = "3:00:00"
        self.getGC_queue = "genD"
        self.getGC_mem = "100G"
        
        #get_seq SPEC PARAMETERS
        self.getseq_qos= "short"
        self.getseq_time = "3:00:00"
        self.getseq_queue = "genD"
        self.getseq_mem = "10G"

        #get_tracks SPEC PARAMETERS
        self.gettracks_qos= "vshort"
        self.gettracks_time = "1:00:00"
        self.gettracks_queue = "genD"
        self.gettracks_mem = "10G"

        #get_bw SPEC PARAMETERS
        self.getbw_qos= "vshort"
        self.getbw_time = "1:00:00"
        self.getbw_queue = "genD"
        self.getbw_mem = "10G"

        #get_tar SPEC PARAMETERS
        self.gettar_qos= "short"
        self.gettar_time = "3:00:00"
        self.gettar_queue = "genD"
        self.gettar_mem = "10G"

        #WILDCARDS
        self.illumina_fastqs = None                    #List with basename of the illumina fastqs
        self.cDNA_fastqs = None                        #List with basename of the cDNA fastqs  
        self.dRNA_fastqs = None                        #List with basename of the dRNA fastqs
        self.pb_fastqs = None                          #List with basename of the pacbio isoseq fastqs
        self.pb_fastas = None                          #List with basename of the pacbio isoseq fasta files 

        ###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.allSpecParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.chunksParameters = {}
        self.chunksSpecParameters = {}
        self.repeatParameters = {}
        self.rmaskSpecParameters = {}
        self.redParameters = {}
        self.redmaskSpecParameters = {}
        self.blastParameters = {}
        self.blastSpecParameters = {}
        self.repgffSpecParameters = {}
        self.illuminaParameters = {}
        self.trimgaloreParameters = {}
        self.starindexSpecParameters = {}
        self.trimgaloreSpecParameters = {}
        self.starSpecParameters = {}
        self.cdnaParameters = {}
        self.drnaParameters = {}
        self.pbParameters = {}
        self.minimapSpecParameters = {}
        self.modelRNAParameters = {}
        self.stringtieSpecParameters = {}
        self.tacoSpecParameters = {}
        self.bam2samSpecParameters = {}
        self.espressoSpecParameters = {}
        self.codingHintsSpecParameters = {}
        self.pasaSpecParameters={}
        self.pasaParameters={}
        self.transdecoderParameters = {}
        self.transdecoderSpecParameters={}
        self.miniprotParameters={}
        self.miniprotSpecParameters={}
        self.getcandSpecParameters={}
        self.trainaugSpecParameters={}
        self.augustusParameters = {}
        self.augustusSpecParameters = {}
        self.augustusArraySpecParameters = {}
        self.geneidParameters = {}
        self.geneidSpecParameters = {}
        self.geneidIntronsParameters = {}  
        self.geneidIntronsSpecParameters = {}
        self.genemarkParameters = {}
        self.genemarkSpecParameters = {}
        self.genemarkETParameters = {}
        self.genemarkETSpecParameters = {}
        self.augustusHintsParameters = {}
        self.augustusHintsSpecParameters = {}
        self.augustusHintsArraySpecParameters = {}
        self.mergegffSpecParameters = {}
        self.evmParameters = {}
        self.pred4evmSpecParameters = {}
        self.prepevmSpecParameters = {}
        self.evmSpecParameters = {}
        self.selectevmSpecParameters = {}
        self.updateParameters = {}
        self.pasaupdateSpecParameters = {}
        self.processupdateSpecParameters = {}
        self.ncRNAannotationParameters = {}
        self.cmsearchSpecParameters = {}
        self.tRNAscanSpecParameters = {}
        self.lncRNASpecParameters = {}
        self.BlastProtSpecParameters = {}
        self.ncRNASpecParameters = {}
        self.getGCSpecParameters = {}
        self.getseqSpecParameters = {}
        self.gettracksSpecParameters = {}
        self.getbwSpecParameters = {}
        self.gettarSpecParameters = {}
        self.wildcardParameters = {}

        ####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_chunks(parser)
        self.register_repeats(parser)
        self.register_red(parser)
        self.register_blast(parser)
        self.register_illumina(parser)
        self.register_trimgalore(parser)
        self.register_cdna(parser)
        self.register_drna(parser)
        self.register_pb(parser)
        self.register_modelRNA(parser)
        self.register_pasa(parser)
        self.register_transdecoder(parser)
        self.register_miniprot(parser)
        self.register_augustus(parser)
        self.register_augustus_hints(parser)
        self.register_geneid(parser)
        self.register_geneid_introns(parser)
        self.register_genemark(parser)
        self.register_genemark_ET(parser)
        self.register_evm(parser)
        self.register_ncRNA_annotation(parser)
        self.register_wildcards(parser)

        ###

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", default=self.configFile, help='Configuration file with the pipeline parameters to be created. Default %s' % self.configFile)
        general_group.add_argument('--specFile', dest="specFile", metavar="specFile", default=self.specFile, help='Cluster specifications JSON  fileto be generated. Default %s' % self.specFile)
        general_group.add_argument('--basename', dest="base_name", metavar="base_name", help='Assembly basename. Default %s' % self.base_name)
        general_group.add_argument('--pipeline-workdir', dest="pipeline_workdir", default=self.pipeline_workdir, help='Base directory for the pipeline run. Default %s' % self.pipeline_workdir)
        general_group.add_argument('--annotation-version', dest="annotation_version", default=self.annotation_version, help='Version of the annotation process. Default %s' % self.annotation_version)
        general_group.add_argument('--project-name', dest="project_name", metavar="project_name", nargs=2, help='Name of the project and version of the annotation space separated, to give the names to the final annotation output.')
        general_group.add_argument('--no-redmask', dest="run_redmask", action="store_false", help='If specified, do not run redmask step.')
        general_group.add_argument('--no-pasa', dest="run_pasa", action="store_false", help='If specified, do not run pasa step.')
        general_group.add_argument('--no-transdecoder', dest="run_transdecoder", action="store_false", help='If specified, do not run transdecoder step.')
        general_group.add_argument('--no-miniprot', dest="run_miniprot", action="store_false", help='If specified, do not run miniprot step.')
        general_group.add_argument('--no-augustus', dest="run_augustus", action="store_false", help='If specified, do not run augustus step.')
        general_group.add_argument('--no-augustus-hints', dest="run_augustus_hints", action="store_false", help='If specified, do not run augustus with hints step.')
        general_group.add_argument('--no-geneid', dest="run_geneid", action="store_false", help='If specified, do not run geneid step.')
        general_group.add_argument('--no-geneid-introns', dest="run_geneid_introns", action="store_false", help='If specified, do not run geneid with introns.')
        general_group.add_argument('--no-genemark', dest="run_genemark", action="store_false", help='If specified, do not run genemark step.')
        general_group.add_argument('--no-genemark-ET', dest="run_genemark_ET", action="store_false", help='If specified, do not run genemark-ET step.')
        general_group.add_argument('--no-EVM', dest="run_EVM", action="store_false", help='If specified, do not run EVM step.')
        general_group.add_argument('--no-update', dest="run_update", action="store_false", help='If specified, do not run the annotation update step.')
        general_group.add_argument('--no-noncoding', dest="run_non_coding", action="store_false", help='If specified, do not run the non_coding annotation step.')
        general_group.add_argument('--no-jbrowse', dest="get_Jbrowse", action="store_false", help='If specified, do not run the get jbrowse tracks step.')
        general_group.add_argument('--keep-intermediate', dest="keep_intermediate", action="store_true", help='If specified, do not delete intermediate files.')

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--scripts-dir', dest="scripts_dir", help='Directory with the different scripts for the pipeline. Default %s' % self.scripts_dir)
        input_group.add_argument('--genome', dest="genome", metavar="genome", help='Path to the genome assembly in fasta format.')
        input_group.add_argument('--glen', dest="glen", metavar="glen", help='Path to the assembly.genome file.')
        input_group.add_argument('--illumina-dir', dest="illumina_dir", help='Directory where the illumina fastqs are stored. Default %s' % self.illumina_dir)
        input_group.add_argument('--cdna-dir', dest="cDNA_dir", help='Directory where the cDNA fastqs are stored. Default %s' % self.cDNA_dir)
        input_group.add_argument('--drna-dir', dest="dRNA_dir", help='Directory where the dRNA fastqs are stored. Default %s' % self.dRNA_dir)
        input_group.add_argument('--pb-dir', dest="pb_dir", help='Directory where the pacbio isoseq fastqs are stored. Default %s' % self.pb_dir)
        input_group.add_argument('--rna-samples-tsv', dest="rna_samples_tsv", metavar="rna_samples_tsv",  help='TSV file describing the samples to use with ESPRESSO. Default %s' % self.rna_samples_tsv) 
        input_group.add_argument('--transcripts', dest="transcripts", metavar="transcripts", help='Path to the fasta with transcript evidence.')
        input_group.add_argument('--pasa-config', dest="pasa_config", metavar="pasa_config", help='Path to the pasa configuration file.')
        input_group.add_argument('--update-config', dest="update_config", metavar="update_config", help='Path to the pasa update configuration file.')
        input_group.add_argument('--proteins', dest="proteins", metavar="proteins", help='Path to the fasta with protein evidence.')
        
    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--genome-masked', dest="genome_masked", metavar="genome_masked", help='Path to the masked genome assembly in fasta format.')
        output_group.add_argument('--repeat-dir', dest="repeat_annotation_dir", help='Directory to the repeat annotation on. Default %s' % self.repeat_annotation_dir)
        output_group.add_argument('--Rep-gff', dest="Repeat_gff", metavar="Repeat_gff", help='Path to the Repeat Annotation gff output.')
        output_group.add_argument('--rna-outdir', dest="rna_outdir", help='Directory to the RNA processing step on. Default %s' % self.rna_outdir)
        output_group.add_argument('--TACO-dir', dest="TACO_dir", metavar="TACO_dir", default=self.TACO_dir, help='Directory to run the final TACO step. Default %s' % self.TACO_dir)
        output_group.add_argument('--junctions', dest="junctions", metavar="junctions", help='Path to the final junctions file. Default %s' % self.junctions)
        output_group.add_argument('--gtf-models', dest="RNAmodels", metavar="RNAmodels",help='Path to the final TACO gtf. Default %s' % self.RNAmodels)
        output_group.add_argument('--dir-masked-chunks', dest="dir_masked_chunks", help='Directory to keep the chunks of the masked genome. Default %s' % self.dir_masked_chunks)
        output_group.add_argument('--annot-dir', dest="annotation_basedir", help='Directory to keep the base annotation steps.  Default %s' % self.annotation_basedir)
        output_group.add_argument('--EVM-dir', dest="EVM_dir", help='Directory to keep the files for the EVM step. Default %s' % self.EVM_dir)
        output_group.add_argument('--miniprot-cds', dest="miniprot_cds", help='Output file for the miniprot output in a cds gff3 format.  Default %s' % self.miniprot_cds)
        output_group.add_argument('--miniprot-gene', dest="miniprot_gene", help='Output file for the miniprot output in a gene gff3 format.  Default %s' % self.miniprot_gene)
        output_group.add_argument('--augustus-prediction', dest="augustus_prediction",  help='Output file for the augustus predictions.  Default %s' % self.augustus_prediction)
        output_group.add_argument('--augustus-preEVM', dest="augustus_preEVM", help='Output file for the augustus predictions converted for EVM.  Default %s' % self.augustus_preEVM)        
        output_group.add_argument('--geneid-prediction', dest="geneid_prediction", help='Output file for the geneid predictions.  Default %s' % self.geneid_prediction)
        output_group.add_argument('--geneid-preEVM', dest="geneid_preEVM", help='Output file for the geneid predictions converted for EVM.  Default %s' % self.geneid_preEVM)
        output_group.add_argument('--genemark-prediction', dest="genemark_prediction", help='Output file for the genemark predictions.  Default %s' % self.genemark_prediction)
        output_group.add_argument('--genemark-preEVM', dest="genemark_preEVM", help='Output file for the genemark predictions converted for EVM.  Default %s' % self.genemark_preEVM)
        output_group.add_argument('--geneid-introns-prediction', dest="geneid_introns_prediction", help='Output file for the geneid with introns predictions.  Default %s' % self.geneid_introns_prediction)
        output_group.add_argument('--geneid-introns-preEVM', dest="geneid_introns_preEVM", help='Output file for the geneid with introns predictions converted for EVM.  Default %s' % self.geneid_introns_preEVM)
        output_group.add_argument('--genemark-ET-prediction', dest="genemark_ET_prediction", help='Output file for the genemark-ET predictions.  Default %s' % self.genemark_ET_prediction)
        output_group.add_argument('--genemark-ET-preEVM', dest="genemark_ET_preEVM", help='Output file for the genemark-ET predictions converted for EVM.  Default %s' % self.genemark_ET_preEVM)
        output_group.add_argument('--augustus-hints-prediction', dest="augustus_hints_prediction",  help='Output file for the augustus with hints predictions.  Default %s' % self.augustus_hints_prediction)
        output_group.add_argument('--augustus-hints-preEVM', dest="augustus_hints_preEVM", help='Output file for the augustus with hints predictions converted for EVM.  Default %s' % self.augustus_hints_preEVM)
        output_group.add_argument('--evm-out', dest="evm_out", help='File with the final EVM models. Default %s' % self.evm_out)
        output_group.add_argument('--update-dir', dest="update_dir", help='Directory to keep the files for annotation update step.   Default %s' % self.update_dir)
        output_group.add_argument('--ncRNA-dir', dest="ncRNA_annotation_dir", help='Directory to keep the files of the ncRNA annotation step. Default %s' % self.ncRNA_annotation_dir)
        
    def register_chunks(self, parser):
        """Register all parameters for making chunks
        with the given argparse parser

        parser -- the argparse parser
        """
        chunks_group = parser.add_argument_group('Chunks')
        chunks_group.add_argument('--masked-chunks', dest="masked_chunks", type=int, default=self.masked_chunks, help='Number of chunks of the masked genome for parallelizing some gene predictors run. Default %s' % self.masked_chunks)
        chunks_group.add_argument('--protein-chunks', dest="protein_chunks", type=int, default=self.protein_chunks, help='Number of chunks to split the protein files for running blast and classify the lncRNAs.  Default %s' % self.protein_chunks)

    def register_repeats(self, parser):
        """Register all repeat annotation parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        repeat_group = parser.add_argument_group('Repeat Annotation')
        repeat_group.add_argument('--species-rdatabase', dest="species_rdatabase", metavar="species_rdatabase", help='Existant database to run a first time Repeat Masker. Default %s' % self.species_rdatabase)
        repeat_group.add_argument('--repeat-library', dest="repeat_library", metavar="repeat_library", help='fasta file containing a pre-existant library of repeats.')
        repeat_group.add_argument('--rmask-cores', dest="rmaskCores", metavar="rmaskCores", type=int, default=self.rmaskCores, help='Default %s' % self.rmaskCores)
        
    def register_red(self, parser):
        """Register all redmask parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        red_group = parser.add_argument_group('RedMask')
        red_group.add_argument('--red-wordlen', dest="red_wordlen", metavar="red_wordlen", type=int, default=self.red_wordlen, help='Redmask wordlen parameter. Default %s' % self.red_wordlen)
        red_group.add_argument('--red-minkmer', dest="red_minkmer", metavar="red_minkmer", default=self.red_minkmer, type=int, help='Redmask minkmer parameter. Default %s' % self.red_minkmer)
        red_group.add_argument('--add-red-option', dest="red_add_opts",  default=self.red_add_opts, help='Option given to add extra options to RedMask')

    def register_blast(self, parser):
        """Register all blast parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        blast_group = parser.add_argument_group('BLAST')
        blast_group.add_argument('--blastdb', dest="blastdb", metavar="blastdb", default= self.blastdb, help='Blast database to check presence of protein families in RedMask library. Default %s' % self.blastdb) 
        blast_group.add_argument('--blast-eval', dest="evalue", metavar="evalue", default= self.evalue, help='Evalue to filter blast hits. Default %s' % self.evalue) 
        blast_group.add_argument('--blast-cores', dest="blastCores", metavar="blastCores", type=int, default=self.blastCores, help='Default %s' % self.blastCores)
                
    def register_illumina(self, parser):
        """Register all the illumina  parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        illumina_group = parser.add_argument_group('Illumina RNA')
        illumina_group.add_argument('--star-genome-dir', dest="star_genome_dir", metavar="star_genome_dir", help='Directory for the star genome index. Default %s' % self.star_genome_dir)
        illumina_group.add_argument('--no-pe', dest="illumina_PE", action="store_false", help='If specified, the input is not paired-end.')
        illumina_group.add_argument('--star-dir', dest="star_dir", metavar="star_dir", help='Directory for mapping the illumina reads. Default %s' % self.star_dir)
        illumina_group.add_argument('--star-cpu', dest="starCores", metavar="starCores", type=int, default = self.starCores, help='Number of threads to run star. Default %s' % self.starCores)
        illumina_group.add_argument('--indexstar-opts', dest="indexstar_additional_options", metavar="indexstar_additional_options", default=self.indexstar_additional_options, help='Additional options to run star index with. Default %s' % self.indexstar_additional_options)
        illumina_group.add_argument('--star-opts', dest="star_additional_options", metavar="star_additional_options", default=self.indexstar_additional_options, help='Additional options to run star with. Default %s' % self.star_additional_options)
        illumina_group.add_argument('--stringtie-illum-opts', dest="stringtie_illumina_opts", metavar="stringtie_illumina_opts", default=self.stringtie_illumina_opts, help='Options to run stringtie in illumina mappings. Default %s' % self.stringtie_illumina_opts)
        illumina_group.add_argument('--TACO-illum-opts', dest="TACO_illumina_opts", metavar="TACO_illumina_opts", default=self.TACO_illumina_opts, help='Options to run TACO in illumina mappings. Default %s' % self.TACO_illumina_opts)
        
    def register_trimgalore(self, parser):
        """Register all trimgalore parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        trimgalore_group = parser.add_argument_group('Trim_Galore')
        trimgalore_group.add_argument('--trim-galore-opts', dest="trim_galore_opts", metavar="trim_galore_opts", default=self.trim_galore_opts, help='Optional parameters for the rule trim_galore. Default %s' % self.trim_galore_opts)
        trimgalore_group.add_argument('--trim-Illumina-cores', type = int, dest="Trim_Illumina_cores", metavar="Trim_Illumina_cores", default=self.Trim_Illumina_cores, help='Number of threads to run the Illumina trimming step. Default %s' % self.Trim_Illumina_cores)

    def register_cdna(self, parser):
        """Register all the cdna parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        cdna_group = parser.add_argument_group('cDNA RNA')
        cdna_group.add_argument('--cdna-mappings', dest="cdna_minimap_dir", metavar="cdna_minimap_dir", help='Directory for the cDNA Minimap2 mappings. Default %s' % self.cdna_minimap_dir)
        cdna_group.add_argument('--minimap2-cdna-opts', dest="minimap2_cDNA_opts", metavar="minimap2_cDNA_opts", default=self.minimap2_cDNA_opts, help='Options to run minimap2 in cDNA mappings. Default %s' % self.minimap2_cDNA_opts)
        cdna_group.add_argument('--stringtie-cdna-opts', dest="stringtie_cDNA_opts", metavar="stringtie_cDNA_opts", default=self.stringtie_cDNA_opts, help='Options to run stringtie in cDNA mappings. Default %s' % self.stringtie_cDNA_opts)
        cdna_group.add_argument('--TACO-cdna-opts', dest="TACO_cDNA_opts", metavar="TACO_cDNA_opts", default=self.TACO_cDNA_opts, help='Options to run TACO in cDNA mappings. Default %s' % self.TACO_cDNA_opts)

    def register_drna(self, parser):
        """Register all the dRNA parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        drna_group = parser.add_argument_group('dRNA RNA')
        drna_group.add_argument('--drna-mappings', dest="drna_minimap_dir", metavar="drna_minimap_dir", help='Directory for the dRNA Minimap2 mappings. Default %s' % self.drna_minimap_dir)
        drna_group.add_argument('--minimap2-drna-opts', dest="minimap2_dRNA_opts", metavar="minimap2_dRNA_opts", default=self.minimap2_dRNA_opts, help='Options to run minimap2 in dRNA mappings. Default %s' % self.minimap2_dRNA_opts)
        drna_group.add_argument('--stringtie-drna-opts', dest="stringtie_dRNA_opts", metavar="stringtie_dRNA_opts", default=self.stringtie_dRNA_opts, help='Options to run stringtie in dRNA mappings. Default %s' % self.stringtie_dRNA_opts)
        drna_group.add_argument('--TACO-drna-opts', dest="TACO_dRNA_opts", metavar="TACO_dRNA_opts", default=self.TACO_dRNA_opts, help='Options to run TACO in dRNA mappings. Default %s' % self.TACO_dRNA_opts)

    def register_pb(self, parser):
        """Register all the pacbio parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        pb_group = parser.add_argument_group('Isoseq')
        pb_group.add_argument('--pb-mappings', dest="pb_minimap_dir", metavar="pb_minimap_dir", help='Directory for the Pacbio Isoseq Minimap2 mappings. Default %s' % self.pb_minimap_dir)
        pb_group.add_argument('--minimap2-pb-opts', dest="minimap2_pb_opts", metavar="minimap2_pb_opts", default=self.minimap2_pb_opts, help='Options to run minimap2 in Pacbio Isoseq mappings. Default %s' % self.minimap2_pb_opts)
        pb_group.add_argument('--stringtie-pb-opts', dest="stringtie_pb_opts", metavar="stringtie_pb_opts", default=self.stringtie_pb_opts, help='Options to run stringtie in Pacbio Isoseq mappings. Default %s' % self.stringtie_pb_opts)
        pb_group.add_argument('--TACO-pb-opts', dest="TACO_pb_opts", metavar="TACO_pb_opts", default=self.TACO_pb_opts, help='Options to run TACO in isoseq mappings. Default %s' % self.TACO_pb_opts)

    def register_modelRNA(self, parser):
        """Register all model RNA parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        modelRNA_group = parser.add_argument_group('Model RNA')
        modelRNA_group.add_argument('--minimap-cpu', dest="minimapCores", metavar="minimapCores", type=int, default = self.minimapCores, help='Number of threads to run Minimap2. Default %s' % self.minimapCores)
        modelRNA_group.add_argument('--taco-cpu', dest="tacoCores", metavar="tacoCores", type=int, default = self.tacoCores, help='Number of threads to run TACO. Default %s' % self.tacoCores)
        modelRNA_group.add_argument('--TACO-all-opts', dest="TACO_global_opts", metavar="TACO_global_opts", default=self.TACO_global_opts, help='Options to run TACO when merging all the datasets. Default %s' % self.TACO_global_opts)
        modelRNA_group.add_argument('--espresso-cpu', dest="espressoCores", metavar="espressoCores", type=int, default = self.espressoCores, help='Number of threads to run espresso. Default %s' % self.espressoCores)
        modelRNA_group.add_argument('--espresso-path', dest="espresso_path", metavar="espresso_path", default=self.espresso_path, help='Path to the ESPRESSO scripts. Default %s' % self.espresso_path)
        modelRNA_group.add_argument('--espresso-outdir', dest="espresso_outdir", metavar="espresso_outdir", help='Directory for runnning ESPRESSO. Default %s' % self.espresso_outdir)

    def register_pasa(self, parser):
        """Register all pasa parameters with given
        argparse parser

        parser -- the argparse parser
        """
        pasa_group = parser.add_argument_group('Pasa parameters')
        pasa_group.add_argument('--pasa-weights', dest="pasa_weights", nargs="+", type=int, default=self.pasa_weights, help='Weights given to pasa mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  8 10 8 ')
        pasa_group.add_argument('--create-database', dest="create_database", action="store_true", help='If specified, create pasa database.')
        pasa_group.add_argument('--aligners', dest="aligners", default=self.aligners, help='Program to map the transcripts with.')
        pasa_group.add_argument('--add-pasa-option', dest="add_option",  default=self.add_option, help='Option given to add extra options to PASA')
        pasa_group.add_argument('--update-rounds', dest="update_rounds", default = self.update_rounds, type=int, help='Number of rounds to run PASA updates. Default %s''' % str(self.update_rounds))

    def register_transdecoder(self, parser):
        """Register all transdecoder parameters with given
        argparse parser

        parser -- the argparse parser
        """
        transdecoder_group = parser.add_argument_group('Transdecoder parameters')
        transdecoder_group.add_argument('--transdecoder-weights', dest="transdecoder_weights", nargs="+", type=int, default=self.transdecoder_weights, help='Weights given to pasa transdecodergff3 output file  when running EVM. Specify the weight for each EVM run separated by a space. Example 3 2 3')

    def register_miniprot(self, parser):
        """Register all miniprot parameters with given
        argparse parser

        parser -- the argparse parser
        """
        miniprot_group = parser.add_argument_group('Miniprot parameters')
        miniprot_group.add_argument('--miniprot-path', dest="miniprot_path", default = self.miniprot_path,  help='Path to Miniprot installation. Default %s' % self.miniprot_path)
        miniprot_group.add_argument('--miniprot-cores', dest="miniprot_cores", type=int, default = self.miniprot_cores, help='''Number of threads. Default %s''' % str(self.miniprot_cores))
        miniprot_group.add_argument('--miniprot-weights', dest="miniprot_weights", nargs="+", type=int, default=self.miniprot_weights, help='Weights given to miniprot mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  10 8 10 ') 
        miniprot_group.add_argument('--additional-miniprot-options', dest="additional_miniprot_options", default=self.additional_miniprot_options, help='Additional miniprot options to run it, see miniprot help for more information about the possible options.')

    def register_augustus(self, parser):
        """Register all augustus parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_group = parser.add_argument_group('Augustus parameters')
        augustus_group.add_argument('--aug-species', dest="aug_species", metavar="aug_species", help='Species name to run augustus with its trained parameters. For augustus and augustus with hints steps.')
        augustus_group.add_argument('--aug-conf', dest="aug_config_path", metavar="aug_config_path", default = self.aug_config_path, help='''Path to the augustus config path. Default %s''' %str(self.aug_config_path))
        augustus_group.add_argument('--optimize-aug-cores', dest="aug_optimize_threads", type=int, default = self.aug_optimize_threads, help='''Number of threads to run the parameter optimization during the augustus training step. Default %s''' % str(self.aug_optimize_threads))
        augustus_group.add_argument('--aug-alternatives-from-sampling', dest="aug_alternatives_from_sampling", default=self.aug_alternatives_from_sampling, choices=['true', 'false'], help='''Report alternative transcripts generated through probabilistic sampling. Default %s''' % str(self.aug_alternatives_from_sampling))
        augustus_group.add_argument('--aug-alternatives-from-evidence', dest="aug_alternatives_from_evidence", default=self.aug_alternatives_from_evidence, choices=['true', 'false'], help='''Report alternative transcripts when they are suggested by hints. Default %s''' % str(self.aug_alternatives_from_evidence))
        augustus_group.add_argument('--aug-uniqueGeneId', dest="aug_uniqueGeneId", default=self.aug_uniqueGeneId, choices = ['true', 'false'], help='''If true, output gene identifyers like this: seqname.gN. For augustus and augustus with hints. Default %s''' % str(self.aug_uniqueGeneId))
        augustus_group.add_argument('--aug-gff3', dest="aug_gff3", default=self.aug_gff3, choices = ['ON', 'OFF', 'on', 'off'], help='''Output in gff3 format. For augustus and augustus with hints. Default %s''' % str(self.aug_gff3))
        augustus_group.add_argument('--aug-sample', dest="aug_sample", type=int, default=self.aug_sample, help='''For augustus and augustus with hints. Default %s''' % str(self.aug_sample))
        augustus_group.add_argument('--aug-noInFrameStop', dest="aug_noInFrameStop", default=self.aug_noInFrameStop, choices = ['true', 'false'], help='''Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. For augustus and augustus with hints. Default %s''' % str(self.aug_noInFrameStop))
        augustus_group.add_argument('--aug-maxtracks', dest="aug_maxtracks", type=int, default=self.aug_maxtracks, help='''Maximum number of tracks allowed. For augustus and augustus with hints. Default %s''' % str(self.aug_maxtracks))
        augustus_group.add_argument('--aug-singlestrand', dest="aug_singlestrand", default=self.aug_singlestrand, choices = ['true', 'false'], help='''Predict genes independently on each strand, allow overlapping genes on opposite strands. For augustus and augustus with hints. Default %s''' % str(self.aug_singlestrand))
        augustus_group.add_argument('--aug-strand', dest="aug_strand", default=self.aug_strand, choices=['both', 'forward', 'backward'], help='''For augustus and augustus with hints. Default %s''' % str(self.aug_strand))
        augustus_group.add_argument('--aug-min-intron-len', dest="aug_min_intron_len", type=int, default=self.aug_min_intron_len, help='''Minimum predicted intron length. For augustus and augustus with hints. Default %s''' % str(self.aug_min_intron_len))
        augustus_group.add_argument('--augustus-weights', dest="augustus_weights", nargs="+", type=int, default=self.augustus_weights, help='Weights given to augustus predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 2 1')
        augustus_group.add_argument('--additional-augustus-options', dest="additional_augustus_options", default=self.additional_augustus_options, help='Additional augustus options to run it, see augustus help for more information about the possible options.')

    def register_augustus_hints(self, parser):
        """Register all augustus with hints parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_hints_group = parser.add_argument_group('Augustus hints parameters')
        augustus_hints_group.add_argument('--extrinsic-file-augustus-hints', dest="extrinsic_file_augustus_hints", metavar="extrinsic_file_augustus_hints", default = self.extrinsic_file_augustus_hints, help='''Path to the Extrinsic file to use when running augustus with hints. Default %s''' % str(self.extrinsic_file_augustus_hints))
        augustus_hints_group.add_argument('--augustus-hints-weights', dest="augustus_hints_weights", nargs="+", type=int, default=self.augustus_hints_weights, help='Weights given to augustus with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        augustus_hints_group.add_argument('--additional-augustus-hints-options', dest="additional_augustus_hints_options", default=self.additional_augustus_hints_options, help='Desired augustus with intron options to run it, see augustus documentation for more information.''')

    def register_geneid(self, parser):
        """Register all geneid parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_group = parser.add_argument_group('Geneid parameters')
        geneid_group.add_argument('--geneid-path', dest="geneid_path", default=self.geneid_path, help='Path to the installation of geneid. Default %s''' % str(self.geneid_path))
        geneid_group.add_argument('--geneid-weights', dest="geneid_weights", nargs="+", type=int, default=self.geneid_weights, help='Weights given to geneid predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 1 2 ')
        geneid_group.add_argument('--geneid-options', dest="geneid_options", default=self.geneid_options, help='Desired geneid options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))
        geneid_group.add_argument('--geneid-parameters', dest="geneid_parameters", metavar="geneid_parameters", help='Path to the geneid parameters file. For geneid, geneid with introns and framefixing (part of annotation update) steps.')

    def register_geneid_introns(self, parser):
        """Register all geneid with introns parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_introns_group = parser.add_argument_group('Geneid Introns parameters')
        geneid_introns_group.add_argument('--geneid-introns-weights', dest="geneid_introns_weights", nargs="+", type=int, default=self.geneid_introns_weights, help='Weights given to geneid with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        geneid_introns_group.add_argument('--geneid-introns-options', dest="geneid_introns_options", default=self.geneid_introns_options,  help='Desired geneid with intron options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))

    def register_genemark(self, parser):
        """Register all genemark parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_group = parser.add_argument_group('Genemark parameters')
        genemark_group.add_argument('--gmk-min-contig', dest="gmk_min_contig", type=int, default=self.gmk_min_contig, help='''Will ignore contigs shorter then min_contig in training. Default %s''' % str(self.gmk_min_contig))
        genemark_group.add_argument('--gmk-max-contig', dest="gmk_max_contig", type=int, default=self.gmk_max_contig, help='''Will split input genomic sequence into contigs shorter than max_contig. Default %s''' % str(self.gmk_max_contig))
        genemark_group.add_argument('--gmk-max-gap', dest="gmk_max_gap", type=int, default=self.gmk_max_gap, help='''Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are interpreted as standing within gaps. Default %s''' % str(self.gmk_max_gap))
        genemark_group.add_argument('--gmk-cores', dest="gmk_cores", type=int, default=self.gmk_cores, help='''Number of threads for running genemark. Default %s''' % str(self.gmk_cores))
        genemark_group.add_argument('--additional-genemark-options', dest="additional_genemark_options", default=self.additional_genemark_options, help='Additional genemark options to run it, see genemark documentation for more information.')
        genemark_group.add_argument('--genemark-weights', dest="genemark_weights", nargs="+", type=int, default=self.genemark_weights, help='Weights given to genemark predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  1 1 1 ') 

    def register_genemark_ET(self, parser):
        """Register all genemark-ET parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_ET_group = parser.add_argument_group('Genemark-ET parameters')
        genemark_ET_group.add_argument('--additional-genemark-ET-options', dest="additional_genemark_ET_options", default=self.additional_genemark_options, help='Additional genemark-ET options to run it, see genemark documentation for more information.')
        genemark_ET_group.add_argument('--genemark-ET-weights', dest="genemark_ET_weights", nargs="+", type=int, default=self.genemark_ET_weights, help='Weights given to genemark-ET predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  3 3 3 ')

    def register_evm(self, parser):
        """Register all evm parameters with given
        argparse parser

        parser -- the argparse parser
        """
        evm_group = parser.add_argument_group('Evm parameters')
        evm_group.add_argument('--evm-path', dest="evm_path", default = self.evm_path, help='Path to the EVM software installation. Default %s''' % str(self.evm_path))
        evm_group.add_argument('--evm-segmentsize', dest="evm_segmentsize", default = self.evm_segmentsize,type = int, help='Size of the genome partitions for EVM. Default %s''' % str(self.evm_segmentsize))
        evm_group.add_argument('--evm-overlapsize', dest="evm_overlapsize", default = self.evm_overlapsize, type=int, help='Size of the overlap between the different EVM partitions. Default %s''' % str(self.evm_overlapsize))
        evm_group.add_argument('--evm-cores', dest="evm_cores", default = self.evm_cores, type=int, help='Number of threads to run EVM. Default %s''' % str(self.evm_cores))
        evm_group.add_argument('--additional-evm-options', dest="additional_evm_options", default = self.additional_evm_options, help='Additional evm options to run it, see evm help for more information about the possible options.')

    def register_ncRNA_annotation(self, parser):
        """Register all ncRNA annotation parameters with given
        argparse parser

        parser -- the argparse parser
        """
        ncRNA_group = parser.add_argument_group('ncRNA Annotation parameters')
        ncRNA_group.add_argument('--ncRNA-version', dest="ncRNA_version", default = self.ncRNA_version, help='Version for the ncRNA annotation. Default %s''' % str(self.ncRNA_version))
        ncRNA_group.add_argument('--Rfam', dest="Rfam", default = self.Rfam, help='CM file with the Rfam library. Default %s''' % str(self.Rfam))
        ncRNA_group.add_argument('--cmsearch-CPUs', dest="cmsearch_CPUs", type=int, default = self.cmsearch_CPUs, help='''Number of CPUs to run cmsearch Default %s''' % str(self.cmsearch_CPUs))

    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--illumina-reads', dest="illumina_fastqs", metavar="illumina_fastqs", help='List with basename of the illumina fastqs. Default %s' % self.illumina_fastqs)
        wildcards_group.add_argument('--cdna-reads', dest="cDNA_fastqs", metavar="cDNA_fastqs", help='List with basename of the cDNA fastqs. Default %s' % self.cDNA_fastqs)
        wildcards_group.add_argument('--drna-reads', dest="dRNA_fastqs", metavar="dRNA_fastqs", help='List with basename of the dRNA fastqs. Default %s' % self.dRNA_fastqs)
        wildcards_group.add_argument('--pb-reads', dest="pb_fastqs", metavar="pb_fastqs", help='List with basename of the Pacbio Isoseq fastqs. Default %s' % self.pb_fastqs)
        wildcards_group.add_argument('--pb-reads-fasta', dest="pb_fastas", metavar="pb_fastas", help='List with basename of the Pacbio Isoseq fastas. Default %s' % self.pb_fastas)

###
    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        if args.pipeline_workdir != None:
          args.pipeline_workdir = os.path.abspath(args.pipeline_workdir) + "/"
        else:
          args.pipeline_workdir = os.getcwd() + "/"

        if args.configFile != None:
          args.configFile = os.path.abspath(args.configFile)
        else:
          args.configFile = os.path.abspath(self.configFile)

        if args.specFile != None:
          args.specFile = os.path.abspath(args.specFile)
        else:
          args.specFile = os.path.abspath(self.specFile)   

        if args.project_name == None:
          parser.print_help()
          print ("ERROR Please specify project name and version") 
          sys.exit(-1)

        ##Checking inputs
        if args.scripts_dir:
          args.scripts_dir = os.path.abspath(args.scripts_dir) + "/"
        else:
          args.scripts_dir = os.path.abspath(self.scripts_dir) + "/"
        if not os.path.exists(args.scripts_dir):
          print (args.scripts_dir + " not found")

        if args.genome == None:
          parser.print_help()
          print ("ERROR Sorry! No genome fasta file defined")
          sys.exit(-1)
        else:
          args.genome = os.path.abspath(args.genome) 
          if not os.path.exists(args.genome):
            print (args.genome + " not found")

          if args.base_name == None:
            args.base_name = os.path.splitext(os.path.basename(args.genome))[0]

        if args.glen:
          if not os.path.exists(args.glen):
            print (args.glen) + " not found"
          else:
            args.glen = os.path.abspath(args.glen)

        if args.repeat_library:
          args.repeat_library = os.path.abspath(args.repeat_library) 
          if not os.path.exists(args.repeat_library):
            print (args.repeat_library + " not found")	
        
        if args.illumina_dir != None:
            args.illumina_dir = os.path.abspath(args.illumina_dir) + "/" 

        if args.cDNA_dir != None:
            args.cDNA_dir = os.path.abspath(args.cDNA_dir) + "/" 

        if args.dRNA_dir != None:
            args.dRNA_dir = os.path.abspath(args.dRNA_dir) + "/" 

        if args.pb_dir != None:
            args.pb_dir = os.path.abspath(args.pb_dir) + "/" 

        if args.proteins != None: 
          args.proteins = os.path.abspath(args.proteins) 
          if not os.path.exists(args.proteins):
            print (args.proteins + " not found")
        else:
          args.proteins = ""
          if args.run_miniprot == True: 
            print ("ERROR: Protein evidence file has not been specified, required to run Miniprot.")
        
        ##Getting outputs
        if args.repeat_annotation_dir:
            args.repeat_annotation_dir = os.path.abspath(args.repeat_annotation_dir) + "/"
        else:
            args.repeat_annotation_dir = args.pipeline_workdir + self.repeat_annotation_dir 

        if args.Repeat_gff:
          args.Repeat_gff = os.path.abspat(args.Repeat_gff)
        else:
          args.Repeat_gff = args.repeat_annotation_dir + "Repeats.4jb.gff3"

        if args.EVM_dir:
          args.EVM_dir = os.path.abspath(args.EVM_dir) + "/"
        else:
          args.EVM_dir = args.pipeline_workdir + "step04_EVM.V" + str(args.annotation_version)  + "/"

        if args.genome_masked:
          args.genome_masked = os.path.abspath(args.genome_masked)
        else:
          args.genome_masked = args.repeat_annotation_dir + args.base_name + ".masked.fa"

        if args.rna_outdir:
            args.rna_outdir = os.path.abspath(args.rna_outdir) + "/"
        else:
            args.rna_outdir = args.pipeline_workdir + self.rna_outdir 
        if not os.path.exists(args.rna_outdir):
          os.makedirs(args.rna_outdir)

        if args.star_genome_dir:
            args.star_genome_dir = os.path.abspath(args.star_genome_dir) + "/"
        else:
            args.star_genome_dir =args.rna_outdir+ "star_genome/"

        if args.star_dir:
            args.star_dir = os.path.abspath(args.star_dir) + "/"
        else:
            args.star_dir = args.rna_outdir + "star/"

        if args.cdna_minimap_dir:
            args.cdna_minimap_dir = os.path.abspath(args.cdna_minimap_dir) + "/"
        else:
            args.cdna_minimap_dir = args.rna_outdir + "cDNA/"

        if args.drna_minimap_dir:
            args.drna_minimap_dir = os.path.abspath(args.drna_minimap_dir) + "/"
        else:
            args.drna_minimap_dir = args.rna_outdir + "dRNA/"

        if args.pb_minimap_dir:
            args.pb_minimap_dir = os.path.abspath(args.pb_minimap_dir) + "/"
        else:
            args.pb_minimap_dir = args.rna_outdir + "isoseq/"

        if args.TACO_dir:
            args.TACO_dir = os.path.abspath(args.TACO_dir) + "/"
        else:
            args.TACO_dir = args.rna_outdir   

        if args.junctions:
          args.junctions = os.path.abspath(args.junctions) 
        else:
          args.junctions = args.rna_outdir + self.junctions    

        if args.RNAmodels:
          args.RNAmodels = os.path.abspath(args.RNAmodels) 
        elif args.illumina_dir or args.cDNA_dir or args.dRNA_dir or args.pb_dir:
          args.RNAmodels = args.rna_outdir + self.RNAmodels 
        elif args.run_pasa:
          print ("WARNING: no RNAmodels given, PASA will be run only with the given transcripts fasta file")

        args.espresso_path = os.path.abspath(args.espresso_path) + "/"

        if args.espresso_outdir == None:
          args.espresso_outdir = args.rna_outdir + self.espresso_outdir
        else:
          args.espresso_outdir = os.path.abspath(args.espresso_outdir)

        if args.dir_masked_chunks:
          args.dir_masked_chunks =os.path.abspath(args.dir_masked_chunks) + "/"
        else:
          args.dir_masked_chunks = args.pipeline_workdir + self.dir_masked_chunks

        if args.annotation_basedir:
          args.annotation_basedir = os.path.abspath(args.annotation_basedir) + "/"
        else:
          args.annotation_basedir =  args.pipeline_workdir + self.annotation_basedir

        if args.miniprot_cds:
          args.miniprot_cds = os.path.abspath(args.miniprot_cds)
        else: 
          args.miniprot_cds =  args.annotation_basedir + "protein_and_transcript_mappings/miniprot/proteins_miniprot_cds.gff3" 

        if args.miniprot_gene:
            args.miniprot_gene = os.path.abspath(args.miniprot_gene)
        else: 
            args.miniprot_gene =  args.annotation_basedir + "protein_and_transcript_mappings/miniprot/proteins_miniprot_gene.gff3" 

        if args.augustus_prediction:
            args.augustus_prediction = os.path.abspath(args.augustus_prediction)
        else:
            args.augustus_prediction =  args.annotation_basedir + "gene_predictions/augustus/augustus_gene_prediction.gff3"

        if args.augustus_preEVM:
            args.augustus_preEVM = os.path.abspath(args.augustus_preEVM)
        else:
            args.augustus_preEVM =  args.annotation_basedir + "gene_predictions/augustus/augustus_preEVM.gff3" 

        if args.geneid_prediction:
            args.geneid_prediction = os.path.abspath(args.geneid_prediction)
        else: 
            args.geneid_prediction = args.annotation_basedir + "gene_predictions/geneid_gene_prediction.gff3"

        if args.geneid_preEVM:
            args.geneid_preEVM = os.path.abspath(args.geneid_preEVM)
        else: 
            args.geneid_preEVM = args.annotation_basedir + "gene_predictions/geneid_preEVM.gff3"

        if args.genemark_prediction:
            args.genemark_prediction = os.path.abspath(args.genemark_prediction)
        else:
            args.genemark_prediction = args.annotation_basedir + "gene_predictions/genemark.gtf"

        if args.genemark_preEVM:
            args.genemark_preEVM = os.path.abspath(args.genemark_preEVM)
        else:
            args.genemark_preEVM = args.annotation_basedir + "gene_predictions/genemark_preEVM.gff3"

        if args.geneid_introns_prediction:
            args.geneid_introns_prediction = os.path.abspath(args.geneid_introns_prediction)
        else: 
            args.geneid_introns_prediction = args.annotation_basedir + "gene_predictions/geneid_introns_gene_prediction.gff3"

        if args.geneid_introns_preEVM:
            args.geneid_introns_preEVM = os.path.abspath(args.geneid_introns_preEVM)
        else: 
            args.geneid_introns_preEVM = args.annotation_basedir + "gene_predictions/geneid_introns_preEVM.gff3"

        if args.genemark_ET_prediction:
            args.genemark_ET_prediction = os.path.abspath(args.genemark_ET_prediction)
        else:
            args.genemark_ET_prediction = args.annotation_basedir + "gene_predictions/genemark-ET.gtf"

        if args.genemark_ET_preEVM:
            args.genemark_ET_preEVM = os.path.abspath(args.genemark_ET_preEVM)
        else:
            args.genemark_ET_preEVM = args.annotation_basedir + "gene_predictions/genemark-ET_preEVM.gff3"

        if args.augustus_hints_prediction:
            args.augustus_hints_prediction = os.path.abspath(args.augustus_hints_prediction)
        else:
            args.augustus_hints_prediction =  args.annotation_basedir + "gene_predictions/augustus_hints/augustus_hints_gene_prediction.gff3"

        if args.augustus_hints_preEVM:
            args.augustus_hints_preEVM = os.path.abspath(args.augustus_hints_preEVM)
        else:
            args.augustus_hints_preEVM =  args.annotation_basedir + "gene_predictions/augustus_hints/augustus_hints_preEVM.gff3" 

        if args.evm_out:
            args.evm_out = os.path.abspath(args.evm_out)
        else:
            args.evm_out = args.EVM_dir + "evm.best.gff3"

        if args.update_dir:
            args.update_dir = os.path.abspath(args.update_dir) + "/"
        else: 
            args.update_dir = args.pipeline_workdir + "step05_annotation_update.V" +str(args.annotation_version)  + "/" 

        if args.ncRNA_annotation_dir:
            args.ncRNA_annotation_dir = os.path.abspath(args.ncRNA_annotation_dir) + "/"
        else:
            args.ncRNA_annotation_dir = args.pipeline_workdir + "step06_ncRNA_annotation.V" + str(args.annotation_version)  + "/"

        ##Get parameters
        args.pasa_dir = self.pasa_dir
        args.pasadb = self.pasadb

        if args.run_pasa == True:
          if args.transcripts != None:
            args.transcripts = os.path.abspath(args.transcripts) 
            if not os.path.exists(args.transcripts):
              print (args.transcripts + " not found")

          if args.pasa_config == None:
            print ("ERROR No pasa configuration file found.")
          else:
            if not os.path.exists(args.pasa_config):
              print (args.pasa_config + " not found")
            else:
              args.pasa_config = os.path.abspath(args.pasa_config)

              with open(args.pasa_config,'r') as f:
                for line in f:
                  if re.search(r'DATABASE=',line):
                    m = line.replace('DATABASE=', '')
                    args.pasa_dir = os.path.dirname(m) + "/"
                    args.pasadb = os.path.basename(m).rstrip()
          if args.create_database:
            print ("WARNING: Remember to turn the create-database parameter to false if it's not the first time that you run PASA!")
          else:
            print ("WARNING: Remember to set the create-database parameter to true the first time you run PASA!")
        else:
          args.pasa_config = ""
          if args.run_transdecoder:
            print ("ERROR: in order to run transdecoder you first need to run PASA, please turn that option ON or do not run transdecoder!")

        if args.run_augustus or args.run_augustus_hints:
          args.aug_config_path = os.path.abspath(args.aug_config_path) + "/"
          if args.aug_species == None:
            print ("ERROR Augustus trained species has not been specified. Please, specify it, remember that if no parameters are available, the pipeline will run the training step.")
          elif not os.path.exists(args.aug_config_path + "species/" + args.aug_species):
            print ("WARNING no augustus parameters found for " + args.aug_species + " in " + args.aug_config_path + "species/, the pipeline will train augustus")
          elif not os.path.exists(os.path.dirname(args.augustus_prediction) + "/" + args.aug_species + ".trained"):
            if not os.path.exists(os.path.dirname(args.augustus_prediction)):
              os.makedirs(os.path.dirname(args.augustus_prediction))
            Path(os.path.dirname(args.augustus_prediction) + "/" + args.aug_species + ".trained").touch()

        if args.run_augustus_hints:
          args.extrinsic_file_augustus_hints = os.path.abspath(args.extrinsic_file_augustus_hints)
          if not os.path.exists(args.extrinsic_file_augustus_hints):
            print (args.extrinsic_file_augustus_hints +  " not found")  

        if args.run_geneid or args.run_geneid_introns:
          if args.geneid_parameters != None:
            args.geneid_parameters = os.path.abspath(args.geneid_parameters) 
            if not os.path.exists(args.geneid_parameters):
              print (args.geneid_parameters + " not found")
          else:
            print ("ERROR Geneid parameter file has not been specified")  

        if args.evm_path:  
          args.evm_path = os.path.abspath(args.evm_path) + "/"

        if args.run_update:
          if args.update_config == None:
            parser.print_help()
            print ("ERROR No pasa update configuration file found.")
            sys.exit(-1)
          else:
            if not os.path.exists(args.update_config):
              print (args.update_config + " not found")
            else:
              args.update_config = os.path.abspath(args.update_config)
        
        ##Assign wildcards
        samples=[]
        if args.illumina_fastqs == None and args.illumina_dir != None:
          illumina_samples = []
          if args.illumina_PE == True:
            args.illumina_fastqs, illumina_samples = get_wildcards(args.illumina_dir, args.illumina_fastqs, '.1.fastq.gz', args.star_dir, "Aligned.sortedByCoord.out.sam")
          else:
            args.illumina_fastqs, illumina_samples = get_wildcards(args.illumina_dir, args.illumina_fastqs, '.fastq.gz', args.star_dir, "Aligned.sortedByCoord.out.sam")
          for i in illumina_samples:
            samples.append(i)

        if args.cDNA_fastqs == None and args.cDNA_dir != None:
          cDNA_samples = []
          args.cDNA_fastqs, cDNA_samples = get_wildcards(args.cDNA_dir, args.cDNA_fastqs, '.fastq.gz', args.cdna_minimap_dir, ".sorted.sam")
          for i in cDNA_samples:
            samples.append(i)

        if args.dRNA_fastqs == None and args.dRNA_dir != None:
          dRNA_samples = []
          args.dRNA_fastqs, dRNA_samples = get_wildcards(args.dRNA_dir, args.dRNA_fastqs, '.fastq.gz', args.drna_minimap_dir, ".sorted.sam")
          for i in dRNA_samples:
            samples.append(i)

        if args.pb_fastqs == None and args.pb_fastqs == None and args.pb_dir != None:
          pb_samples = []
          args.pb_fastqs, pb_samples = get_wildcards(args.pb_dir, args.pb_fastqs, '.fastq.gz', args.pb_minimap_dir, ".sorted.sam")
          args.pb_fastas, pb_samples = get_wildcards(args.pb_dir, args.pb_fastas, '.fasta', args.pb_minimap_dir, ".sorted.sam")
          for i in pb_samples:
            samples.append(i)

        if args.rna_samples_tsv == None:
          args.rna_samples_tsv = args.rna_outdir + self.rna_samples_tsv
        else:
          args.rna_samples_tsv = os.path.abspath(args.rna_samples_tsv)
        
        if not samples:
          if args.run_geneid_introns or args.run_genemark_ET or args.run_augustus_hints:
            parser.print_help()
            print ("ERROR No RNAseq data given to extract junctions, please consider switching off the options for running gene predictors with evidence or provide junctions and/or RNA evidence")
            sys.exit(-1)
        elif not os.path.exists(args.rna_samples_tsv):
          f = open(args.rna_samples_tsv, "w")
          for line in samples:
            f.write(line + "\n")
          f.close()
        else:
          print ("WARNING: " + args.rna_samples_tsv + " already exists, please remove it if you want to update it")

###      
    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["specFile"] = args.specFile
        self.generalParameters["base_name"] = args.base_name
        self.generalParameters["base_dir"] = args.pipeline_workdir
        self.generalParameters["annotation_version"] = args.annotation_version
        self.generalParameters["project_name"] = args.project_name
        self.generalParameters["run_redmask"] = args.run_redmask
        self.generalParameters["run_pasa"] = args.run_pasa
        self.generalParameters["run_transdecoder"] = args.run_transdecoder
        self.generalParameters["run_miniprot"] = args.run_miniprot
        self.generalParameters["run_augustus"] = args.run_augustus
        self.generalParameters["run_augustus_hints"] = args.run_augustus_hints
        self.generalParameters["run_geneid"] = args.run_geneid
        self.generalParameters["run_geneid_introns"] = args.run_geneid_introns
        self.generalParameters["run_genemark"] = args.run_genemark
        self.generalParameters["run_genemark-ET"] = args.run_genemark_ET
        self.generalParameters["run_EVM"] = args.run_EVM
        self.generalParameters["run_update"] = args.run_update
        self.generalParameters["run_non_coding"] = args.run_non_coding
        self.generalParameters["get_JBrowse"] = args.get_Jbrowse
        self.generalParameters["keep_intermediate"] = args.keep_intermediate
        self.allParameters  ["Parameters"] = self.generalParameters

    def storeallSpecParameters(self,args):
        """Updates rule all cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.allSpecParameters["name"] = "{rule}_{base}_annotation_pipeline"
        self.allSpecParameters["qos"] = self.all_qos
        self.allSpecParameters["time"] = self.all_time
        self.allSpecParameters["queue"] = self.all_queue
        self.allParameters ["all"] = self.allSpecParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["scripts_dir"] = args.scripts_dir        
        self.inputParameters["genome"] = args.genome
        self.inputParameters["genome_lengths"] = args.glen
        self.inputParameters["illumina_dir"] = args.illumina_dir
        self.inputParameters["cDNA_dir"] = args.cDNA_dir
        self.inputParameters["dRNA_dir"] = args.dRNA_dir
        self.inputParameters["PB_dir"] = args.pb_dir
        self.inputParameters["RNA_Samples_tsv"] = args.rna_samples_tsv
        self.inputParameters["transcripts"] = args.transcripts
        self.inputParameters["pasa_config"] = args.pasa_config
        self.inputParameters["update_config"] = args.update_config
        self.inputParameters["proteins"] = args.proteins
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["genome_masked"] = args.genome_masked
        self.outputParameters["Repeat_Annotation_dir"] = args.repeat_annotation_dir
        self.outputParameters["Repeat_gff"] = args.Repeat_gff
        self.outputParameters["RNA_outdir"] = args.rna_outdir
        self.outputParameters["TACO_dir"] = args.TACO_dir
        self.outputParameters["junctions"] = args.junctions
        self.outputParameters["GTF models"] = args.RNAmodels
        self.outputParameters["dir_masked_chunks"] = args.dir_masked_chunks
        self.outputParameters["annotation_basedir"] = args.annotation_basedir
        self.outputParameters["EVM_dir"] = args.EVM_dir
        self.outputParameters["pasa_dir"] = args.pasa_dir
        self.outputParameters["miniprot_cds"] = args.miniprot_cds
        self.outputParameters["miniprot_gene"] = args.miniprot_gene
        self.outputParameters["augustus_prediction"] = args.augustus_prediction
        self.outputParameters["augustus_preEVM"] = args.augustus_preEVM
        self.outputParameters["geneid_prediction"] = args.geneid_prediction
        self.outputParameters["geneid_preEVM"] = args.geneid_preEVM        
        self.outputParameters["genemark_prediction"] = args.genemark_prediction
        self.outputParameters["genemark_preEVM"] = args.genemark_preEVM
        self.outputParameters["geneid_introns_prediction"] = args.geneid_introns_prediction
        self.outputParameters["geneid_introns_preEVM"] = args.geneid_introns_preEVM
        self.outputParameters["genemark_ET_prediction"] = args.genemark_ET_prediction
        self.outputParameters["genemark_ET_preEVM"] = args.genemark_ET_preEVM
        self.outputParameters["augustus_hints_prediction"] = args.augustus_hints_prediction
        self.outputParameters["augustus_hints_preEVM"] = args.augustus_hints_preEVM
        self.outputParameters["evm_out"] = args.evm_out
        self.outputParameters["update_dir"] = args.update_dir
        self.outputParameters["ncRNA_annotation_dir"] = args.ncRNA_annotation_dir
        self.allParameters["Outputs"] = self.outputParameters

    def storeChunksParameters(self,args):
        """Updates chunks parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.chunksParameters["masked_chunks"] = args.masked_chunks
        self.chunksParameters["protein_chunks"] = args.protein_chunks
        self.allParameters["Chunks"] = self.chunksParameters

    def storechunksSpecParameters(self,args):
        """Updates get_chunks_fasta cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.chunksSpecParameters["name"] = "{rule}_{base}_{input.fasta}"
        self.chunksSpecParameters["qos"] = self.chunks_qos
        self.chunksSpecParameters["time"] = self.chunks_time
        self.chunksSpecParameters["queue"] = self.chunks_queue
        self.chunksSpecParameters["mem"] = self.chunks_mem
        self.allParameters ["get_chunks_fasta"] = self.chunksSpecParameters

    def storeRepeatParameters(self,args):
        """Updates repeat annotation parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.repeatParameters["species_repeat_database"] = args.species_rdatabase
        self.repeatParameters["Repeat_library"] = args.repeat_library
        self.repeatParameters["rmaskCores"] = args.rmaskCores
        self.allParameters["Repeat Annotation"] = self.repeatParameters

    def storermaskSpecParameters(self,args):
        """Updates repeat masker cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.rmaskSpecParameters["name"] = "{rule}_{wildcards.masked}"
        self.rmaskSpecParameters["qos"] = self.rmask_qos
        self.rmaskSpecParameters["time"] = self.rmask_time
        self.rmaskSpecParameters["queue"] = self.rmask_queue
        self.rmaskSpecParameters["mem"] = self.rmask_mem
        self.allParameters ["repeat_masker"] = self.rmaskSpecParameters

    def storeRedParameters(self,args):
        """Updates redmask parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.redParameters["minkmer"] = args.red_minkmer
        self.redParameters["wordlen"] = args.red_wordlen
        self.redParameters["add_option"] = args.red_add_opts
        self.allParameters["RedMask"] = self.redParameters

    def storeredmaskSpecParameters(self,args):
        """Updates RedMask cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.redmaskSpecParameters["name"] = "{rule}_{base}"
        self.redmaskSpecParameters["qos"] = self.redmask_qos
        self.redmaskSpecParameters["time"] = self.redmask_time
        self.redmaskSpecParameters["queue"] = self.redmask_queue
        self.redmaskSpecParameters["mem"] = self.redmask_mem
        self.allParameters ["redmask"] = self.redmaskSpecParameters

    def storeBlastParameters(self,args):
        """Updates blast parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.blastParameters["blastdb"] = args.blastdb
        self.blastParameters["evalue"] = args.evalue
        self.blastParameters["blastCores"] = args.blastCores
        self.allParameters["BLAST"] = self.blastParameters

    def storeblastSpecParameters(self,args):
        """Updates BLAST cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.blastSpecParameters["name"] = "{rule}_{base}"
        self.blastSpecParameters["qos"] = self.blast_qos
        self.blastSpecParameters["time"] = self.blast_time
        self.blastSpecParameters["queue"] = self.blast_queue
        self.blastSpecParameters["mem"] = self.blast_mem
        self.allParameters ["filter_prot"] = self.blastSpecParameters

    def storerepgffSpecParameters(self,args):
        """Updates get repeats gff cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.repgffSpecParameters["name"] = "{rule}_{base}"
        self.repgffSpecParameters["qos"] = self.repgff_qos
        self.repgffSpecParameters["time"] = self.repgff_time
        self.repgffSpecParameters["queue"] = self.repgff_queue
        self.repgffSpecParameters["mem"] = self.repgff_mem
        self.allParameters ["get_repeats_gff"] = self.repgffSpecParameters

    def storeIlluminaParameters(self,args):
        """Updates illumina parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.illuminaParameters["star_genome_dir"] = args.star_genome_dir
        self.illuminaParameters["PE"] = args.illumina_PE
        self.illuminaParameters["star_dir"] = args.star_dir
        self.illuminaParameters["starCores"] = args.starCores
        self.illuminaParameters["star_index_additional_opts"] = args.indexstar_additional_options
        self.illuminaParameters["star_additional_opts"] = args.star_additional_options
        self.illuminaParameters["stringtie_opts"] = args.stringtie_illumina_opts
        self.illuminaParameters["TACO_opts"] = args.TACO_illumina_opts
        self.allParameters ["Illumina RNA"] = self.illuminaParameters

    def storeTrimgaloreParameters(self,args):
        """Updates the Trim_Galore parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreParameters["options"] = args.trim_galore_opts
        self.trimgaloreParameters["Trim_Illumina_cores"] = args.Trim_Illumina_cores
        self.allParameters ["Trim_Galore"] = self.trimgaloreParameters

    def storestarindexSpecParameters(self,args):
        """Updates star index cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.starindexSpecParameters["name"] = "{rule}_{base}"
        self.starindexSpecParameters["qos"] = self.star_index_qos
        self.starindexSpecParameters["time"] = self.star_index_time
        self.starindexSpecParameters["queue"] = self.star_index_queue
        self.starindexSpecParameters["mem"] = self.star_index_mem   
        self.allParameters ["star_index"] = self.starindexSpecParameters

    def storetrimgaloreSpecParameters(self,args):
        """Updates trimgalore cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.file}"
        self.trimgaloreSpecParameters["qos"] = self.trimgalore_qos
        self.trimgaloreSpecParameters["time"] = self.trimgalore_time
        self.trimgaloreSpecParameters["queue"] = self.trimgalore_queue
        self.trimgaloreSpecParameters["mem"] = self.trimgalore_mem
        self.allParameters ["trim_galore"] = self.trimgaloreSpecParameters

    def storestarSpecParameters(self,args):
        """Updates star cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.starSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.file}"
        self.starSpecParameters["qos"] = self.star_qos
        self.starSpecParameters["time"] = self.star_time
        self.starSpecParameters["queue"] = self.star_queue
        self.starSpecParameters["mem"] = self.star_mem   
        self.allParameters ["star"] = self.starSpecParameters

    def storeCdnaParameters(self,args):
        """Updates cDNA parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.cdnaParameters["minimap_dir"] = args.cdna_minimap_dir
        self.cdnaParameters["minimap2_opts"] = args.minimap2_cDNA_opts
        self.cdnaParameters["stringtie_opts"] = args.stringtie_cDNA_opts
        self.cdnaParameters["TACO_opts"] = args.TACO_cDNA_opts
        self.allParameters ["cDNA RNA"] = self.cdnaParameters

    def storeDrnaParameters(self,args):
        """Updates dRNA parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.drnaParameters["minimap_dir"] = args.drna_minimap_dir
        self.drnaParameters["minimap2_opts"] = args.minimap2_dRNA_opts
        self.drnaParameters["stringtie_opts"] = args.stringtie_dRNA_opts
        self.drnaParameters["TACO_opts"] = args.TACO_dRNA_opts
        self.allParameters ["dRNA RNA"] = self.drnaParameters

    def storePbParameters(self,args):
        """Updates pacbio isoseq parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pbParameters["minimap_dir"] = args.pb_minimap_dir
        self.pbParameters["minimap2_opts"] = args.minimap2_pb_opts
        self.pbParameters["stringtie_opts"] = args.stringtie_pb_opts
        self.pbParameters["TACO_opts"] = args.TACO_pb_opts
        self.allParameters ["Isoseq"] = self.pbParameters

    def storeminimapSpecParameters(self,args):
        """Updates minimap cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.minimapSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.file}"
        self.minimapSpecParameters["qos"] = self.minimap_qos
        self.minimapSpecParameters["time"] = self.minimap_time
        self.minimapSpecParameters["queue"] = self.minimap_queue
        self.minimapSpecParameters["mem"] = self.minimap_mem
        self.allParameters ["minimap2"] = self.minimapSpecParameters

    def storeModelRNAParameters(self,args):
        """Updates the parameters to model RNA to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.modelRNAParameters["minimapCores"] = args.minimapCores
        self.modelRNAParameters["TACO global options"] = args.TACO_global_opts
        self.modelRNAParameters["TACOCores"] = args.tacoCores
        self.modelRNAParameters["EspressoCores"] = args.espressoCores
        self.modelRNAParameters["ESPRESSO path"] = args.espresso_path
        self.modelRNAParameters["ESPRESSO outdir"] = args.espresso_outdir
        self.allParameters ["Model RNA"] = self.modelRNAParameters

    def storestringtieSpecParameters(self,args):
        """Updates stringtie cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.stringtieSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.file}"
        self.stringtieSpecParameters["qos"] = self.stringtie_qos
        self.stringtieSpecParameters["time"] = self.stringtie_time
        self.stringtieSpecParameters["queue"] = self.stringtie_queue
        self.stringtieSpecParameters["mem"] = self.stringtie_mem
        self.allParameters ["stringtie"] = self.stringtieSpecParameters

    def storetacoSpecParameters(self,args):
        """Updates TACO cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.tacoSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.dir}"
        self.tacoSpecParameters["qos"] = self.taco_qos
        self.tacoSpecParameters["time"] = self.taco_time
        self.tacoSpecParameters["queue"] = self.taco_queue
        self.tacoSpecParameters["mem"] = self.taco_mem
        self.allParameters ["join_models"] = self.tacoSpecParameters

    def storebam2samSpecParameters(self,args):
        """Updates BAM2SAM cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.bam2samSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.dir}"
        self.bam2samSpecParameters["qos"] = self.bam2sam_qos
        self.bam2samSpecParameters["time"] = self.bam2sam_time
        self.bam2samSpecParameters["queue"] = self.bam2sam_queue
        self.bam2samSpecParameters["mem"] = self.bam2sam_mem
        self.allParameters ["bam2sam"] = self.bam2samSpecParameters

    def storeespressoSpecParameters(self,args):
        """Updates ESPRESSO cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.espressoSpecParameters["name"] = "{rule}_{base}"
        self.espressoSpecParameters["qos"] = self.espresso_qos
        self.espressoSpecParameters["time"] = self.espresso_time
        self.espressoSpecParameters["queue"] = self.espresso_queue
        self.espressoSpecParameters["mem"] = self.espresso_mem
        self.allParameters ["ESPRESSO"] = self.espressoSpecParameters

    def storecodingHintsSpecParameters(self,args):
        """Updates get_coding_hints cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.codingHintsSpecParameters["name"] = "{rule}_{base}"
        self.codingHintsSpecParameters["qos"] = self.codingHints_qos
        self.codingHintsSpecParameters["time"] = self.codingHints_time
        self.codingHintsSpecParameters["queue"] = self.codingHints_queue
        self.codingHintsSpecParameters["mem"] = self.codingHints_mem
        self.allParameters ["get_coding_junctions"] = self.codingHintsSpecParameters

    def storePasaParameters(self,args):
        """Updates pasa parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pasaParameters["create_database"] = args.create_database
        self.pasaParameters["pasadb"] = args.pasadb
        self.pasaParameters["pasa_weights"] = args.pasa_weights
        self.pasaParameters["aligners"] = args.aligners
        self.pasaParameters["add_option"] = args.add_option
        self.pasaParameters["update_rounds"] = args.update_rounds
        self.allParameters["pasa"] = self.pasaParameters

    def storePasaSpecParameters(self,args):
        """Updates PASA cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pasaSpecParameters["name"] = "{rule}_{base}"
        self.pasaSpecParameters["qos"] = self.PASA_qos
        self.pasaSpecParameters["time"] = self.PASA_time
        self.pasaSpecParameters["queue"] = self.PASA_queue
        self.pasaSpecParameters["mem"] = self.PASA_mem
        self.allParameters ["pasa"] = self.pasaSpecParameters

    def storeTransdecoderParameters(self,args):
        """Updates transdecoder parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.transdecoderParameters["transdecoder_weights"] = args.transdecoder_weights
        self.allParameters["transdecoder"] = self.transdecoderParameters

    def storeTransdecoderSpecParameters(self,args):
        """Updates transdecoder cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.transdecoderSpecParameters["name"] = "{rule}_{base}"
        self.transdecoderSpecParameters["qos"] = self.transdecoder_qos
        self.transdecoderSpecParameters["time"] = self.transdecoder_time
        self.transdecoderSpecParameters["queue"] = self.transdecoder_queue
        self.transdecoderSpecParameters["mem"] = self.transdecoder_mem
        self.allParameters ["Transdecoder"] = self.transdecoderSpecParameters

    def storeMiniprotParameters(self,args):
        """Updates miniprot parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.miniprotParameters["miniprot_path"] = args.miniprot_path
        self.miniprotParameters["miniprot_cores"] = args.miniprot_cores
        self.miniprotParameters["miniprot_weights"] = args.miniprot_weights
        self.miniprotParameters["additional_miniprot_options"] = args.additional_miniprot_options      
        self.allParameters ["miniprot"] = self.miniprotParameters

    def storeMiniprotSpecParameters(self,args):
        """Updates Miniprot cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.miniprotSpecParameters["name"] = "{rule}_{base}"
        self.miniprotSpecParameters["qos"] = self.miniprot_qos
        self.miniprotSpecParameters["time"] = self.miniprot_time
        self.miniprotSpecParameters["queue"] = self.miniprot_queue
        self.miniprotSpecParameters["mem"] = self.miniprot_mem
        self.allParameters ["miniprot"] = self.miniprotSpecParameters

    def storegetcandSpecParameters(self,args):
        """Updates get training candidates cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.getcandSpecParameters["name"] = "{rule}_{base}"
        self.getcandSpecParameters["qos"] = self.getcand_qos
        self.getcandSpecParameters["time"] = self.getcand_time
        self.getcandSpecParameters["queue"] = self.getcand_queue
        self.getcandSpecParameters["mem"] = self.getcand_mem
        self.allParameters ["get_training_candidates"] = self.getcandSpecParameters
  
    def storetrainaugSpecParameters(self,args):
        """Updates training augustus cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trainaugSpecParameters["name"] = "{rule}_{base}"
        self.trainaugSpecParameters["qos"] = self.trainaug_qos
        self.trainaugSpecParameters["time"] = self.trainaug_time
        self.trainaugSpecParameters["queue"] = self.trainaug_queue
        self.trainaugSpecParameters["mem"] = self.trainaug_mem
        self.allParameters ["train_augustus"] = self.trainaugSpecParameters

    def storeAugustusParameters(self,args):
        """Updates augustus parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusParameters["species"] = args.aug_species
        self.augustusParameters["optimize_threads"] = args.aug_optimize_threads
        self.augustusParameters["config_path"] = args.aug_config_path
        self.augustusParameters["alternatives_from_sampling"] = args.aug_alternatives_from_sampling
        self.augustusParameters["alternatives_from_evidence"] = args.aug_alternatives_from_evidence
        self.augustusParameters["uniqueGeneId"] = args.aug_uniqueGeneId
        self.augustusParameters["gff3"] = args.aug_gff3
        self.augustusParameters["sample"] = args.aug_sample
        self.augustusParameters["noInFrameStop"] = args.aug_noInFrameStop
        self.augustusParameters["maxtracks"] = args.aug_maxtracks
        self.augustusParameters["singlestrand"] = args.aug_singlestrand
        self.augustusParameters["strand"] = args.aug_strand
        self.augustusParameters["min_intron_len"] = args.aug_min_intron_len
        self.augustusParameters["weights"] = args.augustus_weights
        self.augustusParameters["additional_augustus_options"] = args.additional_augustus_options
        self.allParameters ["augustus"] = self.augustusParameters

    def storeAugustusSpecParameters(self,args):
        """Updates augustus cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusSpecParameters["name"] = "{rule}_{base}"
        self.augustusSpecParameters["qos"] = self.augustus_qos
        self.augustusSpecParameters["time"] = self.augustus_time
        self.augustusSpecParameters["queue"] = self.augustus_queue
        self.augustusSpecParameters["mem"] = self.augustus_mem
        self.allParameters ["augustus"] = self.augustusSpecParameters

    def storeAugustusArraySpecParameters(self,args):
        """Updates augustus cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusArraySpecParameters["name"] = "{rule}_{base}"
        self.augustusArraySpecParameters["qos"] = self.augustus_qos
        self.augustusArraySpecParameters["time"] = self.augustus_time
        self.augustusArraySpecParameters["queue"] = self.augustus_queue
        self.augustusArraySpecParameters["mem"] = self.augustus_mem
        self.augustusArraySpecParameters["array"] = "1-{masked_chunks}%15"
        self.allParameters ["augustus_jobarray"] = self.augustusArraySpecParameters

    def storeAugustusHintsParameters(self,args):
        """Updates augustus with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusHintsParameters["extrinsic_file_augustus_hints"] = args.extrinsic_file_augustus_hints
        self.augustusHintsParameters["augustus_hints_weights"] = args.augustus_hints_weights
        self.augustusHintsParameters["additional_augustus_hints_options"] = args.additional_augustus_hints_options       
        self.allParameters ["augustus_hints"] = self.augustusHintsParameters

    def storeAugustusHintsSpecParameters(self,args):
        """Updates augustus with hintscluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusHintsSpecParameters["name"] = "{rule}_{base}"
        self.augustusHintsSpecParameters["qos"] = self.augustus_hints_qos
        self.augustusHintsSpecParameters["time"] = self.augustus_hints_time
        self.augustusHintsSpecParameters["queue"] = self.augustus_hints_queue
        self.augustusHintsSpecParameters["mem"] = self.augustus_hints_mem
        self.allParameters ["augustus_hints"] = self.augustusHintsSpecParameters

    def storeAugustusHintsArraySpecParameters(self,args):
        """Updates augustus with hints cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusHintsArraySpecParameters["name"] = "{rule}_{base}"
        self.augustusHintsArraySpecParameters["qos"] = self.augustus_hints_qos
        self.augustusHintsArraySpecParameters["time"] = self.augustus_hints_time
        self.augustusHintsArraySpecParameters["queue"] = self.augustus_hints_queue
        self.augustusHintsArraySpecParameters["mem"] = self.augustus_hints_mem
        self.augustusHintsArraySpecParameters["array"] = "1-{masked_chunks}%15"
        self.allParameters ["augustus_hints_jobarray"] = self.augustusHintsArraySpecParameters

    def storeGeneidParameters(self,args):
        """Updates geneid parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidParameters["parameters"] = args.geneid_parameters
        self.geneidParameters["path"] = args.geneid_path
        self.geneidParameters["weights"] = args.geneid_weights
        self.geneidParameters["options"] = args.geneid_options       
        self.allParameters ["geneid"] = self.geneidParameters

    def storeGeneidSpecParameters(self,args):
        """Updates geneid cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidSpecParameters["name"] = "{rule}_{base}"
        self.geneidSpecParameters["qos"] = self.geneid_qos
        self.geneidSpecParameters["time"] = self.geneid_time
        self.geneidSpecParameters["queue"] = self.geneid_queue
        self.geneidSpecParameters["mem"] = self.geneid_mem
        self.allParameters ["geneid"] = self.geneidSpecParameters

    def storeGeneidIntronsParameters(self,args):
        """Updates geneid with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidIntronsParameters["weights"] = args.geneid_introns_weights
        self.geneidIntronsParameters["options"] = args.geneid_introns_options       
        self.allParameters ["geneid_introns"] = self.geneidIntronsParameters

    def storeGeneidIntronsSpecParameters(self,args):
        """Updates geneid introns cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidIntronsSpecParameters["name"] = "{rule}_{base}"
        self.geneidIntronsSpecParameters["qos"] = self.geneid_introns_qos
        self.geneidIntronsSpecParameters["time"] = self.geneid_introns_time
        self.geneidIntronsSpecParameters["queue"] = self.geneid_introns_queue
        self.geneidIntronsSpecParameters["mem"] = self.geneid_introns_mem
        self.allParameters ["geneid_introns"] = self.geneidIntronsSpecParameters

    def storeGenemarkParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkParameters["max_contig"] = args.gmk_max_contig
        self.genemarkParameters["min_contig"] = args.gmk_min_contig
        self.genemarkParameters["max_gap"] = args.gmk_max_gap
        self.genemarkParameters["cores"] = args.gmk_cores
        self.genemarkParameters["weights"] = args.genemark_weights
        self.genemarkParameters["additional_genemark_options"] = args.additional_genemark_options   
        self.allParameters ["genemark"] = self.genemarkParameters

    def storeGenemarkSpecParameters(self,args):
        """Updates genemark cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkSpecParameters["name"] = "{rule}_{base}"
        self.genemarkSpecParameters["qos"] = self.genemark_qos
        self.genemarkSpecParameters["time"] = self.genemark_time
        self.genemarkSpecParameters["queue"] = self.genemark_queue
        self.genemarkSpecParameters["mem"] = self.genemark_mem
        self.allParameters ["genemark"] = self.genemarkSpecParameters

    def storeGenemarkETParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkETParameters["genemark_ET_weights"] = args.genemark_ET_weights
        self.genemarkETParameters["additional_genemark_ET_options"] = args.additional_genemark_ET_options 
        self.allParameters ["genemark_ET"] = self.genemarkETParameters

    def storeGenemarkETSpecParameters(self,args):
        """Updates genemark-ET cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkETSpecParameters["name"] = "{rule}_{base}"
        self.genemarkETSpecParameters["qos"] = self.genemark_et_qos
        self.genemarkETSpecParameters["time"] = self.genemark_et_time
        self.genemarkETSpecParameters["queue"] = self.genemark_et_queue
        self.genemarkETSpecParameters["mem"] = self.genemark_et_mem
        self.allParameters ["genemark_ET"] = self.genemarkETSpecParameters

    def storemergegffSpecParameters(self,args):
        """Updates merge_gffs cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.mergegffSpecParameters["name"] = "{rule}_{base}_{wildcards.name}"
        self.mergegffSpecParameters["qos"] = self.mergegff_qos
        self.mergegffSpecParameters["time"] = self.mergegff_time
        self.mergegffSpecParameters["queue"] = self.mergegff_queue
        self.mergegffSpecParameters["mem"] = self.mergegff_mem
        self.allParameters ["merge_gffs"] = self.mergegffSpecParameters

    def storepred4evmSpecParameters(self,args):
        """Updates predictions4EVM cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pred4evmSpecParameters["name"] = "{rule}_{base}_{wildcards.name}"
        self.pred4evmSpecParameters["qos"] = self.pred4evm_qos
        self.pred4evmSpecParameters["time"] = self.pred4evm_time
        self.pred4evmSpecParameters["queue"] = self.pred4evm_queue
        self.pred4evmSpecParameters["mem"] = self.pred4evm_mem
        self.allParameters ["predictions4EVM"] = self.pred4evmSpecParameters

    def storeEvmParameters(self,args):
        """Updates evm parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.evmParameters["evm_path"] = args.evm_path
        self.evmParameters["segmentsize"] = args.evm_segmentsize
        self.evmParameters["overlapsize"] = args.evm_overlapsize
        self.evmParameters["cores"] = args.evm_cores
        self.evmParameters["additional_evm_options"] = args.additional_evm_options
        self.allParameters["EVM"] = self.evmParameters

    def storeprepevmSpecParameters(self,args):
        """Updates prepare EVM cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.prepevmSpecParameters["name"] = "{rule}_{base}"
        self.prepevmSpecParameters["qos"] = self.prepevm_qos
        self.prepevmSpecParameters["time"] = self.prepevm_time
        self.prepevmSpecParameters["queue"] = self.prepevm_queue
        self.prepevmSpecParameters["mem"] = self.prepevm_mem
        self.allParameters ["prepare_evm"] = self.prepevmSpecParameters

    def storeevmSpecParameters(self,args):
        """Updates EVM cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.evmSpecParameters["name"] = "{rule}_{base}_{wildcards.w}"
        self.evmSpecParameters["qos"] = self.evm_qos
        self.evmSpecParameters["time"] = self.evm_time
        self.evmSpecParameters["queue"] = self.evm_queue
        self.evmSpecParameters["mem"] = self.evm_mem
        self.allParameters ["EVM2"] = self.evmSpecParameters

    def storeselectevmSpecParameters(self,args):
        """Updates select EVM cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.selectevmSpecParameters["name"] = "{rule}_{base}"
        self.selectevmSpecParameters["qos"] = self.selectevm_qos
        self.selectevmSpecParameters["time"] = self.selectevm_time
        self.selectevmSpecParameters["queue"] = self.selectevm_queue
        self.selectevmSpecParameters["mem"] = self.selectevm_mem
        self.allParameters ["select_EVM"] = self.selectevmSpecParameters

    def storepasaupdateSpecParameters(self,args):
        """Updates PASA update cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pasaupdateSpecParameters["name"] = "{rule}_{wildcards.round}_{base}"
        self.pasaupdateSpecParameters["qos"] = self.pasaupdate_qos
        self.pasaupdateSpecParameters["time"] = self.pasaupdate_time
        self.pasaupdateSpecParameters["queue"] = self.pasaupdate_queue
        self.pasaupdateSpecParameters["mem"] = self.pasaupdate_mem
        self.allParameters ["PASA_update"] = self.pasaupdateSpecParameters

    def storeprocessupdateSpecParameters(self,args):
        """Updates process update cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.processupdateSpecParameters["name"] = "{rule}_{base}_{params.project}{params.version}"
        self.processupdateSpecParameters["qos"] = self.processupdate_qos
        self.processupdateSpecParameters["time"] = self.processupdate_time
        self.processupdateSpecParameters["queue"] = self.processupdate_queue
        self.processupdateSpecParameters["mem"] = self.processupdate_mem
        self.allParameters ["process_update"] = self.processupdateSpecParameters

    def storencRNAannotationParameters(self,args):
        """Updates ncRNA Annotation parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ncRNAannotationParameters["ncRNA_version"] = args.ncRNA_version
        self.ncRNAannotationParameters["cmsearch_CPUs"] = args.cmsearch_CPUs
        self.ncRNAannotationParameters["Rfam"] = args.Rfam
        self.allParameters["ncRNA_annotation"] = self.ncRNAannotationParameters

    def storecmsearchSpecParameters(self,args):
        """Updates cmsearch cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.cmsearchSpecParameters["name"] = "{rule}_{base}"
        self.cmsearchSpecParameters["qos"] = self.cmsearch_qos
        self.cmsearchSpecParameters["time"] = self.cmsearch_time
        self.cmsearchSpecParameters["queue"] = self.cmsearch_queue
        self.cmsearchSpecParameters["mem"] = self.cmsearch_mem
        self.allParameters ["cmsearch"] = self.cmsearchSpecParameters

    def storetRNAscanSpecParameters(self,args):
        """Updates tRNAscan cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.tRNAscanSpecParameters["name"] = "{rule}_{base}"
        self.tRNAscanSpecParameters["qos"] = self.tRNAscan_qos
        self.tRNAscanSpecParameters["time"] = self.tRNAscan_time
        self.tRNAscanSpecParameters["queue"] = self.tRNAscan_queue
        self.tRNAscanSpecParameters["mem"] = self.tRNAscan_mem
        self.allParameters ["tRNAscan"] = self.tRNAscanSpecParameters

    def storelncRNASpecParameters(self,args):
        """Updates lncRNA cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.lncRNASpecParameters["name"] = "{rule}_{base}"
        self.lncRNASpecParameters["qos"] = self.lncRNA_qos
        self.lncRNASpecParameters["time"] = self.lncRNA_time
        self.lncRNASpecParameters["queue"] = self.lncRNA_queue
        self.lncRNASpecParameters["mem"] = self.lncRNA_mem
        self.allParameters ["lncRNAannotation"] = self.lncRNASpecParameters

    def storeBlastProtSpecParameters(self,args):
        """Updates BLAST prot cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.BlastProtSpecParameters["name"] = "{rule}_{base}_{wildcards.dirs}{wildcards.i}"
        self.BlastProtSpecParameters["qos"] = self.blast_prot_qos
        self.BlastProtSpecParameters["time"] = self.blast_prot_time
        self.BlastProtSpecParameters["queue"] = self.blast_prot_queue
        self.BlastProtSpecParameters["mem"] = self.blast_prot_mem
        self.allParameters ["Blast_prot"] = self.BlastProtSpecParameters

    def storencRNASpecParameters(self,args):
        """Updates ncRNA cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ncRNASpecParameters["name"] = "{rule}_{base}"
        self.ncRNASpecParameters["qos"] = self.ncRNA_qos
        self.ncRNASpecParameters["time"] = self.ncRNA_time
        self.ncRNASpecParameters["queue"] = self.ncRNA_queue
        self.ncRNASpecParameters["mem"] = self.ncRNA_mem
        self.allParameters ["ncAnnotation"] = self.ncRNASpecParameters

    def storegetGCSpecParameters(self,args):
        """Updates getGC for jbrowse track cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.getGCSpecParameters["name"] = "{rule}_{base}"
        self.getGCSpecParameters["qos"] = self.getGC_qos
        self.getGCSpecParameters["time"] = self.getGC_time
        self.getGCSpecParameters["queue"] = self.getGC_queue
        self.getGCSpecParameters["mem"] = self.getGC_mem
        self.allParameters ["get_GCcontent"] = self.getGCSpecParameters

    def storegetseqSpecParameters(self,args):
        """Updates browse seq for jbrowse track cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.getseqSpecParameters["name"] = "{rule}_{base}"
        self.getseqSpecParameters["qos"] = self.getseq_qos
        self.getseqSpecParameters["time"] = self.getseq_time
        self.getseqSpecParameters["queue"] = self.getseq_queue
        self.getseqSpecParameters["mem"] = self.getseq_mem
        self.allParameters ["browse_seq"] = self.getseqSpecParameters

    def storegettracksSpecParameters(self,args):
        """Updates browse tracks for jbrowse track cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.gettracksSpecParameters["name"] = "{rule}_{base}"
        self.gettracksSpecParameters["qos"] = self.gettracks_qos
        self.gettracksSpecParameters["time"] = self.gettracks_time
        self.gettracksSpecParameters["queue"] = self.gettracks_queue
        self.gettracksSpecParameters["mem"] = self.gettracks_mem
        self.allParameters ["browse_tracks"] = self.gettracksSpecParameters

    def storegetbwSpecParameters(self,args):
        """Updates browse bw for jbrowse track cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.getbwSpecParameters["name"] = "{rule}_{base}"
        self.getbwSpecParameters["qos"] = self.getbw_qos
        self.getbwSpecParameters["time"] = self.getbw_time
        self.getbwSpecParameters["queue"] = self.getbw_queue
        self.getbwSpecParameters["mem"] = self.getbw_mem
        self.allParameters ["get_bw"] = self.getbwSpecParameters

    def storegettarSpecParameters(self,args):
        """Updates get tar for jbrowse track cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.gettarSpecParameters["name"] = "{rule}_{base}"
        self.gettarSpecParameters["qos"] = self.gettar_qos
        self.gettarSpecParameters["time"] = self.gettar_time
        self.gettarSpecParameters["queue"] = self.gettar_queue
        self.gettarSpecParameters["mem"] = self.gettar_mem
        self.allParameters ["get_tar"] = self.gettarSpecParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["illumina_fastqs"] = args.illumina_fastqs
        self.wildcardParameters["cDNA_fastqs"] = args.cDNA_fastqs
        self.wildcardParameters["dRNA_fastqs"] = args.dRNA_fastqs
        self.wildcardParameters["isoseq_fastqs"] = args.pb_fastqs
        self.wildcardParameters["isoseq_fastas"] = args.pb_fastas
        self.allParameters ["Wildcards"] = self.wildcardParameters

#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()
specManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the repeat annotation pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeChunksParameters(args)
configManager.storeRepeatParameters(args)
configManager.storeRedParameters(args)
configManager.storeBlastParameters(args)
configManager.storeTrimgaloreParameters(args)
configManager.storeIlluminaParameters(args)
configManager.storeCdnaParameters(args)
configManager.storeDrnaParameters(args)
configManager.storePbParameters(args)
configManager.storeModelRNAParameters(args)
configManager.storePasaParameters(args)
configManager.storeTransdecoderParameters(args)
configManager.storeMiniprotParameters(args)
configManager.storeAugustusParameters(args)
configManager.storeAugustusHintsParameters(args)
configManager.storeGeneidParameters(args)
configManager.storeGeneidIntronsParameters(args)
configManager.storeGenemarkParameters(args)
configManager.storeGenemarkETParameters(args)
configManager.storeEvmParameters(args)
configManager.storencRNAannotationParameters(args)
configManager.storeWildcardParameters(args)

specManager.storeallSpecParameters(args)
specManager.storechunksSpecParameters(args)
specManager.storermaskSpecParameters(args)
if args.run_redmask:
  specManager.storeredmaskSpecParameters(args)
  specManager.storeblastSpecParameters(args)
specManager.storerepgffSpecParameters(args)

if args.illumina_dir != None:
  specManager.storetrimgaloreSpecParameters(args)
  specManager.storestarindexSpecParameters(args)
  specManager.storestarSpecParameters(args)

if args.cDNA_dir != None or args.dRNA_dir != None or args.pb_dir != None:
  specManager.storeminimapSpecParameters(args)

if args.illumina_dir != None or args.cDNA_dir != None or args.dRNA_dir != None or args.pb_dir != None:
  specManager.storetacoSpecParameters(args)
  specManager.storebam2samSpecParameters(args)
  specManager.storeespressoSpecParameters(args)
  specManager.storestringtieSpecParameters(args)

specManager.storecodingHintsSpecParameters(args)

if args.run_pasa:
  specManager.storePasaSpecParameters(args)

if args.run_transdecoder:
  specManager.storeTransdecoderSpecParameters(args)

if args.run_miniprot:
  specManager.storeMiniprotSpecParameters(args)

specManager.storegetcandSpecParameters(args)

if args.run_augustus:
  if args.masked_chunks > 1:   
    specManager.storeAugustusArraySpecParameters(args)
    
  else:
    specManager.storeAugustusSpecParameters(args)
  
if args.run_augustus_hints:
  if args.masked_chunks > 1:   
    specManager.storeAugustusHintsArraySpecParameters(args)
  else:
    specManager.storeAugustusHintsSpecParameters(args)

if args.run_augustus or args.run_augustus_hints:
  specManager.storetrainaugSpecParameters(args)
  specManager.storepred4evmSpecParameters(args)
  if args.masked_chunks > 1: 
    specManager.storemergegffSpecParameters(args)

if args.run_geneid:
  specManager.storeGeneidSpecParameters(args)

if args.run_geneid_introns:
  specManager.storeGeneidIntronsSpecParameters(args)

if args.run_genemark:
  specManager.storeGenemarkSpecParameters(args)

if args.run_genemark_ET:
  specManager.storeGenemarkETSpecParameters(args)

if args.run_EVM:  
  specManager.storeprepevmSpecParameters(args)
  specManager.storeevmSpecParameters(args)
  specManager.storeselectevmSpecParameters(args)

if args.run_update:  
  specManager.storepasaupdateSpecParameters(args)
  specManager.storeprocessupdateSpecParameters(args)

if args.run_non_coding:  
  specManager.storecmsearchSpecParameters(args)
  specManager.storetRNAscanSpecParameters(args)
  specManager.storelncRNASpecParameters(args)  
  specManager.storeBlastProtSpecParameters(args)
  specManager.storencRNASpecParameters(args)

if args.get_Jbrowse:  
  specManager.storegetGCSpecParameters(args)
  specManager.storegetseqSpecParameters(args)
  specManager.storegettracksSpecParameters(args)
  specManager.storegetbwSpecParameters(args)
  specManager.storegettarSpecParameters(args)

#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)
with open(args.specFile, 'w') as of:
    json.dump(specManager.allParameters, of, indent=2)
