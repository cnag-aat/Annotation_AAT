# Annotation_AAT
Structural genome annotation pipeline used in the CNAG AA Team for eukaryote genomes. 

It is based on Evidence-Modeler + PASA approach. It can annotate a genome using transcript evidence (optional but highly recommended), protein evidence (compulsory) and ab initio gene predictions. 

It needs a config file and a spec file (json file with instructions on which resources should slurm use for each of the jobs). Both files are created by the script "create_config_annotation.py" that is located in the bin directory. To check all the options accepted by the script, do:

```
bin/create_config_annotation.py -h
```

Besides the help, if you are just planning to run it with default parameters, there is also a way to make it  request all the compulsory options:

```
bin/create_config_file.annotation.py
ERROR Sorry! No masked genome fasta file defined
ERROR Sorry! No genome fasta file defined
Sorry! No junctions gff file given.
ERROR Geneid parameter file has not been specified
ERROR Augustus trained species has not been specified
ERROR Protein evidence file has not been specified
Sorry! No pasadb name given
ERROR Transcript evidence file has not been specified
ERROR No pasa configuration file found.
WARNING: Remember to set the create-database parameter to true the first time you run PASA!
WARNING:no transcripts gtf given, PASA will be run only with the given fasta file
No annotCompare.config file given!
Sorry no project name and version have been given
```

The final command would look like this: 

```
bin/create_config_file.annotation.py --genome evidence/your_assembly.fa --genome-masked evidence/your_masked_assembly.fa --proteins evidence/uniprot_sprot.fasta --geneid-param evidence/pipeline_inputs/M.galloprovincialis.geneid.optimized.U12.param --species human --pasadb human.sqlite --pasa-conf evidence/pipeline_inputs/alignAssembly.v2.5.2.pasa2.config --update-conf evidence/pipeline_inputs/annotCompare_v2.5.2.pasa2.config --transcripts evidence/pipeline_inputs/transcript_evidence.fa
```

* Tricks to obtain the different input options: *

If you do not know which **geneid parameters** to use, you can download the most appropriate ones from: https://genome.crg.es/software/geneid/index.html#parameters. Also, for Augustus species (--species option), check the available species  already trained for augustus or train your own).

Modify the annotCompare_v2.5.2.pasa2.config and alignAssembly.v2.5.2.pasa2.config config files with the name and path to tha pasa database that the pipeline will create. Examples of this files can be found under the "config_examples/" directory.


Once the 2 config files are produced, the pipeline can be launched using snakemake (tested and developed for Snakemake v6.3.0) like this:

``snakemake --notemp -j 999 --snakefile annotation_AAT.smk --configfile annotation.config --is --cluster-conf annotation.spec --use-conda --use-envmodules``

If you are using an HPC cluster, please check how should you run snakemake to launch the jobs to the cluster. 

Most of the tools used will be installed via conda using the environments of the "envs" directory after providing the "--use-conda" option to snakemake. However, a few tools cannot be installed via conda and will have to be available in your PATH, or as a module in the cluster. Those tools are:

- GeneMark 
- GeneID

# Description of all options
```
bin/create_config_file.annotation.py -h
usage: create_configuration_file [-h] [--configFile configFile]
                                 [--specFile specFile] [--basename base_name]
                                 [--no-augustus] [--no-augustus-hints]
                                 [--no-geneid] [--no-geneid-introns]
                                 [--no-genemark] [--no-genemark-ET]
                                 [--no-miniprot] [--no-pasa]
                                 [--no-transdecoder] [--no-EVM] [--no-update]
                                 [--no-noncoding] [--no-jbrowse]
                                 [--scripts-dir SCRIPTS_DIR] [--genome genome]
                                 [--glen glen] [--genome-masked genome_masked]
                                 [--junctions junctions]
                                 [--geneid-parameters geneid_parameters]
                                 [--proteins proteins]
                                 [--transcripts transcripts]
                                 [--pasa-config pasa_config]
                                 [--trans-gtf trans_gtf]
                                 [--update-config update_config]
                                 [--RM-gff RM_gff]
                                 [--annotation-step ANNOTATION_STEP]
                                 [--annotation-version ANNOTATION_VERSION]
                                 [--pipeline-workdir PIPELINE_WORKDIR]
                                 [--output-dir OUTPUT_DIR] [--EVM-dir EVM_DIR]
                                 [--dir-masked-chunks DIR_MASKED_CHUNKS]
                                 [--dir-genome-chunks DIR_GENOME_CHUNKS]
                                 [--augustus-prediction AUGUSTUS_PREDICTION]
                                 [--augustus-preEVM AUGUSTUS_PREEVM]
                                 [--augustus-hints-prediction AUGUSTUS_HINTS_PREDICTION]
                                 [--augustus-hints-preEVM AUGUSTUS_HINTS_PREEVM]
                                 [--geneid-prediction GENEID_PREDICTION]
                                 [--geneid-preEVM GENEID_PREEVM]
                                 [--geneid-introns-prediction GENEID_INTRONS_PREDICTION]
                                 [--geneid-introns-preEVM GENEID_INTRONS_PREEVM]
                                 [--genemark-prediction GENEMARK_PREDICTION]
                                 [--genemark-preEVM GENEMARK_PREEVM]
                                 [--genemark-ET-prediction GENEMARK_ET_PREDICTION]
                                 [--genemark-ET-preEVM GENEMARK_ET_PREEVM]
                                 [--miniprot-cds MINIPROT_CDS]
                                 [--miniprot-gene MINIPROT_GENE]
                                 [--pasa-dir PASA_DIR] [--evm-out EVM_OUT]
                                 [--update-dir UPDATE_DIR]
                                 [--ncRNA-dir NCRNA_ANNOTATION_DIR]
                                 [--out-cmsearch OUT_CMSEARCH]
                                 [--out-tRNAscan OUT_TRNASCAN]
                                 [--masked-chunks MASKED_CHUNKS]
                                 [--genome-chunks GENOME_CHUNKS]
                                 [--protein-chunks PROTEIN_CHUNKS]
                                 [--species aug_species]
                                 [--aug-alternatives-from-sampling {true,false}]
                                 [--aug-alternatives-from-evidence {true,false}]
                                 [--aug-uniqueGeneId {true,false}]
                                 [--aug-gff3 {ON,OFF,on,off}]
                                 [--aug-sample AUG_SAMPLE]
                                 [--aug-noInFrameStop {true,false}]
                                 [--aug-maxtracks AUG_MAXTRACKS]
                                 [--aug-singlestrand {true,false}]
                                 [--aug-strand {both,forward,backward}]
                                 [--aug-min-intron-len AUG_MIN_INTRON_LEN]
                                 [--augustus-weights AUGUSTUS_WEIGHTS [AUGUSTUS_WEIGHTS ...]]
                                 [--additional-augustus-options ADDITIONAL_AUGUSTUS_OPTIONS]
                                 [--extrinsic-file-augustus-hints extrinsic_file_augustus_hints]
                                 [--augustus-hints-weights AUGUSTUS_HINTS_WEIGHTS [AUGUSTUS_HINTS_WEIGHTS ...]]
                                 [--additional-augustus-hints-options ADDITIONAL_AUGUSTUS_HINTS_OPTIONS]
                                 [--geneid-path GENEID_PATH]
                                 [--geneid-weights GENEID_WEIGHTS [GENEID_WEIGHTS ...]]
                                 [--geneid-options GENEID_OPTIONS]
                                 [--geneid-introns-weights GENEID_INTRONS_WEIGHTS [GENEID_INTRONS_WEIGHTS ...]]
                                 [--geneid-introns-options GENEID_INTRONS_OPTIONS]
                                 [--gmk-min-contig GMK_MIN_CONTIG]
                                 [--gmk-max-contig GMK_MAX_CONTIG]
                                 [--gmk-max-gap GMK_MAX_GAP]
                                 [--gmk-cores GMK_CORES]
                                 [--additional-genemark-options ADDITIONAL_GENEMARK_OPTIONS]
                                 [--genemark-weights GENEMARK_WEIGHTS [GENEMARK_WEIGHTS ...]]
                                 [--additional-genemark-ET-options ADDITIONAL_GENEMARK_ET_OPTIONS]
                                 [--genemark-ET-weights GENEMARK_ET_WEIGHTS [GENEMARK_ET_WEIGHTS ...]]
                                 [--miniprot-path MINIPROT_PATH]
                                 [--miniprot-cores MINIPROT_CORES]
                                 [--miniprot-weights MINIPROT_WEIGHTS [MINIPROT_WEIGHTS ...]]
                                 [--additional-miniprot-options ADDITIONAL_MINIPROT_OPTIONS]
                                 [--pasa-cores PASA_CORES] [--pasadb pasadb]
                                 [--pasa-weights PASA_WEIGHTS [PASA_WEIGHTS ...]]
                                 [--create-database] [--aligners ALIGNERS]
                                 [--add-pasa-option ADD_OPTION]
                                 [--transdecoder-weights TRANSDECODER_WEIGHTS [TRANSDECODER_WEIGHTS ...]]
                                 [--evm-path EVM_PATH]
                                 [--evm-segmentsize EVM_SEGMENTSIZE]
                                 [--evm-overlapsize EVM_OVERLAPSIZE]
                                 [--evm-cores EVM_CORES]
                                 [--additional-evm-options ADDITIONAL_EVM_OPTIONS]
                                 [--project-name project_name project_name]
                                 [--update-rounds UPDATE_ROUNDS]
                                 [--ncRNA-version NCRNA_VERSION] [--Rfam RFAM]
                                 [--cmsearch-CPUs CMSEARCH_CPUS]
                                 [--blast-threads BLAST_THREADS]

Create a configuration json file for the repeat annotation pipeline.

optional arguments:
  -h, --help            show this help message and exit

General Parameters:
  --configFile configFile
                        Configuration file with the pipeline parameters to be
                        created. Default Annotation.config
  --specFile specFile   Cluster specifications JSON fileto be generated.
                        Default Annotation.spec
  --basename base_name  Assembly basename. Default None
  --no-augustus         If specified, do not run augustus step.
  --no-augustus-hints   If specified, do not run augustus with hints step.
  --no-geneid           If specified, do not run geneid step.
  --no-geneid-introns   If specified, do not run geneid with introns.
  --no-genemark         If specified, do not run genemark step.
  --no-genemark-ET      If specified, do not run genemark-ET step.
  --no-miniprot         If specified, do not run miniprot step.
  --no-pasa             If specified, do not run pasa step.
  --no-transdecoder     If specified, do not run transdecoder step.
  --no-EVM              If specified, do not run EVM step.
  --no-update           If specified, do not run the annotation update step.
  --no-noncoding        If specified, do not run the non_coding annotation
                        step.
  --no-jbrowse          If specified, do not run the get jbrowse tracks step.

Inputs:
  --scripts-dir SCRIPTS_DIR
                        Directory with the different scripts for the pipeline.
                        Default bin/../scripts/
  --genome genome       Path to the genome assembly in fasta format.
  --glen glen           Path to the assembly.genome file.
  --genome-masked genome_masked
                        Path to the masked genome assembly in fasta format.
  --junctions junctions
                        Path to the junctions gff file to run gene predictors
                        with introns.
  --geneid-parameters geneid_parameters
                        Path to the geneid parameters file. For geneid, geneid
                        with introns and framefixing (part of annotation
                        update) steps.
  --proteins proteins   Path to the fasta with protein evidence.
  --transcripts transcripts
                        Path to the fasta with transcript evidence.
  --pasa-config pasa_config
                        Path to the pasa configuration file.
  --trans-gtf trans_gtf
                        Path to the fasta with transcript evidence.
  --update-config update_config
                        Path to the Pasa configuration file.
  --RM-gff RM_gff       Path to the Repeat Masker gff output.

Outputs:
  --annotation-step ANNOTATION_STEP
                        Step of the annotation pipeline in the annotation
                        process. Default 3
  --annotation-version ANNOTATION_VERSION
                        Version of the annotation process. Default 01
  --pipeline-workdir PIPELINE_WORKDIR
                        Base directory for the pipeline run. Default
                        /software/assembly/pipelines/Annotation_AAT/
  --output-dir OUTPUT_DIR
                        Directory to keep the outputs of the first annotation
                        steps. Default /software/assembly/pipelines/Annotation
                        _AAT/step03_annotation_pipeline.V01
  --EVM-dir EVM_DIR     Directory to keep the files for the EVM step. Default
                        step04_EVM.V01
  --dir-masked-chunks DIR_MASKED_CHUNKS
                        Directory to keep the chunks of the masked genome.
                        Default /software/assembly/pipelines/Annotation_AAT/st
                        ep03_annotation_pipeline.V01/chunks_masked_reference/
  --dir-genome-chunks DIR_GENOME_CHUNKS
                        Directory to keep the chunks of the genome. Default /s
                        oftware/assembly/pipelines/Annotation_AAT/step03_annot
                        ation_pipeline.V01/chunks_genome_reference/
  --augustus-prediction AUGUSTUS_PREDICTION
                        Output file for the augustus predictions. Default /sof
                        tware/assembly/pipelines/Annotation_AAT/step03_annotat
                        ion_pipeline.V01/gene_predictions/augustus/augustus_ge
                        ne_prediction.gff3
  --augustus-preEVM AUGUSTUS_PREEVM
                        Output file for the augustus predictions converted for
                        EVM. Default /software/assembly/pipelines/Annotation_A
                        AT/step03_annotation_pipeline.V01/gene_predictions/aug
                        ustus/augustus_preEVM.gff3
  --augustus-hints-prediction AUGUSTUS_HINTS_PREDICTION
                        Output file for the augustus with hints predictions.
                        Default /software/assembly/pipelines/Annotation_AAT/st
                        ep03_annotation_pipeline.V01/gene_predictions/augustus
                        _hints/augustus_hints_gene_prediction.gff3
  --augustus-hints-preEVM AUGUSTUS_HINTS_PREEVM
                        Output file for the augustus with hints predictions
                        converted for EVM. Default /software/assembly/pipeline
                        s/Annotation_AAT/step03_annotation_pipeline.V01/gene_p
                        redictions/augustus_hints/augustus_hints_preEVM.gff3
  --geneid-prediction GENEID_PREDICTION
                        Output file for the geneid predictions. Default /softw
                        are/assembly/pipelines/Annotation_AAT/step03_annotatio
                        n_pipeline.V01/gene_predictions/geneid/geneid_gene_pre
                        diction.gff3
  --geneid-preEVM GENEID_PREEVM
                        Output file for the geneid predictions converted for
                        EVM. Default /software/assembly/pipelines/Annotation_A
                        AT/step03_annotation_pipeline.V01/gene_predictions/gen
                        eid/geneid_preEVM.gff3
  --geneid-introns-prediction GENEID_INTRONS_PREDICTION
                        Output file for the geneid with introns predictions.
                        Default /software/assembly/pipelines/Annotation_AAT/st
                        ep03_annotation_pipeline.V01/gene_predictions/geneid_w
                        ith_introns/geneid_introns_gene_prediction.gff3
  --geneid-introns-preEVM GENEID_INTRONS_PREEVM
                        Output file for the geneid with introns predictions
                        converted for EVM. Default /software/assembly/pipeline
                        s/Annotation_AAT/step03_annotation_pipeline.V01/gene_p
                        redictions/geneid_with_introns/geneid_introns_preEVM.g
                        ff3
  --genemark-prediction GENEMARK_PREDICTION
                        Output file for the genemark predictions. Default /sof
                        tware/assembly/pipelines/Annotation_AAT/step03_annotat
                        ion_pipeline.V01/gene_predictions/genemark.gtf
  --genemark-preEVM GENEMARK_PREEVM
                        Output file for the genemark predictions converted for
                        EVM. Default /software/assembly/pipelines/Annotation_A
                        AT/step03_annotation_pipeline.V01/gene_predictions/gen
                        emark_preEVM.gff3
  --genemark-ET-prediction GENEMARK_ET_PREDICTION
                        Output file for the genemark-ET predictions. Default /
                        software/assembly/pipelines/Annotation_AAT/step03_anno
                        tation_pipeline.V01/gene_predictions/genemark-ET.gtf
  --genemark-ET-preEVM GENEMARK_ET_PREEVM
                        Output file for the genemark-ET predictions converted
                        for EVM. Default /software/assembly/pipelines/Annotati
                        on_AAT/step03_annotation_pipeline.V01/gene_predictions
                        /genemark-ET_preEVM.gff3
  --miniprot-cds MINIPROT_CDS
                        Output file for the miniprot output in a cds gff3
                        format. Default /software/assembly/pipelines/Annotatio
                        n_AAT/step03_annotation_pipeline.V01/protein_and_trans
                        cript_mappings/miniprot/proteins_miniprot_cds.gff3
  --miniprot-gene MINIPROT_GENE
                        Output file for the miniprot output in a gene gff3
                        format. Default /software/assembly/pipelines/Annotatio
                        n_AAT/step03_annotation_pipeline.V01/protein_and_trans
                        cript_mappings/miniprot/proteins_miniprot_gene.gff3
  --pasa-dir PASA_DIR   Directory to keep all the pasa outputs. Default /softw
                        are/assembly/pipelines/Annotation_AAT/step03_annotatio
                        n_pipeline.V01/protein_and_transcript_mappings/pasa/
  --evm-out EVM_OUT     File with the final EVM models. Default
                        step04_EVM.V01/evm.best.gff3
  --update-dir UPDATE_DIR
                        Directory to keep the files for annotation update
                        step. Default step05_annotation_update.V01/
  --ncRNA-dir NCRNA_ANNOTATION_DIR
                        Directory to keep the files of the ncRNA annotation
                        step. Default step06_ncRNA_annotation.V01/
  --out-cmsearch OUT_CMSEARCH
                        Output file to keep the cmsearch results. Default
                        step06_ncRNA_annotation.V01//cmsearch.tbl
  --out-tRNAscan OUT_TRNASCAN
                        Output file to keep the tRNAscan-SE results. Default
                        step06_ncRNA_annotation.V01//tRNAscan-SE/tRNAscan.out

Chunks:
  --masked-chunks MASKED_CHUNKS
                        Number of chunks of the masked genome for
                        parallelizing some gene predictors run. Default 50
  --genome-chunks GENOME_CHUNKS
                        Number of chunks of the genome for parallelizing
                        cmsearch. Default 1
  --protein-chunks PROTEIN_CHUNKS
                        Number of chunks to split the protein files for
                        running blast and classify the lncRNAs. Default 20

Augustus parameters:
  --species aug_species
                        Species name to run augustus with its trained
                        parameters. For augustus and augustus with hints
                        steps.
  --aug-alternatives-from-sampling {true,false}
                        Report alternative transcripts generated through
                        probabilistic sampling. Default true
  --aug-alternatives-from-evidence {true,false}
                        Report alternative transcripts when they are suggested
                        by hints. Default true
  --aug-uniqueGeneId {true,false}
                        If true, output gene identifyers like this:
                        seqname.gN. For augustus and augustus with hints.
                        Default true
  --aug-gff3 {ON,OFF,on,off}
                        Output in gff3 format. For augustus and augustus with
                        hints. Default on
  --aug-sample AUG_SAMPLE
                        For augustus and augustus with hints. Default 60
  --aug-noInFrameStop {true,false}
                        Do not report transcripts with in-frame stop codons.
                        Otherwise, intron-spanning stop codons could occur.
                        For augustus and augustus with hints. Default true
  --aug-maxtracks AUG_MAXTRACKS
                        Maximum number of tracks allowed. For augustus and
                        augustus with hints. Default 2
  --aug-singlestrand {true,false}
                        Predict genes independently on each strand, allow
                        overlapping genes on opposite strands. For augustus
                        and augustus with hints. Default false
  --aug-strand {both,forward,backward}
                        For augustus and augustus with hints. Default both
  --aug-min-intron-len AUG_MIN_INTRON_LEN
                        Minimum predicted intron length. For augustus and
                        augustus with hints. Default 30
  --augustus-weights AUGUSTUS_WEIGHTS [AUGUSTUS_WEIGHTS ...]
                        Weights given to augustus predictions when running
                        EVM. Specify the weight for each EVM run separated by
                        a space. Example 2 2 1
  --additional-augustus-options ADDITIONAL_AUGUSTUS_OPTIONS
                        Additional augustus options to run it, see augustus
                        help for more information about the possible options.

Augustus hints parameters:
  --extrinsic-file-augustus-hints extrinsic_file_augustus_hints
                        Path to the Extrinsic file to use when running
                        augustus with hints. Default /software/assembly/conda/
                        augustus3.5.0/config/extrinsic/extrinsic.E.cfg
  --augustus-hints-weights AUGUSTUS_HINTS_WEIGHTS [AUGUSTUS_HINTS_WEIGHTS ...]
                        Weights given to augustus with intron predictions when
                        running EVM. Specify the weight for each EVM run
                        separated by a space. Example 3 3 3
  --additional-augustus-hints-options ADDITIONAL_AUGUSTUS_HINTS_OPTIONS
                        Desired augustus with intron options to run it, see
                        augustus documentation for more information.

Geneid parameters:
  --geneid-path GENEID_PATH
                        Path to the installation of geneid. Default
                        /software/assembly/src/geneid/
  --geneid-weights GENEID_WEIGHTS [GENEID_WEIGHTS ...]
                        Weights given to geneid predictions when running EVM.
                        Specify the weight for each EVM run separated by a
                        space. Example 2 1 2
  --geneid-options GENEID_OPTIONS
                        Desired geneid options to run it, see geneid
                        documentation for more information. Default -3U

Geneid Introns parameters:
  --geneid-introns-weights GENEID_INTRONS_WEIGHTS [GENEID_INTRONS_WEIGHTS ...]
                        Weights given to geneid with intron predictions when
                        running EVM. Specify the weight for each EVM run
                        separated by a space. Example 3 3 3
  --geneid-introns-options GENEID_INTRONS_OPTIONS
                        Desired geneid with intron options to run it, see
                        geneid documentation for more information. Default -3U

Genemark parameters:
  --gmk-min-contig GMK_MIN_CONTIG
                        Will ignore contigs shorter then min_contig in
                        training. Default 50000
  --gmk-max-contig GMK_MAX_CONTIG
                        Will split input genomic sequence into contigs shorter
                        than max_contig. Default 5000000
  --gmk-max-gap GMK_MAX_GAP
                        Will split sequence at gaps longer than max_gap.
                        Letters 'n' and 'N' are interpreted as standing within
                        gaps. Default 5000
  --gmk-cores GMK_CORES
                        Number of threads for running genemark. Default 24
  --additional-genemark-options ADDITIONAL_GENEMARK_OPTIONS
                        Additional genemark options to run it, see genemark
                        documentation for more information.
  --genemark-weights GENEMARK_WEIGHTS [GENEMARK_WEIGHTS ...]
                        Weights given to genemark predictions when running
                        EVM. Specify the weight for each EVM run separated by
                        a space. Example 1 1 1

Genemark-ET parameters:
  --additional-genemark-ET-options ADDITIONAL_GENEMARK_ET_OPTIONS
                        Additional genemark-ET options to run it, see genemark
                        documentation for more information.
  --genemark-ET-weights GENEMARK_ET_WEIGHTS [GENEMARK_ET_WEIGHTS ...]
                        Weights given to genemark-ET predictions when running
                        EVM. Specify the weight for each EVM run separated by
                        a space. Example 3 3 3

Miniprot parameters:
  --miniprot-path MINIPROT_PATH
                        Path to Miniprot installation. Default
                        /software/assembly/src/miniprot/miniprot/
  --miniprot-cores MINIPROT_CORES
                        Number of threads. Default 8
  --miniprot-weights MINIPROT_WEIGHTS [MINIPROT_WEIGHTS ...]
                        Weights given to miniprot mappings when running EVM.
                        Specify the weight for each EVM run separated by a
                        space. Example 10 8 10
  --additional-miniprot-options ADDITIONAL_MINIPROT_OPTIONS
                        Additional miniprot options to run it, see miniprot
                        help for more information about the possible options.

Pasa parameters:
  --pasa-cores PASA_CORES
                        Number of CPUs to run Pasa. Default 1
  --pasadb pasadb       Name of the pasa database, it must coincide with the
                        name in pasa_config.
  --pasa-weights PASA_WEIGHTS [PASA_WEIGHTS ...]
                        Weights given to pasa mappings when running EVM.
                        Specify the weight for each EVM run separated by a
                        space. Example 8 10 8
  --create-database     If specified, create pasa database.
  --aligners ALIGNERS   Program to map the transcripts with.
  --add-pasa-option ADD_OPTION
                        Option given to add extra options to PASA

Transdecoder parameters:
  --transdecoder-weights TRANSDECODER_WEIGHTS [TRANSDECODER_WEIGHTS ...]
                        Weights given to pasa transdecodergff3 output file
                        when running EVM. Specify the weight for each EVM run
                        separated by a space. Example 3 2 3

Evm parameters:
  --evm-path EVM_PATH   Path to the EVM software installation. Default /scratc
                        h/project/devel/aateam/src/EVidenceModeler-1.1.1/
  --evm-segmentsize EVM_SEGMENTSIZE
                        Size of the genome partitions for EVM. Default 2000000
  --evm-overlapsize EVM_OVERLAPSIZE
                        Size of the overlap between the different EVM
                        partitions. Default 1000000
  --evm-cores EVM_CORES
                        Number of threads to run EVM. Default 24
  --additional-evm-options ADDITIONAL_EVM_OPTIONS
                        Additional augusevmtus options to run it, see evm help
                        for more information about the possible options.

Update parameters:
  --project-name project_name project_name
                        Name of the project and version of the annotation
                        space separated, to give the names to the final
                        annotation output.
  --update-rounds UPDATE_ROUNDS
                        Number of rounds to run PASA updates. Default 2

ncRNA Annotation parameters:
  --ncRNA-version NCRNA_VERSION
                        Version for the ncRNA annotation. Default A
  --Rfam RFAM           CM file with the Rfam library. Default
                        /scratch/devel/jgomez/RFAM_db/290721/Rfam.cm
  --cmsearch-CPUs CMSEARCH_CPUS
                        Number of CPUs to run cmsearch Defau lt 16
  --blast-threads BLAST_THREADS
                        Number of CPUs to run blast Defau lt 4
```
