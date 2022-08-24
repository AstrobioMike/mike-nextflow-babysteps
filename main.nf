#!/usr/bin/env nextflow
/*
==========================================================================================
Largely modified from the nf-core/methylseq workflow: https://github.com/nf-core/methylseq
==========================================================================================
*/

// Declare syntax version
nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

if (params.help) {
    log.info "\n┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅"
    log.info "┇          GeneLab Methyl-seq Workflow: $workflow.manifest.version            ┇"
    log.info "┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅\n"
    log.info "    Usage example (after setting parameters in the 'nextflow.config' file):"
    log.info "        `nextflow run main.nf`\n"
    exit 0
}

////////////////////////////////////////////////////
/* --                PROCESSES                 -- */
////////////////////////////////////////////////////

include { FASTQC } from './modules/QC.nf'

workflow {

    ch_input_reads = Channel.fromFilePairs( params.input_reads, size: params.single_end ? 1 : 2 ) { file -> file.name.replaceAll( /.fastq.gz|.fq.gz/,'' ) }
    

    reads_list_obj = ch_input_reads | map { it -> it[0] } 
                                    | collectFile( name: 'samples.txt', newLine: true, storeDir: "./" )

    FASTQC(ch_input_reads)
    
}
