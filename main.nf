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
    // ch_input_reads | view

    // ch_input_reads_var = ch_input_reads | flatten
    // ch_input_reads_var | view


    // ch_input_reads.subscribe { $it[0] }

    // reads_list_obj = ch_input_reads.toList().val[0]
    // println reads_list_obj

    // ch_input_reads | collect() | view

    // ch_input_reads | view

    // mylist = [1,2,3]

    // println mylist.getClass()


    // println reads_list_obj

    // println reads_list_obj.getClass()

    // for ( entry in reads_list_obj ) {

    //     println(entry[0])

    // }

    // reads_list_obj.getAt(0)

    // out_file = file("samples.txt")

    // out_file.text = 'stuff?'


    // FASTQC(ch_input_reads)

    // FASTQC.out.zip | view
    
}
