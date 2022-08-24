/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {

    tag "fastqc on: $name"

    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("${name}_fastqc.zip"), path("${name}_fastqc.html"), emit: fastqc

    script:

        """
        fastqc $reads
        """

}

process MULTIQC {

    tag "multiqc on: ${ params.MQCLabel }"

    publishDir params.multiqc_outputs_dir, mode: 'link'

    input:
        // path("samples.txt")
        path("mqc_in/*") // any number of multiqc compatible files
        // path(multiqc_config)

    output:
        path("${ params.MQCLabel }_multiqc.html"), emit: html
        path("${ params.MQCLabel }_multiqc_report.zip"), emit: zipped_report

    script:

        """
        multiqc --interactive -o ${ params.MQCLabel }_multiqc_report \
                -n ${ params.MQCLabel }_multiqc mqc_in \

        mv ${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html .
        
        zip -m -r '${ params.MQCLabel }_multiqc_report.zip' '${ params.MQCLabel }_multiqc_report'
        """

}

process TRIMGALORE {

    debug true

    tag "trimgalore on: $name"

    publishDir params.filtered_reads_dir, mode: 'link'

    input:
        tuple val(name), path(reads)

    script:
    
        """
        # this depends on the lib_type and then if paired-end or not
        if [ ${params.lib_type} == 1 ]; then

            if [ ${params.single_end} == 'true' ]; then

                printf "    first lib type, single-end\n"

            else

                printf "    first lib type, paired-end\n"

            fi
        
        elif [ ${params.lib_type} == 2 ]; then

            printf "    second lib type\n"

        elif [ ${params.lib_type} == 3 ]; then

            printf "    third lib type\n"

        fi
        """
}
