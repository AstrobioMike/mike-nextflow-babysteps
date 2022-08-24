/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
    publishDir params.multiqc_outputs_dir, mode: 'link'

    tag "$name"

    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("${name}_fastqc.zip"), path("${name}_fastqc.html"), emit: zip

    script:

        """
        fastqc $reads
        """

}

// process MULTIQC {

//     input:
//         path("samples.txt")
//         path("mqc_in/*") // any number of multiqc compatible files
//         path(multiqc_config)

//   output:
//     path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html"), emit: html
//     path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc_data"), emit: data
//     path("${ params.MQCLabel }_multiqc_report.zip"), emit: zipped_report
//     path("${ params.MQCLabel }_multiqc_report"), emit: unzipped_report
//     path("versions.txt"), emit: version

//   script:
//     config_arg =  multiqc_config.name != "NO_FILE" ? "--config ${ multiqc_config }" : ""
//     """
//     multiqc --sample-names samples.txt  \
//             --interactive -o ${ params.MQCLabel }_multiqc_report \
//             -n ${ params.MQCLabel }_multiqc mqc_in \
//             ${ config_arg }
//     zip -r '${ params.MQCLabel }_multiqc_report.zip' '${ params.MQCLabel }_multiqc_report'
//     multiqc --version > versions.txt
//     """
// }