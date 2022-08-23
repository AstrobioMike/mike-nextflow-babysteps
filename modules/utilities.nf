/*
 * Processes for general utilities.
 */

 process WRITE_SAMPLE_NAMES_TO_FILE {

    input:
        tuple val(name), path(reads)

    output:
        path("samples.txt")

    exec:

        out_file = path("samples.txt")

}
