#!/usr/bin/env nextflow

params.singleEnd = false
params.pairedEnd = false
params.bam = false
params.reads = "data/*{1,2}.fastq.gz"
params.readPaths = false
params.outdir="$baseDir/results"


if (!params.bam){
     Channel
        .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs " +
            "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nValid input file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_input_for_BamToFastq } 

} else {
     Channel
        .fromPath( params.reads )
        .filter { it =~/.*.bam/ }
        .map { row -> [file( row )]}
        .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nValid input file types: '.bam'" +
            "to be enclosed in quotes!\n" }
        .into {ch_input_for_BamToFastq}

}



/*
 * PREPROCESSING - Convert BAM to FastQ if BAM instead of FastQ file(s)
 *
 */


process BamToFastq {

        publishDir "$params.outdir/BamToFastq", mode: 'copy'

                tag "$bam"

                when: params.bam


                input:
                file bam from ch_input_for_BamToFastq

                output:
                set val("${base}"), file ("*.fastq.gz") into ch_output_BamToFastq

                script:
                base = "${bam.baseName}"
                """
                samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.fastq.gz
                """
}

