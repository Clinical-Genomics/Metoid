#!/usr/bin/env nextflow

params.singleEnd = false
params.pairedEnd = false
params.bam = false
params.reads = "data/*{1,2}.fastq.gz"
params.readPaths = false
params.outdir="$baseDir/results"


/* 
 * Check if we get 4 files for paired-end reads
 *
 */


if (!params.bam){
     Channel
        .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs " +
            "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nValid input file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --single_end on the command line." }
        .into {ch_input_bamtofastq; ch_input_fastqc} 

} else {
     Channel
        .fromPath( params.reads )
        .filter { it =~/.*.bam/ }
        .map { row -> [file( row )]}
        .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nValid input file types: '.bam'" +
            "to be enclosed in quotes!\n" }
        .into {ch_input_bamtofastq; ch_input_fastqc}

}



/*
 * PREPROCESSING - Convert BAM to FastQ if BAM instead of FastQ file(s)
 *
 */


process BamToFastq {

	tag "$bam"
        publishDir "$params.outdir/BamToFastq", mode: 'copy'

   
        when: params.bam


        input:
        file bam from ch_input_bamtofastq

        output:
        set val("${base}"), file ("*.fastq.gz") into ch_output_BamToFastq

        script:
        base = "${bam.baseName}"
        """
        samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.fastq.gz
        """
}


/*
 * QC- fastqc
 * Adaptor trimming- fastp
 * Do we need fastqc control on trimmed data
 */



process fastqc {

	tag "$name"
	publishDir  "${params.outdir}/FastQC", mode: 'copy'


	input:
	set val(name), file(reads) from ch_input_fastqc

	output:
	file "*_fastqc.{zip,html}" into fastqc_results


	script:
	"""
    	fastqc -t "${task.cpus}" -q $reads
    	"""


}

