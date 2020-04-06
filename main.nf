#!/usr/bin/env nextflow

params.singleEnd = false
params.pairedEnd = false
params.bam = false
params.reads = "data/*{1,2}.fastq.gz"
params.outdir="$baseDir/results"


/* 
 * Check if we get 4 files for paired-end reads
 *
 */


if (!params.bam){
	Channel
	.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
	.filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nValid input file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --singleEnd on the command line." }
	.view() 
	.into{ch_input_bamtofastq;ch_input_fastqc;ch_input_fastp} 

} else {
	Channel
	.fromPath( params.reads )
	.filter { it =~/.*.bam/ }
	.map { row -> [file( row )]}
	.ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nValid input file types: '.bam'" }
	.view() 
	.set{ch_input_bamtofastq} 


}


/*
 * PREPROCESSING - Change name and convert bam to fastq
 *
 */


process BamToFastq {

	tag "$bam"
	publishDir "$params.outdir/BamToFastq", mode: 'copy'


	when: params.bam 

	input:
	file bam from ch_input_bamtofastq

	output:
	set val("${base}"), file ("*.fastq.gz") into ch_output_bamtofastq

	script:
	base = "${bam.baseName}"
	"""
	new_bam="\$(samtools view -h ${bam} | grep -P '\tSM:' | head -n 1 | sed 's/.\\+SM:\\(.\\+\\)/\\1/' | sed 's/\t.\\+//' | sed 's/\\s/_/g')"	
	echo "\${new_bam}" 
	samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.fastq.gz
	mv ${base}.fastq.gz "\${new_bam}".fastq.gz
	"""
}


if (params.bam) {
	ch_output_bamtofastq
	.into {ch_input_fastqc; ch_input_fastp}

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

process fastp {

	tag "$name"
	publishDir  "${params.outdir}/fastp", mode: 'copy'

	input:
	set val(name), file(reads) from ch_input_fastp

	output:
	set val(name), file("${name}_*trimmed.fastq.gz") into trimmed_reads	
	file("${name}.html")
	file ("${name}_fastp.json")

	script:
	if(params.singleEnd || params.bam){

	"""
	fastp -i "${reads[0]}" -o "${name}_trimmed.fastq.gz" -j "${name}_fastp.json" -h "${name}.html" 
	"""
	}	
	else {
	"""
	fastp -i "${reads[0]}" -I "${reads[1]}" -o "${name}_1.trimmed.fastq.gz" -O "${name}_2.trimmed.fastq.gz" -j "${name}_fastp.json" -h "${name}.html" 

	"""
	}


}

process multiqc {

	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
	file ('fastqc/*') from fastqc_results.collect()

	output:
	file('multiqc_report.html')

	script:
	"""
	multiqc .
	"""
}


