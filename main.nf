#!/usr/bin/env nextflow

params.singleEnd = false
params.pairedEnd = false
params.bam = false
params.reads = "data/*{1,2}.fastq.gz"
params.outdir="$baseDir/results"
/*params.GRCh38="ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"*/
/*params.GRCh38="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt"*/
params.krakendb="$baseDir/results/databases/HumanViral"
params.kraken2_build = false


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

process renameBam {

        tag "$bam"
        publishDir "${params.outdir}/renameBam", mode: 'copy'

        when: params.bam

        input:
        file bam from ch_input_bamtofastq
      
 	output:
	file ("*.bam") into ch_renamebam
		
 
	
        script:

        new_bam="\$(samtools view -h ${bam} | grep -P '\tSM:' | head -n 1 | sed 's/.\\+SM:\\(.\\+\\)/\\1/' | sed 's/\t.\\+//' | sed 's/\\s/_/g')"
                       
        """
        cp ${bam} "$new_bam".bam
        """

}


process BamToFastq {
	
	tag "$bam"
	publishDir "${params.outdir}/BamToFastq", mode: 'copy'


	when: params.bam 

	input:
        file bam from ch_renamebam

	output:
	set val("${base}"), file("*.fastq.gz") into ch_output_bamtofastq

	script:
        base = "${bam.baseName}"
	"""
	samtools fastq -tn ${bam} | gzip > ${base}.fastq.gz 
	"""
} 


if (params.bam) {
	ch_output_bamtofastq
        .view()
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
	fastqc -t "${task.cpus}" -q $reads --extract
	"""
}

process fastp {

	tag "$name"
	publishDir  "${params.outdir}/fastp", mode: 'copy'

	input:
	set val(name), file(reads) from ch_input_fastp

	output:
	set val(name), file("*_trimmed.fastq.gz") into trimmed_reads	
	file("*.html")
	file ("*_fastp.json")

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

	
	tag "$name"
	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
	set val(name), file ('fastqc/*') from fastqc_results.collect()

	output:
	file('*_multiqc_report.html')

	script:
	"""
	multiqc .
	"""
} 

process krakenBuild {

	publishDir "${params.outdir}/databases", mode: 'copy'

	when: params.kraken2_build


	script:

	"""
        mkdir -p ${params.outdir}/databases
	UpdateKrakenDatabases.py ${params.outdir}/databases

			
	"""
	

} 


process kraken2 {
	
	publishDir "${params.outdir}/kraken2", mode: 'copy'


	input:
        set val(name), file(reads) from trimmed_reads


	output:
        set val(name), file('*.kraken.out') into kraken_out
        set val(name), file('*.kraken.report') into kraken_report
	set val(name), file("results.krona") into ch_krona 


	script:
        out = name+".kraken.out"
        kreport = name+".kraken.report"
        if (params.pairedEnd){
            """
            kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]} 
            """    
        } else {
            """
            kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport ${reads[0]}
            """
        }



}
