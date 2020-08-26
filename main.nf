#!/usr/bin/env nextflow

params.singleEnd = false
params.pairedEnd = false
params.bam = false
params.longreads = false
params.reads = "data/*{1,2}.fastq.gz"
params.outdir="$baseDir/results"
params.krakendb="$baseDir/results/databases/HumanViral"
params.kaiju_db="$baseDir/results/databases/virus_kaiju"
params.db_build = false
params.GRCh38="/srv/rs6/sofia/Metoid/Metoid/results/databases/GCF_000001405.39_GRCh38.p13_genomic.fna"
human_ref=file(params.GRCh38)
params.accession_list="$baseDir/bin/accession_list.txt"
params.contaminants="/srv/rs6/sofia/Metoid/Metoid/results/Contaminants/contaminants.fna"
contaminants_file=file(params.contaminants)
params.porechopParam = "-t 4"
params.porechop = true

/* 
 * Get input data
 *
 */


if (params.bam){
        Channel
        .fromPath( params.reads )
        .filter { it =~/.*.bam/ }
        .map { row -> [file( row )]}
        .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nValid input file types: '.bam'" }
        .set{ch_input_bamtofastq}
}
else if (params.longreads){
        Channel
        .fromFilePairs( params.reads, size: 1 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nValid input file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }
        .into{ch_input_fastqc;ch_input_porechop} 
}
else {
        Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nValid input file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into{ch_input_fastqc;ch_input_fastp}
}


/*
 * PREPROCESSING BAM INPUT - Change name and convert bam to fastq
 *
 */

if (params.bam) {
process renameBam {

        tag "$bam"
        publishDir "${params.outdir}/renameBam", mode: 'copy'

        input:
        file bam from ch_input_bamtofastq
      
 	output:
	file ("*.bam") into ch_renamebam
		
        script:

        new_bam="\$(samtools view -h ${bam} | grep -P '\tSM:' | head -n 1 | sed 's/.\\+SM:\\(.\\+\\)/\\1/' | sed 's/\t.\\+//' | sed 's/\\s/_/g')"
                       
        """
        cp ${bam} "$new_bam".bam
        """

}} 

if (params.bam) {
process BamToFastq {
	
	tag "$bam"
	publishDir "${params.outdir}/BamToFastq", mode: 'copy'

	input:
        file bam from ch_renamebam

	output:
	set val("${base}"), file("*.fastq.gz") into ch_output_bamtofastq

	script:
        base = "${bam.baseName}"
	"""
	samtools fastq -tn ${bam} | gzip > ${base}.fastq.gz 
	"""
}} 


if (params.bam) {
	ch_output_bamtofastq
        .view()
	.into {ch_input_fastqc; ch_input_fastp} 

} 


/*
 * QC- fastqc
 * Adaptor trimming - fastp, porechop
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

if (params.longreads) {
process porechop {
        
        publishDir  "${params.outdir}/Porechop", mode: 'copy'
        
        when: params.porechop        

        input:
        set val(name), file(reads) from ch_input_porechop

        
        output:
        set val(name), file("*_chopped.fastq.gz") into (ch_input_fastp)

        script:
        """
        porechop -i $reads -o ${name}_chopped.fastq.gz --format fastq.gz $params.porechopParam
        """
       
}}

if (params.longreads && !params.porechop) {
        ch_input_porechop
        .view()
        .into {ch_input_fastp}
}


process fastp {

	tag "$name"
	publishDir  "${params.outdir}/fastp", mode: 'copy'

	input:
	set val(name), file(reads) from ch_input_fastp

	output:
	set val(name), file("*_trimmed.fastq.gz") into (trimmed_reads_bowtie2, trimmed_reads_fastqc)	
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


process fastqc_after_trimming {

	tag "$name"
        publishDir  "${params.outdir}/FastQC", mode: 'copy'

        input:
        set val(name), file(reads) from trimmed_reads_fastqc

        output:
        file "*_fastqc.{zip,html}" into fastqc_after_trimming


        script:
        """
        fastqc -t "${task.cpus}" -q $reads --extract
        """

}


process multiqc {

	
	tag "$name"
	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
	file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report		

	script:
	"""
	multiqc .
	"""
}


process retrieve_contaminants {

	publishDir "${params.outdir}/Contaminants", mode: 'copy'

	when: params.db_build

        script:

        """
        mkdir -p ${params.outdir}/Contaminants
        RetrieveContaminants.py ${params.outdir}/Contaminants ${params.accession_list} 
	"""

}


process index_contaminants {

	input:
	file cont_genomes from contaminants_file

	output:
	file 'index*' into ch_index_contaminants

	script:
	"""
	bowtie2-build $cont_genomes index
	"""

}


process index_host {

	when: params.db_build


	input:
	file host_genome from human_ref

	output:
	file 'index*' into ch_index_host

	script:
        """
	bowtie2-build $host_genome index
        """
} 

ch_index_host
	.mix (ch_index_contaminants)
	.set {bowtie2_input}


process bowtie2 {
	
	tag "$name"
	publishDir "${params.outdir}/bowtie2", mode: 'copy'

	input:
	set val(name), file(reads) from trimmed_reads_bowtie2
	file(index) from bowtie2_input.collect()

	output:
        set val(name), file("*_removed.fastq") into (ch_bowtie2_kraken, ch_bowtie2_kaiju)
	
	script:
	samfile = name+".sam"
	bamfile = name+"_mapped_and_unmapped.bam"	
	unmapped_bam = name+ "_only_unmapped.bam"
	mapped_bam = name+"_mapped.bam"
	bam_sorted = name+"_sorted.bam"
	bam_mapped_sorted = name+"_mapped_sorted.bam"

	index_name = index.toString().tokenize(' ')[0].tokenize('.')[0]	
	
	if (params.pairedEnd){
        """
        bowtie2 -x $index_name --local -1 ${reads[0]} -2 ${reads[1]} > $samfile
	samtools view -Sb $samfile > $bamfile
	samtools view -b -f4 $bamfile > $unmapped_bam
	samtools sort -n $unmapped_bam -o $bam_sorted
	bedtools bamtofastq -i $bam_sorted -fq "${name}_1.removed.fastq" -fq2 "${name}_2.removed.fastq"
        """
        }else {
        """
        bowtie2 -x $index_name --local -U $reads > $samfile
	samtools view -Sb $samfile > $bamfile
        samtools view -b -f4 $bamfile > $unmapped_bam
        samtools sort -n $unmapped_bam -o $bam_sorted
	bedtools bamtofastq -i $bam_sorted -fq "${name}_removed.fastq"
        """
        }

}


process krakenBuild {

	publishDir "${params.outdir}/databases", mode: 'copy'

	when: params.db_build


	script:

	"""
        mkdir -p ${params.outdir}/databases
	UpdateKrakenDatabases.py ${params.outdir}/databases 
		
	"""
	
} 


process kraken2 {
	
	tag "$name"
	
	publishDir "${params.outdir}/kraken2", mode: 'copy'


	input:
        set val(name), file(reads) from ch_bowtie2_kraken


	output:
        set val(name), file("*.kraken.out.txt") into kraken_out
        set val(name), file("*.kraken.report") into kraken_report 


	script:
        out = name+".kraken.out.txt"
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

process krona_taxonomy {
    

    output:
    file("taxonomy/taxonomy.tab") into file_krona_taxonomy


    script:
    """
    ktUpdateTaxonomy.sh taxonomy
    """
}




process krona_kraken {

    tag "$name"

    publishDir "${params.outdir}/krona_kraken", mode: 'copy'


    input:
    set val(name), file(report) from kraken_out
    file("taxonomy/taxonomy.tab") from file_krona_taxonomy

    output:
    file("*")

    script:

    """
    ktImportTaxonomy -o ${name}.krona.html -t 3 -s 4 ${name}.kraken.out.txt -tax taxonomy
    """
}

process kaiju {

    tag "$name"
   
    publishDir "${params.outdir}/kaiju", mode: 'copy'

    input:
    set val(name), file(reads) from ch_bowtie2_kaiju

    output:
    set val(name), file("*.kaiju.out") into kaiju_out
    set val(name), file("*.kaiju.out.krona") into ch_kaiju_krona

    script:
    out = name+".kaiju.out"	
    krona_kaiju = name + ".kaiju.out.krona"

    if (params.pairedEnd){
    """
    kaiju -x -v -t ${params.kaiju_db}/nodes.dmp -f ${params.kaiju_db}/kaiju_db_viruses.fmi -i ${reads[0]} -j ${reads[1]} -o $out
    kaiju2krona -t ${params.kaiju_db}/nodes.dmp -n ${params.kaiju_db}/names.dmp -i $out -o $krona_kaiju    
    """
    } else {
    """
    kaiju -x -v -t ${params.kaiju_db}/nodes.dmp -f ${params.kaiju_db}/kaiju_db_viruses.fmi -i ${reads[0]} -o $out
    kaiju2krona -t ${params.kaiju_db}/nodes.dmp -n ${params.kaiju_db}/names.dmp -i $out -o $krona_kaiju
            
    """
    }

}

process krona_kaiju {

    tag "$name"

    publishDir "${params.outdir}/krona_kaiju", mode: 'copy'


    input:
    set val(name), file(report) from ch_kaiju_krona
    file("taxonomy/taxonomy.tab") from file_krona_taxonomy

    output:
    file("*")

    script:

    """
    ktImportText -o ${name}.kaiju.html ${name}.kaiju.out.krona
    """
}


