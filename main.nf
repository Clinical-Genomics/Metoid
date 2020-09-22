#!/usr/bin/env nextflow
def help() {
        log.info"""
        nextflow <nextflow options> run main.nf <options>
        
        IMPORTANT
        -profile Comma-separated list of run profiles to use
        
                 Process executor profiles:
                 local - run processes locally
                 slurm - run processes using a slurm resource manager
                 
                 Container and environments profiles:
                 docker - run processes in docker images
                 singularity - run processes in singularity images
        
        --reads  input files for sample [default:"data/*{1,2}.fastq.gz"]
        --outdir directory for output [default:"results/"]
        
        --singleEnd   Input is single read data [default:false]
        --pairedEnd   Input is paired end data [default:false]
        --bam         Input bam file [default:false]
        --longRead    Input is long read data, e.g. nanopore [default:false]
        
        CONFIG FILE
        Nextflow will by default find configuration files in the following locations:
                  "nextflow.config" in current directory
                  "nextflow.config" in script directory
                  "~/.nextflow/config"

        Specify config file with nextflow options:        
        -c        Path to additional config file
        -C        Path to config file. Other config files are ignored.
        
        OPTIONS
        --build              Build all databases and index all references [default:false]
        --buildKraken2       Build Kraken2 database [default:false]
        --buildHost          Index host reference [default:false]
        --buildContaminants  Index contaminant references [default:false]

        --porechop           Run Porechop trimming (only if long read data) [default:true]
        --fastp              Run Fastp trimming [default:true]
        --fastqcPre          Run fastqc on raw data [default:true]
        --fastqcPost         Run fastqc on trimmed data [default:true]
        --filterBowtie2      Run filtering of host and/or contaminants [default:true]
        --kraken2            Run kraken2 taxonomic classification [default:true]
        --kaiju              Run kaiju taxonomic classification [default:true]
        --multiqc            Run multiqc [default:true]
        --help               Display this help message and exit
        
        PARAMETERS
        Refer to the config file to see or change the default parameters used in the pipeline.
        
        To change a parameter for a single run, use:
        --<param> "value"
        
        --porechopParam          Settings for trimming with Porechop
        --fastpParamLong         Settings for trimming with Fastp for long read data
        --fastpParamShort        Settings for trimming with Fastp for short read data
        --krakenDB               Kraken2 database path
        --kaijuDB                Kaiju database path
        --references             Location of indexed references
        --hostReference          Host reference fasta path
        --hostFTP                Alternatively: FTP for downloading host reference
        --contaminantReferences  Fasta file with contaminant species
        --contaminantList        Alternatively: File with accession numbers for contaminants
        """
}

if (params.help) {
        help()
        exit 0
}

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
        params.fastpParam = params.fastpParamShort
}
else if (params.longRead){
        Channel
        .fromFilePairs( params.reads, size: 1 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nValid input file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }
        .into{ch_input_fastqc;ch_input_porechop}
        params.fastpParam = params.fastpParamLong 
}
else {
        Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nValid input file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into{ch_input_fastqc;ch_input_fastp}
        params.fastpParam = params.fastpParamShort
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

        when params.fastqcPre
	
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

if (params.longRead) {
process porechop {
        
        publishDir  "${params.outdir}/Porechop", mode: 'copy'
        
        when: params.porechop        

        input:
        set val(name), file(reads) from ch_input_porechop

        output:
        set val(name), file("*_chopped.fastq.gz") into ch_input_fastp

        script:
        """
        porechop -i $reads -o ${name}_chopped.fastq.gz --format fastq.gz $params.porechopParam
        """
       
}}

if (params.longRead && !params.porechop) {
        ch_input_porechop
        .view()
        .set {ch_input_fastp}
}

 if (params.fastp) {
process fastp {

	tag "$name"
	publishDir  "${params.outdir}/fastp", mode: 'copy'

        when: params.fastp

	input:
	set val(name), file(reads) from ch_input_fastp

	output:
	set val(name), file("*.trimmed.fastq.gz") into (trimmed_reads_bowtie2, trimmed_reads_fastqc)	
	file("*.html")
	file ("*_fastp.json")

	script:
	if(params.singleEnd || params.bam){
        """
        fastp -i "${reads[0]}" -o "${name}.trimmed.fastq.gz" -j "${name}_fastp.json" -h "${name}.html" $params.fastpParam 
        """
        }       
        else {
        """
        fastp -i "${reads[0]}" -I "${reads[1]}" -o "${name}_1.trimmed.fastq.gz" -O "${name}_2.trimmed.fastq.gz" -j "${name}_fastp.json" -h "${name}.html" $params.fastpParam 
        """
	} 
}}

else {
        ch_input_fastp
        .view()
        .into {trimmed_reads_bowtie2; trimmed_reads_fastqc}
}

process fastqc_after_trimming {

        when params.fastqcPost

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

        when: params.multiqc
	
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

/*
 * BUILD DATABASES
 *
 */

if (params.contaminantReferences == "") {
process retrieve_contaminants {

	publishDir "${params.references}/contaminants/", mode: 'copy'

	when: params.build || params.buildContaminants
        
        input:
        path contaminants from params.contaminantList

        output:
        file "contaminants.fna" into ch_contaminants_reference

        script:
        """
        RetrieveContaminants.py . ${contaminants}
        """
}}

else {
        Channel
        .fromPath(params.contaminantReferences)
        .set {ch_contaminants_reference}
}

process index_contaminants {

        publishDir "${params.references}/contaminants/", mode: 'copy'

        when: params.build || params.buildContaminants

        input:
        file genome from ch_contaminants_reference

        output:
        file 'index*' into ch_contaminants_index

        script:
        """
        bowtie2-build ${genome} index
        """
}

// Get prebuilt contaminant indices
if (params.filterContaminants && !params.build && !params.buildContaminants) {
        Channel
        .fromPath("${params.references}/contaminants/*.bt2")
        .set {ch_contaminants_index}
}

if (params.hostReference == "") {
process retrieve_host {

        publishDir "${params.references}/host", mode: 'copy'

        when: params.build || params.buildHost

        output:
        file "*.fna.gz" into ch_host_reference        

        script:
        """
        wget ${params.hostFTP} 
        """
}}

else {
        Channel
        .fromPath(params.hostReference)
        .set {ch_host_reference}
}

process index_host {

        publishDir "${params.references}/host", mode: 'copy'

        when: params.build || params.buildHost
      
        input:
        file genome from ch_host_reference
 
        output:
        file 'index*' into ch_host_index
 
        script:
        """
        bowtie2-build ${genome} index
        """
}

// Get prebuilt host indices
if (!params.build && !params.buildHost) {
        Channel
        .fromPath("${params.references}/host/*.bt2")
        .set {ch_host_index}
}

// Join host and contaminant indices
if (params.filterContaminants) {
ch_host_index
	.mix (ch_contaminants_index)
	.set {bowtie2_input} 
}
else {
ch_host_index.set {bowtie2_input}
}

process krakenBuild {

        publishDir "${params.outdir}/databases", mode: 'copy'

        when: params.build || params.buildKraken2

        script:

        """
        mkdir -p ${params.outdir}/databases
        UpdateKrakenDatabases.py ${params.outdir}/databases
        """
}

/*
 * FILTER READS
 *
 */

if (params.filterBowtie2) {
process bowtie2 {
	
	tag "$name"
	publishDir "${params.outdir}/bowtie2", mode: 'copy'

	input:
	set val(name), file(reads) from trimmed_reads_bowtie2
	file(index) from bowtie2_input.collect()

	output:
        set val(name), file("*.removed.fastq.gz") into (ch_bowtie2_kraken, ch_bowtie2_kaiju)
	
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
        gzip "*.fastq"
        """
        }else {
        """
        bowtie2 -x $index_name --local -U $reads > $samfile
	samtools view -Sb $samfile > $bamfile
        samtools view -b -f4 $bamfile > $unmapped_bam
        samtools sort -n $unmapped_bam -o $bam_sorted
	bedtools bamtofastq -i $bam_sorted -fq "${name}.removed.fastq"
        gzip "*.fastq"
        """
        }
}}

else  {
        trimmed_reads_bowtie2
        .view()
        .into {ch_bowtie2_kraken; ch_bowtie2_kaiju}
}

/*
 * TAXONOMIC CLASSIFICATION
 *
 */

if (params.kraken2) {
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
            kraken2 --db ${params.krakenDB} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]}
            """    
        } else {
            """
            kraken2 --db ${params.krakenDB} --threads ${task.cpus} --output $out --report $kreport ${reads[0]} 
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
}} 


if (params.kaiju) {
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
    kaiju -x -v -t ${params.kaijuDB}/nodes.dmp -f ${params.kaijuDB}/kaiju_db_viruses.fmi -i ${reads[0]} -j ${reads[1]} -o $out
    kaiju2krona -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $krona_kaiju    
    """
    } else {
    """
    kaiju -x -v -t ${params.kaijuDB}/nodes.dmp -f ${params.kaijuDB}/kaiju_db_viruses.fmi -i ${reads[0]} -o $out
    kaiju2krona -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $krona_kaiju
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
}}
 
