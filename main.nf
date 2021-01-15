#!/usr/bin/env nextflow
//================================================================================================//
// ToDo:                                                                                          //
//Update README how to run the pipeline based on the reads                                        //
//Do not build the index of contaminants and host sequence every times the pipeline is run        //
//use when clause                                                                                 //
//Bacteria and virus for DNA, only virus for RNA                                                  //

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
		--longReads   Input is long read data, e.g. nanopore [default:false]

		CONFIG FILE
		Nextflow will by default find configuration files in the following locations:
		"nextflow.config" in current directory
		"nextflow.config" in script directory
		"~/.nextflow/config"
		Specify config file with nextflow options:        
		-c        Path to additional config file
		-C        Path to config file. Other config files are ignored.

		OPTIONS
		--build       Build all databases [default:false]
		--porechop    Run Porechop trimming (only if long read data) [default:true]
		--fastp       Run Fastp trimming [default:true]
		--help        Display this help message and exit

		PARAMETERS
		Refer to the config file to see or change the default parameters used in the pipeline.

		To change a parameter for a single run, use:
		--<param> "value"

		--porechopParam    Settings for trimming with Porechop
		--fastpParamLong   Settings for trimming with Fastp for long read data
		--fastpParamShort  Settings for trimming with Fastp for short read data
		--krakenDB         Kraken2 database path
		--kaijuDB          Kaiju database path
		--hostReference    Host reference fasta path
		--contaminants     Fasta file with contaminant species
		--contaminantList  File with list of contaminant species
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
        params.fastpParam = params.fastpSingleEnd
}
else if (params.longReads){
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
        params.fastpParam = params.fastpPairedEnd
}



 
/*
 * PREPROCESSING BAM INPUT - Change name and convert bam to fastq
 *
 */

if (params.bam) {
	process renameBam {

		tag "$bam"

		//when: params.bam
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

	}
} 

if (params.bam) {
	process BamToFastq {

		tag "$bam"

		when: params.bam
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
	}
} 

def getCommonName(file,classifier ) {
	file.split('_')[0..1].join('_')+('_')+classifier
		/*file.tokenize('_')[0]*/
		/*file.substring(0, file.lastIndexOf("_"))*/       
}

/*def getCommonName (strs) {

  if (strs.size() == 0) {
  return "";
  }
  String result = strs[0];
  for (int i = 1; i < strs.size(); i++) {
  while (!strs[i].startsWith(result)) {
  result = result.substring(0, result.length() - 1);
  }
  }
  return result;
  }*/

if (params.bam) {
	ch_output_bamtofastq
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

if (params.longReads && !params.bam) {
process porechop {

    publishDir  "${params.outdir}/Porechop", mode: 'copy'

    when: params.porechop        

    input:
    set val(name), file(reads) from ch_input_porechop

    output:
    set val(name), file("*_chopped.fastq.gz") into (trimmed_reads_minimap2, trimmed_reads_fastqc)

    script:
    """
    porechop -i $reads -o ${name}_chopped.fastq.gz --format fastq.gz -t 4 $params.porechopParam
    """
    }
}

if (params.longReads && !params.porechop) {
	ch_input_porechop
		.into {trimmed_reads_minimap2}
}

if (!params.longReads) {
process fastp {

	tag "$name"
		publishDir  "${params.outdir}/fastp", mode: 'copy'

		when: params.fastp

		input:
                set val(name), file(reads) from ch_input_fastp

		output:
		set val(name), file("*.trimmed.fastq.gz") into (trimmed_reads_bowtie2, trimmed_reads_fastqc,trimmed_reads_bwa)	
		file("*.html")
		file ("*_fastp.json")

		script:
		if(params.singleEnd || params.bam ){
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

if (!params.fastp) {
	ch_input_fastp
		.into {trimmed_reads_bowtie2; trimmed_reads_fastqc}
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

/*process multiqc {


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
}*/

/*
 * BUILD DATABASES
 *
 */

/*process retrieve_contaminants {
  publishDir "${params.outdir}/Contaminants", mode: 'copy'
when: params.build
script:
"""
mkdir -p ${params.outdir}/Contaminants
RetrieveContaminants.py ${params.outdir}/Contaminants ${params.contaminantList} 
"""
} */

/*process index_contaminants {
input:
file cont_genomes from file(params.contaminants)
output:
file 'contaminants*' into ch_index_contaminants
script:
"""
bowtie2-build $cont_genomes contaminants
"""
} */

/*
if (!params.longReads) {
process index_host {

//when: params.build

input:
file host_genome from file(params.hostReference)

output:
file 'human*' into bowtie2_input

script:
"""
bowtie2-build $host_genome human
"""
}}*/
/*
// Get prebuilt host indices
if (!params.build) {
Channel
.fromPath("${params.hostIndex}/*.bt2")
.set {bowtie2_input}
}*/

if (!params.longReads) {
process index_host_bwa {
	
	input:
	file host_genome from file(params.hostReference)
	
	output:
	file '*' into bwa_input

	script:
	"""
	bwa index $host_genome 
"""
}}

//Exclude this part
/*ch_index_host
  .mix (ch_index_contaminants)
  .set {bowtie2_input} */
//

process krakenBuild {

	publishDir "${params.outdir}/databases", mode: 'copy'
	when: params.build

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

if (!params.longReads) {
process bwa {

	tag "$name"

	publishDir "${params.outdir}/bwa", mode: 'copy'

	input:
	set val(name), file(reads) from trimmed_reads_bwa
	file(index) from bwa_input.collect()
	file host_genome from file(params.hostReference)

	output:
	set val(name), file("*_unmapped.fastq") into (ch_bowtie2_kraken, ch_bowtie2_kaiju,ch_bowtie2_centrifuge)

	script:
	samfile = name + ".sam"
	bamfile = name + "_mapped_and_unmapped.bam"
	unmapped_bam = name + "_only_unmapped.bam"
	mapped_bam = name + "_mapped.bam"
	bam_sorted = name + "_sorted.bam"
	bam_mapped_sorted = name + "_mapped_sorted.bam"
	index_name = index.toString().tokenize(' ')[0].tokenize('.')[0]

	if (params.bam || params.singleEnd || params.longReads){
	"""
	bwa mem ${host_genome} ${reads} -t 4 > $samfile
	samtools view -Sb $samfile > $bamfile
	samtools view -b -f4 $bamfile > $unmapped_bam 
	samtools sort -n $unmapped_bam -o $bam_sorted
	bedtools bamtofastq -i $bam_sorted -fq "${name}_unmapped.fastq"
	"""
	}	
	else {
	"""
	bwa mem ${host_genome} ${reads[0]} ${reads[1]} -t 4 > $samfile
	samtools view -Sb $samfile > $bamfile
	samtools view -b -f4 $bamfile > $unmapped_bam
	samtools sort -n $unmapped_bam -o $bam_sorted
	bedtools bamtofastq -i $bam_sorted -fq "${name}_1.removed.fastq" -fq2 "${name}_2.removed.fastq"
	"""
	}
}} 

/*if (!params.longReads) {
process bowtie2 {

  tag "$name"
  publishDir "${params.outdir}/bowtie2", mode: 'copy'

input:
set val(name), file(reads) from trimmed_reads_bowtie2
file(index) from bowtie2_input.collect()

output:
set val(name), file("*_unmapped.fastq") into (ch_bowtie2_kraken, ch_bowtie2_kaiju,ch_bowtie2_centrifuge)

script:
samfile = name + ".sam"
bamfile = name + "_mapped_and_unmapped.bam"	
unmapped_bam = name + "_only_unmapped.bam"
mapped_bam = name + "_mapped.bam"
bam_sorted = name + "_sorted.bam"
bam_mapped_sorted = name + "_mapped_sorted.bam"

index_name = index.toString().tokenize(' ')[0].tokenize('.')[0]	

if (params.bam || params.singleEnd || params.longReads ){
"""
bowtie2 -x $index_name --local -U $reads > $samfile
samtools view -Sb $samfile > $bamfile
samtools view -b -f4 $bamfile > $unmapped_bam
samtools sort -n $unmapped_bam -o $bam_sorted
bedtools bamtofastq -i $bam_sorted -fq "${name}_unmapped.fastq"
"""
}else {
"""
bowtie2 -x $index_name --local -1 ${reads[0]} -2 ${reads[1]} > $samfile
samtools view -Sb $samfile > $bamfile
samtools view -b -f4 $bamfile > $unmapped_bam
samtools sort -n $unmapped_bam -o $bam_sorted
bedtools bamtofastq -i $bam_sorted -fq "${name}_1.removed.fastq" -fq2 "${name}_2.removed.fastq"
"""
}

}}*/

if (params.longReads) {
process minimap2 {
	
        cpus 4
        time "3h"
        tag "$name"      
        
	publishDir "${params.outdir}/minimap2", mode: 'copy'
	
        input:
        set val(name), file(reads) from trimmed_reads_minimap2
        file host_genome from file(params.hostReference)

	output:
	file "*.bam*"
        set val(name), file("${name}_unmapped.fastq") into (ch_bowtie2_kraken, ch_bowtie2_kaiju, ch_bowtie2_centrifuge)

	script:
        samfile = name + ".sam"
        bamfile = name + ".bam"
        unmapped_bam = name + "_unmapped.bam"
        unmapped_sorted_bam = name + "_unmapped_sorted.bam"

        """
        minimap2 -ax map-ont -t $task.cpus -a $host_genome $params.minimap2Param $reads > $samfile
        samtools view -bh $samfile > $bamfile
        samtools view -bh -f4 $bamfile > $unmapped_bam
        samtools sort -n $unmapped_bam -o $unmapped_sorted_bam
        bedtools bamtofastq -i $unmapped_sorted_bam -fq "${name}_unmapped.fastq" 
        """
}}

/*
 * TAXONOMIC CLASSIFICATION
 *
 */

process kraken2 {

	tag "$name"
        memory '80 GB'
	publishDir "${params.outdir}/kraken2", mode: 'copy'


	input:
	set val(name), file(reads) from ch_bowtie2_kraken


	output:
	set val("kraken2"), val(name), file("*.kraken.out.txt") into kraken2_krona
	set val("kraken2"), val(name), file("*.kraken.report") into kraken_report 


	script:
	out = name+".kraken.out.txt"
	kreport = name+".kraken.report"
	if (params.bam || params.singleEnd || params.longReads){
	"""
	kraken2 --db ${params.krakenDB} --threads ${task.cpus} --output $out --report $kreport ${reads[0]}
	"""    
	} else {
	"""
	kraken2 --db ${params.krakenDB} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]}
	"""
	}	
}

process kaiju {

	tag "$name"
        memory '130 GB'

	publishDir "${params.outdir}/kaiju", mode: 'copy'
        database = Channel.fromPath("${params.kaijuDB}/*.fmi") 

	input:
	set val(name), file(reads) from ch_bowtie2_kaiju
        path fmi from database 

	output:
	set val("kaiju"), val(name), file("*.kaiju.out") into kaiju_out
	set val("kaiju"), val(name), file("*.kaiju_summary.tsv") into kaiju_summary   
	set val(name), file("*.kaiju_names.out")

	script:
	out = name+".kaiju.out"	
	krona_kaiju = name + ".kaiju.out.krona"
	summary = name+".kaiju_summary.tsv"
	taxon = name+".kaiju_names.out"


	if (params.pairedEnd){
	"""
	kaiju -t ${params.kaijuDB}/nodes.dmp -f $fmi -i ${reads[0]} -j ${reads[1]} -o $out
	kaiju2krona -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $krona_kaiju    
	kaiju2table -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -r species -o $summary $out
	kaiju-addTaxonNames -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $taxon
	"""
	} else {
	"""
	kaiju -t ${params.kaijuDB}/nodes.dmp -f $fmi -i ${reads[0]} -o $out
	kaiju2krona -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $krona_kaiju
	kaiju2table -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -r species -o $summary $out
	kaiju-addTaxonNames -t ${params.kaijuDB}/nodes.dmp -n ${params.kaijuDB}/names.dmp -i $out -o $taxon
	"""
	}
} 

process centrifuge {

	tag "$name"
        memory '80 GB'
	cpus 4

	publishDir "${params.outdir}/centrifuge", mode: 'copy'

	input:
	set val(name), file(reads) from ch_bowtie2_centrifuge     

	output:
	set val("centrifuge"), val(name), file("*.kreport.txt") into centrifuge_out
	set val("centrifuge"), val(name), file("*.results.txt") into centrifuge_krona
	file("*.report.txt")

	script:
	kreport = name+".kreport.txt"
	report = name + ".report.txt"
	results = name + ".results.txt"    

	if (params.bam || params.singleEnd || params.longReads) {
	"""
	centrifuge -x ${params.centrifugeDB} -p ${task.cpus} --report-file $report -S $results -U ${reads}
	centrifuge-kreport -x  ${params.centrifugeDB} $results > $kreport
	#cat results.txt | cut -f 1,3 > results.krona
	"""
	}	

	else {
	"""
	centrifuge -x ${params.centrifugeDB} -p ${task.cpus} --report-file $report -S $results -1 ${reads[0]} -2 ${reads[1]}
	centrifuge-kreport -x  ${params.centrifugeDB} $results > $kreport
	#cat results.txt | cut -f 1,3 > results.krona
	"""
	}
}

kraken2_krona
.mix(kaiju_out, centrifuge_krona)
.set { ch_krona }
 

process krona_taxonomy {

	output:
	file("taxonomy/taxonomy.tab") into file_krona_taxonomy


	script:
	"""
	ktUpdateTaxonomy.sh taxonomy
	"""
} 

process krona {

	tag "${name}_${tool}"

	publishDir "${params.outdir}/${tool}", mode: 'copy'

	input:
	set val(tool), val(name), file(report) from ch_krona
	file("taxonomy/taxonomy.tab") from file_krona_taxonomy

	output:
	file("*") 

	script:

	"""
	ktImportTaxonomy -o ${name}_${tool}_krona.html -t 3 -s 4 ${report} -tax taxonomy
	"""
} 

/*kraken_report
  .concat(kaiju_summary,centrifuge_out)
  .into{ ch_controls;ch_all_samples }*/

if (params.bam) {
	kraken_report
		.concat(kaiju_summary,centrifuge_out)
		.map { classifier,file1,file2 -> tuple(getCommonName(file1,classifier),classifier,file1,file2) } 
		.transpose()
		.set {parser_control}

	parser_control
		.branch { 
		control: it[2] =~ /NK/
		sample: it[2] !=~ /NK/
	}
	.set { ch_parser_control }

	ch_parser_control.sample.set {parser_sample}
	ch_parser_control.control.set {parser_control}

	parser_control.cross(parser_sample)
	.map { input -> tuple(input[0][1], input[0][2], input[0][3], input[1][2], input[1][3])}.view()
	.set { ch_sample_control} 
}

else {
	kraken_report
	.concat(kaiju_summary,centrifuge_out)
	.set {ch_parser}
}
 
if (!params.bam) {
	process parse_classifier {

	publishDir "${params.outdir}/filtering", mode: 'copy'

	tag "${name}_${tool}"

	input:
	set val(tool), val(name), file(report) from ch_parser

	output:
	set val(tool), val(name), file("*filtered_taxonomic_hits.tsv") into ch_accessions 
	file("*.tsv")

	script:
	"""
	extractReferences.py ${name} ${tool} 100 ${report} 
	"""

	}	
}

if (params.bam) {

	process parse_classifier_control {

	publishDir "${params.outdir}/filtering", mode: 'copy'

	tag "${name}_${tool}"

	input:
	set val(tool),val(control),file(report_control),val(sample),file(report_sample) from ch_sample_control

	output:
        set val(tool), val(name), file("*filtered_taxonomic_hits.tsv") into ch_accessions
        file("*.tsv")

	script:
	"""
	extractReferences.py ${name} ${tool} 200,50,0.05,20 ${report_sample} ${report_control}
	"""
	}
} 
