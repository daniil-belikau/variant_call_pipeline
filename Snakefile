configfile: "/proj/workflow_config.yml"
ruleorder: annotate_vars > multiqc

def get_bwa_mem_input(wildcards):
	files = [
		f"data/input/ref_genomes/{config['genome']}",
		config["reads"][wildcards.sample]["mate1"],
		config["reads"][wildcards.sample]["mate2"]
	]
	return files
	

def get_reqs(wildcards):
	extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
	return [f"data/input/ref_genomes/{config['genome']}{ext}" for ext in extensions]


rule all:
	input:
		expand("data/results/vcf/{sample}/{sample}.filtered.annot.vcf", sample = config["sample"]),
		expand("data/results/fastqc/{sample}/{sample}{ext}_fastqc.html", sample = config["sample"], ext = config["ext"]),
		expand("data/results/alignments/{sample}/{sample}.aligned.sorted.markdup.flagstat.txt", sample = config["sample"]),
		expand("data/results/multiqc/{sample}/{sample}.done", sample = config["sample"])
		

rule fastqc:
	input: 
		lambda wildcards: [mate for mate in config["reads"][wildcards.sample].values()]
	output:
		multiext("data/results/fastqc/{sample}/{sample}{ext}_fastqc", ".html", ".zip") 
	params:
		out_dir = "data/results/fastqc/{sample}",
		log_dir = "logs/fastqc"
	log:
		"logs/fastqc/{sample}.{ext}.log"
	threads: 
		2
	shell:
		"mkdir -p {params} && "
		"fastqc -t {threads} -o {params.out_dir} {input} > {log} 2>&1"


rule bwa_index:
	input:
		"data/input/ref_genomes/{genome}"
	output:
		multiext("data/input/ref_genomes/{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa") 
	params:
		log_dir = "logs/bwa_index"
	log:
		"logs/bwa_index/{genome}.log"
	shell:
		"mkdir -p {params.log_dir} && "
		"bwa index {input} 2> {log}"	


rule bwa_mem:
	input:
		input = get_bwa_mem_input,			
		required = get_reqs
	output: 
		sam = "data/results/alignments/{sample}/{sample}.aligned.sam",
		bam = protected("data/results/alignments/{sample}/{sample}.aligned.bam"),
	params:
		read_group=r"@RG\tID:{sample}\tSM:{sample}",
		out_dir = "data/results/alignments/{sample}",
		log_dir = "logs/bwa_mem"
	log:
		"logs/bwa_mem/{sample}.log"
	threads: 
		4
	shell:
		"mkdir -p {params.log_dir} {params.out_dir} && "
		"bwa mem -t {threads} -R '{params.read_group}' {input.input} 1> {output.sam} 2> {log} && "
		"sambamba view -t {threads} -S -h -f bam {output.sam} > {output.bam}" 


rule bam_index:
	input:
		"data/results/alignments/{sample}/{sample}.aligned.bam"
	output: 
		"data/results/alignments/{sample}/{sample}.aligned.bam.bai"
	log:
		"logs/sambamba/{sample}.log"
	shell:
		"sambamba index {input} > {log}"


rule sort_bam:
	input:
		rules.bwa_mem.output.bam
	output: 
		bam = "data/results/alignments/{sample}/{sample}.aligned.sorted.bam",
		bai = "data/results/alignments/{sample}/{sample}.aligned.sorted.bam.bai"
	shell:
		"sambamba sort {input}"


rule markdup:
	input:
		rules.sort_bam.output.bam
	output: 
		"data/results/alignments/{sample}/{sample}.aligned.sorted.markdup.bam"
	log:
		"logs/sambamba/{sample}.log"
	params:
		log_dir = "logs/sambamba"
	threads: 
		4
	shell:
		"mkdir -p {params.log_dir} && "
		"sambamba markdup -t {threads} {input} {output} > {log}"


rule flagstat:
	input:
		rules.markdup.output
	output: 
		"data/results/alignments/{sample}/{sample}.aligned.sorted.markdup.flagstat.txt"
	log:
		"logs/sambamba/{sample}.log"
	threads: 
		4
	shell:
		"sambamba flagstat -t {threads} {input} 1> {output} 2>> {log}"


rule call_vars:
	input:
		bam = rules.markdup.output,
		ref = expand("data/input/ref_genomes/{genome}", genome=config["genome"])
	output: 
		"data/results/vcf/{sample}/{sample}.vcf"
	params:
		log_dir = "logs/freebayes"
	log:
		"logs/freebayes/{sample}.log"
	shell:
		"mkdir -p {params.log_dir} && "
		"freebayes -f {input.ref} --max-complex-gap 35 -j -p 1 {input.bam} 1> {output} 2> {log}"


rule filter_vars:
	input:
		rules.call_vars.output
	output: 
		"data/results/vcf/{sample}/{sample}.filtered.vcf"
	params:
		log_dir = "logs/bcftools"
	log:
		"logs/bcftools/{sample}.log"
	threads: 
		4
	shell:
		"mkdir -p {params.log_dir} && "
		"$BCF/bcftools view --threads {threads} -e 'QUAL < 20 || INFO/DP < 10 || INFO/DP > 100' {input} 1> {output} 2>> {log}"


rule annotate_vars:
	input:
		rules.filter_vars.output
	output: 
		vcf = "data/results/vcf/{sample}/{sample}.filtered.annot.vcf",
		tmp = touch("data/results/vcf/{sample}/{sample}.done")
	params:
		log_dir = "logs/snpeff",
		html_path = "data/results/vcf/{sample}/{sample}.annot.summary.html",
		csv_path = "data/results/vcf/{sample}/{sample}.annot.summary.csv",
		genome = ".".join(config["genome"].split(".")[:-1])
	log:
		"logs/snpeff/{sample}.log"
	shell:
		"mkdir -p {params.log_dir} && "
		"/opt/conda/pkgs/snpeff-5.0-hdfd78af_1/bin/snpEff ann "
		"-htmlStats {params.html_path} -csvStats {params.csv_path} "
		"{params.genome} {input} 1> {output.vcf} 2> {log}"


rule multiqc:
	input:
		rules.annotate_vars.output.tmp
	output:
		touch("data/results/multiqc/{sample}/{sample}.done")
	params:
		out_dir = "data/results/multiqc/{sample}",
		log_dir = "logs/multiqc"
	log:
		"logs/multiqc/{sample}.log"
	shell:
		"mkdir -p {params} && "
		"multiqc -o {params.out_dir} data/ > {log}"


