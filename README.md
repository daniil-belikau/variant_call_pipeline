# Variant Calling Containerized Workflow

## Tool Description
This is a variant calling workflow defined using ```snakemake``` and packaged inside a Docker image. The steps defined as part of the workflow are:
* ```fastqc```:
  * QC of sequencing reads
* ```bwa```:
  * Reference index construction
  * Paired end read alignment
* ```sambamba```:
  * SAM --> BAM conversion
  * BAM indexing
  * BAM sorting 
  * Marking technical duplicate reads
  * Summarizing flag statistics
* ```freebayes```:
  * Calling variants
* ```bcftools```:
  * Filtering called variants
* ```snpeff```:
  * Annotating called variants 
* ```multiqc```:
  * QC of multiple pipeline steps across samples

The entire pipeline can be invoked simultaneously, or each step could be called individually. 

## Using The Pipeline
The pipeline can be run by executing the Docker container directly. Here is the simple setup:

### Directory structure
The container is going to interact with the host file system using a volume, so we will need to provide our input files in structure that the pipeline expects:
* proj/
  * Snakemake
  * workflow_config.yml
  * logs/
  * data/
    * input
      * reads
      * ref_genomes
    * results

### Setup
To set up we simply need to get the Docker image and we're ready to run the pipeline.
```bash
docker pull daniilbelikau/var_call:latest
```

### Running the pipeline
A minimal sample dataset is included in this repository. To run the pipeline as a container we use the following:
```bash
docker run -v {host_proj_dir_path}:/proj var_call
```
This invokes the ```snakemake``` workflow specified in ```workflow_config.yml``` by automatically running this inside the container:
```bash
snakemake -c 4 -p
```
You can make the container run something else by providing your arguments like this:
```bash
docker run -v {host_proj_dir_path}:/proj var_call [your_commands_and_args]
```
### Producing workflow DAG
If you want to visualize the workflow that ```snakemake``` is about to perform, run:
```bash
snakemake --forceall --dag | dot -Tpdf > workflow_dag.pdf
```
## Unit Tests
You can generate pipeline unit tests by providing a sample dataset and running:
```bash
snakemake -c {num_threads} --forceall --notemp
snakemake --generate-unit-tests
```
SnakeMake will generate a directory called ```.tests``` and populate it with python functions testing each step of the workflow. It also stores the sample input and output data for each step. These data and functions can be incorporated into a larger testing framework. Note that when running the tests, the default behaviour is to compare the current and expected output byte by byte.
A suite of example unit tests for the minimal dataset can be found in this repository. 

## Larger Example
The minimal example is convenient for testing and developing, but we can try to process a more demanding input. Let's try to characterize the genetic variation of the SARS-CoV-2 Delta strain compared to the MN908947.3 reference genome. We need to procure some sequencing data:
```bash
mkdir -p data/input/reads/SRS10100888
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR15829420/BCSIR-NILMRC757_S19_L001_R1_001.fastq.gz.1 &
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR15829420/BCSIR-NILMRC757_S19_L001_R2_001.fastq.gz.1
mv BCSIR-NILMRC757_S19_L001_R1_001.fastq.gz.1 data/input/reads/SRS10100888/SRS10100888_R1_001.fastq.gz
mv BCSIR-NILMRC757_S19_L001_R2_001.fastq.gz.1 data/input/reads/SRS10100888/SRS10100888_R2_001.fastq.gz
```
Now we have to define a new pipeline configuration file ```delta_config.yml```, and all that's left is to run the container:
```bash
docker run -v {host_proj_dir_path}:/proj var_call snakemake --configfile delta_config.yml -c 4 -p
```
Now we can examine our findings and hypothesize about the effects of S furin cleavage site mutations on the virulence of SARS-CoV-2.

## Thoughts On Heterozygosity Ratio
In an organism with a non-haploid genome, the heterozygosity ratio can serve as a proxy for evaluating the extent of unexpected genotype being found in the sample. For instance, in humans most populations have a heterozygosity ratio between 1.5 and 2. However, in a sample with genetic mosaicism loci are less likely to be homozygous for the alternative allele, thus you would expect a higher ratio with more extensive mosaicism. Similar logic can be used in looking at contamination. 
