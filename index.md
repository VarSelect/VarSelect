# Before you start
The following tools are required for running VarSelect, please read carefully and install all the following tools for full benefits. Please note that the current version of VarSelect supports only human genome reference hg19/GRCh37, and stay in touch for GRCH38 update. 

## Gemini

Gemini is a framework for human genetic variations analysis by the SQLite database engine, and can be downloaded and installed by following the instruction at https://gemini.readthedocs.io/en/latest/content/installation.html . Once installed, please follow the commands at https://gemini.readthedocs.io/en/latest/content/installation.html#updating-the-gemini-executables-and-annotations to download the data file of the additional annotations, such as the GERP and CADD scores.
 
## VEP
VEP is part of the Ensembl project and is a comprehensive variants annotation tool. To install VEP, please follow the instruction at http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html. Please note that the pre-computed cache files are required to speeding up the annotation process. Please follow the instruction step by step at http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html to download the VEP of GRCh37 version.

## VEP plugins
VEP supports plugin modules to incorporate annotations from the external datasets. VarSelect requires the dbNSFP plugin, which are available at https://github.com/Ensembl/VEP_plugins/blob/release/86/dbNSFP.pm, respectively. 

## dbNSFP
dbNSFP annotates the non-synonymous single nucleotide variations (nsSNVs). The data file of dbNSFP is available at https://sites.google.com/site/jpopgen/dbNSFP .

## ANNOVAR
ANNOVAR  is a variant annotation tool with high efficiency to a variety of annotation databases, and is available at http://annovar.openbioinformatics.org/en/latest/user-guide/download/. Please note that a license is required. Please follow the instruction at http://annovar.openbioinformatics.org/en/latest/user-guide/startup/ to install scripts into proper directories when all the required packages are downloaded. Databases will be automatically installed by VarSelect installation script.

## snpEff
snpEff annotates and predicts the impact of genetic variants, and is available at http://snpeff.sourceforge.net/download.html. After downloading the software, the pre-built snpEff database for human genome of GRCh37 is needed. Please download it with following command:

```
java –jar /path/to/your/snpEff.jar download –v GRCh37.75
```
## vcftools
vcftools are a set of tools for manipulating genetic variants in the VCF-formatted files, and are available at https://vcftools.github.io/index.html. Please follow the instruction to install vcftools at https://vcftools.github.io/examples.html.

## bcftools, bgzip, tabix
bcftools, bgzip and tabix are tools to compress, index and manipulate VCF files. bcftools are available at http://www.htslib.org/download/, and include the bgzip and the tabix tools in the release.

# Download VarSelect
The latest version of VarSelect package can be downloaded from https://github.com/VarSelect/VarSelect

# Install VarSelect
Please make sure that you have already downloaded all the required packages and resources for running VarSelect. When you are all set, please run the following command to decompress the VarSelect file.

```
tar zxvf VarSelect-latest.tar.gz
```

After extracting the package, run the VarSelect installation script

```
/path/to/your/VarSelect/install_VarSelect.pl
```

Add the VarSelect path to your system's $PATH settings

```
export PATH=/path/to/your/VarSelect/dir:$PATH
```


# Quick start

VarSelect script is executable only on command line. Please use -h flag for the basic usage information. 

```
varselect.pl –h
```

VarSelect annotates and analyzes sequence variants, and compares the results from different primary analyses. To start using VarSelect, please use "annotate" to process your vcf file(s) of interests.

```
varselect.pl annotate -v /path/to/vcf_files_list \
                      -p /path/to/ped/file \
                      -m workflow_mode
```

The annotation function combines the variants from different samples and comprehensively annotates all reported variants. The -v option specifies the file recording the links between samples and corresponding variant files. The links is specified by a comma separator, one file per line as the following format:

```
sample1,/path/to/vcf/file1
sample1,/path/to/vcf/file2
sample2,/path/to/vcf/file3 
```

The gender and phenotype information is also required, and is specify by -p option followed by a PED file containing the pedigree information. An example of the PED file is available at http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

Once annotated, VarSelect creates a SQLite database file by using the GEMINI framework. The database can be queried and filtered by using the GEMINI commands to extract the annotation information.

```
gemini query —header -q 'select * from variants limit 10' \
/path/to/varselect.db
```

Besides annotation, VarSelect also integrates copy number information to the variants. The option –c specifies the relationship between sample and copy number variation call from CNVkit. The relationship is specified in the file by comma separated pair, one file per line as following format:

```
sample1,/path/to/cns/file1
sample2,/path/to/cns/file2
```

VarSelect will annotate with tags: cnv_samples, cnv_log2. Quantitative change by log base 2 between the paired samples is computed and annotated on the cnv_fc_samples and cnv_foldchange_log2 tags. 

VarSelect also incorporate user's gene expression profile to annotate variants and supports the output files from DEXseq, sailfish or featureCounts. VarSelect annotates the variants on tags: xpr_samples, xpr_readcounts, xpr_fpkm, xpr_tpm.

```
sample1,readcount,/path/to/file/of/reads_counts
sample1,exoncount,/path/to/file/of/dexseq_generated_exoncounts
sample1,fpkm,/path/to/combined.fpkm
sample1,tpm,/path/to/combined.gene.sf.tpm
sample2,readcount,/path/to/file/of/reads_counts
sample2,exoncount,/path/to/file/of/dexseq_generated_exoncounts
sample2,fpkm,/path/to/combined.fpkm
sample2,tpm,/path/to/combined.gene.sf.tpm
```

Quantitative change between the paired sample is pairwise computed and annotated on tags: xpr_fc_samples, xpr_foldchange_readcount, xpr_foldchange_fpkm , xpr_foldchange_tpm

## Built-in analytic workflow for family and paired sample analysis

VarSelect provides two basic analytic workflows: 1) family analysis or 2) paired case-control analysis, by specifying the option “-m family” or “-m paired”, respectively.

For family analysis, VarSelect analyzes the genetic variants for five genetic models, including autosomal recessive (AR), compound heterozygosity (CH), de novo recessive (DNR), two-hit recessive (THR) and X-link recessive (XLR). The variants follow each of the model will be labelled on the tags: is_AR, is_CH, is_DR, is_TH and is_XL, respectively.

For paired case-control analysis, VarSelect labels the changes of nucleotides between the paired samples with 'loss of heterozygous (LOH)' or 'de novo somatic mutations'. The changes are labelled on tags: is_LOH and is_de_novo. 
The result of each analysis is recorded by a tag in_analysis_jobid, while the 'jobid' is the time stamp when the analysis is performed by VarSelect.

## Start from vcf files by multiple callers
Different variant callers often result in inconsistent variant calling reports while mostly are consistent with some inconsistent variants. VarSelect is engineered to deal with such situation by processing multiple VCF files by different variants callers in two ways: 1) unify all reported variants or 2) intersect variants reported by all callers. The option -k flag triggers this function, followed by option -u (union) or option -i (intersection), depending on the analytic purpose.

```
varselect.pl annotate -v /path/to/vcf/files/list \
                      -p /path/to/ped/file \
                      -m workflow_mode  \
                      -k -i

```

The format of the list of VCF files is different from single-caller mode since now a sample would have multiple VCF files. User must specify sample, variant caller, and the associated VCF file. An extra comma separated field with the name of variant caller is added to the end of each line as follows:

```
sample1,/path/to/vcf/file1,caller1
sample1,/path/to/vcf/file2,caller2
sample2,/path/to/vcf/file3,caller1
sample2,/path/to/vcf/file4,caller2
```

Please note that regardless filtering for the union or intersection, inconsistent calls among different callers are regarded ambiguity and are marked and removed from further analysis. The list of removed variants will be stored in the result directory.

## Re-analysis and update the analysis database 
Samples in a complex study design can be annotated together and analyzed according to various analytic purpose. The label of the phenotypic information (e.g. tumor and normal; affected and unaffected) of the samples can be rewritten and thereafter analyzed accordingly, according to the labels specified in the PED file. Samples begun with a hash character '#' in PED file are excluded from the downstream analysis. An example PED file is shown here, and the samples “uncle” and “aunt” are excluded from the downstream analysis

```
#family_id  sample_id   paternal_id maternal_id    sex     phenotype
family1     father       0           0         1       1
family1     mother      0           0         2       1
family1     daughter    father      mother      2       2
#family1    uncle       0           0         2       2
#family1    aunt        0           0         2       2
```

An initial full primary generates a VarSelect database, which is required for recording the results of the subsequent re-analysis. The -d option specifies the location of where the database file is kept. 

```
varselect.pl analysis -d /path/to/gemini/db \
                      -p /path/to/modified/ped/file \
                      -m workflow mode \
                      -k -u
```
After re-analysis, a new analysis directory will be created with new jobid. The results and logs will be stored in new directory. Filtered variants will be assigned new tag in_analysis_new_jobbid in the VarSelect database.

## Secondary analysis: compare results from any two primary analysis 

Analysis in full or in part of the annotation, single or multiple variant callers, and analytic workflow is termed ‘primary analysis’, while the results from different primary analysis can be compared for specific purposes and is termed ‘secondary analysis’. Comparison between any two primary analyses is specified by the options “–a” and “–b” as follows.

```
varselect.pl compare  -a /dir/to/analysis/of/datasetA \
                   -b /dir/to/analysis/of/datasetB \
                   -c [1-4] 
```

The secondary analysis includes four comparisons including: 1) the union of analysis A and B. 2) intersection of analysis A and B. 3) variants presented in the analysis A but not in the analysis B. 4) variants present in the analysis B but not in the analysis A. 

Results of new secondary analysis will be stored in a new directory. The filtered variants will be assigned a tag in_analysis with the Job ID.

# Family analysis
The family analysis workflow includes three steps: 1) analyzing the core family members including father, mother and sick child either male or female, the parents could be both healthy or either one is sick; 2) excluding the variants appeared in the other healthy members in the pedigree; 3) keeping the variants appeared in the other sick members in the pedigree. The inheritance models include: autosomal recessive (AR), compound heterozygosity (CH), de novo recessive (DNR), two-hits recessive (THR), and X chromosome-linked recessive (XLR). In humans, the pseudoautosomal region (PAR) on the X chromosome is diploid while the rest of the nucleotide sequence is haploid in male. Variants found on the X chromosome that fall in the CH, DNR, and THR models are reported on each of the category, and excluded from the XLR model. 

## Autosomal recessive inheritance 
In the AR model, both the pathogenic alleles came from each of the parents while the genotype of both of the parents are heterogeneous and includes exactly one pathogenic allele. When one of the parents is also sick, the sick parent shall carry homozygous variants, while the other parent shall be heterozygous and carries only one pathogenic variant.

## Compound heterozygosity inheritance
For discovering the pathogenic alleles of compound heterozygosity (CH), two candidate alleles shall locate in the same gene but on different positions yet chromosomes from each of the patients. The genotypes of the CH variants in both parents have to be heterozygous and includes exactly one pathogenic allele. If one of the parents is sick, the candidate gene shall be also compound heterozygous as the sick child. For the genes on the X chromosome, only the PAR region are scanned to the sick boys. Because of healthy father, only the PAR region are scanned for sick girls. 

## De novo recessive mutation 
In the de novo recessive (DNS) model, the variants of the sick child shall be homozygous and have to be different to that of the parents. For that on the X chromosome, the variants on the non-PAR region in the male and female patients shall be hemizygous and homozygous, respectively, while the variants on the PAR shall be homozygous for both the male and female patients. When one of the parents is sick, the homozygous variant in the sick child shall be different to that in the sick parent. 

## Two-hit recessive mutation
In the two-hit model, the sick child inherits one pathogenic allele from one of the parents while the other allele is somatically acquired. When healthy, both parents are homozygous wild-type or heterozygous carriers. For the variants on the X chromosome, the entire chromosome is scanned for sick girls, and is only the PAR region for sick boys. When one of the parents is sick, he/she shall carry two homozygous pathogenic alleles, while the other is either homozygous wild-type or not carrying the pathogenic allele. 

## X-linked recessive inheritance
For detecting X-linked recessive (XLR) inheritance, the variants on the non-PAR region in the male and female patients are hemizygous and homozygous, respectively. The variants on the PAR shall be homozygous for both the male and female patients. Typically, the sick boy shall inherit the pathogenic allele on the non-PAR region from his healthy or sick mother who carries heterozygous or homozygous variants, respectively. When the father is healthy, the girl received healthy hemizygous sequences on the non-PAR region from her father. Her disease is likely not related to the variants on non-PAR region, and the non-PAR region is therefore skipped unless the father is also sick. Depending on which parent is sick, either one or two pathogenic alleles shall be seen on the non-PAR region of the father or mother, respectively, while are always two pathogenic alleles on the PAR region. 

# Paired case/control analysis
In the case/control workflow, the variants of the matched samples are classified as either ‘loss of heterozygosity (LOH)’ or ‘de novo mutations’. All the different variants in the vcf files of both the cases and control, regardless single nucleotide, multi-nucleotide replacement, or small INDEL, are compared for LOH otherwise de novo mutations. VarSelect supports the comparison between one or more than one cases versus exactly one control sample. Comparison of using other sample as a control can be performed by setting up new secondary analysis. 

## Loss of heterozygosity (LOH) events
Variants that are heterozygous in the vcf file of control sample lose it heterozygosity while are homozygous for one of the two heterozygous alleles is so called LOH. When multiple cases are analyzed, each case is individually compared against the same control sample and annotated on the “is_LOH” tag. 

## De novo mutation
De novo mutation can occur on the alleles of heterozygosity and homozygosity, as long as the alternative nucleotide in the case isn’t present in the control sample. When multiple cases are analyzed, each case is individually compared against the same control sample and annotated on the “is_de novo” tag.


# Description of VarSelect scripts
VarSelect is composed of many individual scripts. Bellows are the description of each VarSelect script. 

* varselect.pl is the main script of VarSelect. It provides three commands included in VarSelect: annotate (initial annotation), analysis (primary and re-analysis) and compare (secondary analysis).
  * Command “annotate” triggers the script vs_annotate.pl to annotate VCF files from scratch and the script vs_analysis.pl to analyze variants through workflow of choice. There are three required options: -v sample-vcf file list, -p PED file, -m workflow mode.
  * Command “analysis” triggers the script vs_analysis.pl. There are three required options: -d gemini db file, -p PED file, -m workflow mode. 
  * Command “compare” triggers the script vs_compare.pl to compare results between two primary analysis. There are three required options: “-a” and “-b” to specify the path of the analysis A and B. “-c” specify the method of comparison (1. union, 2. intersection, 3. A only, and 4. B only)

* vs_annotate.pl is called by varselect.pl for three functions. 
  * Firstly, it prepares VCF file for annotation. VCF files from same sample will be joined together by vcf-concat of VCFtools. VCF files of different samples are then merged into a VCF file by vcf-merge in the VCFtools. The genotypes of sex chromosome on merged VCF file are curated by vcf-fix-ploidy in the VCFtools.
  * Secondly, the script triggers VEP, snpEff and ANNOVAR for annotation
  * Thirdly, the script triggers gemini framework to generate a sqlite db for downstream analysis.

* vs_analysis.pl is triggered by varselect.pl to classify and label variants by the genotypes and affected status.
  * If -k option supplied, multi-caller mode is turned on with union (option -u) or intersection (option -i) set of variants from different callers.
  * If -c option supplied, the script calls up cnvkit_parse.pl to annotate copy number information.
  * If -x option supplied, the script calls xprprofile_parse.pl to annotate gene expression.
  * The script supports two basic analytical workflows by specifying -m option. For the paired case-control workflow, it triggers loh_detector.pl and denovo_detector.pl scripts to classify variants of loss of heterozygosity or of de novo changes. For family workflow, the script triggers a total of five scripts including Autosomal-recessive.py, Compound-het.py, Denovo-recessive.py, Two-hits.py and X-linked.py to classify the variants into the related genetic models. All variants selected by the workflow will be annotated on the tag is_analysis_JobID in the variants table of the database.

* vs_compare.pl is triggered by varselect.pl to extract the intersection/union/subtraction of variants between two primary analyses.
* run_vep.pl is triggered by vs_annotate.pl to annotate the VCF file through Ensembl VEP. The script also enables the dbNSFP plugin to annotate non-synonymous variants.
* run_snpeff.pl is triggered by vs_annotate.pl to annotate the VCF file through snpEff. 
* run_annovar.pl is triggered by vs_annotate.pl to annotate the VCF file through ANNOVAR.
* GO_parse.pl is triggered by vs_annotate.pl to annotate variants on Gene Ontology (GO) terms. The GO terms are mapped by the corresponding ensembl id of transcripts where variants reside.
* pathway_parse.pl is triggered by vs_annotate.pl to annotate variants on KEGG pathways. Correlated pathways will be extracted by the built-in function in the GEMINI framework “gemini pathways” and updated to the variants table.
* cnvkit_parser.pl is triggered by vs_analysis.pl to annotate variants with the information of copy number variation. It extracts the log2 values where each variant reside and calculates fold-change between the paired samples.
* xprprofile_parser.pl is triggered by vs_analysis.pl to annotate variants with gene express profiles. The script supports the annotation of four types of gene expression profile: gene-based read counts, exon-based read counts, gene-based tpm (Transcripts Per Kilobase Million), gene-based fpkm (Fragments Per Kilobase Million). Please specify the -x option followed by the file of gene expression profiles.
* loh_detector.pl is triggered by vs_analysis.pl for the paired case/control workflow. The script checks variant genotypes between affected and unaffected samples and determine whether the variant is LOH. Variants classified as LOH are labeled on the tag “is_LOH”. Further description of the paired case/control workflow is described in the next section. 
* denovo_detector.pl is triggered by vs_analysis.pl for the paired case/control workflow. The script compares the variants in the control and in the case sample for de novo somatic mutations, and labels on the tag “is_de_novo” Further description of the paired case/control workflow is described in the next section.
* Autosomal-recessive.py is triggered by vs_analysis.pl for the family workflow. The script compares the variants between family samples to classify whether the variant is caused by autosomal recessive inheritance and label on the tag “is_AR”. Further description of the family workflow is described in the next section.
* Compound-het.py is triggered by vs_analysis.pl for the family workflow. The script compares the variants between family samples to classify whether the variant is caused by compound heterozygosity inheritance and label on the tag “is_CH”. Further description of the family workflow is described in the next section.
* Denovo-recessive.py is triggered by vs_analysis.pl for the family workflow. The script compares the variants between family samples to classify whether the variant is caused by De novo recessive mutation and label on the tag “is_DR”. Further description of the family workflow is described in the next section.
* Two-hits.py is triggered by vs_analysis.pl for the family workflow. The script compares the variants between family samples to classify whether the variant is caused by two-hit recessive mutation and label on the tag “is_TH”. Further description of the family workflow is described in the next section.
* X-linked.py is triggered by vs_analysis.pl for the family workflow. The script compares the variants between family samples to classify whether the variant is caused by X-linked recessive inheritance and label on the tag “is_XL”. Further description of the family workflow is described in the next section


# Examples
The following examples present some common scenario to start with.


## Examples 1 - samples of a family study

In this example, the variants of the chromosome 22 of a family trio, including NA12878, NA12891, NA12892 are demonstrated as an example for family-based analysis.

All required files are stored in directory "varselect/examples/example1/" :

* NA12878_chr22.vcf.gz 
* NA12891_chr22.vcf.gz 
* NA12892_chr22.vcf.gz 
* example1.txt : A comma separated file that describe the relationship between sample and corresponding VCF file.  
* example1.ped : A tab separated file that describes the pedigree information of this family, including the gender and the disease status. In this example, we assume NA12878 is a sick child while NA12891 and NA12892 are unaffected parents. 

For detail of ped format, please check the page (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) 

To start the analysis, please type in following command:

```
varselect.pl annotate -v example1.txt -p example1.ped -m family
```

VarSelect records the timestamp as job id when analysis is submitted (ex: 20170105162548). Directory “VarSelectAnalysisResult_20170105162548/” is created when VarSelect begins. Log files and intermediate results are stored in this directory.

When the primary analysis is finished, example1_varselect.db is generated. User can type in the gemini commands to query the database for further analysis. 


```
gemini db_info example1_varselect.db |grep '^variants'
```
Through the family workflow, five specific columns come with your job id: is_AR_20170105162548 , is_CH_20170105162548, is_DR_20170105162548, is_TH_20170105162548 and is_XL_20170105162548 corresponding the following genetic models: autosomal recessive, compound heterozygosity, de novo recessive, two-hit recessive and X-linked recessive, respectively.


The following command extracts the variants following compound heterozygosity inheritance with the information of chromosome, position, ref allele, alt allele, gene and genotype of all sample 


```
gemini query --header
       –q ‘select chrom,start,ref,alt,gene,gts from variants
           where is_CH_20170105162548 = 1’
       example1_varselect.db
```

Column “in_analysis_20170105162548” records variants qualifying the workflow analysis, and, in this case, is a union of the five genetic models. Variants falling in any of the five genetic models can be filtered out by following command:


```
gemini query --header
       –q ‘select chrom,start,ref,alt,is_AR_20170105162548 ,
           is_CH_20170105162548, is_DR_20170105162548, is_TH_20170105162548,
           is_XL_20170105162548 from variants where in_analysis_20170105162548 = 1'
       example1_varselect.db   
```

## Examples 2 - paired case/control samples
In this case, we use the variants of chromosome 22 of three samples, including blood , primary tumor ccRCC (clear cell Renal Cell Carcinoma), and a metastasis lung cancer available at https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063388.
All files are stored in subdirectory examples/example2/

* blood-gatk-chr22.vcf.gz 
* ccRCC-gatk-chr22.vcf.gz
* meta-lung-gatk-chr22.vcf.gz
* example2.txt: A comma separated file that describes the links between sample and corresponding VCF file.
* example2.ped: A tab separated file that describes the relationship between each sample, the gender and the disease status. In this case, blood is assigned as an unaffected (control) sample, while ccRCC and meta-lung are assigned as affected samples.

To start the analysis, please type in the following command. 

```
varselect.pl annotate -v example2.txt -p example2.ped -m paired
```

When the analysis is finished, a file “example2_varselect.db” is created at the same directory.
To check which column is in the variations table first.

```
gemini db_info example2_varselect.db |grep '^variants'
```

There are several columns includes is_LOH (loss of heterozygosity ), LOH_samples, is_denovo and denovo_gt denovo_samples specific in paired mode. VarSelect compares the genotype difference between control (unaffected) and case (affected) samples, and assign LOH or de_novo tag.

The following commands filters out the de_novo somatic variants.

```
gemini query --header 
             –q ‘select chrom, start, ref, alt, gts, denovo_samples_20170106175124 
                 from variants where is_denovo_20170106175124=1 ’  
             example2_varselect.db
```

## Examples 3 - Re-analysis and comparison (secondary analysis)
Multiple primary analysis (re-analysis) can be performed for various study purposes. Comparison between multiple primary analysis provides the flexibility of hierarchical comparison, namely secondary analysis. User can repeat the analysis by changing the label recorded in the ped file. For example, in the case of ccRCC, you can filter the common de novo mutations presented in both ccRCC and the lung meta samples by performing “blood vs. ccRCC” and “blood vs. meta-lung’ analyses. The third sample inn the ped file marked with '#' will be excluded in this analysis.

Firstly, replicate the ped file in the sample 2 as follows.

```
cp example2.ped example3_mark1.ped
```

Edit the example3_mark1.ped file and exclude meta-lung sample by ‘#” sign.

```
example3        blood   0       0       1       1
example3        ccRCC   0       0       1       2
#example3        meta-lung       0       0       1       2
```

Perform re-analysis by using this new ped file. User can use “varselect.pl analysis” to perform analysis and update the existing db file as follows.

```
varselect.pl analysis -d example2_varselect.db -p example3_mark1.ped -m paired
```

A new job id is generated at the time of beginning the re-analysis (ex: 20170109183415), and a new directory “VarSelectAnalysisResult_20170109183415” is created to store the log files and intermediate result files of this new analysis job.
When the re-analysis is finished, new columns tagged with new job id are added into the variants table. The command is as follows.

```
gemini db_info example2_varselect.db | grep '20170109183415'
```

User can use following command to filter out the variants from de_novo mutation present in the ccRCC.

```
gemini query --header -q 'select chrom,start,ref,alt,gts from variants where is_de novo_20170109183415 = 1' example2_varselect.db
```

Replicate the example3_mark1.ped file and save as the example3_mark2.ped. Mark a ‘#’ sign on the sample meta-lung . 

```
example3        blood   0       0       1       1
#example3       ccRCC   0       0       1       2
example3        meta-lung       0       0       1       2

```

Perform re-analysis again with the newly edited ped file to identify the de novo variants present in the metastasis lung tumor.

```
varselect.pl analysis -d example2_varselect.db -p example3_mark2.ped -m paired
```
New directory “./VarSelectAnalysisResult_20170110174312” is created with corresponding tags updated in the database.
You can extract the de novo mutation presented in the metastasis tumor by typing in the following command.

```
gemini query --header 
    -q 'select chrom,start,ref,alt,gts from variants where 
        is_denovo_20170110174312 = 1' 
        example2_varselect.db
```

To identify the de novo mutations that present in both ccRCC and lung meta samples, we compare the results from the two (primary) analysis (jobs) by intersecting the results. Please note that only analyses stored in the same database (db) file can be compared (namely secondary analysis). In this example, the common de novo mutations are extracted by typing in the following command.

```
varselect.pl compare  
    -a VarSelectAnalysisResult_20170109183415  
    -b VarSelectAnalysisResult_20170110174312 
    -c 2
```

A new job id is assigned (say, 20170111082651) and the database (db) will be updated with the results of the secondary analysis. The following command lists the results stored in the database (db). 

```
gemini query --header 
    -q 'select chrom, start, ref, alt, gts from variants
    where in_analysis_20170111082651 = 1'
    example2_varselect.db
```

If you want to extract novel mutation only present in the metastasis tumor, change the option –c to ‘4’ (unique to the analysis in the ‘-b’ option). 

```
varselect.pl compare  
    -a VarSelectAnalysisResult_20170109183415 
    -b VarSelectAnalysisResult_20170110174312 
    -c 4
```

## Examples 4 – same sample with results from multiple variants callers
For various reasons, union and/or intersection of variant calls from different callers could be of interests to some researchers. This example demonstrates the integrative results of two popular variant callers: GATK-HaplotypeCaller and freebayes.

ll files are stored in directory examples/example4/:

* NA12878-gatk-chr22.vcf.gz
* NA12891-gatk-chr22.vcf.gz
* NA12892-gatk-chr22.vcf.gz
* NA12892-freebayes-chr22.vcf.gz
* NA12878-freebayes-chr22.vcf.gz
* NA12891-freebayes-chr22.vcf.gz
* example4.txt: A comma separated file that describes the links between samples and corresponding VCF files. For multi-caller function, an extra column is required to describing the corresponding variant caller.

```
    NA12878,NA12878-gatk-chr22.vcf.gz,gatk-haplotype,
    NA12891,NA12891-gatk-chr22.vcf.gz,gatk-haplotype
    NA12892,NA12892-gatk-chr22.vcf.gz,gatk-haplotype
    NA12878,NA12878-freebayes-chr22.vcf.g,freebayes z
    NA12891,NA12891-freebayes-chr22.vcf.gz,freebayes
    NA12892,NA12892-freebayes-chr22.vcf.gz,freebayes
```

* example4.ped: A tab separated file that describes the relationship between each sample, the gender and the disease status.

The ‘-k’ option triggers the multi-caller function, followed by ‘-u’ option for preparing the union of all variants from the two variant callers. 

```
varselect.pl annotate 
    -v example4.txt 
    -p example4.ped 
    -m family 
    -k -u 
```

A new job id (ex:20170111162057) is tagged on a new directory “VarSelectAnalysisResult_20170111162057”. Please note that the variants with inconsistent calls by different callers will be removed for the downstream analysis.

The full list of the removed variants could be found in the file “multicaller_union_inconsistant_20170111162057.txt”.

The following command lists the variants in this analysis.

```
gemini query --header -q 'select chrom, start, ref, alt, gts 
    from variants where in_analysis_20170111162057 = 1' 
    example4_varselect.db
```

To analyze on the intersection of variants identified by both variant callers, replace the option ‘-u’ to ‘-i’ in the original command as following.

```
    varselect.pl analysis 
        -d example4_varselect.db 
        -p example4.ped 
        -m family 
        -k -i 
```

Please note that this two options ‘-u’ and ‘-i’ are mutually exclusive. Inconsistent calls by different callers are removed for the downstream analysis and listed in file “multicaller_intersect_inconsistant_20170111190712.txt” in the same directory of analysis output.

## Known issue
* To make VCF files be merged smoothly, decomposed multiple alt-allele variants with same position or wrong number of format value are excluded from analysis. 


