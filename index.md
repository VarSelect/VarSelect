# VarSelect Manual 20170413

## Before you start

The following tools are required for running VarSelect, please follow the instruction and install all the following tools for full benefits. Please note that the current version of VarSelect supports only the human genome reference hg19/GRCh37. The support for GRCh38 will be available shortly. 


**Gemini**

Gemini[1] is an analytic framework for human genetic variations by taking advantage of the SQLite database engine. It can be downloaded and installed by following the instruction at [`https://gemini.readthedocs.io/en/latest/content/installation.html`](https://gemini.readthedocs.io/en/latest/content/installation.html). Once installed, please follow the commands to download the data files of the additional annotations, such as the GERP[2] and CADD[3] scores at [`https://gemini.readthedocs.io/en/latest/content/installation.html#updating-the-gemini-executables-and-annotations`](https://gemini.readthedocs.io/en/latest/content/installation.html#updating-the-gemini-executables-and-annotations).


**VEP**

 VEP[4] is a comprehensive variants annotation tool, and is part of the Ensembl project. To install VEP, please follow the instruction at [`http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html`](http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html). Please note that the pre-built cache files are required for speeding up the annotation process. Please follow the instruction step by step at [`http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html`](http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html) to download VEP of GRCh37 version. 


**VEP plugins**

 VEP supports plugin modules to incorporate annotations from the external datasets. VarSelect requires the dbNSFP plugin, which are available at [`https://github.com/Ensembl/VEP_plugins/blob/release/86/dbNSFP.pm`](https://github.com/Ensembl/VEP_plugins/blob/release/86/dbNSFP.pm), respectively. 


**dbNSFP**

 dbNSFP[5, 6] annotates genome-wide non-synonymous single nucleotide variations (nsSNVs). The data file of dbNSFP is available at [`https://sites.google.com/site/jpopgen/dbNSFP`](https://sites.google.com/site/jpopgen/dbNSFP). 


**ANNOVAR**

 ANNOVAR[7] is a variant annotation tool with high efficiency to a variety of annotation databases, and is available at [`http://annovar.openbioinformatics.org/en/latest/user-guide/download/`](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). Please note that a license is required. Please follow the instruction at [`http://annovar.openbioinformatics.org/en/latest/user-guide/startup/`](http://annovar.openbioinformatics.org/en/latest/user-guide/startup/) to install scripts into proper directories when all the required packages are downloaded. Databases will be automatically installed by the VarSelect installation script.


**snpEff**

 snpEff[8] annotates and predicts the impact of genetic variants, and is available at [`http://snpeff.sourceforge.net/download.html`](http://snpeff.sourceforge.net/download.html). After downloading, please also download the pre-built snpEff annotation database. Please note that the current version of VarSelect is based on human genome reference hg19/GRCh37. Please download it with following command:
java –jar /path/to/your/snpEff.jar download –v GRCh37.75
 

**vcftools**

[`vcftools`](https://vcftools.github.io/)[9] is a set of tools for manipulating genetic variants in the VCF-formatted files, and is available at [`https://vcftools.github.io/index.html`](https://vcftools.github.io/index.html). Please follow the instruction to install vcftools at [`https://vcftools.github.io/examples.html`](https://vcftools.github.io/examples.html).


**bcftools, bgzip, tabix**

 bcftools, bgzip and tabix are tools to compress, index and manipulate VCF files. bcftools is available at [`http://www.htslib.org/download/`](http://www.htslib.org/download/), and includes the bgzip and the tabix tools in the software release.


## **Download VarSelect**

The latest version of VarSelect is available at [`https://github.com/VarSelect/VarSelect`](https://github.com/VarSelect/VarSelect).


## **Install VarSelect**

Please make sure that you have downloaded all the required packages and resources needed by VarSelect. When you are all set, please run the following command to decompress the VarSelect files.

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

## **Quick Start**

VarSelect script is executable only on command line. Please use -h flag for the basic usage information. 

```
varselect.pl -h
```

VarSelect annotates and analyzes sequence variants, and compares the results from different primary analyses. To start using VarSelect, please use "annotate" to process your vcf file(s) of interests.

```
varselect.pl annotate -v /path/to/vcf_files_list
                      -p /path/to/ped/file 
                      -m workflow_mode
```

The annotation function combines the variants from different samples and annotates all reported variants of 406 features. The -v option specifies the file describing the links between samples and corresponding variant files. The link is specified by a comma separator, one file per line as the following format:

```
**sample1,/path/to/vcf/file1
 sample1,/path/to/vcf/file2
 sample2,/path/to/vcf/file3 **
```

The gender and phenotype information is also required, and is specified by the -p option followed by a PED file containing the pedigree information. An example of the PED file is available at http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped .
 
When the variants are annotated, VarSelect creates a SQLite database file by using the `GEMINI framework.``[1]` The database can be queried and filtered by using the GEMINI and/or the SQL commands to extract the annotation information.

```
**gemini query --header -q 'select * from variants limit 10' /path/to/varselect.db**
```

VarSelect also annotates copy number variation (CNV) produced by CNVkit[10] to the variants. The option -c specifies the file describing the list of CNV calls as following format:

```
**sample1,/path/to/cns/file1
 sample2,/path/to/cns/file2**
```

VarSelect annotates CNV to the database with tags ‘cnv_samples’ and ‘cnv_log2’, and computes quantitative change by log base 2 between the paired samples. The results are tagged ‘cnv_fc_samples’ and ‘cnv_foldchange_log2’. 
 
VarSelect also integrates the gene expression profiles, if provided, to annotate sequence variants. VarSelect currently supports the results from DEXseq,[11, 12] sailfish[13] and featureCounts.[14] VarSelect annotates the variants with the following tags ‘xpr_samples’, ‘xpr_readcounts’, ‘xpr_fpkm’, and ‘xpr_tpm’. The option -x specifies the file describing the list of gene expression profiles as following format:

```
**sample1,readcount,/path/to/file/of/reads_counts
 sample1,exoncount,/path/to/file/of/dexseq_generated_exoncounts
 sample1,fpkm,/path/to/combined.fpkm
 sample1,tpm,/path/to/combined.gene.sf.tpm
 sample2,readcount,/path/to/file/of/reads_counts
 sample2,exoncount,/path/to/file/of/dexseq_generated_exoncounts
 sample2,fpkm,/path/to/combined.fpkm
 sample2,tpm,/path/to/combined.gene.sf.tpm**
```

 
Quantitative change between the paired samples is pairwise computed and annotated with the tags ‘xpr_fc_samples’, ‘xpr_foldchange_readcount’, ‘xpr_foldchange_fpkm’, and ‘xpr_foldchange_tpm’.


## ****Built-in analytic workflow for family and paired sample analysis****

In the primary analysis, VarSelect provides two common analytic workflows: 1) family and 2) paired case-control analysis, by specifying the option “-m family” or “-m paired”, respectively.

For family analysis, VarSelect analyzes the genetic variants for five genetic models, including autosomal recessive (AR), compound heterozygosity (CH), de novo recessive (DNR), two-hit recessive (THR) and X-link recessive (XLR). The variants satisfied the criteria of each of the models will be labelled with the tags: is_AR, is_CH, is_DR, is_TH and is_XL, respectively.
 
For paired case-control analysis, VarSelect labels the changes of nucleotides between the paired samples with either 'loss of heterozygous (LOH)' or 'de novo somatic mutations'. The changes are labelled with the tags: is_LOH and is_denovo. The result of each analysis is recorded by a tag ‘in_analysis_{jobid}’, while the 'jobid' is the time when the analysis begins.


## ****Start from vcf files by multiple callers ****

Different variant callers deliver inconsistent variant calling reports while the majority are consistent.[15] VarSelect deals with such situation by processing multiple VCF files by different variants callers in two ways: 1) unify all reported variants or 2) intersect variants reported by all callers. The option -k flag triggers this function, followed by the options -u (union) or -i (intersection), depending on the analytic purpose.

```
**varselect.pl annotate* *-v /path/to/vcf/files/list 
* **                     *-p /path/to/ped/file 
* **                     *-m workflow_mode 
* **                     *-k * *-i**
```

The list specifying VCF files is slightly different from the single-caller mode since a sample would now have multiple VCF files. User must specify sample, variant callers, and the associated VCF files. An extra comma separated field with the name of variant caller is added to the end of each line as follows:

```
**sample1,/path/to/vcf/file1,**caller1**
 sample1,/path/to/vcf/file2,**caller2**
 sample2,/path/to/vcf/file3,**caller1**
 sample2,/path/to/vcf/file4,**caller2****
```


Please note that regardless selection for the union or intersection, inconsistent calls among different callers are regarded ambiguity and are marked and removed from further analysis. The list of removed variants will also be stored in the result directory.


## ****Primary analysis and updates of VarSelect database ****

Samples can go through either the family or the paired case-control workflow (namely primary analysis) depending on the study design. Primary analysis can be performed repeatedly by inclusion and/or exclusion of samples. The label of the phenotypic information (e.g. tumor and normal; affected and unaffected) of the samples can be changed and thereafter analyzed, according to the labels specified in the PED file. Samples begun with a hash character '#' in the PED file are excluded from the downstream analysis. In the example of the PED file shown here, the samples “uncle” and “aunt” are herein excluded from the downstream analysis. 

```
#family_id   sample_id    paternal_id     maternal_id   sex   phenotype
 family1     father          0               0           1         1
 family1     mother          0               0           2         1
 family1     daughter      father          mother        2         2
 #family1     uncle           0               0           2         2
 #family1     aunt            0               0           2         2
```

A first full run of primary analysis creates a VarSelect database, which is required for recording all the results of subsequent re-analyses (repetitive primary analyses). The -d option specifies the location of where the database file is stored. 

```
**varselect.pl analysis* *-d /path/to/gemini/db 
* **                     *-p /path/to/modified/ped/file 
* **                     *-m workflow mode 
       *               *-k * *-u**
```

After re-analysis, a new analysis directory will be created named with a new jobid. The results and logs will also be stored in new directory. Filtered variants will be assigned new tag in_analysis_{new_jobid} in the VarSelect database.

## ****Secondary analysis: compare results from any two primary analyses ****

The results from different primary analyses can be compared for specific purposes and is termed ‘secondary analysis’. Comparison between any two primary analyses is specified by the options “-a” and “-b” as follows. 

```
**varselect.pl compare*  *-a jobid_of_analysisA 
* **                     *-b jobid_of_analysisB 
* **                     *-c [1-4]
        **              *-d /path/to/db_file*
```

The secondary analysis includes four comparisons including: 1) the union of analysis A and B. 2) intersection of analysis A and B. 3) variants presented in the analysis A but not in the analysis B. 4) variants present in the analysis B but not in the analysis A. Results of new secondary analysis will be stored in a new directory. The filtered variants will be assigned a tag in_analysis with the Job ID.

## **Description of VarSelect scripts**

Bellows are the short description of each VarSelect script. 

* **varselect.pl** is the main script of VarSelect, and includes three commands: annotate (initial annotation), analysis (primary and re-analysis) and compare (secondary analysis). 
    * Command “annotate” triggers the script vs_annotate.pl to annotate VCF files from scratch. It also triggers the script vs_analysis.pl to specify workflow of choice. There are three required options: -v sample-vcf file list, -p PED file, -m workflow mode. 
    * Command “analysis” triggers the script vs_analysis.pl. There are three required options: -d gemini db file, -p PED file, -m workflow mode. 
    * Command “compare” triggers the script vs_compare.pl to compare results between two primary analysis. There are three required options: “-a” and “-b” to specify the Job Id of the analysis A and B. Option “-c” specifies the method of comparison (1. union, 2. intersection, 3. A only, and 4. B only). Option “-d” specifies the VarSelect database file.
* **vs_annotate.pl** is triggered by varselect.pl for three functions. 
    * Firstly, it processes the VCF files, of which from same sample will be joined together by vcf-concat of VCFtools.[9] VCF files of different samples are then merged into a single VCF file by vcf-merge included in the VCFtools. The variants in the sex chromosome in merged VCF file will be fixed for ploidy by vcf-fix-ploidy in the VCFtools.
    * Secondly, the script triggers VEP,[4] snpEff[8] and ANNOVAR[7] for annotation.
    * Thirdly, the script triggers the Gemini framework to generate a SQLite database for downstream analysis. 
* **vs_analysis.pl** is triggered by varselect.pl to analyze variants by taking into account the genotypes, disease status, and the provided copy number and gene expression information. 
    * If option -k specified, multi-caller mode is turned on with union (option -u) or intersection (option -i) set of variants from different callers.
    * If option -c specified, the script triggers cnvkit_parse.pl to annotate copy number variation to the variants.
    * If option -x specified, the script triggers xprprofile_parse.pl to annotate gene expression information to the variants.
    * The script supports two built-in analytic workflows by specifying -m option. For the paired case-control workflow, it triggers loh_detector.pl and denovo_detector.pl scripts to classify variants of loss of heterozygosity or of de novo changes. For the family workflow, the script triggers five scripts including Autosomal-recessive.py, Compound-het.py, Denovo-recessive.py, Two-hits.py and X-linked.py to classify the variants into the related genetic models. All variants selected by the workflow will be annotated on the tag is_analysis_{JobID} in the varselect_variants table of the database.
* **vs_compare.pl** is triggered by varselect.pl to extract the intersection, union, or subtraction of variants between any two primary analyses.
* **vsl_annot.pl** is triggered by vs_annotation.pl, vs_analysis.pl and vs_compare.pl to add new columns and that attaching annotation in the VarSelect database.
* **create_view_for_varselect.pl** is triggered by vs_annotation.pl, vs_analysis.pl and vs_compare.pl to create a view by merging the updated tables in the VarSelect database file.
* **run_vep.pl** is triggered by vs_annotate.pl to annotate VCF file through the Ensembl VEP.[4] The script also enables the dbNSFP plugin[5, 6] to annotate the non-synonymous variants.
* **run_snpeff.pl** is triggered by vs_annotate.pl to annotate VCF file through snpEff.[8] 
* **run_annovar.pl** is triggered by vs_annotate.pl to annotate VCF file through ANNOVAR.[7]
* **GO_parse.pl** is triggered by vs_annotate.pl to annotate variants with the Gene Ontology (GO) terms.[16] The GO terms are assigned based on the transcripts where variants reside.
* **pathway_parse.pl** is triggered by vs_annotate.pl to annotate variants on KEGG pathways.[17-19] Associated pathways will be extracted by the built-in “gemini pathways” in the GEMINI framework.
* **cnvkit_parser.pl** is triggered by vs_analysis.pl to annotate variants with the information of copy number variation. It extracts the log2 values where each variant resides and compute the CNV changes between the paired samples.
* **xprprofile_parser.pl** is triggered by vs_analysis.pl to annotate variants with gene express profiles. The script supports the annotation of four gene expression profiles, including gene-based read counts, exon-based read counts, gene-based tpm (Transcripts Per Kilobase Million), and gene-based fpkm (Fragments Per Kilobase Million). Please specify the -x option followed by the file of gene expression profiles.
* **loh_detector.pl** is triggered by vs_analysis.pl for the paired case/control analytic workflow. The script checks variant genotype between the affected and unaffected samples, and classifies the loss of heterozygosity (LOH) variants. Variants classified as LOH are labeled with the tag “is_LOH”. Further description of the paired case/control workflow is described in the next section. 
* **denovo_detector.pl** is triggered by vs_analysis.pl for the paired case/control analytic workflow. The script compares variants both in the control and the case sample for *de novo* somatic mutations, and labels with the tag “is_denovo”. Further description of the paired case/control workflow is described in the next section.
* **Autosomal-recessive.py** is triggered by vs_analysis.pl for the family analytic workflow. The script filters variants of autosomal recessive inheritance. Filtered variants are labelled with the tag “is_AR”. 
* **Compound-het.py** is triggered by vs_analysis.pl for the family analytic workflow. The script filters variants of compound heterozygosity inheritance. Filtered variants are labelled with the tag “is_CH”. 
* **Denovo-recessive.py** is triggered by vs_analysis.pl for the family analytic workflow. The script filters variants of *de novo* recessive mutations. Filtered variants are labelled with the tag “is_DR”. 
* **Two-hits.py** is triggered by vs_analysis.pl for the family analytic workflow. The script filters variants of two-hit recessive mutations. Filtered variants are labelled with the tag “is_TH”. 
* **X-linked.py **is triggered by vs_analysis.pl for the family analytic workflow. The script filters variants of X-chromosome linked recessive inheritance. Filtered variants are labelled with the tag “is_XL”. 

## ****Examples****

The following examples are some common scenario.

## ****Example 1 - samples of a family study****

In this example, the variants on the chromosome 22 of a family trio, including NA12878, NA12891, NA12892 are used to demonstrate for family-based analysis.
 All required files are stored in the directory “varselect/examples/example1/”. 

    * NA12878_chr22.vcf.gz
    * NA12891_chr22.vcf.gz
    * NA12892_chr22.vcf.gz
    * example1.txt: A comma separated file that describes the links between the samples and the associated VCF files. 
    * example1.ped: A tab separated file that describes the pedigree information of this family, including gender and disease status. In this example, we hypothesize NA12878 is a sick child and that NA12891 and NA12892 are unaffected parents. The instruction of creating a ped file is available at the following website: [`http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped`](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped)

To run the analysis, please type in following command:

```
**varselect.pl annotate -v example1.txt***** ** ** 
*                      -p example1.ped***** 
                      *-m family**
```

VarSelect writes down the timestamp as Job Id when analysis begins (ex: 20170105162548), and, therefore, directory “VarSelectAnalysisResult_20170105162548/” is created. The log files and the intermediate results are stored in this directory. When the primary analysis is finished, example1_varselect.db is created. 
Through the family workflow, five specific columns come with the Job Id: is_AR_20170105162548, is_CH_20170105162548, is_DR_20170105162548, is_TH_20170105162548 and is_XL_20170105162548 corresponding the following genetic models: autosomal recessive, compound heterozygosity, de novo recessive, two-hit recessive and X-linked recessive, respectively.


The following command filters the variants of compound heterozygosity inheritance and that shows the information of chromosome, position, ref allele, alt allele, gene and genotype of the selective variants.

```
**gemini query --header 
* **            *-q 'select chrom,start,ref,alt,gene,gts from variants 
                 ****where is_CH_20170105162548 = 1'***** 
             **example1_varselect.db**                   ***
```


 Column “in_analysis_20170105162548” records variants qualifying the criteria of workflow analysis, and, in this case, is a union of the five genetic models. Variants of the five genetic models can be filtered out by following command: 

```
**gemini query --header* *-q 'select chrom,start,ref,alt,is_AR_20170105162548,***** 
*                          is_CH_20170105162548,is_DR_20170105162548, 
                          ****is_TH_20170105162548,is_XL_20170105162548 ****from variants 
                          where in_analysis_20170105162548=1'*****   
              *example1_varselect.db*   ***
```

 

## **Example 2 - paired case/control samples**

In this case, we use the variants of chromosome 22 of three samples, including blood, primary tumor ccRCC (clear cell Renal Cell Carcinoma), and a metastasis lung cancer available at [`https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063388`](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063388).
 All files are stored in subdirectory examples/example2/

* blood-gatk-chr22.vcf.gz 
* ccRCC-gatk-chr22.vcf.gz
* meta-lung-gatk-chr22.vcf.gz
* example2.txt: A comma separated file that describes the  links between the sample and the corresponding VCF file.
* example2.ped: A tab separated file that describes the  relationship between each sample, gender and disease status. In this case,  blood is assigned as an unaffected (control) sample, while ccRCC and  meta-lung samples are assigned as affected samples.

To start the analysis, please type in the following command. 

```
**varselect.pl annotate -v example2.txt ***** 
                      *-p example2.ped***** 
                      *-m paired**
```

When the analysis is finished, the file “example2_varselect.db” is created at the same directory.VarSelect compares the genotypic difference between the control (unaffected) and the case (affected) samples, and assign ‘LOH’ or ‘de_novo’ tags for loss of heterozygosity or *de novo* somatic mutations. There are several columns including ‘is_LOH’, ‘LOH_samples’, ‘is_denovo’ and** ‘**denovo_samples’ specific in the paired case/control workflow. 

 The following commands filters out the *de novo* somatic variants.

```
**gemini query --header ****-q 'select chrom, start, ref, alt, gts, 
****                          denovo_samples_20170106175124 from ****variants 
                          where is_denovo_20170106175124=1'
             ****example2_varselect.db**
```

## ****Example 3 – Re-analysis and comparison (secondary analysis)****

Multiple primary analysis (re-analysis) can be performed for various study purposes. Comparison between any two primary analyses provides the flexibility of hierarchical comparison, namely secondary analysis. User can repeat analysis by marking labels written in the ped file. For example, in the case of ccRCC, you can filter the common *de novo* mutations presented in both ccRCC and the metastatic lung samples by performing “blood vs. ccRCC” and “blood vs. meta-lung’ analyses. The third sample in the ped file marked with '#' will be excluded in this analysis.Firstly, replicate the ped file in the sample 2 as follows.

```
**cp example2.ped example3_mark1.ped**
```


 Edit the example3_mark1.ped file and exclude meta-lung sample by marking the ‘#” sign.

```
**example3*        *blood*              *0*             *0*             *1*             *1
 example3*        *ccRCC*              *0*             *0*             *1*             *2
 #example3 *    *meta-lung*             *0*             *0*             *1*             *2**
```

Perform primary analysis (re-analysis) by using this new ped file. User can use “varselect.pl analysis” to perform analysis and update the existing database file as the following command.

```
**varselect.pl analysis -d example2_varselect.db 
                      ****-p example3_mark1.ped 
                      ****-m paired**
```

A new Job Id is generated at the time of beginning the re-analysis (ex: 20170109183415), and a new directory “VarSelectAnalysisResult_20170109183415” is created to store the log files and intermediate result files of this new analysis job. When the re-analysis is finished, new columns tagged with new Job Id are added into the varselect_variants table. 
 
User can use following command to filter out the *de novo* somatic mutations in the ccRCC.

```
**gemini query --header ****-q 'select chrom,start,ref,alt,gts from variants 
                          ****where is_denovo_20170109183415 = 1'
             ****example2_varselect.db**
```

 
Replicate the example3_mark1.ped file and save as the example3_mark2.ped. Mark the sample meta-lung with the ‘#’ sign. 

```
**example3 *       *blood *          *0*             *0*             *1*             *1
 #example3*  **      *ccRCC *          *0*             *0*             *1*             *2
 example3 *      *meta-lun *        *0*             *0*             *1*             *2**
```

 
Perform primary analysis again with the new edited ped file to identify the *de novo* variants present in the metastatic lung tumor.

```
*** *varselect.pl analysis -d example2_varselect.db
                       **-p example3_mark2.ped
                       **-m paired******
```

A new directory “./VarSelectAnalysisResult_20170110174312” is created with the corresponding tags updated in the database.
You can extract the *de novo* mutation presented in the metastasis tumor by the following command.

```
*** *gemini query --header ****-q 'select chrom,start,ref,alt,gts from variants 
                           where ****is_denovo_20170110174312 = 1'
              ****example2_varselect.db**
```

 
To identify the de novo mutations that present in both ccRCC and lung meta samples, we can compare the results from the two (primary) analysis (jobs) by intersecting the results. Please note that only analyses stored in the same database (db) file can be compared (namely secondary analysis). In this example, common *de novo *mutations can be selected by typing in the following command.

```
**varselect.pl compare -a 20170109183415
                     ****-b 20170110174312
                     ****-c 2 *****              
     *                -d example2.db**
```

 
A new Job Id is assigned (say, 20170111082651) and the database (db) will be updated with the results of the secondary analysis. The following command lists the results stored in the database (db). 

```
**gemini query --header ****-q 'select chrom, start, ref, alt, gts from variants***** 
*                          where in_analysis_20170111082651 = 1'***** 
             *example2_varselect.db**
```

 
If you want to select the novel mutations present only in the metastasis tumor, change the option –c to ‘4’ (unique to the analysis in the ‘-b’ option). 

```
**varselect.pl compare -a 20170109183415 ***** 
*                     -b 20170110174315*****  **  
                     *-c 4
                     ****-d example2_varselect.db**
```

## ****Example 4 – comparison of multiple variants callers****

Union and/or intersection of variant calls from different variants callers could be of interests. This example demonstrates manipulation of results of two popular variant callers: ‘GATK-HaplotypeCaller’ and ‘freebayes’.
 
All files are stored in directory examples/example4/:

    * NA12878-gatk-chr22.vcf.gz
    * NA12891-gatk-chr22.vcf.gz
    * NA12892-gatk-chr22.vcf.gz
    * NA12892-freebayes-chr22.vcf.gz: Variant calls from freebayes of chromosome 22 in these genomes
    * NA12878-freebayes-chr22.vcf.gz
    * NA12891-freebayes-chr22.vcf.gz
    * example4.txt: A comma separated file that describes the links between the samples and the corresponding VCF files. For multi-caller function, an extra column is required to describing the corresponding variant caller.

```
**NA12878,NA12878-gatk-chr22.vcf.gz,gatk-haplotype,
 NA12891,NA12891-gatk-chr22.vcf.gz,gatk-haplotype
 NA12892,NA12892-gatk-chr22.vcf.gz,gatk-haplotype
 NA12878,NA12878-freebayes-chr22.vcf.gz,freebayes
 NA12891,NA12891-freebayes-chr22.vcf.gz,freebayes
 NA12892,NA12892-freebayes-chr22.vcf.gz,freebayes**
```

    * example4.ped: A tab separated file that describes the relationship between each sample, gender and disease status.

 
The ‘-k’ option triggers the multi-caller function, followed by ‘-u’ option for preparing the union of all variants from the two variant callers. 

```
*** *varselect.pl annotate ****-v example4.txt 
                       ****-p example4.ped 
                       ****-m family 
                       ****-k -u **
```

 
A new Job Id (ex:20170111162057) is tagged on a new directory “VarSelectAnalysisResult_20170111162057”. Please note that the variants with inconsistent calls by different callers are removed from the downstream analysis.
 
The full list of the removed variants could be found in the file “multicaller_union_inconsistant_20170111162057.txt”.
 
The following command lists the variants in this analysis.

```
**gemini query --header ****-q 'select chrom, start, ref, alt, gts from variants
                         ****where in_analysis_20170111162057 = 1'
             ****example4_varselect.db**
```

 
To analyze the intersection of variants identified by both variant callers, replace the option ‘-u’ to ‘-i’ in the original command as following.

```
varselect.pl analysis -d example4_varselect.db 
                      -p example4.ped 
                      -m family 
                      -k -i 
```

Please note that the two options ‘-u’ and ‘-i’ are mutually exclusive. Inconsistent calls by different callers are removed for the downstream analysis and listed in file “multicaller_intersect_inconsistant_20170111190712.txt” in the same directory of analysis output.


## **Known issues**

1. Decomposition of multiple alternative alleles on the same position results in inconsistent number of allelic quality value and that causing errors on merging VCF files. These variants are currently excluded from downstream analysis.
2. Variants of multiple-nucleotide substitution are incorrectly annotated as “indel” by the Gemini framework. This problem only appeared in the scenario when the VCF files come from the Ion Torrent platform or converted by the CGA (Complete Genomics Analysis) tools. 



## **Family analysis**

The family analysis workflow includes three steps: 1) analyzing the core family members including father, mother and sick child either male or female, the parents could be both healthy or either one is sick; 2) excluding the variants appeared in the other healthy members in the pedigree; 3) keeping the variants appeared in the other sick members in the pedigree. The inheritance models include: autosomal recessive (AR), compound heterozygosity (CH), *de novo* recessive (DNR), two-hits recessive (THR), and X chromosome-linked recessive (XLR). In humans, the pseudoautosomal region (PAR) on the X chromosome is diploid in male while the rest of the nucleotide sequence is haploid. Variants found on the X chromosome that fall in the CH, DNR, and THR models are reported on each of the category, and excluded from the result of the XLR model. 
 
**Autosomal recessive inheritance**
 
In the AR model, both pathogenic alleles came from each of the parents while the genotype of both parents are heterogeneous and includes exactly one pathogenic allele. When one of the parents is also sick, the sick parent shall carry homozygous variants, while the other parent shall be heterozygous and carries only one pathogenic variant. 
 
**Compound heterozygosity inheritance**
 
For discovering the pathogenic alleles of compound heterozygosity (CH), two candidate alleles shall locate in the same gene but on different positions yet chromosomes from each of the patients. The genotypes of the CH variants in both parents have to be heterozygous and includes exactly one pathogenic allele. If one of the parents is sick, the candidate gene shall be also compound heterozygous as the sick child. For the genes on the X chromosome, only the PAR region is scanned on the sick boys. Because of the healthy father, only the PAR region is scanned for the sick girls. 
 
***De novo*** **recessive mutation**
 
In the *de novo* recessive (DNS) model, the variants of the sick child shall be homozygous and have to be different to that of the parents. For that on the X chromosome, the variants on the non-PAR region in the male and female patients shall be hemizygous and homozygous, respectively, while the variants on the PAR shall be homozygous for both the male and the female patients. When one of the parents is sick, the homozygous variant in the sick child shall be different to that in the sick parent. 
 
**Two-hit recessive mutation**
 
In the two-hit model, the sick child inherits one pathogenic allele from one of the parents while the other allele is somatically acquired. When healthy, one of the parents has to be a heterozygous carrier. For the variants on the X chromosome, the entire chromosome is scanned for sick girls, and is only the PAR region for sick boys. When one of the parents is sick, he/she shall carry homozygous alleles of pathogenicity, while the other is either homozygous wild-type or not carrying the pathogenic allele. 
 
**X-linked recessive inheritance**
 
For detecting X-linked recessive (XLR) inheritance, the variants on the non-PAR region in the male and female patients are hemizygous and homozygous, respectively. The variants on the PAR shall be homozygous for both the male and female patients. Typically, the sick boy shall inherit the pathogenic allele on the non-PAR region from his healthy or sick mother who carries heterozygous or homozygous variants, respectively. When the father is healthy, the girl received healthy hemizygous sequences on the non-PAR region from her father. Her disease is likely not related to the variants on non-PAR region, and the non-PAR region is therefore skipped unless the father is also sick. Depending on which parent is sick, either one or two pathogenic alleles shall be seen on the non-PAR region of the father or mother, respectively, while are always two pathogenic alleles on the PAR region. 
 
## **Paired case/control analysis**
 
In the case/control workflow, the variants of the matched samples are classified as either ‘loss of heterozygosity (LOH)’ or ‘*de novo* mutations’. All the different variants in the vcf files of both the case and control, regardless single and multi-nucleotide substitutions, or small INDEL, are compared for LOH otherwise classified as *de novo* mutations. VarSelect supports the comparison between one or more cases versus one control sample.
 
**Loss of heterozygosity (LOH) events**
 
Variants that are heterozygous in the control sample lose it heterozygosity is so called LOH. When multiple cases are analyzed, each case is individually compared against the same control sample and annotated on the “is_LOH” tag. 
 
***De novo*** **mutation**

*De novo* mutation can occur on the alleles of heterozygosity and homozygosity, as long as the alternative nucleotide in the case isn’t present in the control sample. When multiple cases are analyzed, each case is individually compared against the same control sample and annotated on the “is_denovo” tag.

