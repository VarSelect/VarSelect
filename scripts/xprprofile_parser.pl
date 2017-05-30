#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Basename ;
use File::Path qw(make_path remove_tree);
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;


my $file_gff_dexseq = $Setting->{gff_dexseq} ;
my $file_gtf = $Setting->{gtf_ensembl_transcript} ;
my $threads_active = 16 ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
# -p $file_ped -x $file_xpr_list -i $file_vcfgz_for_workflow -o $file_vcfgz_xpr -j $jobid -l $dir_log 
my $opts = {} ;
getopts("x:i:j:l:p:o:" , $opts) ;

my $jobid	    = $opts->{j} ;
my $file_xpr_list   = $opts->{x} ;
my $file_vcf_input  = $opts->{i} ;
my $file_ped	    = $opts->{p} ;
my $file_output	    = $opts->{o} ;
my ($input_name,$input_path,$input_ext) = fileparse($file_vcf_input, qw/.vcf .vcf.gz/) ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $file_log = "$dir_log/log_xprprofile_parser_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;


#step 0 read files
# step 0-1 dexseq gff3
# zzz
$log->write("Loading $file_gff_dexseq start") ;
my $agenepos = {} ;	# aggregate_gene
my $ageneexonpos = {} ; # aggregate_gene
open (my $SRC_gff , "$file_gff_dexseq") ;
while (<$SRC_gff>) {
    chomp ;
    my @data = split /\t/ ;
    my $anno = $data[8] ;

    if ($data[2] eq 'aggregate_gene') {
	$anno =~ s/\"//g ;
	$anno =~ s/gene_id // ;
	$agenepos->{$anno} = "$data[0]\t$data[3]\t$data[4]" ;

    } elsif ($data[2] eq 'exonic_part') {
	$anno =~ s/\"//g ;
	my ($transcript,$exon,$gene) = split (/; / , $anno) ;
	$transcript =~ s/transcript // ;
	$exon =~ s/exonic_part_number // ;
	$gene =~ s/gene_id // ;
	my $exontag = "$gene:$exon" ;
	$ageneexonpos->{"$exontag"}  =  "$data[0]\t$data[3]\t$data[4]" ;
    }
}
close $SRC_gff ;
$log->write("Loading $file_gff_dexseq start") ;

# step 0-2 gtf
# zzz
$log->write("Loading $file_gtf start") ;
my $genepos = {} ;
open (my $SRC_gtf, "$file_gtf") ;
while(<$SRC_gtf>) {
    chomp ;
    my ($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attr) = split /\t/ ;
    my @data =  split(/; /,$attr) ;

    foreach my $attrdata (@data) {
	my ($idx,$val) = split / / , $attrdata ;

	if ($idx eq 'gene_id' && $feature eq 'gene') {
	    $val =~ s/\"//g ;
	    $genepos->{$val} = "$chr\t$start\t$end" ;
	}
   }
}
close $SRC_gtf ;
$log->write("Loading $file_gtf finish") ;


# step 0-3 xpr files
my $sample2file = {} ;
open (my $SRC,"$file_xpr_list") ;
while (my $line = <$SRC>) {
    chomp $line ;
    my ($sample,$type,$file) = split /\,/ , $line ;
    if ($type eq 'readcount') { 
	$type = 'rc' ;
    }

    if ($type eq 'exoncount') {
	$type = 'exon' ;
    }

    $sample2file->{$sample}->{$type} = $file ;
}
close $SRC ;
my $sample_list = [keys %$sample2file] ;

# step 0-4 
my $affected_samples = [] ;
my $unaffected_samples = [] ;
open (my $SRC_ped , "$file_ped") ;
while (<$SRC_ped>) {
    my ($family,$sample,$fatherid,$motherid,$sexid,$phenotype,$other) = split /\t/ ;
    push @$affected_samples   , $sample if ($phenotype == 2) ;
    push @$unaffected_samples , $sample if ($phenotype == 1) ;
}
close $SRC_ped ;


#step 1-0 convert into bed files
my $files_to_combine_in_vannot2 = {} ;
foreach my $sample (keys %$sample2file) {

    # dexseq files convert to bed then vcf-annotate
    #
    # readcount
    # step 1-1 vannot *.counts.bed -> *.rctmp.vcf.gz
    my ($filename,$filepath,$fileext) = fileparse($sample2file->{$sample}->{rc}) ;
    my $file_bed = "$dir_work/$sample\_$filename.bed" ;
    my $file_vcfgz_rc_tmp = "$dir_work/$sample\_$filename.rctmp.vcf.gz" ;

    run_vannot1 ($sample2file->{$sample}->{rc} , $file_bed , $file_vcf_input , $file_vcfgz_rc_tmp , "readcount" , "Integer") ;
    $files_to_combine_in_vannot2->{$sample}->{rc} = $file_vcfgz_rc_tmp ;

    # tpm 
    # step 1-2 vannot *.fpkm -> bed -> *tpmtmp.vcf.gz
    my ($tpmfilename,$tpmfilepath,$tpmfileext) = fileparse($sample2file->{$sample}->{tpm}) ;
    my $file_bed_tpm = "$dir_work/$sample\_$tpmfilename.bed" ;
    my $file_vcfgz_tpm_tmp = "$dir_work/$sample\_$tpmfilename.tpmtmp.vcf.gz" ;

    run_vannot1 ($sample2file->{$sample}->{tpm} , $file_bed_tpm , $file_vcf_input , $file_vcfgz_tpm_tmp , "tpm" , "Float") ;
    $files_to_combine_in_vannot2->{$sample}->{tpm} = $file_vcfgz_tpm_tmp ;

    # fpkm
    # step 1-3 vannot *.tpm -> bed -> *.fpkmtmp.vcf.gz
    my ($fpkmfilename,$fpkmfilepath,$fpkmfileext) = fileparse($sample2file->{$sample}->{fpkm}) ;
    my $file_bed_fpkm = "$dir_work/$sample\_$fpkmfilename.bed" ;
    my $file_vcfgz_fpkm_tmp = "$dir_work/$sample\_$fpkmfilename.fpkmtmp.vcf.gz" ;

    run_vannot1 ($sample2file->{$sample}->{fpkm} , $file_bed_fpkm , $file_vcf_input , $file_vcfgz_fpkm_tmp , "fpkm" , "Float") ;
    $files_to_combine_in_vannot2->{$sample}->{fpkm} = $file_vcfgz_fpkm_tmp ;

    # exon
    # step 1-4 vannot *.dexseq -> bed -> *.exontmp.vcf.gz
    my ($exonfilename,$exonfilepath,$exonfileext) = fileparse($sample2file->{$sample}->{exon}) ;
    my $file_bed_exon = "$dir_work/$sample\_$exonfilename.bed" ;
    my $file_vcfgz_exon_tmp = "$dir_work/$sample\_$exonfilename.exontmp.vcf.gz" ;

    run_vannot1_agene ($sample2file->{$sample}->{exon} , $file_bed_exon , $file_vcf_input , $file_vcfgz_exon_tmp , "exonreadcount" , "Integer") ;
    $files_to_combine_in_vannot2->{$sample}->{exon} = $file_vcfgz_exon_tmp ;
}

#step 2 combine *.xpr_xxx.vcf -> combined.xpr_rc.vcf    etc....
my $rc_vcf_parsers = {} ;
my $tpm_vcf_parsers = {} ;
my $fpkm_vcf_parsers = {} ;
my $exon_vcf_parsers = {} ;
foreach my $sample (keys %$sample2file) {
    #rc
    $rc_vcf_parsers->{$sample} = Vcf_parser->new (file=>$files_to_combine_in_vannot2->{$sample}->{rc})  ;

    #tpm
    $tpm_vcf_parsers->{$sample} = Vcf_parser->new (file=>$files_to_combine_in_vannot2->{$sample}->{tpm}) ;

    #fpkm
    $fpkm_vcf_parsers->{$sample} = Vcf_parser->new (file=>$files_to_combine_in_vannot2->{$sample}->{fpkm}) ;

    #exon
    $exon_vcf_parsers->{$sample} = Vcf_parser->new (file=>$files_to_combine_in_vannot2->{$sample}->{exon}) ;
}

# go through each line
my $file_combine_xpr = "$dir_results_output/$input_name.xpr_combine.$jobid.txt" ;
open (my $TGT_combine_output , ">$file_combine_xpr") ;
my $first_sample = shift @$sample_list ;

while (my $line_sample1 = $rc_vcf_parsers->{$first_sample}->next_var) {
    my $sample2rc = {
	$first_sample => $line_sample1->{INFO}->{readcount} +0 ,
    } ;

    my $sample2tpm = {
	$first_sample => $tpm_vcf_parsers->{$first_sample}->next_var->{INFO}->{tpm}+0 ,
    } ;

    my $sample2fpkm = {
	$first_sample => $fpkm_vcf_parsers->{$first_sample}->next_var->{INFO}->{fpkm}+0 ,
    } ;

    my $line_exon = $exon_vcf_parsers->{$first_sample}->next_var ;
    my $sample2exonrc = {
	$first_sample => $line_exon->{INFO}->{exonreadcount}+0 ,
    } ;


    my $output_sample = [$first_sample] ;

    foreach my $sample (@$sample_list) {
	$sample2rc->{$sample} = $rc_vcf_parsers->{$sample}->next_var->{INFO}->{readcount} + 0 ;
	$sample2tpm->{$sample} = $tpm_vcf_parsers->{$sample}->next_var->{INFO}->{tpm}  + 0;
	$sample2fpkm->{$sample} = $fpkm_vcf_parsers->{$sample}->next_var->{INFO}->{fpkm}  + 0;
	$sample2exonrc->{$sample} = $exon_vcf_parsers->{$sample}->next_var->{INFO}->{exonreadcount}  + 0;

	push @$output_sample , $sample ;
    }

    my $output_rc = [map {$sample2rc->{$_} } @$output_sample] ;
    my $output_tpm = [map {$sample2tpm->{$_} } @$output_sample] ;
    my $output_fpkm = [map {$sample2fpkm->{$_} } @$output_sample] ;
    my $output_exonrc = [map{$sample2exonrc->{$_} } @$output_sample] ;
    my $output_exonid = $line_exon->{INFO}->{exonid} ;


    # Calculate fold change for any aff/unaff sample pair
    my $output_pair = [] ;
    my $output_fc_rc = [] ;
    my $output_fc_tpm = [] ;
    my $output_fc_fpkm = [] ;
    my $output_fc_exonrc = [] ;

    foreach my $sample_aff (@$affected_samples) {
        foreach my $sample_unaff (@$unaffected_samples) {
	    my $sample_pair = "$sample_aff/$sample_unaff" ;
	    my $fc_rc = ($sample2rc->{$sample_unaff})? $sample2rc->{$sample_aff} / $sample2rc->{$sample_unaff} : '-' ;
	    my $fc_tpm = ($sample2tpm->{$sample_unaff})? $sample2tpm->{$sample_aff} / $sample2tpm->{$sample_unaff} : '-' ;
	    my $fc_fpkm = ($sample2fpkm->{$sample_unaff})? $sample2fpkm->{$sample_aff} / $sample2fpkm->{$sample_unaff} : '-' ;
	    my $fc_exonrc = ($sample2exonrc->{$sample_unaff})? $sample2exonrc->{$sample_aff} / $sample2exonrc->{$sample_unaff} : '-' ;

	    push @$output_pair , $sample_pair ; 
	    push @$output_fc_rc , $fc_rc ;
	    push @$output_fc_tpm , $fc_tpm ;
	    push @$output_fc_fpkm , $fc_fpkm ;
	    push @$output_fc_exonrc , $fc_exonrc ;
	}
    }

    my $output_vcf_line = join("\t" , map {$line_sample1->{$_}} qw/CHROM POS REF ALT/) ;
    $output_vcf_line .=  "\t" . join("," , @$output_sample) ;
    $output_vcf_line .=  "\t" . join("," , @$output_rc) ;
    $output_vcf_line .=  "\t" . join("," , @$output_tpm) ;
    $output_vcf_line .=  "\t" . join("," , @$output_fpkm) ;
    $output_vcf_line .=  "\t" . join("," , @$output_exonrc) ;
    $output_vcf_line .=  "\t$output_exonid" ;
    $output_vcf_line .=  "\t" . join("," , @$output_pair) ;
    $output_vcf_line .=  "\t" . join("," , @$output_fc_rc) ;
    $output_vcf_line .=  "\t" . join("," , @$output_fc_tpm) ;
    $output_vcf_line .=  "\t" . join("," , @$output_fc_fpkm) ;
    $output_vcf_line .=  "\t" . join("," , @$output_fc_exonrc) ;

    print $TGT_combine_output "$output_vcf_line\n" ;
    
}
close $TGT_combine_output ;
#`bgzip $file_combine_xpr` ;
tabix_vcf("$file_combine_xpr") ;

#step 3 vcf-annotate final 
my $cmd_vannot2 = "zcat $file_vcf_input | vcf-annotate -a $file_combine_xpr.gz " ;
$cmd_vannot2 .= " -c " . join("," , qw"CHROM POS REF ALT INFO/xpr_samples INFO/xpr_readcount INFO/xpr_tpm INFO/xpr_fpkm INFO/xpr_exonreadcount INFO/xpr_exonid INFO/xpr_fc_samples INFO/xpr_foldchange_readcount INFO/xpr_foldchange_tpm INFO/xpr_foldchange_fpkm INFO/xpr_foldchange_exonreadcount") ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_samples,Number=1,Type=String,Description='xpr_samples provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_readcount,Number=1,Type=String,Description='xpr_readcount provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_tpm,Number=1,Type=String,Description='xpr_tpm provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_fpkm,Number=1,Type=String,Description='xpr_fpkm provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_exonreadcount,Number=1,Type=String,Description='xpr_exon read count provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_exonid,Number=1,Type=String,Description='xpr_exonid provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_fc_samples,Number=.,Type=String,Description='xpr_foldchange_samples provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_foldchange_readcount,Number=.,Type=String,Description='xpr_foldchange_readcount provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_foldchange_tpm,Number=.,Type=String,Description='xpr_foldchange_tpm provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_foldchange_fpkm,Number=.,Type=String,Description='xpr_foldchange_fpkm provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=xpr_foldchange_exonreadcount,Number=.,Type=String,Description='xpr_foldchange_exonreadcount provided by VarSelect'" ;
$cmd_vannot2 .= " 2> $dir_log/stderr_xpr_vannot_$jobid\_combine.log " ;
$cmd_vannot2 .= " | bgzip -c > $file_output " ;

$log->andRun($cmd_vannot2) ;


sub run_vannot1_agene {
    # dexseq files -> bed -> vannot -> vcf.gz

    my $dexseq_file = shift ;
    my $file_bed = shift ;
    my $file_vcf_input = shift ;
    my $file_vcfgz_output_tmp = shift ;
    my $tag = shift ;
    my $type = shift ;

    open (my $TGT_bed , ">$file_bed") ;
    open (my $SRC,"$dexseq_file") ;
    while (<$SRC>) {
        chomp ;
        my ($idx,$val) = split /\t/ ;
	next unless (exists $ageneexonpos->{$idx}) ;

        print $TGT_bed "$ageneexonpos->{$idx}\t$val\t$idx\n" ; # output extra field exonid
    }
    close $SRC ;
    close $TGT_bed ;
    my $cmd_sort_bed = "sort -k1,1V -k2,2n -k3,3n $file_bed > $file_bed.sorted ; mv $file_bed.sorted $file_bed" ;
    $log->andRun($cmd_sort_bed) ;
    tabix_bed ($file_bed) ;

    my $cmd_vannot1 = "zcat $file_vcf_input | vcf-annotate -a $file_bed.gz " ;
    $cmd_vannot1 .= " -c " . join ("," , (qw/CHROM FROM TO/, "INFO/$tag" , "INFO/exonid") ) ;
    $cmd_vannot1 .= " -d key=INFO,ID=$tag,Number=1,Type=$type,Description='dexseq $tag' " ;
    $cmd_vannot1 .= " -d key=INFO,ID=exonid,Number=1,Type=String,Description='dexseq exonid' " ;
    $cmd_vannot1 .= " 2> $dir_log/stderr_xpr_vannot_$jobid\_aggreated_bed_gz.log " ;
    $cmd_vannot1 .= " |bgzip -c > $file_vcfgz_output_tmp" ;

    $log->andRun($cmd_vannot1) ;
}

sub run_vannot1 {
    # dexseq files -> bed -> vannot -> vcf.gz

    my $dexseq_file = shift ;
    my $file_bed = shift ;
    my $file_vcf_input = shift ;
    my $file_vcfgz_output_tmp = shift ;
    my $tag = shift ;
    my $type = shift ;

    open (my $TGT_bed , ">$file_bed") ;
    open (my $SRC,"$dexseq_file") ;
    while (<$SRC>) {
        chomp ;
        my ($idx,$val) = split /\t/ ;
	next unless (exists $genepos->{$idx}) ;

        print $TGT_bed "$genepos->{$idx}\t$val\n" ;
    }
    close $SRC ;
    close $TGT_bed ;

    my $cmd_sort_bed = "sort -k1,1V -k2,2n -k3,3n $file_bed > $file_bed.sorted ; mv $file_bed.sorted $file_bed" ;
    $log->andRun($cmd_sort_bed) ;
    tabix_bed ($file_bed) ;

    my $cmd_vannot1 = "zcat $file_vcf_input | vcf-annotate -a $file_bed.gz " ;
    $cmd_vannot1 .= " -c " . join ("," , (qw/CHROM FROM TO/ , "INFO/$tag")) ;
    $cmd_vannot1 .= " -d key=INFO,ID=$tag,Number=1,Type=$type,Description='dexseq $tag' " ;
    $cmd_vannot1 .= " |bgzip -c > $file_vcfgz_output_tmp" ;

    $log->andRun($cmd_vannot1) ;
}

sub tabix_bed {
    my $file = shift ;

    if ($file !~ /gz$/) {
        my $file_orig = $file ;
        $file = "$file_orig.gz" ;
        bgzip($file_orig,$file) ;
    }

    my $cmd_tabix_vcf = "tabix -p bed -f $file " ;
    $log->andRun($cmd_tabix_vcf) ;
}

