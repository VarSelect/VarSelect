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
use Storable ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
my $opts = {} ;
getopts("c:i:j:p:o:d:" , $opts) ;
my $jobid		= $opts->{j} || getts() ;
my $file_vcfgz_input	= $opts->{i} || die "Error: -i VCF file is required!\n\n"	. usage() ;
my $file_list_cns	= $opts->{c} || die "Error: -c cns file list is required!\n\n"	. usage() ;
my $file_ped		= $opts->{p} || die "Error: -p ped file is required!\n\n"	. usage() ;
my $file_db		= $opts->{d} || die "Error: -d db file is required!\n\n"	. usage() ;

my ($input_name,$input_path,$input_ext) = fileparse($file_vcfgz_input, qw/.vcf .vcf.gz/) ;
my $file_vcfgz_output	= $opts->{o} || $input_path . $input_name . "\_cnv.vcf.gz" ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
make_path($dir_log , $dir_work ,{ chmod => 0755 ,} ) unless (-e $dir_work && -e $dir_log) ;

my $file_output_cnv = "$dir_work/$input_name\_$jobid\_cnv.txt" ;

# start logging
my $file_log = "$dir_log/log_cnvkitparser_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("CNVkit parser start") ;


# Step 0. Loading PED, get sample infomation
my $affected_samples = [] ;
my $unaffected_samples = [] ;
open (my $SRC_ped , "$file_ped") ;
while (<$SRC_ped>) {
    my ($family,$sample,$fatherid,$motherid,$sexid,$affectid,$other) = split /\t/ ;
    push @$affected_samples , $sample if ($affectid == 2) ;
    push @$unaffected_samples , $sample if ($affectid == 1) ;
} 
close $SRC_ped ;

# Step 1. v-annotate *.cns to *.cnstmp.vcf
my $files_tocombine = [] ;
my $unique_sample_hash = {} ;
open (my $SRC_list_cns , "$file_list_cns") ;
while (my $cns_line = <$SRC_list_cns>) {
    chomp $cns_line ;
    my ($sample,$file_cns) = split /\,/ , $cns_line ;
    $unique_sample_hash->{$sample} ++ ;

    # step 1.1 bgzip tabix cns
    my ($fn_cns,$fpath_cns,$fext_cns) = fileparse($file_cns,qw/.cns/) ;
    my $file_gz_cns_for_analysis = "$dir_work/$sample\_$fn_cns.cns.gz" ;

    my $cmd_bgzip = "bgzip -c $file_cns > $file_gz_cns_for_analysis " ;
    $log->andRun($cmd_bgzip) ;

    my $cmd_tabix = "tabix -s 1 -b 2 -e 3 -S 1 -f $file_gz_cns_for_analysis" ;
    $log->andRun($cmd_tabix) ;

    my $file_vcfgz_cnstmp = "$dir_work/$input_name.cnstmp.$sample.vcf.gz" ;
    push @$files_tocombine , $file_vcfgz_cnstmp ;
    
    # step 1.2 v-annotate 
    my $cmd_vannot1 = "zcat $file_vcfgz_input | vcf-annotate -a $file_gz_cns_for_analysis " ;
    $cmd_vannot1 .= " -c " . join("," , qw"CHROM FROM TO - INFO/cnslog2 - ") ;
    $cmd_vannot1 .= " -d key=INFO,ID=cnslog2,Number=1,Type=Float,Description='cnslog2' " ;
    $cmd_vannot1 .= " 2> $dir_log/stderr_cnv_vannot_$jobid\_$sample.log " ;
    $cmd_vannot1 .= " |bgzip -c > $file_vcfgz_cnstmp" ;

    $log->andRun($cmd_vannot1) ;
    
}
close $SRC_list_cns ;

# Step 2. combine *.cnstmp.vcf into combine.cnvlog2
# step 2.1 create vcf parsers
my $sample_list = [sort keys %$unique_sample_hash] ;
my $vcf_parsers = {} ;
foreach my $sample (@$sample_list) {
    my $file_vcf_cnstmp = "$dir_work/$input_name.cnstmp.$sample.vcf.gz" ;
    $vcf_parsers->{$sample} = Vcf_parser->new (file => $file_vcf_cnstmp) ;
}

# step 2.2 go through each line
my $file_combinecns = "$dir_work/$input_name.cns_combine.$jobid.txt" ;
open (my $TGT_combine_output,">$file_combinecns") ;
my $first_sample = shift @$sample_list ;
while (my $line_sample1 = $vcf_parsers->{$first_sample}->next_var) {

    my $sample2log2 = {
	$first_sample => $line_sample1->{INFO}->{cnslog2} ,
    } ;
    my $output_sample = [$first_sample] ;

    foreach my $sample (@$sample_list) {
	$sample2log2->{$sample} = $vcf_parsers->{$sample}->next_var->{INFO}->{cnslog2} ;
	push @$output_sample , $sample ;
    }

    my $output_log2 = [ map {$sample2log2->{$_} } @$output_sample ] ;

    
    # Calculate fold change for any aff/unaff sample pair
    my $output_pair = [] ;
    my $output_fc = [] ;
    foreach my $sample_aff (@$affected_samples) {
	foreach my $sample_unaff (@$unaffected_samples) {
	    my $sample_pair = "$sample_aff/$sample_unaff" ;
	    my $fc = ($sample2log2->{$sample_unaff}) ? $sample2log2->{$sample_aff}/$sample2log2->{$sample_unaff} : '-' ;

	    push @$output_pair , $sample_pair ;
	    push @$output_fc , $fc ;
	}
    }

    print $TGT_combine_output join("\t" , map {$line_sample1->{$_}}  qw/CHROM POS REF ALT/) ."\t" . join (",", @$output_sample) ."\t".  join("," , @$output_log2) . "\t" .join(",", @$output_pair) . "\t" . join("," , @$output_fc) . "\n" ;
}
close $TGT_combine_output ;

tabix_vcf($file_combinecns) ;	# include bgzip

# Step 3. vcf-annotate again
my $cmd_vannot2 = "zcat $file_vcfgz_input | vcf-annotate -a $file_combinecns.gz " ;
$cmd_vannot2 .= " -c " . join("," , qw"CHROM POS REF ALT INFO/cnv_samples INFO/cnv_log2 INFO/cnv_fc_samples INFO/cnv_foldchange_log2") ;
$cmd_vannot2 .= " -d key=INFO,ID=cnv_samples,Number=1,Type=text,Description='cnv_samples provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=cnv_log2,Number=1,Type=text,Description='cnv_log2 provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=cnv_fc_samples,Number=.,Type=text,Description='cnv_foldchange_samples provided by VarSelect'" ;
$cmd_vannot2 .= " -d key=INFO,ID=cnv_foldchange_log2,Number=.,Type=text,Description='foldchange of cnv log2 provided by VarSelect'" ;
$cmd_vannot2 .= " 2> $dir_log/stderr_cnv_vannot_$jobid\_combine.log " ;
$cmd_vannot2 .= " |bgzip -c >$file_vcfgz_output " ;

$log->andRun($cmd_vannot2) ;

tabix_vcf($file_vcfgz_output) ;

$log->write("CNVkit parser finish") ;

sub usage {
return <<EOUsage
$0 
Usage:
    $0 -i file -o file -p file

    -i  VCF file input                      (Required)
    -o  VCF file output                     (Required)
    -c  cns file list                       (Required)

    -h  this page

EOUsage
}


