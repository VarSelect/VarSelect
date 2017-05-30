#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Path qw(make_path remove_tree);
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;

my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $jar_snpeff = $Setting->{jar_snpeff} ;
my $dir_snpeff = $Setting->{dir_snpeff} ;
my $genome_build = $Setting->{genome_snpeff} ;

# Options handling
my $opts = {} ;
getopts("i:o:j:n:d:", $opts) ;

my $jobid             =	$opts->{j} ;
my $file_vcfgz_input  = $opts->{i} ;
my $file_vcfgz_output = $opts->{o} ;
my $file_geminidb     = $opts->{d} ;
my $threads_active    =	$opts->{n} || 16 ;

my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $dir_temp_for_gemini = "$dir_work/temp_for_gemini_snpEff" ;
make_path( $dir_temp_for_gemini , { chmod => 0755 ,} ) unless (-e $dir_temp_for_gemini) ;

my $file_log = "$dir_log/running_snpeff_$jobid.log" ;

my $file_vcf_snpeff_tmp = "$dir_work/snpeff_tmp.vcf" ;
my $file_db_snpeff_tmp = "$dir_work/snpeff_tmp.db" ;
my $file_gemini_dump  = "$dir_work/snpeff_geminidump.txt" ;
my $file_gemini_dump_tmp  = "$dir_work/snpeff_geminidump_tmp.txt" ;

my $log = VarSelect::Log->new(file => $file_log) ;

# Step 1. run snpeff
open (my $LOG, ">$file_log") ;
my $file_csvstats = "$dir_work/snpeff_csvstats.csv" ;
my $file_stat = "$dir_work/snpeff_stats.html" ;

my $cmd_snpeff =  "java -jar $jar_snpeff -dataDir $dir_snpeff -csvStats $file_csvstats -s $file_stat $genome_build $file_vcfgz_input > $file_vcf_snpeff_tmp 2> $dir_log/stderr_snpEff_jar_$jobid.log" ;
$log->andRun($cmd_snpeff) ;

bgzip($file_vcf_snpeff_tmp, "$file_vcf_snpeff_tmp.gz") ;
tabix_vcf("$file_vcf_snpeff_tmp.gz") ;

# Step 2. gemini db load
my $cmd_gemini_load = "gemini load " ;
$cmd_gemini_load .= " -t snpEff " ;
$cmd_gemini_load .= " -v $file_vcf_snpeff_tmp.gz " ;
$cmd_gemini_load .= " --cores $threads_active " ;
$cmd_gemini_load .= " --tempdir $dir_temp_for_gemini " ;
$cmd_gemini_load .= " $file_geminidb 2> $dir_log/stderr.geminiload.snpeff.log" ;

$log->andRun($cmd_gemini_load) ;

