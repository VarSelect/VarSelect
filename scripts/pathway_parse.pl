#!/usr/bin/perl
use strict ;
use Getopt::Std ;
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



my $opts = {} ;
getopts("d:v:o:l:j:n:", $opts) ;

my $file_geminidb = $opts->{d} ;
my $file_vcfgz_output = $opts->{o} ;
my $file_vcfgz_input = $opts->{v} ;
my $jobid = $opts->{j} ;
my $threads_active = $opts->{n} || 16 ;

my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;
my $file_config = "$dir_results_output/varselect.config" ;
my $file_log = "$dir_log/log_pathwayparse_$jobid.log" ;
my $file_pathway_output = "$dir_results_output/gemini_pathways.txt" ;
my $file_pathway_header = "$dir_work/header_pathway.txt" ;

my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("Generate VCF header for pathways data") ;
open (my $TGT_pheader , ">$file_pathway_header") ;
print $TGT_pheader "##INFO=<ID=pathway,Number=.,Type=String,Description='gemini pathways'>" ;
close $TGT_pheader ;


my $cmd_pathway = "gemini pathways $file_geminidb " ;
$log->write($cmd_pathway) ;
open (my $SRC_pathway , "$cmd_pathway |") ;

#chrom   start   end     ref     alt     impact  sample  genotype        gene    transcript      pathway
my $header = <$SRC_pathway> ;
chomp $header ;
my $header_col = [ split(/\t/,$header)] ;

my $record = {} ;
while (my $line = <$SRC_pathway>) {
    chomp $line ;

    my $i = 0 ;
    my $data = { map {$header_col->[$i++] => $_ } split (/\t/,$line) } ;

    $data->{start} ++ ; # gemini output is 0 based, convert to 1 based to match back vcf

    my $tag = join("_" , map {$data->{$_}} qw/chrom start ref alt/) ;

    my @pathways = split(/\,/ ,  $data->{pathway}) ;
    foreach my $pathway_of_transcript (@pathways) {
	$record->{$tag}->{$pathway_of_transcript} ++ ;
    }
}

close $SRC_pathway ;


my $pos_count = 0 ;
open (my $TGT_output_temp ,">$file_pathway_output") ;
foreach my $tag (sort keys %$record) {
    
    my $pathways_of_this_tag = join("," , sort keys %{$record->{$tag}}) ;
    $pos_count++ ;

    # CHROM POS REF ALT INFO/pathway
    my ($chrom,$pos,$ref,$alt) = split(/\_/,$tag) ;
    print $TGT_output_temp join("\t" , ($chrom,$pos,$ref,$alt,$pathways_of_this_tag) ) ."\n" ;

    $chrom =~ s/chr// ; # dirty hack for chromsome id
    print $TGT_output_temp join("\t" , ($chrom,$pos,$ref,$alt,$pathways_of_this_tag) ) ."\n" ;
}
$log->write ("$pos_count variants with pathways annotation",1) ;
close $TGT_output_temp ;

my $cmd_sort_bgzip = "vcf-sort $file_pathway_output |bgzip -@ $threads_active -c > $file_pathway_output.gz" ;
$log->andRun($cmd_sort_bgzip) ;

my $cmd_tabix_pahtway = "tabix -s 1 -b 2 -e 2 -f $file_pathway_output.gz " ;
$log->andRun($cmd_tabix_pahtway) ;

#my $cmd_bannot = "bcftools annotate " ;
#$cmd_bannot .= " -a $file_pathway_output.gz " ;
#$cmd_bannot .= " -c " . join ("," , qw"CHROM POS REF ALT INFO/pathway") ;
#$cmd_bannot .= " -h $file_pathway_header " ;
##'##INFO=<ID=pathway,Number=.,Type=String,Description=\"gemini pathways\">' " ;
#$cmd_bannot .= " -O z " ;
#$cmd_bannot .= " -o $file_vcfgz_output " ;
#$cmd_bannot .= " $file_vcfgz_input " ;
#$cmd_bannot .= " 2> $dir_log/stderr_bannot_pathway_$jobid.log" ; 
##$log->andRun($cmd_bannot) ;

my $cmd_vannot = "zcat $file_vcfgz_input | vcf-annotate " ;
$cmd_vannot .= " -a $file_pathway_output.gz " ;
$cmd_vannot .= " -c " . join ("," , qw"CHROM POS REF ALT INFO/pathway") ;
$cmd_vannot .= " -d key=INFO,ID=pathway,Number=.,Type=String,Description='gemini pathways' " ;
#$cmd_vannot .= " $file_vcfgz_input " ;
$cmd_vannot .= " 2> $dir_log/stderr_vannot_pathway_$jobid.log" ;
$cmd_vannot .= " |bgzip -@ $threads_active -c > $file_vcfgz_output " ;
$log->andRun($cmd_vannot) ;

my $cmd_tabix_vcfoutput = "tabix -p vcf -f $file_vcfgz_output " ;
$log->andRun($cmd_tabix_vcfoutput) ;

my $cmd_gannot = "gemini annotate " ;
$cmd_gannot .= " -f $file_vcfgz_output " ;
$cmd_gannot .= " -a extract " ;
$cmd_gannot .= " -c pathway " ;
$cmd_gannot .= " -e pathway " ;
$cmd_gannot .= " -t text " ;
$cmd_gannot .= " -o first " ;
$cmd_gannot .= " 2>$dir_log/stderr_gannot_pathway_$jobid.log  $file_geminidb " ;
#$log->andRun($cmd_gannot) ;


sub bgzip {
    my $file_src = shift ;
    my $file_tgt = shift ;

    my $cmd_bgzip = "bgzip -@ $threads_active -c $file_src > $file_tgt" ;
    $log->andRun($cmd_bgzip) ;
}



