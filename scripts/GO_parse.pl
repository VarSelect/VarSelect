#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Basename ;
use Storable ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;

my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;

my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $file_hashref_enst2go = "$dir_db/ENST2GO.hashref" ;

my $opts = {} ;
getopts("d:v:o:l:j:n:", $opts) ;
my $file_db = $opts->{d} ;
my $file_vcfgz_output = $opts->{o} ;
my $file_vcfgz_input = $opts->{v} ;
my $jobid = $opts->{j} ;
my $threads_active = $opts->{n} || 16 ;

my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;


my $file_log = "$dir_log/log_goparse_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("Running GO_parse.pl start",1) ;

my $file_go_output = "$dir_results_output/transcript_go_$jobid.txt" ;
my $file_go_header = "$dir_work/go_header_$jobid.txt" ;

$log->write("Output header to $file_go_header",1) ;

open (my $TGT_header,">$file_go_header") ;
print $TGT_header "##INFO=<ID=GO_Biological_Process,Number=.,Type=String,Description='GO_Biological_Process provided by VarSelect'>\n" ;
print $TGT_header "##INFO=<ID=GO_Cellular_Component,Number=.,Type=String,Description='GO_Cellular_Component provided by VarSelect'>\n" ;
print $TGT_header "##INFO=<ID=GO_Molecular_Function,Number=.,Type=String,Description='GO_Molecular_Function provided by VarSelect'>\n" ;
close $TGT_header ;

$log->write("Loading mapping of Ensembl transcript to database $file_hashref_enst2go start",1) ;
print "$file_hashref_enst2go\n" ;
my $enst2go = retrieve($file_hashref_enst2go) ;
$log->write("Loading mapping of Ensembl transcript to database $file_hashref_enst2go finish",1) ;


#my $cmd_gemini_transcript = "gemini query  -q 'select chrom,start+1,ref,alt,transcript from vep_variants' $file_db " ;
my $cmd_get_transcript = "sqlite3 $file_db 'select chrom,start+1,ref,alt,transcript from vep_variants' " ;
$log->write("$cmd_get_transcript") ;

open (my $SRC_transcript , "$cmd_get_transcript |") ;
open (my $TGT_go , ">$file_go_output") ;

$log->write("Extracting transcript infomation from $file_db", 1) ;
while (my $line = <$SRC_transcript>) {
    chomp $line ;
#    my ($chrom,$pos,$ref,$alt,$transcript) = split /\t/, $line ;
    my ($chrom,$pos,$ref,$alt,$transcript) = split /\|/, $line ;

    print $TGT_go "$chrom\t$pos\t$ref\t$alt\t" . join("\t" , map {$enst2go->{$transcript}->{$_}} qw/BP CC MF/) . "\n" ;

    $chrom =~ s/chr// ; # dirty hack for chromsome id
    print $TGT_go "$chrom\t$pos\t$ref\t$alt\t" . join("\t" , map {$enst2go->{$transcript}->{$_}} qw/BP CC MF/) . "\n" ;
}

close $TGT_go ;
close $SRC_transcript ;

my $cmd_sort_bgzip = "vcf-sort $file_go_output |bgzip -@ $threads_active -c > $file_go_output.gz" ;
$log->andRun($cmd_sort_bgzip) ;

my $cmd_tabix_go_output = "tabix -s 1 -b 2 -e 2 -f $file_go_output.gz " ;
$log->andRun($cmd_tabix_go_output) ;

my $cmd_bcfannotate = "bcftools annotate " ;
$cmd_bcfannotate .= " -a $file_go_output.gz " ;
$cmd_bcfannotate .= " -c " . join("," , qw"CHROM POS REF ALT INFO/GO_Biological_Process INFO/GO_Cellular_Component INFO/GO_Molecular_Function" )  ;
$cmd_bcfannotate .= " -h $file_go_header " ;
$cmd_bcfannotate .= " -O z" ; # output as vcf.gz
$cmd_bcfannotate .= " -o $file_vcfgz_output " ;
$cmd_bcfannotate .= " $file_vcfgz_input " ;
$cmd_bcfannotate .= " 2> $dir_log/stderr_bannot_go_$jobid.log " ;


#$log->andRun($cmd_bcfannotate) ; 

my $cmd_vannot = "zcat $file_vcfgz_input | vcf-annotate " ;
$cmd_vannot .= " -a $file_go_output.gz " ;
$cmd_vannot .= " -c " . join("," , qw"CHROM POS REF ALT INFO/GO_Biological_Process INFO/GO_Cellular_Component INFO/GO_Molecular_Function" )  ;
$cmd_vannot .= " -d key=INFO,ID=GO_Biological_Process,Number=.Type=String,Description='GO_Biological_Process' "  ;
$cmd_vannot .= " -d key=INFO,ID=GO_Cellular_Component,Number=.Type=String,Description='GO_Cellular_Component' "  ;
$cmd_vannot .= " -d key=INFO,ID=GO_Molecular_Function,Number=.Type=String,Description='GO_Molecular_Function' "  ;
$cmd_vannot .= " 2> $dir_log/stderr_vannot_go_$jobid.log " ;
$cmd_vannot .= " | bgzip -c > $file_vcfgz_output " ;

$log->andRun($cmd_vannot) ;

my $cmd_tabix_vcfoutput = "tabix -p vcf -f $file_vcfgz_output " ;
$log->andRun($cmd_tabix_vcfoutput) ;

#my $cmd_gannot = "gemini annotate " ;
#$cmd_gannot .= " -f $file_vcfgz_output " ;
#$cmd_gannot .= " -a extract " ;
#$cmd_gannot .= " -c GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function " ;
#$cmd_gannot .= " -e GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function " ;
#$cmd_gannot .= " -t text,text,text " ;
#$cmd_gannot .= " -o first,first,first " ;
#$cmd_gannot .= " 2> $dir_log/stderr_gannot_go_$jobid.log $file_db " ;

#$log->andRun($cmd_gannot) ;


$log->write("Running GO_parse.pl finish",1) ;


