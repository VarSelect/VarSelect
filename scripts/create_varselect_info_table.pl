#!/usr/bin/perl
#
# a script to create varselect_info Table
use strict ;
use Getopt::Std ;
use DBI qw(:sql_types) ;

my $table_varselectinfo = 'varselect_info' ;


my $opts={} ;
getopts("d:r:v:m:",$opts) ;

my $file_db		    = $opts->{d} || die "-d file_db is required!" ;
my $file_vcfgz_result	    = $opts->{r} || die "-r result VCF of VarSelect for analysis is required" ;
my $file_list_vcfinput	    = $opts->{v} || die "-v vcf file list for vs_annotate is required" ;
my $flag_multiple_caller    = $opts->{m} || 0 ;

# Load bgzipped VCF blob from file
my ($vcf_blob, $buff) ;
open (my $SRC_vcfinput,"$file_vcfgz_result") ;
while (read $SRC_vcfinput, $buff , 1024) {
    $vcf_blob .=  $buff ;
}
close $SRC_vcfinput ;

# Load VCF file list
my $sample_list = [] ;
my $caller_list = [] ;
open (my $SRC_filelist , "$file_list_vcfinput") ;
while (my $line = <$SRC_filelist>) {
    chomp $line ;
    my ($sample, $file_vcf, $caller) = split (/\,/ , $line) ;

    push @$sample_list , $sample ;
    push @$caller_list , $caller ;

    print "Load sample list: $sample ......\n" if $sample ;
    print "Load caller list: $caller ......\n" if $caller ;
}
close $SRC_filelist ;

# Connect SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","") ;

# CREATE TABLE varselect_info
my $sql_varsel_info_create = "CREATE TABLE IF NOT EXISTS $table_varselectinfo ( " ;
$sql_varsel_info_create .= " sn INTEGER " ;
$sql_varsel_info_create .= " , vcfgz_for_analysis BLOB " ; #vcf_input
$sql_varsel_info_create .= " , sample_list TEXT " ;
$sql_varsel_info_create .= " , caller_list TEXT " ;
$sql_varsel_info_create .= " , multicaller BOOLEAN " ;
$sql_varsel_info_create .= ")" ;

$dbh->do($sql_varsel_info_create) ;

print "$sql_varsel_info_create\n" ;


# INSERT
my $sql_varsel_info_insert = "INSERT INTO $table_varselectinfo " ;
$sql_varsel_info_insert .= " (sn , vcfgz_for_analysis , sample_list, caller_list, multicaller) values " ;
$sql_varsel_info_insert .= " (0, ?, ?, ?, ?) " ;

my $sth_varsel_info_insert = $dbh->prepare($sql_varsel_info_insert) ;

$sth_varsel_info_insert->bind_param(1, $vcf_blob, SQL_BLOB) ;
$sth_varsel_info_insert->bind_param(2, join(",",@$sample_list) ) ;
$sth_varsel_info_insert->bind_param(3, join(",",@$caller_list) ) ;
$sth_varsel_info_insert->bind_param(4, $flag_multiple_caller) ;

$sth_varsel_info_insert->execute() ;

print "$sql_varsel_info_insert\n" ;
