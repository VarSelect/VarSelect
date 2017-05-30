#!/usr/bin/perl
use strict ;

use Getopt::Std ;
use File::Basename ;
use Storable ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;

use DBI ;

# Default Global Setting
#=======================================
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $script_vslannot = "$dir_script/vsl_annot.pl" ;

my $method_compare_avail= {
    1 => 'A or B' ,
    2 => 'A and B' ,
    3 => 'A only' ,
    4 => 'B only' ,
} ;

# Options handle
#=======================================
my $opts = {} ;
getopts("a:b:c:d:j:n:q:",$opts) ;

my $jobid	    = $opts->{j} ;
my $jobidA	    = $opts->{a} ;
#|| die "\n-a is required\n\n" ;
my $jobidB	    = $opts->{b} ;
#|| die "\n-b is required\n\n" ;
my $jobids_query    = $opts->{q} ;
my $method_compare  = $opts->{c} ;
my $file_db	    = $opts->{d} || die "\n-d is required\n\n" ;
my $threads_enable  = $opts->{n} || 16 ;

die "-c should be 1-4!\n"  unless (exists $method_compare_avail->{$method_compare}) ;

# Get all jobinfo of db file
#=======================================
my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","") ;
my $table_varselectinfo = 'varselect_info' ;
my $table_jobinfo = 'jobinfo' ;

my $sql_selall_jobs = "select * from jobinfo" ;
my $sth_selall_jobs = $dbh->prepare("$sql_selall_jobs") ;
$sth_selall_jobs->execute() ;

my $joblist = [] ;
my $jobs = {} ;

while (my $job = $sth_selall_jobs->fetchrow_hashref) {
    # jobid ,  ped , multicaller , workflow

    my $jobid = $job->{jobid} ;
    $jobs->{$jobid} = $job ;

#    display_jobs($job) ;
}

#die "ERROR: Job ID $jobidA doesn't exist in $file_db!" unless (exists $jobs->{$jobidA}) ;
#die "ERROR: Job ID $jobidB doesn't exist in $file_db!" unless (exists $jobs->{$jobidB}) ;


# Jobs specific path
#=======================================
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $file_log = "$dir_log/log_vscompare_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

my ($dbfn,$dbpath,$dbext) = fileparse($file_db , qw/.db/) ;
my $file_vcf_finaloutput = "VarSelect_compare_$jobid.vcf" ;
my $file_txt_finaloutput = "VarSelect_compare_$jobid.txt" ;

my $iaA = "in_analysis_$jobidA" ;
my $iaB = "in_analysis_$jobidB" ;
my $jobids = [] ;

if ($jobids_query) {
    $jobids = [split(/\,/,$jobids_query)] ;
} else {
    push @$jobids, $jobidA ;
    push @$jobids, $jobidB ;
}
    

my $compare_mode = '' ;
my $query_sql = "select chrom,start+1,ref,alt,1 from variants where " ;
my $query_sql_dump = "select * from variants where " ;

if ($method_compare == 1) { 
    # A or B
    my $query_where_text = join (" or " , map {"in_analysis_$_ = \"1\""} @$jobids) ;

    $query_sql .= $query_where_text ;
    $query_sql_dump .= $query_where_text ;
    $compare_mode = 'union' ;

} elsif ($method_compare == 2) { 
    # A and B
    my $query_where_text = join (" and " , map {"in_analysis_$_ = \"1\""} @$jobids) ;
    $query_sql .= $query_where_text ;
    $query_sql_dump .= $query_where_text ;
    $compare_mode = 'intersection' ;

} elsif ($method_compare == 3) { 
    # A only
    $query_sql .= " $iaA = \"1\" and ($iaB = \"0\" or $iaB is NULL)" ;
    $query_sql_dump .= " $iaA = \"1\" and ($iaB = \"0\" or $iaB is NULL)" ;
    $compare_mode = 'Aonly' ;

} elsif ($method_compare == 4) { 
    # B only
    $query_sql .= " ($iaA = \"0\" or $iaA is NULL) and $iaB = \"1\"" ;
    $query_sql_dump .= " ($iaA = \"0\" or $iaA is NULL) and $iaB = \"1\"" ;
    $compare_mode = 'Bonly' ;

} else {  # should not be here
}
    
print "Job id : $jobid\n" ;

my $result_compare_tab_annot = "$dir_results_output/vs_compare_output_annot_$jobid.txt" ;
#my $result_compare_tab = "$dir_results_output/vs_compare_output_$jobid.txt" ;
my $result_compare_vcf = "$dir_results_output/vs_compare_output_$jobid.vcf" ;

gemini_query     ("$query_sql"	, $file_db , "$result_compare_tab_annot" ,0) ;
#gemini_query     ("$query_sql_dump" , $file_db , "$result_compare_tab" ,1) ;
gemini_query_vcf ("$query_sql_dump" , $file_db , $result_compare_vcf , 1) ;

#    my $cmd_bgzip = "bgzip $result_compare_vcf -c > $result_compare_vcf.gz"  ;
#    tabix_vcf($result_compare_tab ,$log) ;
tabix_vcf($result_compare_tab_annot ,$log) ;
tabix_vcf($result_compare_vcf ,$log) ;

my $file_vcfgz_outputtmp = "$dir_work/compare_output_tmp.vcf.gz" ;


my $cmd_vannot = "cat $result_compare_vcf | vcf-annotate " ;
$cmd_vannot .= " -a $result_compare_tab_annot.gz " ;
$cmd_vannot .= " -c " . join ("," , (qw/CHROM POS REF ALT/ , "INFO/in_analysis_$jobid") ) ;
$cmd_vannot .= " -d  key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='Secondary analysis $jobid' " ;
$cmd_vannot .= " 2> $dir_log/stderr_vannot_compare_$jobid.log " ;
$cmd_vannot .= " | bgzip -@ $threads_enable -c > $file_vcfgz_outputtmp " ;

$log->andRun($cmd_vannot , 1) ;
tabix_vcf ($file_vcfgz_outputtmp ,$log) ;

my $cmd_vsl_annot =  "$script_vslannot" ;
$cmd_vsl_annot .= " -f $file_vcfgz_outputtmp " ;
$cmd_vsl_annot .= " -e in_analysis_$jobid " ;
$cmd_vsl_annot .= " -t boolean " ;
$cmd_vsl_annot .= " -c in_analysis_$jobid " ;
$cmd_vsl_annot .= " -d $file_db " ;
$cmd_vsl_annot .= " -b varselect_variants " ;
$cmd_vsl_annot .= " 2> $dir_log/stderr_vslannot_compare_$jobid.log " ;

$log->andRun($cmd_vsl_annot,1) ;

# CREATE VIEW
#=======================================
my $cmd_createview = "$dir_script/create_view_for_varselect.pl -j $jobid -d $file_db" ;
$log->andRun($cmd_createview, 1) ;

$log->write("VarSelect compare finish") ;


sub gemini_query_vcf {
    my $sql = shift ;
    my $db = shift ;
    my $file_output = shift ;
    my $header = shift || 0 ;

    my $cmd_gemini_query = "gemini query " ;
    $cmd_gemini_query .= " --format vcf " ;
    $cmd_gemini_query .= " --header " if ($header) ;
    $cmd_gemini_query .= " -q '$sql' $db |sed 's/^chrM/MT/' |sed 's/^chr//' > $file_output" ;
    $log->andRun($cmd_gemini_query,1) ;
}


sub gemini_query {
    my $sql = shift ;
    my $db = shift ;
    my $file_output = shift ;
    my $header = shift || 0 ;

    my $cmd_gemini_query = "gemini query " ;
    $cmd_gemini_query .= " --header " if ($header) ;
    $cmd_gemini_query .= " -q '$sql' $db |sed 's/^chrM/MT/' |sed 's/^chr//'> $file_output" ;
    $log->andRun($cmd_gemini_query,1) ;
}


sub load_config {
    my $file = shift ;
    die "$file isn't exist! please check it's a varselect_analysis_dir or not" unless (-e $file) ;

    my $output = {} ;
    open (my $SRC , "$file") ;
    while (<$SRC>) {
        chomp ;
        my ($key,$val) = split (/\t/ , $_ , 2) ;
        if (exists $output->{$key}) {
            $output->{$key} .= "\n$val" ;
        } else {
            $output->{$key} = $val ;
        }
    }
    close $SRC ;

    return $output ;
}
    

sub display_jobs {
    my $job = shift ;

    $job->{ped} =~ s/\n/\n\t/g ;

    print "Jobid: $job->{jobid}\n" ;
    print "\tWorkflow: $job->{workflow}\n" ;
    print "\tMultiCaller: $job->{multicaller}\n" ;
    print "\tPED:\n\t$job->{ped}\n" ;
}

sub CHROM_nochr {
    my $chr_name = shift ;

    if ($chr_name eq 'chrM') {
        return 'MT' ;
    } elsif ($chr_name =~ /^chr/) {
        return substr($chr_name,3) ;
    } else {
        return $chr_name ;
    }
}

