#!/usr/bin/perl
use strict ;
use Getopt::Std ;
use File::Basename ;
use File::Path qw(make_path remove_tree);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use DBI ;

my $VERSION = "0.9.6 (20170524)-dev" ;
my $program_ver = "VarSelect v$VERSION" ;

# should be put in global config
my $DEFAULT_threads = 16 ;

my $jobid = getts() ;

# Commands handle
my $commands_avail = {map {$_ => 1} qw/annotate analysis compare list/}  ;
my $commands_desc = {
    annotate => 'Annotate VCFs from scratch and create VarSelect db' ,
    analysis => 'Original analysis from VarSelect db' ,
    compare  => 'Comparative analysis between two or multiple original analyses' ,
    list  => 'Output analysis list in VarSelect db' ,
} ;

my $command = shift @ARGV ;
die usage_all() unless $command ;
die usage_all() . "\"$command\" is not a valid VarSelect command!\n\n" unless exists $commands_avail->{$command} ;

my $workflow_avail = {map {$_ => 1} qw/family paired none/} ;

# VarSelect Path settings
my $dir_vsHOME = dirname(abs_path($0)) ;
my $dir_db	  = "$dir_vsHOME/db" ;
my $dir_lib	  = "$dir_vsHOME/lib" ;
my $dir_scripts   = "$dir_vsHOME/scripts" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $script_vsannotate = "$dir_scripts/vs_annotate.pl" ;
my $script_vsanalysis = "$dir_scripts/vs_analysis.pl" ;
my $script_vscompare  = "$dir_scripts/vs_compare.pl" ;

my $opts = {} ;
if ($command eq 'annotate') {

    # Options for annotate
    getopts("v:p:m:c:x:hkuin:" , $opts) ;

    die usage_annotate()  if ($opts->{h}) ;
    my $file_vcflist	    = $opts->{v} || die "-v vcf_file_list is required!\n\n" . usage_annotate() ;
    my $file_ped	    = $opts->{p} || die "-p ped file is required!\n\n" . usage_annotate() ;
    my $mode_workflows	    = $opts->{m} || die "-m workflow mode is required, [paired, family, none] !\n\n" . usage_annotate()  ;
    my $flag_multicallers   = $opts->{k} ;
    my $threads_enable	    = $opts->{n} || $DEFAULT_threads ;
    my $file_cnv	    = $opts->{c} ;
    my $file_xpr	    = $opts->{x} ;

    my $flag_combine_union  = $opts->{u} ;
    my $flag_combine_intersect = $opts->{i} ;

    my $type_combine = "union" ;
    if ($flag_combine_union && $flag_combine_intersect) {
	die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
    } elsif ($flag_combine_intersect) {
	$type_combine = "intersection" ;
    }

    die "Only [ " . join (",", sort keys %$workflow_avail) . " ] are available for option -m \n" unless (exists $workflow_avail->{$mode_workflows}) ;

    # Output setting
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    my $file_stderr  = "$dir_log/stderr_vs_annotate.$jobid.log" ;
    my $file_stderr2 = "$dir_log/stderr_vs_analysis.$jobid.log" ;
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("VarSelect Start",1) ;
    $log->write("Version: $VERSION",1) ;
    $log->write("Jobid: $jobid",1) ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    my ($fname_vfile,$fpath_vfile,$fext_vfile) = fileparse($file_vcflist , qw/.txt/) ;
    my $file_db = "./$fname_vfile\_$jobid\_varselect.db" ;


    print "VarSelect Job id: $jobid\n" ;

    # Step 1. Annotate from vcf
    my $cmd_vs_annotate = "$script_vsannotate " ;
    $cmd_vs_annotate .= " -v $file_vcflist " ;
    $cmd_vs_annotate .= " -p $file_ped " ;
    $cmd_vs_annotate .= " -j $jobid " ;
    $cmd_vs_annotate .= " -n $threads_enable " ;
    $cmd_vs_annotate .= " -k " if ($flag_multicallers) ;
    $cmd_vs_annotate .= " -d $file_db " ; # make sure where to put db file
    $cmd_vs_annotate .= " 2>&1 | tee $file_stderr " ;

    $log->write("VarSelect annotate start" , 1) ;
    $log->andRun($cmd_vs_annotate) ;
    $log->loadlog("$dir_log/log_vsannotate_$jobid.log") ;
    $log->write("VarSelect annotate finish", 1) ;

    # Step 2. Analysis according to PED info
    my $cmd_vs_analysis = "$script_vsanalysis " ; 
    $cmd_vs_analysis .= " -j $jobid " ;
    $cmd_vs_analysis .= " -d $file_db " ;
    $cmd_vs_analysis .= " -p $file_ped " ;
    $cmd_vs_analysis .= " -m $mode_workflows " ;
    $cmd_vs_analysis .= " -k " if ($flag_multicallers) ;
    $cmd_vs_analysis .= " -i " if ($flag_combine_intersect) ;
    $cmd_vs_analysis .= " -u " if ($flag_combine_union) ;
    $cmd_vs_analysis .= " -c $file_cnv " if ($file_cnv) ;
    $cmd_vs_analysis .= " -x $file_xpr " if ($file_xpr) ;
    $cmd_vs_analysis .= " 2>&1 |tee $file_stderr2 " ;

    $log->write("VarSelect analysis start" , 1) ;

    $log->andRun($cmd_vs_analysis) ;
    $log->loadlog("$dir_log/log_vsanalysis_$jobid.log") ;
    $log->write("VarSelect analysis finish", 1) ;
    print "VarSelect job $jobid finish\n" ;

} elsif ($command eq 'analysis') {
    # Options for analysis
    getopts("d:p:m:hkuic:x:n:",$opts) ;

    die usage_analysis() if $opts->{h} ;

    my $file_db		= $opts->{d} || die "\n\n-d is required!\n" . usage_analysis() ;
    my $file_ped	= $opts->{p} || die "\n\n-p is required!\n" . usage_analysis() ;
    my $mode_workflows	= $opts->{m} || die "\n\n-m is required!\n" . usage_analysis() ;

    my $flag_multicallers   = $opts->{k} ;
    my $threads_enable	    = $opts->{n} || $DEFAULT_threads ;
    my $file_cnv	    = $opts->{c} ;
    my $file_xpr	    = $opts->{x} ;
    my $freq_vgt_filter	    = $opts->{f} || 1;

    my $flag_combine_union  = $opts->{u} ;
    my $flag_combine_intersect = $opts->{i} ;

    my $type_combine = "union" ;
    if ($flag_combine_union && $flag_combine_intersect) {
	die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
    } elsif ($flag_combine_intersect) {
	$type_combine = "intersection" ;
    }

    die "\ndb $file_db doesn't exist!\n" unless (-e $file_db) ;
    die "\ndb $file_db is not writable!\n" unless (-w $file_db) ;

    die "Only [ " . join ("," , sort keys %$workflow_avail) . " ] are available for option -m \n" unless (exists $workflow_avail->{$mode_workflows}) ;

    print "VarSelect Job id: $jobid\n" ;

    # Output setting, create output directory
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    print "$dir_results_output is created.\n" ;

    my $file_stderr = "$dir_log/stderr.vs_analysis.$jobid.log" ;

    # Start logging
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("Version: $VERSION",1) ;
    $log->write("Jobid: $jobid",1) ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    # pipeline setting
    my $cmd_vs_analysis = "$script_vsanalysis " ; 
    $cmd_vs_analysis .= " -j $jobid " ;
    $cmd_vs_analysis .= " -d $file_db " ;
    $cmd_vs_analysis .= " -p $file_ped " ;
    $cmd_vs_analysis .= " -m $mode_workflows " ;
    $cmd_vs_analysis .= " -k " if ($flag_multicallers) ;
    $cmd_vs_analysis .= " -i " if ($flag_combine_intersect) ;
    $cmd_vs_analysis .= " -u " if ($flag_combine_union) ;
    $cmd_vs_analysis .= " -c $file_cnv " if ($file_cnv) ;
    $cmd_vs_analysis .= " -x $file_xpr " if ($file_xpr) ;

    $log->write("VS_analysis start" , 1) ;
    $log->andRun($cmd_vs_analysis) ;
    $log->loadlog("$dir_log/log_vsanalysis_$jobid.log") ;
    $log->write("VS_analysis finish", 1) ;

    print "VarSelect job $jobid finish\n" ;
    print localtime() ."\n" ;

} elsif ($command eq 'compare') {
    # options for compare
    getopts("a:b:c:d:q:h",$opts) ;

    die usage_compare() if ($opts->{h}) ;
    my $datasetA	= $opts->{a} ; # || die "\n-a is required!\n" . usage_compare() ;
    my $datasetB	= $opts->{b} ; # || die "\n-b is required!\n" . usage_compare() ;
    my $method_compare	= $opts->{c} || die "\n-c is required!\n" . usage_compare() ;
    my $file_db		= $opts->{d} || die "\n-d is required!\n" . usage_compare() ;
    my $jobids_query	= $opts->{q} ; # comma seprated jobs id 

    die "\ndb $file_db doesn't exist!\n" unless (-e $file_db) ;
    die "\ndb $file_db is not writable!\n" unless (-w $file_db) ;

    if ($jobids_query) {
	# -q enabled
	die "\n-c 3 or -c 4 is for two jobs comparison only.\n" . usage_compare() if ($method_compare eq '3' || $method_compare eq '4') ;

    } else {
	# check if -a -b
	die "\n-a is required!\n" . usage_compare() unless($datasetA) ;
	die "\n-b is required!\n" . usage_compare() unless($datasetB) ;
    }

    # Output setting
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    my $file_stderr = "$dir_log/stderr_vs_compare.$jobid.log" ;
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("VarSelect Start") ;
    $log->write("Jobid: $jobid") ;
    $log->write("VarSelect compare start",1) ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    my $cmd_vs_compare = "$script_vscompare " ;
    $cmd_vs_compare .= " -a $datasetA " if ($datasetA) ;
    $cmd_vs_compare .= " -b $datasetB " if ($datasetB) ;
    $cmd_vs_compare .= " -c $method_compare " ;
    $cmd_vs_compare .= " -d $file_db " ;
    $cmd_vs_compare .= " -q $jobids_query" if ($jobids_query) ;
    $cmd_vs_compare .= " -j $jobid " ;
    
    $log->andRun($cmd_vs_compare) ;
    $log->write("VarSelect compare finish",1 ) ;

} elsif ($command eq 'list') {
    getopts("d:h",$opts) ;
    die usage_joblist() if ($opts->{h}) ;
    my $file_db = $opts->{d} || die "-d file_db is required!\n" ;
    die "\ndb $file_db doesn't exist!\n" unless (-e $file_db) ;

    my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","") ;
    my $table_jobinfo = 'jobinfo' ;

    my $sql_selall_jobs = "select * from $table_jobinfo" ;
    my $sth_selall_jobs = $dbh->prepare("$sql_selall_jobs") ;
    $sth_selall_jobs->execute() ;

    my $joblist = [] ;
    my $jobs = {} ;

    while (my $job = $sth_selall_jobs->fetchrow_hashref) {
	# jobid ,  ped , multicaller , workflow

	my $jobid = $job->{jobid} ;
	$jobs->{$jobid} = $job ;

        display_jobs($job) ;
    }

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

    # PED
    # ID
    # columns
    return $output ;
}
    

sub usage_all {
return <<EOusage

$program_ver
Usage: $0 <command> [options]

Commands:
    annotate - $commands_desc->{annotate}
    analysis - $commands_desc->{analysis}
    compare  - $commands_desc->{compare} 
    list     - $commands_desc->{list}

use $0 <command> -h to for more help information

EOusage
}

sub usage_annotate {
    my $text_available_mode = join (", " , keys %$workflow_avail) ;
    return <<EOusage1

$program_ver 
Command: annotate - $commands_desc->{annotate} 

Usage: $0 annotate -v <VCF files list> -p <ped file> -m <workflow mode> [-k] [-u] [-i] [-c <cns files list>] [-x <expression files list>]

Options:
    -v VCF files list					 (required)
    -p PED file						 (required)
    -m workflow mode: [$text_available_mode]		 (required)

    -k multiple callers mode				 (optional)
    -u get union of variants between callers		 (work only with -k option)
    -i get intersection of variants between callers	 (work only with -k option)

    -c CNVkit log2 calls cns file list			 (optional)
    -x expression profile list				 (optional)

    -n number of threads (default:$DEFAULT_threads)			 (optional)
    -h this help page

EOusage1
}

sub usage_analysis {
    return <<EOusage2

$program_ver
Command: analysis - $commands_desc->{analysis} 

Usage: $0 analysis -d <VarSelect db file> -p <ped file> -m <workflow mode> [-k] [-u] [-i] [-c <cns files list>] [-x <expression files list>]

Options:
    -d VarSelect db file		    		 (required)
    -p PED file						 (required)
    -m workflow mode: [paired, family, none]	    	 (required)

    -k multiple callers mode				 (optional)
    -u get union of variants between callers		 (work only with -k option)
    -i get intersection of variants between callers	 (work only with -k option)

    -c CNVkit log2 calls cns file list			 (optional)
    -x expression profile list				 (optional)

    -n number of threads (default:$DEFAULT_threads)			 (optional)
    -h this help page

EOusage2
}

sub usage_compare {
    return <<EOusage3

$program_ver
Command: compare - $commands_desc->{compare}

Usage: $0 compare -a <analysis id of dataset A> -b <analysis id of dataset B> -c <compare mode> -d <VarSelect db file>
       $0 compare -q <analysis id,analysis id,analysis id...> -c <compare mode> -d <VarSelect db file>

Options:
    -a Analysis id of dataset A		    (for comparison between two original analyses)
    -b Analysis id of dataset B		    (for comparison between two original analyses)    
    -q Comma separated multiple analysis id (for comparison between multiple original analyses)

    -c compare mode [1-4]		    (required)
       1: A or B
       2: A and B
       3: A only			    (valid only for comparison between two original analyses) 
       4: B only			    (valid only for comparison between two original analyses)

    -d VarSelect db file		    (required)
    -h this help page

EOusage3
}

sub usage_joblist {
    return <<EOusage4

$program_ver
Command: list - $commands_desc->{list}

Usage: $0 list -d <VarSelect db file>

Options:
    -d VarSelect db file
    -h this help page

EOusage4
}


sub display_jobs {
    my $job = shift ;

    $job->{ped} =~ s/\n/\n\t/g ;

    print "Analysis id: $job->{jobid}\n" ;
    print "\tWorkflow: $job->{workflow}\n" ;
    print "\tMultiCaller: $job->{multicaller}\n" ;
    print "\tPED:\n\t$job->{ped}\n" ;
}

