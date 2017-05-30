#!/usr/bin/perl

# the script is to drop and re-create view variants in VarSelect.db

use strict ;
use Getopt::Std ;
use File::Basename ;
use Cwd qw(abs_path) ;
use lib dirname(dirname( abs_path($0) )) . '/lib' ;
use VarSelect ;
use VarSelect::Log ;

my $debug_mode = 0 ;

# Default Global Setting
#=======================================
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db        = "$dir_vsHOME/db" ;
my $dir_lib       = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;


# Options handle
#=======================================
my $opts = {} ;
getopts("j:d:",$opts) ;

my $jobid      = $opts->{j} ;
my $file_db    = $opts->{d} ;


# Jobs specific path
#=======================================
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log  = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $ts = time ;
my $file_log = "$dir_log/log_create_view_$jobid.log.$ts" ;
my $log = VarSelect::Log->new(file => $file_log) ;

# VarSelect specific settings:
#=======================================
my $table_vep       = 'vep_variants' ;
my $table_snpeff    = 'snpeff_variants' ;
my $table_annovar   = 'annovar_variants' ;
my $table_varselect = 'varselect_variants' ;


# CREATE VIEW
#=======================================
my $vep_collist         = get_col_list($table_vep,$file_db) ;
#my $snpeff_collist      = get_col_list($table_snpeff,$file_db) ;
my $annovar_collist     = get_col_list($table_annovar,$file_db,[qw/variant_id chrom start end ref alt/]) ;
my $varselect_collist   = get_col_list($table_varselect,$file_db,[qw/variant_id chrom start end ref alt/]) ;

my $sql_create_view = "CREATE VIEW variants AS SELECT " ;
$sql_create_view .= join (" , " , map {"v.'$_' as '$_'"} @$vep_collist) ;
#$sql_create_view .= ", " . join (" , " , map {(/^gt/)?"s.'$_' as 'gt_snpeff_$_'":"s.'$_' as 'snpeff_$_'"} @$snpeff_collist) ;
$sql_create_view .= ", " . join (" , " , map {"a.'$_' as 'anv_$_'" } @$annovar_collist) ;
$sql_create_view .= ", " . join (" , " , map {"l.'$_' as '$_'"} @$varselect_collist ) ;
$sql_create_view .= " FROM $table_vep v " ;
#$sql_create_view .= " INNER JOIN $table_snpeff s ON s.variant_id = v.variant_id " ;
$sql_create_view .= " INNER JOIN $table_annovar a ON a.variant_id = v.variant_id " ;
$sql_create_view .= " INNER JOIN $table_varselect l ON l.variant_id = v.variant_id " ;

my $cmd_drop_view_first = "sqlite3 $file_db 'DROP VIEW variants' 2> $dir_log/stderr.dropview.log_$ts" ;
$log->andRun($cmd_drop_view_first, $debug_mode) ;

my $cmd_create_view = "sqlite3 $file_db \"$sql_create_view\" 2> $dir_log/stderr.createview.log_$ts"  ;
$log->andRun($cmd_create_view, $debug_mode) ;

sub get_col_list {
    my $table = shift ;
    my $db = shift ;
    my $drop_list = shift ;
    my $output = [] ;

    my $in_droplist = {map {$_ => 1} @$drop_list} ;

    my $sql_PRAGMA = "PRAGMA table_info($table) " ;

    open (my $SRC_tableinfo, "sqlite3 $db '$sql_PRAGMA' |") ;
    while (<$SRC_tableinfo>) {
        chomp ;
        my ($sn,$col_name,$col_type,$null_state,$default,$pk) = split(/\|/,$_) ;
        push @$output , $col_name unless (exists $in_droplist->{$col_name}) ; # ignore column listed in $drop_list
    }

    close $SRC_tableinfo ;

    return $output ;
}

