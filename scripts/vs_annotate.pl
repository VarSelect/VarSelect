#!/usr/bin/perl
use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;

use IO::Handle ;
autoflush STDERR 1 ;
autoflush STDOUT 1 ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db	  = "$dir_vsHOME/db" ;
my $dir_lib	  = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
my $opts = {} ;
getopts("v:p:j:n:kd:",$opts) ;

my $file_vcf_list	= $opts->{v} ;
my $file_ped		= $opts->{p} ;
my $jobid		= $opts->{j} ;
my $flag_multicallers   = $opts->{k} || 0 ;
my $file_db		= $opts->{d} ;
my $threads_enable      = $opts->{n} ;
my $threads_active = $threads_enable ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log  = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $file_config = "$dir_results_output/varselect.config" ;
my $file_log = "$dir_log/log_vsannotate_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("VS annot start",1) ;
$log->write("Jobid: $jobid",1) ;
$log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ,1) ;


# zz: should be stored in table:varselect_jobs in the furure
open (my $CONFIG , ">$file_config") ;
print $CONFIG "ID\t$jobid\n" ;

# zz: shold be remove
my $file_hashref_vannot_tmp = "$dir_work/vannot_tmp.hashref" ;

my ($fname_db,$fpath_db,$fext_db) = fileparse($file_db , qw/.db/) ;
my $file_db_prefix = "$fpath_db/$fname_db" ;

# zz: should be stored in table:varselect_jobs in the furure
my $file_vcfgz_foranalysis = "$file_db_prefix.vcf.gz" ;
my $file_config_db = "$file_db_prefix.config" ;

my $file_prefix = "$dir_results_output/$fname_db" ;
my $work_prefix = "$dir_work/$fname_db" ;
my ($file_prefix0 , $file_path , $file_ext) = fileparse($file_vcf_list , qw/.txt .csv/) ;

# PreProcess 01: Load vcf files of samples and get sample/caller list
#=======================================
my $sample2file = {} ;
my $caller_list = {} ;
my $sample_list = {} ;
$log->write ("Load and check VCF file list start", 1) ;
open (my $SRC_filelist , "$file_vcf_list") ;
while (my $line = <$SRC_filelist>) {
    chomp $line ;

    my ($sample, $caller, $file_vcf) ;
    my $sample_tag ;

    if ($flag_multicallers) {	# multi-caller mode
        ($sample, $file_vcf, $caller) = split (/\,/ , $line) ;
        die "It is multi-caller mode\n Please provide \"Sample,FileLocation,VariantCaller\" for each vcf file." unless (defined $caller) ;

        $caller_list->{$caller} ++ ;
        $sample_tag = "$sample\_$caller" ;

    } else {			# single-caller mode
        ($sample, $file_vcf , $caller) = split (/\,/ , $line) ;
	die "It is single-caller mode\nWe detect 3rd column in input file $file_vcf_list. Please provide ONLY 2 columns -- \"Sample,FileLocation\" -- in input file." if ($caller) ;

        $sample_tag = $sample ;
    }

    die ("VCF file $file_vcf is not exists.\n") unless (-e $file_vcf) ;
    $sample_list->{$sample} ++ ;

    my ($filevcf_fn , $filevcf_path, $filevcf_ext) = fileparse($file_vcf , qw/.vcf.gz .vcf/) ;

    # make sure every input vcf is bgziped and tabixed
    my $file_vcfbgz_input = "$dir_work/$filevcf_fn.vcf.gz" ;
    $log->write("Clean variants with multiple line or wrong sample format values: $file_vcf" , 1) ;
    my $cmd_force_bgzipped_input_vcf = "$dir_script/drop_multiline.pl -j $jobid -i $file_vcf |bgzip -@ $threads_enable -c > $file_vcfbgz_input " ;
    $log->andRun($cmd_force_bgzipped_input_vcf) ;
    tabix_vcf($file_vcfbgz_input) ;

    push @{$sample2file->{$sample_tag}} , $file_vcfbgz_input ;
}
close $SRC_filelist ;
$log->write ("Load and check VCF file list finish" , 1) ;

print $CONFIG "multicaller\t$flag_multicallers\n" ;
print $CONFIG "caller_list\t" . join("," , sort keys %$caller_list) . "\n" if ($flag_multicallers) ;
print $CONFIG "sample_list\t" . join("," , sort keys %$sample_list) . "\n" ;

# PreProcess 02: Load PED file, output Sample sex, and  get sample list
#=======================================
my $file_sample_sex = $file_prefix . "_sample_sex.txt" ;
my $file_sample_sex_mcaller = $file_prefix . "_sample_sex_mcaller.txt" ;
my $list_samples = {} ;

my ($ped_fn,$ped_path,$ped_ext) = fileparse($file_ped,qw/.ped/) ;
my $file_ped_mcaller =  "$dir_work/$ped_fn\_mcaller.ped" ;
open (my $TGT_pedmcaller,">$file_ped_mcaller") if ($flag_multicallers) ;

$log->write ("Load ped file start",1) ;
open (my $TGT_samplesex , ">$file_sample_sex") ;
open (my $TGT_samplesex_mcaller , ">$file_sample_sex_mcaller") if ($flag_multicallers) ;

open (my $SRC_ped , "$file_ped") ;
while (<$SRC_ped>) {
    next if /^#/ ;  # skip line start with #
    chomp ;
    my ($family,$sample,$fatherid,$motherid,$sexid,$affectid,$other) = split /\t/ ;

    my $sex = '' ;
    if ($sexid == 1) {
        $sex = 'M' ;
    } elsif ($sexid == 2) {
        $sex = 'F' ;
    } else {
        die "PED file $file_ped with wrong sex id $sexid\n\n" ;
    }

    print $TGT_samplesex "$sample $sex\n" ;
    if ($flag_multicallers) {
	# alternate ped file
        foreach my $caller (sort keys %$caller_list) {
            my $sample_tag = $sample . "_" . $caller ;
            die "$sample_tag doesn't have related vcf file, please check." unless (exists $sample2file->{$sample_tag}) ;
            print $TGT_pedmcaller join("\t" , ($family,$sample_tag,($fatherid)?"$fatherid\_$caller":$fatherid , ($motherid)?"$motherid\_$caller":$motherid , $sexid, $affectid, $other)) ;
            print $TGT_pedmcaller "\n" ;
	    # alternate sex file
	    print $TGT_samplesex_mcaller "$sample_tag $sex\n" ;
            $list_samples->{$sample_tag} ++ ;

            # reheader by bcftools sample => sample_caller
	    my $file_new_header = "$dir_work/newheader_$sample_tag.txt" ;
            open (my $TGT_nh,">$file_new_header") ;
            print $TGT_nh "$sample $sample_tag\n" ;
            close $TGT_nh ;

            foreach my $old_vcfgz (@{$sample2file->{$sample_tag}}) {
                my $file_new_vcf = "$dir_work/$sample_tag\_nh.vcf" ;
                my $cmd_bcf_reheader = "bcftools reheader --sample $file_new_header $old_vcfgz |bcftools view > $file_new_vcf" ;
                $log->andRun($cmd_bcf_reheader) ;

                my $file_orig = $old_vcfgz . ".orig" ;
                my $cmd_mvorig = "mv $old_vcfgz $file_orig" ;
                $log->andRun($cmd_mvorig) ;
                bgzip($file_new_vcf , "$file_new_vcf.gz") ;
                my $cmd_mvnewone = "mv $file_new_vcf.gz $old_vcfgz " ;
                $log->andRun($cmd_mvnewone) ;
            }
        }

    } else {
        $list_samples->{$sample} ++ ;
    }
}
close $SRC_ped ;
close $TGT_samplesex ;

my $file_ped_orig = $file_ped if($flag_multicallers) ;
my $file_sample_sex_orig = $file_sample_sex  if ($flag_multicallers) ;

if ($flag_multicallers) {
    close $TGT_pedmcaller ;
    close $TGT_samplesex_mcaller ;

    $file_ped = $file_ped_mcaller ;
    $file_sample_sex = $file_sample_sex_mcaller ;
}
$log->write ("Output sex infomation of Sample to $file_sample_sex",1) ;
$log->write ("Load ped file finish",1) ;

# PreProcess 03: Input files concat and merge
#=======================================

my $log_stderr_vcf_concat = "$dir_log/stderr_vcfconcat.log" ;
my $log_stderr_vcf_sort   = "$dir_log/stderr_vcfsort.log" ;
my $file_merged_vcfgz     = "$work_prefix\_merged.vcf.gz" ;

my $files_to_process = [] ;

# PreProcess 03-1: Concat vcf files from same sample 
#=======================================
$log->write("Check and concat vcf files from same sample start.",1) ;

foreach my $sample (keys %$sample2file) {
    my $filenum_in_one_sample = scalar @{$sample2file->{$sample}} ;
    if ($filenum_in_one_sample > 1) {
	my $all_files = join (" ", @{$sample2file->{$sample}}) ;
	my $file_concated_vcfgz = "$sample\_concated.vcf.gz" ;

	my $cmd_vcf_concat = "vcf-concat $all_files 2>$log_stderr_vcf_concat | vcf-sort -p $threads_active 2>$log_stderr_vcf_sort | bgzip -@ $threads_active -c > $file_concated_vcfgz" ;
	$log->andRun($cmd_vcf_concat) ;

	push @$files_to_process , $file_concated_vcfgz ;

    } else {
	push @$files_to_process , $sample2file->{$sample}->[0] ;
    }
}

foreach my $file_to_process (@$files_to_process) {
    tabix_vcf($file_to_process) ;
}
$log->write("Check and concat vcf files from same sample finish.",1) ;


# PreProcess 03-2: Merge vcf files from different sample
#=======================================
$log->write ("Merge vcf files from different samples start.",1) ;
my $log_stderr_vcf_merge = "$dir_log/stderr_vcfmerge.log" ;

my $cmd_merge = "vcf-merge -R '0/0' " . join (" ", @$files_to_process) . " 2> $log_stderr_vcf_merge |bgzip -@ $threads_active -c > $file_merged_vcfgz" ;
#    my $cmd_merge = "bcftools merge " . join (" ", @$files_to_process) . " 2> $log_stderr_vcf_merge |bgzip -@ $threads_active -c > $file_merged_vcfgz" ;
$log->andRun($cmd_merge) ;
$log->write ("Merge vcf files from different samples finish.",1) ;

# PreProcess 04: fix vcf ploidy
#=======================================
my $file_vcf_cleanXY = $file_prefix . "_cleanXY.vcf" ;
my $file_vcfgz_cleanXY = $file_prefix . "_cleanXY.vcf.gz" ;
my $log_stderr_fixploidy = "$dir_log/stderr_fixploidy.log" ;
$log->write ("Fixing ploidy of sex chromosome start",1) ;

my $cmd_fix_ploidy = "zcat $file_merged_vcfgz |vcf-fix-ploidy --samples $file_sample_sex >$file_vcf_cleanXY 2> $log_stderr_fixploidy" ;
$log->andRun($cmd_fix_ploidy) ;

my $cmd_bgzip_vcf_cleanXY = "bgzip -@ $threads_active -c $file_vcf_cleanXY > $file_vcfgz_cleanXY 2>> $log_stderr_fixploidy" ;
$log->andRun("$cmd_bgzip_vcf_cleanXY") ;

tabix_vcf($file_vcfgz_cleanXY) ;

$log->write ("Fixing ploidy of sex chromosome finish",1) ;


# for old config style info transfer
close $CONFIG ;
#my $cmd_cp_config = "cp $file_config $file_config_db" ;
#$log->andRun($cmd_cp_config) ;
#my $cmd_cp_foranalysis = "cp $file_vcfgz_cleanXY $file_vcfgz_foranalysis" ;
#$log->andRun($cmd_cp_foranalysis) ;



# Preprocess done. Do Annotation
#=======================================
my $cmd_vsann2 = "$dir_script/vsann2.pl -v $file_vcfgz_cleanXY -p $file_ped -j $jobid -d $file_db -n $threads_active 2>$dir_log/stderr_vsann2_$jobid.log" ;
print "\n$cmd_vsann2\n" ;
$log->andRun($cmd_vsann2) ;
	   
$log->write("VSannot finish",1) ;


sub tabix_vcf {
    my $file = shift ;

    if ($file !~ /gz$/) {
        my $file_orig = $file ;
        $file = "$file_orig.gz" ;
        bgzip($file_orig,$file) ;
    }

    my $cmd_tabix_vcf = "tabix -p vcf -f $file " ;
    $log->andRun($cmd_tabix_vcf) ;
}

sub bgzip {
    my $file_src = shift ;
    my $file_tgt = shift ;

    my $cmd_bgzip = "bgzip -@ $threads_active -c $file_src > $file_tgt" ;
    $log->andRun($cmd_bgzip) ;
}


