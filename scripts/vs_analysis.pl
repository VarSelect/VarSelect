#!/usr/bin/perl
use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use File::Path qw(make_path) ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;
use DBI ;

use threads ;
use threads::shared ;
use Thread::Queue ;


# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;
my $script_vslannot = "$dir_script/vsl_annot.pl" ;

# options handle
my $opts = {} ;
getopts("j:d:p:m:n:kuix:c:",$opts) ;
my $jobid               = $opts->{j} ;
my $file_db		= $opts->{d} ;
my $file_ped            = $opts->{p} ;
my $mode_workflow	= $opts->{m} ;
my $flag_multicallers   = $opts->{k} || 0 ;
my $threads_enable      = $opts->{n} || 16;
my $file_cnv_list	= $opts->{c} ;
my $file_xpr_list	= $opts->{x} ;

my $multicaller_insert = 0 ;
my $type_combine = "union" ;
$multicaller_insert = 'u' if ($flag_multicallers) ;

if ($opts->{u} && $opts->{i}) {
    die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
} elsif ($opts->{i}) {
    $type_combine = "intersection" ;
    $multicaller_insert = 'i' ;
}

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
make_path($dir_log , $dir_work ,{ chmod => 0755 ,} ) unless (-e $dir_work && -e $dir_log) ;

#my $file_config = "$dir_results_output/varselect.config" ;

my ($fname_db, $fpath_db , $file_ext) = fileparse($file_db, qw/.db/) ;
my $file_prefix = "$dir_results_output/$fname_db" ;
#my $file_config_db = "$fpath_db/$fname_db.config" ;

# Start logging
my $file_log = "$dir_log/log_vsanalysis_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("VSanlz start") ;
$log->write("Jobid: $jobid") ;
$log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ; # list of options


my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","") ;
$dbh->{RaiseError} = 1 ;
$dbh->{PrintError} = 1 ;

my $table_varselectinfo = 'varselect_info' ;
my $sth = $dbh->prepare("select * from $table_varselectinfo") ;
$sth->execute() ;
my $vslinfo = $sth->fetchrow_hashref ;

my $caller_list = {map {$_ => 1} split(/\,/, $vslinfo->{caller_list}) } ;
my $sample_list = {map {$_ => 1} split(/\,/, $vslinfo->{sample_list}) } ;

my $file_vcfgz_foranalysis = "$dir_results_output/for_analysis_$jobid.vcf.gz" ;
open (my $TGT_for_analysis,">$file_vcfgz_foranalysis") ;
print $TGT_for_analysis $vslinfo->{vcfgz_for_analysis} ;

while (my $vslinfo = $sth->fetchrow_hashref) {
    print $TGT_for_analysis $vslinfo->{vcfgz_for_analysis} ;
}

close $TGT_for_analysis ;

#my $cmd_gunzip_and_bzip_vcf = "gunzip -c $file_vcfgz_foranalysis |bgzip -c -@ $threads_enable > $file_vcfgz_foranalysis.tmp" ;
#$log->andRun($cmd_gunzip_and_bzip_vcf,1) ;
#`mv $file_vcfgz_foranalysis.tmp $file_vcfgz_foranalysis` ;
tabix_vcf($file_vcfgz_foranalysis,$log) ; 


# jobinfo table
my $sql_create_jobinfo = "CREATE TABLE IF NOT EXISTS jobinfo ( " ;
$sql_create_jobinfo .= " jobid vahrchar(20) PRIMARY KEY" ;
$sql_create_jobinfo .= " , ped TEXT DEFAULT NULL" ;
$sql_create_jobinfo .= " , multicaller TEXT DEFAULT NULL" ;
$sql_create_jobinfo .= " , workflow varchar(20) DEFAULT NULL" ;
$sql_create_jobinfo .= ")"  ;

my $sth_create_jobinfo = $dbh->prepare($sql_create_jobinfo) ;
$sth_create_jobinfo->execute() ;



my $file_sample_sex = $file_prefix ."_sample_sex.txt" ;
my $file_sample_sex_mcaller = $file_prefix ."_sample_sex_mcaller.txt" ;

## Update PED
#my $cmd_gemini_amend = "gemini amend --sample $file_ped $file_db 2> $dir_log/stderr_ped_amend_$jobid.log" ;
#$log->andRun ($cmd_gemini_amend) ;

my ($ped_fn,$ped_path,$ped_ext) = fileparse($file_ped,qw/.ped/) ;
my $cmd_cp_ped_to_anadir = "cp $file_ped $dir_results_output/$ped_fn\_$jobid.ped" ;

# Loading PED
# check affected status of sample
my $affect_samples = [] ;
my $unaffect_samples = [] ;
my $samples = [] ;
my $sex_of_sample = {} ;
my $ped_text ;
open (my $SRC_ped,"$file_ped") ;
while (my $line = <$SRC_ped>) {
    chomp $line;
    $ped_text .= "$line\n" ;
    next if ($line =~ /^#/) ;

    my ($family,$sample,$father,$mother,$sex,$affect_status) = split(/\t/,$line) ;

    $sex_of_sample->{$sample} = $sex ;

    if ($affect_status == 1) {
	push @$unaffect_samples , $sample ;

    } elsif ($affect_status == 2) {
	push @$affect_samples , $sample ;
    }

    push @$samples , $sample ;
}
close $SRC_ped ;

my $sql_insert_jobinfo = "INSERT INTO jobinfo (jobid,ped,workflow,multicaller) values (?,?,?,?) " ; 
my $sth_insert_jobinfo = $dbh->prepare($sql_insert_jobinfo) ;

$sth_insert_jobinfo->bind_param(1,$jobid) ;
$sth_insert_jobinfo->bind_param(2,$ped_text) ;
$sth_insert_jobinfo->bind_param(3,$mode_workflow) ;
$sth_insert_jobinfo->bind_param(4,$multicaller_insert) ;
$sth_insert_jobinfo->execute() ;


#if (-e $file_vcfgz_foranalysis) {
#
#} else {
#    die("$file_vcfgz_foranalysis not found! please check your db create directory and put the $fname_db.vcf.gz to the current location of your db!\n") ;
#}

# Start thread
$log->write("Threads running start",1) ;

# Thread I: create thread and queue for vsl_annot
my $vslannotQueue = Thread::Queue->new ;
my $vslannot_thread = async {
    $log->write("thread for vsl_annot start",1) ;
    while ( my $cmd_vslannot = $vslannotQueue->dequeue) {
	$log->andRun($cmd_vslannot,1) ;
	$log->write("$cmd_vslannot Finish", 1) ;
    }
    $log->write("thread for vsl_annot finish",1) ;
} ;



# Check multi-caller sample
if ($flag_multicallers) {
    $log->write("Multicaller processing start" , 1) ;
    if ($type_combine eq "union") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;

        my $file_output_union_conflict = "$dir_results_output/multicaller_union_inconsistant_$jobid.txt" ;
        my $file_vcf_union_set = "$dir_results_output/multicaller_unionset_$jobid.vcf";

        open (my $TGT_union_conflict , ">$file_output_union_conflict" ) ;
        open (my $TGT_union_vcf, ">$file_vcf_union_set") ;

#        print $TGT_union_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , sort keys %$sample_list)) . "\n" ;
        print $TGT_union_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , @$samples)) . "\n" ;

        while (my $var = $vcf_parser->next_var) {
            my $flag_inanalysis = 1 ;
            my $reason_drop = '' ;
            my $gt_all_in_one_locus = [] ;
            my $list_gt_by_sample = [] ;

	    my $chr = $var->{CHROM} ;
            # conflict if more than one gt type in gt_got
#            foreach my $sample (sort keys %$sample_list) {
            foreach my $sample (@$samples) {
                my $gt_got = {} ;
		my $sex = $sex_of_sample->{$sample} ;

                foreach my $caller (sort keys %$caller_list) {
                    my $sample_tag = $sample . "_" . $caller ;
                    my $gt = $var->{sample_val}->{$sample_tag}->{GT} ;
                    $gt =~ s/\|/\// ;   # ignore phased, treat as an un-phased genotype
                    next if ($gt eq '0/0') ;
                    $gt_got->{$gt} ++ ;

                    push @$gt_all_in_one_locus , "$sample_tag:$gt" ;
                }

                my $list_of_got_gt = [keys %$gt_got] ;
                my $gt_type_count = scalar @$list_of_got_gt ;

                if ($gt_type_count > 1) {
		    #union con-flict
                   $flag_inanalysis = 0 ;
                    $reason_drop .= "Conflict genotype in 1+ callers of sample $sample,"

                } elsif ($gt_type_count == 1) {

		    if ($chr eq 'Y' && $sex eq '2') { # female chrY should be '.'
			push @$list_gt_by_sample , '.' ;

		    } else {
			push @$list_gt_by_sample , $list_of_got_gt->[0] ;
		    }

                } else {
                    push @$list_gt_by_sample , '0/0' ;
                }
            }

            if ($flag_inanalysis) {
		my $i = 0 ;
#		my $info_text = join(";" , map {"genotype_union_$_=" . $list_gt_by_sample->[$i++]} sort keys %$sample_list) ;
		my $info_text = join(";" , map {"genotype_union_$_=" . $list_gt_by_sample->[$i++]} @$samples) ;
                # output vcf with sample only (no caller now) for workflow analysis
                print $TGT_union_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\t$info_text\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;

            } else { #inconsistence call
                print $TGT_union_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) ."\t$reason_drop\n"  ;
            }
        }
        close $TGT_union_conflict ;
        close $TGT_union_vcf ;

        tabix_vcf($file_vcf_union_set) ;


	my $cmd_vsl_annot = "$dir_script/vsl_annot.pl " ;
	$cmd_vsl_annot .= " -f $file_vcf_union_set.gz " ;
	$cmd_vsl_annot .= " -e " . join("," , map {"genotype_union_$_"} @$samples) ;
	$cmd_vsl_annot .= " -t " . join("," , map {"text"} @$samples ) ;
	$cmd_vsl_annot .= " -c " . join("," , map {"genotype_union_$_\_$jobid"} @$samples) ;
	$cmd_vsl_annot .= " -d $file_db " ;
	$cmd_vsl_annot .= " -b varselect_variants " ;
	$cmd_vsl_annot .= " 2> $dir_log/stderr_vslannot_uniongt_$jobid.log" ;
	#$cmd_vsl_annot .= " -e " . join("," , map {"genotype_union_$_"} sort keys %$sample_list) ;
	#$cmd_vsl_annot .= " -c " . join("," , map {"genotype_union_$_\_$jobid"} sort keys %$sample_list) ;

	#$log->andRun($cmd_vsl_annot) ;
	$log->write("enqueue $cmd_vsl_annot" , 1) ;
	$vslannotQueue->enqueue($cmd_vsl_annot) ;


        $file_vcfgz_foranalysis = "$file_vcf_union_set.gz" ;

    } elsif ($type_combine eq "intersection") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;
        my $file_output_isect_conflict = "$dir_results_output/multicaller_intersect_inconsistant_$jobid.txt" ;
        my $file_vcf_isect_set = "$dir_results_output/multicaller_insectset_$jobid.vcf" ;

        open (my $TGT_isect_conflict , ">$file_output_isect_conflict") ;
        open (my $TGT_isect_vcf , ">$file_vcf_isect_set") ;

#        print $TGT_isect_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , sort keys %$sample_list)) . "\n" ;
        print $TGT_isect_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , @$samples)) . "\n" ;

        while(my $var=$vcf_parser->next_var) {
            my $flag_inanalysis = 1 ;
            my $reason_drop = '' ;
            my $gt_all_in_one_locus = [] ;
            my $list_gt_by_sample = [] ;
	    my $chr = $var->{CHROM} ;

#            foreach my $sample (sort keys %$sample_list) {
            foreach my $sample (@$samples) {
                my $gt_got = {} ;
		my $sex = $sex_of_sample->{$sample} ;
                foreach my $caller (sort keys %$caller_list) {
                    my $sample_tag = $sample . "_" . $caller ;
                    my $gt = $var->{sample_val}->{$sample_tag}->{GT} ;
                    $gt =~ s/\|/\// ;   # ignore phased, treat as an un-phased genotype
                    $gt_got->{$gt} ++ ;
                    push @$gt_all_in_one_locus, "$sample_tag:$gt" ;
                }

                my $list_of_got_gt = [keys %$gt_got] ;
                my $gt_type_count = scalar @$list_of_got_gt ;

                if ($gt_type_count == 1) {

                    $flag_inanalysis *= 1 ;

		    if ($chr eq 'Y' && $sex eq '2') { # female chY should be '.'
                        push @$list_gt_by_sample , '.' ;

                    } else {
                        push @$list_gt_by_sample , $list_of_got_gt->[0] ;
                    }
		    
                } else {
		    # Conflict genotype in 1+ callers of sample $sample  in insect
                    $flag_inanalysis *= 0 ;

                    $reason_drop .= "Conflict genotype in 1+ callers of sample $sample," ;
                }
            }

            if ($flag_inanalysis) {
		my $i = 0 ;
		my $info_text = join(";" , map {"genotype_isec_$_=" . $list_gt_by_sample->[$i++]} @$samples) ;
                print $TGT_isect_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\t$info_text\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;

            } else {
                print $TGT_isect_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) . "\t$reason_drop\n"  ;
            }
        }
        close $TGT_isect_conflict ;
        close $TGT_isect_vcf ;

        tabix_vcf($file_vcf_isect_set) ;

	my $cmd_vsl_annot = "$dir_script/vsl_annot.pl " ;
	$cmd_vsl_annot .= " -f $file_vcf_isect_set.gz " ;
	$cmd_vsl_annot .= " -e " . join("," , map {"genotype_isec_$_"} @$samples) ;
	$cmd_vsl_annot .= " -t " . join("," , map {"text"} @$samples ) ;
	$cmd_vsl_annot .= " -c " . join("," , map {"genotype_isec_$_\_$jobid"} @$samples) ;
	$cmd_vsl_annot .= " -d $file_db " ;
	$cmd_vsl_annot .= " -b varselect_variants " ;
	$cmd_vsl_annot .= " 2> $dir_log/stderr_vslannot_osecgt_$jobid.log" ;

	#$log->andRun($cmd_vsl_annot) ;
	$log->write("enqueue $cmd_vsl_annot" , 1) ;
	$vslannotQueue->enqueue($cmd_vsl_annot) ;


        $file_vcfgz_foranalysis = "$file_vcf_isect_set.gz" ;
    }
    $log->write("Multicaller processing finish" , 1) ;

    # the end of multi_caller mode
}



# Check GT frequency
$log->write("Genotype frequency calculating start" , 1) ;
my $file_output_gtfreq = "$dir_results_output/aff_gt_freq_$jobid.txt" ;
my $file_vcfgz_gtfreq  = "$dir_results_output/aff_gt_freq_$jobid.vcf.gz" ;
my $file_output_gtfreq_unaff = "$dir_results_output/unaff_gt_freq_$jobid.txt" ;
my $file_vcfgz_gtfreq_unaff  = "$dir_results_output/unaff_gt_freq_$jobid.vcf.gz" ;

open (my $TGT_gtfreq ,">$file_output_gtfreq") ;
open (my $TGT_gtfreq_unaff ,">$file_output_gtfreq_unaff") ;

$log->write("Loading genotype frequency of each sample, start." , 1) ;

my $vcf_parser = Vcf_parser->new(file=>$file_vcfgz_foranalysis) ;
while (my $var = $vcf_parser->next_var) {
    my $affect_sample_num = scalar @$affect_samples ;
    my $unaffect_sample_num = scalar @$unaffect_samples ;

    if ($affect_sample_num) {
	my $gt_count = {} ;
	my $list_output_gt = [] ;
	my $list_output_freq = [] ;

	foreach my $sample (@$affect_samples) {
	    my $genotype = $var->{sample_val}->{$sample}->{GT} ;
	    $gt_count->{$genotype} ++ unless ($genotype eq '0/0'|| $genotype eq '0|0' || $genotype eq '0') ; # record , but skip ref-gt
	}
	my $gt_freq = {map {$_ => $gt_count->{$_} / $affect_sample_num } keys %$gt_count } ;
	my $alt_sum = 0 ;
	foreach my $gt (sort keys %$gt_freq) {
	    push @$list_output_gt   , $gt ;
	    push @$list_output_freq , $gt_freq->{$gt} ;
	    $alt_sum += $gt_freq->{$gt} ;
	}

	my $output_count = scalar @$list_output_gt ;
	if ($output_count) {
	    print $TGT_gtfreq join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
	    print $TGT_gtfreq "\t" . join("," , @$list_output_gt) ;
	    print $TGT_gtfreq "\t" . join("," , @$list_output_freq) ;
	    print $TGT_gtfreq "\t$alt_sum" ;
	    print $TGT_gtfreq "\n" ;
	}
    }

    if ($unaffect_sample_num) {
	my $gt_count = {} ;
	my $list_output_gt_unaff = [] ;
	my $list_output_freq_unaff = [] ;

	foreach my $sample (@$unaffect_samples) {
            my $genotype = $var->{sample_val}->{$sample}->{GT} ;
            $gt_count->{$genotype} ++ unless ($genotype eq '0/0'|| $genotype eq '0|0' || $genotype eq '0') ; # record , but skip ref-gt
        }
	my $gt_freq = {map {$_ => $gt_count->{$_} / $unaffect_sample_num } keys %$gt_count } ;

	my $alt_sum =0 ;
	foreach my $gt (sort keys %$gt_freq) {
                push @$list_output_gt_unaff   , $gt ;
                push @$list_output_freq_unaff , $gt_freq->{$gt} ;
		$alt_sum += $gt_freq->{$gt} ;
	}

	my $output_count = scalar @$list_output_gt_unaff ;
	if ($output_count) {
	    print $TGT_gtfreq_unaff join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
	    print $TGT_gtfreq_unaff "\t" . join("," , @$list_output_gt_unaff) ;
	    print $TGT_gtfreq_unaff "\t" . join("," , @$list_output_freq_unaff) ;
	    print $TGT_gtfreq_unaff "\t$alt_sum" ;
	    print $TGT_gtfreq_unaff "\n" ;
	}

    }
}
$log->write("Loading genotype frequency of each sample, finish." , 1) ;

close $TGT_gtfreq ;
close $TGT_gtfreq_unaff ;

tabix_vcf($file_output_gtfreq,$log) ;
tabix_vcf($file_output_gtfreq_unaff,$log) ;


my $cmd_vannot_var_freq = "zcat $file_vcfgz_foranalysis| vcf-annotate " ;
$cmd_vannot_var_freq .= " -a $file_output_gtfreq.gz " ;
$cmd_vannot_var_freq .= " -c " . join("," , qw"CHROM POS REF ALT INFO/affected_vgt INFO/affected_vgfreq INFO/affected_altfreqsum") ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_vgt,Number=.,Type=String,Description='Variant genotype in affected samples' " ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_vgfreq,Number=.,Type=String,Description='Variant genotype frequency in affected samples' " ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_altfreqsum,Number=.,Type=Float,Description='Summary of affected-alt-frequency variant genotype in affected samples' " ;
$cmd_vannot_var_freq .= " 2> $dir_log/stderr_vannot_gtfreq_$jobid.log " ;
$cmd_vannot_var_freq .= " | bgzip -c > $file_vcfgz_gtfreq " ;

$log->write("vcf-annot affected vgt start" , 1) ;
$log->andRun($cmd_vannot_var_freq) ;
tabix_vcf($file_vcfgz_gtfreq) ;
$log->write("vcf-annot affected vgt finish" , 1) ;

my $cmd_vsl_gtfreq = "$script_vslannot" ;
$cmd_vsl_gtfreq   .= " -a extract " ;
$cmd_vsl_gtfreq   .= " -f $file_vcfgz_gtfreq " ;
$cmd_vsl_gtfreq	  .= " -e " . join ("," , qw/affected_vgt affected_vgfreq affected_altfreqsum/) ;
$cmd_vsl_gtfreq   .= " -t " . join ("," , qw/text text float/) ;
$cmd_vsl_gtfreq	  .= " -c " . join ("," , map {"$_\_$jobid"} qw/affected_vgt affected_vgfreq affected_altfreqsum/) ;
$cmd_vsl_gtfreq	  .= " -d $file_db " ;
$cmd_vsl_gtfreq	  .= " -b varselect_variants " ;

#$log->andRun($cmd_vsl_gtfreq) ;
$log->write("enqueue $cmd_vsl_gtfreq" , 1) ;
$vslannotQueue->enqueue($cmd_vsl_gtfreq) ;

my $cmd_vannot_var_freq_unaff = "zcat $file_vcfgz_foranalysis| vcf-annotate " ;
$cmd_vannot_var_freq_unaff .= " -a $file_output_gtfreq_unaff.gz " ;
$cmd_vannot_var_freq_unaff .= " -c " . join("," , qw"CHROM POS REF ALT INFO/unaffected_vgt INFO/unaffected_vgfreq INFO/unaffected_altfreqsum") ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_vgt,Number=.,Type=String,Description='variant genotype in unaffected samples' " ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_vgfreq,Number=.,Type=String,Description='Variant genotype frequency in unaffected samples' " ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_altfreqsum,Number=.,Type=Float,Description='Summary of unaffected-alt-frequency variant genotype in affected samples' " ;
$cmd_vannot_var_freq_unaff .= " 2> $dir_log/stderr_vannot_gtfreq_$jobid.log " ;
$cmd_vannot_var_freq_unaff .= " | bgzip -c > $file_vcfgz_gtfreq_unaff " ;

$log->write("vcf-annot un-affected vgt start" , 1) ;
$log->andRun($cmd_vannot_var_freq_unaff) ;
tabix_vcf($file_vcfgz_gtfreq_unaff) ;
$log->write("vcf-annot un-affected vgt finish" , 1) ;

my $cmd_vsl_gtfreq_unaff = "$script_vslannot" ;
$cmd_vsl_gtfreq_unaff   .= " -a extract" ;
$cmd_vsl_gtfreq_unaff   .= " -f $file_vcfgz_gtfreq_unaff" ;
$cmd_vsl_gtfreq_unaff	.= " -e " . join ("," , qw/unaffected_vgt unaffected_vgfreq unaffected_altfreqsum/) ;
$cmd_vsl_gtfreq_unaff   .= " -t " . join ("," , qw/text text float/) ;
$cmd_vsl_gtfreq_unaff	.= " -c " . join ("," , map {"$_\_$jobid"} qw/unaffected_vgt unaffected_vgfreq unaffected_altfreqsum/)  ;
$cmd_vsl_gtfreq_unaff   .= " -d $file_db " ;
$cmd_vsl_gtfreq_unaff   .= " -b varselect_variants" ;

#$log->andRun($cmd_vsl_gtfreq_unaff) ;
$log->write("enqueue $cmd_vsl_gtfreq_unaff" , 1) ;
$vslannotQueue->enqueue($cmd_vsl_gtfreq_unaff) ;

$log->write("Genotype frequency calculating finish" , 1) ;


my $vcf_input_for_next_step = $file_vcfgz_gtfreq_unaff ;

# if -c is set 
if ($file_cnv_list) {
    $log->write("CNV annotate start" , 1) ;

    my $file_vcfgz_cnv = "$dir_results_output/$jobid\_cnv.vcf.gz" ;
    my $cmd_cnv = "$dir_script/cnvkit_parser.pl " ;
    $cmd_cnv .= " -d $file_db " ;
    $cmd_cnv .= " -p $file_ped " ;
    $cmd_cnv .= " -c $file_cnv_list " ;
    $cmd_cnv .= " -i $file_vcfgz_foranalysis " ;
    $cmd_cnv .= " -o $file_vcfgz_cnv " ;
    $cmd_cnv .= " -j $jobid" ;

    $log->andRun($cmd_cnv) ;

    my $collist_cnvkit = [qw/cnv_samples cnv_log2 cnv_fc_samples cnv_foldchange_log2/] ;

    my $cmd_vslannot_cnv = "$script_vslannot " ;
    $cmd_vslannot_cnv .= " -a extract " ;
    $cmd_vslannot_cnv .= " -f $file_vcfgz_cnv " ;
    $cmd_vslannot_cnv .= " -e " . join ("," , @$collist_cnvkit)  ;
    $cmd_vslannot_cnv .= " -t " . join ("," , map {'text'} @$collist_cnvkit) ;
    $cmd_vslannot_cnv .= " -c " . join ("," , map {"$_\_$jobid"} @$collist_cnvkit );
    $cmd_vslannot_cnv .= " -d $file_db " ;
    $cmd_vslannot_cnv .= " -b varselect_variants" ;

#    $log->andRun($cmd_vslannot_cnv) ;
    $log->write("enqueue $cmd_vslannot_cnv" , 1) ;
    $vslannotQueue->enqueue($cmd_vslannot_cnv) ;

    $log->write("CNV annotate finish" , 1) ;
}

# if -x is set
if ($file_xpr_list) {
    $log->write("Expression profile parse and annotate start" , 1) ;
    my $script_xprparer = "$dir_script/xprprofile_parser.pl" ;
    my $file_vcfgz_xpr = "$file_prefix\_xpr.vcf.gz" ;

    my $cmd_xpr = "$script_xprparer " ;
    $cmd_xpr .= " -p $file_ped " ;
    $cmd_xpr .= " -x $file_xpr_list " ;
    $cmd_xpr .= " -i $file_vcfgz_foranalysis " ;
#    $cmd_xpr .= " -i $vcf_input_for_next_step " ;
    $cmd_xpr .= " -o $file_vcfgz_xpr " ;
    $cmd_xpr .= " -j $jobid " ;
    $cmd_xpr .= " -l $dir_log " ;
    $log->andRun($cmd_xpr) ;

    # prepare for gemini annotate
    my $collist_xpr = [qw/xpr_samples xpr_readcount xpr_tpm xpr_fpkm xpr_exonreadcount xpr_exonid xpr_fc_samples xpr_foldchange_readcount xpr_foldchange_tpm xpr_foldchange_fpkm xpr_foldchange_exonreadcount/] ;

    my $cmd_vslannot_xpr = "$script_vslannot" ;
    $cmd_vslannot_xpr .= " -f $file_vcfgz_xpr " ;
    $cmd_vslannot_xpr .= " -a extract " ;
    $cmd_vslannot_xpr .= " -e " . join("," , @$collist_xpr)  ;
    $cmd_vslannot_xpr .= " -t " . join("," , map {'text'} @$collist_xpr ) ;
    $cmd_vslannot_xpr .= " -c " . join("," , map {"$_\_$jobid"} @$collist_xpr) ;
    $cmd_vslannot_xpr .= " -d $file_db" ;
    $cmd_vslannot_xpr .= " -b varselect_variants" ;

#    $log->andRun($cmd_vslannot_xpr) ;
    $log->write("enqueue $cmd_vslannot_xpr" , 1) ;
    $vslannotQueue->enqueue($cmd_vslannot_xpr) ;

    $log->write("Expression profile parse and annotate finish" , 1) ;
}

# -m = paired or family
my $col_list_for_wf_union_set = [] ;
if ($mode_workflow eq 'family') {
    $log->write("Family workflow start" , 1) ;
    my $file_vcf_for_workflow = "$dir_work/for_genetic_workflow_$jobid.vcf" ;
    my $cmd_extract_for_workflow = "zcat $file_vcfgz_foranalysis > $file_vcf_for_workflow" ;
    $log->andRun($cmd_extract_for_workflow) ;

    my $file_stderr_wf_ad = "$dir_log/stderr_runworkflow_AD_$jobid.log" ;
    my $file_stderr_wf_ar = "$dir_log/stderr_runworkflow_AR_$jobid.log" ;
    my $file_stderr_wf_ch = "$dir_log/stderr_runworkflow_CH_$jobid.log" ;
    my $file_stderr_wf_dd = "$dir_log/stderr_runworkflow_DD_$jobid.log" ;
    my $file_stderr_wf_dr = "$dir_log/stderr_runworkflow_DR_$jobid.log" ;
    my $file_stderr_wf_sh = "$dir_log/stderr_runworkflow_SH_$jobid.log" ;
    my $file_stderr_wf_xl = "$dir_log/stderr_runworkflow_XL_$jobid.log" ;

    my $file_gz_output_ad = "$dir_results_output/workflow_AD_output_$jobid.txt.gz" ;
    my $file_gz_output_ar = "$dir_results_output/workflow_AR_output_$jobid.txt.gz" ;
    my $file_gz_output_ch = "$dir_results_output/workflow_CH_output_$jobid.txt.gz" ;
    my $file_gz_output_dd = "$dir_results_output/workflow_DD_output_$jobid.txt.gz" ;
    my $file_gz_output_dr = "$dir_results_output/workflow_DR_output_$jobid.txt.gz" ;
    my $file_gz_output_sh = "$dir_results_output/workflow_SH_output_$jobid.txt.gz" ;
    my $file_gz_output_xl = "$dir_results_output/workflow_XL_output_$jobid.txt.gz" ;

    my $col_ad = "is_AD_$jobid" ;
    my $col_ar = "is_AR_$jobid" ;
    my $col_ch = "is_CH_$jobid" ;
    my $col_dd = "is_DD_$jobid" ;
    my $col_dr = "is_DR_$jobid" ;
    my $col_sh = "is_SH_$jobid" ;
    my $col_xl = "is_XL_$jobid" ;

    my $file_vcfgz_ad_output = "$dir_work/AutosomalDominant_$jobid.vcf.gz" ;
    my $file_vcfgz_ar_output = "$dir_work/AutosomalRecessive_wAd_$jobid.vcf.gz" ;
    my $file_vcfgz_ch_output = "$dir_work/CompoundHet_wAdAr_$jobid.vcf.gz" ;
    my $file_vcfgz_dd_output = "$dir_work/DenovoDominant_wAdArCh_$jobid.vcf.gz" ;
    my $file_vcfgz_dr_output = "$dir_work/DenovoRecessive_wAdArChDd_$jobid.vcf.gz" ;
    my $file_vcfgz_sh_output = "$dir_work/SecondHit_wAdArChDdDr_$jobid.vcf.gz" ;
    my $file_vcfgz_xl_output = "$dir_work/XLinked_wAdArChDdDrTh_$jobid.vcf.gz" ;

    push @$col_list_for_wf_union_set, map { $_.$jobid } qw/is_AD_ is_AR_ is_CH_ is_DD_ is_DR_ is_SH_ is_XL_/ ;

    # Autosomal-dominant
    my $exe_wf_ad  = "$dir_workflows/autosomal-dominant/Autosomal-dominant.py" ;
    my $exe_wf2_ad = "$dir_workflows/autosomal-dominant/Compare.py" ;
    running_workflow($exe_wf_ad , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ad , $file_gz_output_ad , $exe_wf2_ad) ;
    wf_vannot($file_gz_output_ad , $col_ad , $vcf_input_for_next_step , $file_vcfgz_ad_output ) ;

    # Autosomal-recessive
    my $exe_wf_ar  = "$dir_workflows/autosomal-recessive/Autosomal-recessive.py" ;
    my $exe_wf2_ar = "$dir_workflows/autosomal-recessive/Compare.py" ;
    running_workflow($exe_wf_ar , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ar , $file_gz_output_ar , $exe_wf2_ar) ;
    wf_vannot($file_gz_output_ar , $col_ar , $file_vcfgz_ad_output , $file_vcfgz_ar_output ) ;

    # Compound-het
    my $exe_wf_ch  = "$dir_workflows/compound-het/Compound-het.py" ;
    my $exe_wf2_ch = "$dir_workflows/compound-het/Compare.py" ;
    running_workflow($exe_wf_ch , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ch , $file_gz_output_ch , $exe_wf2_ch , $file_db) ;
    wf_vannot($file_gz_output_ch , $col_ch , $file_vcfgz_ar_output , $file_vcfgz_ch_output) ;

    # De Novo-dominant
    my $exe_wf_dd  =  "$dir_workflows/denovo-dominant/Denovo-dominant.py" ; 
    my $exe_wf2_dd =  "$dir_workflows/denovo-dominant/Compare.py" ;
    running_workflow($exe_wf_dd , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_dd , $file_gz_output_dd, $exe_wf2_dd) ;
    wf_vannot($file_gz_output_dd , $col_dd , $file_vcfgz_ch_output , $file_vcfgz_dd_output) ;

    # De Novo-recessive
    my $exe_wf_dr  =  "$dir_workflows/denovo-recessive/Denovo-recessive.py" ;
    my $exe_wf2_dr =  "$dir_workflows/denovo-recessive/Compare.py" ;
    running_workflow($exe_wf_dr , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_dr , $file_gz_output_dr, $exe_wf2_dr) ;
    wf_vannot($file_gz_output_dr , $col_dr , $file_vcfgz_dd_output , $file_vcfgz_dr_output) ;

    # Second-hit
    my $exe_wf_sh  = "$dir_workflows/second-hit/Second-hit.py" ;
    my $exe_wf2_sh = "$dir_workflows/second-hit/Compare.py" ;
    running_workflow($exe_wf_sh , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_sh , $file_gz_output_sh , $exe_wf2_sh) ;
    wf_vannot($file_gz_output_sh , $col_sh , $file_vcfgz_dr_output , $file_vcfgz_sh_output) ;

    # X-linked
    my $exe_wf_xl  = "$dir_workflows/x-linked/X-linked.py" ;
    my $exe_wf2_xl = "$dir_workflows/x-linked/Compare.py" ;
    running_workflow($exe_wf_xl , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_xl , $file_gz_output_xl , $exe_wf2_xl) ;
    wf_vannot($file_gz_output_xl , $col_xl , $file_vcfgz_sh_output , $file_vcfgz_xl_output) ;

    $vcf_input_for_next_step = $file_vcfgz_xl_output ;

    my $in_analysis_file_list = {
	$file_gz_output_ad => 5,
	$file_gz_output_ar => 5,
	$file_gz_output_ch => 5,
	$file_gz_output_dd => 5,
	$file_gz_output_dr => 5,
	$file_gz_output_sh => 5,
	$file_gz_output_xl => 5,
    } ;

    my $file_gz_inanalysis = "$dir_results_output/inanalysis_$jobid.txt.gz" ;
    my $file_vcfgz_inanalysis = "$dir_results_output/inanalysis_$jobid.vcf.gz" ;

    get_union($in_analysis_file_list , $file_gz_inanalysis ) ;
    tabix_vcf($file_gz_inanalysis) ;

    my $cmd_vannot_inanalysis = "zcat $vcf_input_for_next_step | vcf-annotate " ;
    $cmd_vannot_inanalysis .= " -a $file_gz_inanalysis " ;
    $cmd_vannot_inanalysis .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/in_analysis_$jobid")) ;
    $cmd_vannot_inanalysis .= " -d key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='in_analysis_$jobid' " ;
    $cmd_vannot_inanalysis .= "2> $dir_log/stderr_vannot_inanalysis_$jobid.log " ;
    $cmd_vannot_inanalysis .= " |bgzip -c > $file_vcfgz_inanalysis " ;

    $log->andRun($cmd_vannot_inanalysis) ;
    tabix_vcf($file_vcfgz_inanalysis) ;
    
    $vcf_input_for_next_step = $file_vcfgz_inanalysis ;

    my $collist_family_model = [($col_ad, $col_ar , $col_ch , $col_dd , $col_dr , $col_sh , $col_xl , "in_analysis_$jobid")] ;

    my $cmd_vslannot_family = "$script_vslannot" ;
    $cmd_vslannot_family .= " -a extract " ;
    $cmd_vslannot_family .= " -f $file_vcfgz_inanalysis " ;
    $cmd_vslannot_family .= " -e " . join("," , @$collist_family_model) ;
    $cmd_vslannot_family .= " -t " . join("," , map {"boolean"} @$collist_family_model) ;
    $cmd_vslannot_family .= " -c " . join("," , @$collist_family_model) ;
    $cmd_vslannot_family .= " -d $file_db " ;
    $cmd_vslannot_family .= " -b varselect_variants" ;

#    $log->andRun($cmd_vslannot_family) ;
    $log->write("enqueue $cmd_vslannot_family" , 1) ;
    $vslannotQueue->enqueue($cmd_vslannot_family) ;

} elsif ($mode_workflow eq 'paired') {
    $log->write("Paired case/control workflow start" , 1) ;

    # LOH detect
    $log->write("LOH detect start") ;
    my $script_lohdetect = "$dir_script/loh_detector.pl" ;
    my $file_vcfgz_isLOH = "$file_prefix\_loh.vcf.gz" ;
    my $file_output_isLOH = "$file_prefix\_loh.txt" ;

    my $cmd_loh = "$script_lohdetect" ;
    $cmd_loh .= " -j $jobid " ;
#    $cmd_loh .= " -i $file_vcfgz_foranalysis " ;
    $cmd_loh .= " -i $vcf_input_for_next_step " ;
    $cmd_loh .= " -v $file_vcfgz_isLOH " ;
    $cmd_loh .= " -o $file_output_isLOH " ;
    $cmd_loh .= " -p $file_ped " ;
    $cmd_loh .= " -d $file_db " ;
    $cmd_loh .= " -n $threads_enable " ;

    $log->andRun($cmd_loh) ;
    $log->loadlog("$dir_log/log_lohdetect_$jobid.log") ;
    $log->write("LOH detect finish") ;
    push @$col_list_for_wf_union_set, "is_LOH_$jobid" ;

    # Peired sample de novo detect
    $log->write("Paired sample denovo variants detect start") ;
    my $script_denovo_detect = "$dir_script/paired_denovo_detector.pl" ;
    my $file_output_is_denovo = "$dir_results_output/is_denovo_$jobid.txt" ;
    my $file_vcfgz_is_denovo = "$file_prefix\_denovo.vcf.gz" ;

    my $cmd_denovo = "$script_denovo_detect " ;
#    $cmd_denovo .= " -i $file_vcfgz_foranalysis " ;
    $cmd_denovo .= " -i $file_vcfgz_isLOH " ;
    $cmd_denovo .= " -o $file_output_is_denovo " ;
    $cmd_denovo .= " -v $file_vcfgz_is_denovo " ;
    $cmd_denovo .= " -p $file_ped " ;
    $cmd_denovo .= " -j $jobid " ;
    $cmd_denovo .= " -d $file_db " ;
    $cmd_denovo .= " -n $threads_enable " ;

    $log->andRun($cmd_denovo) ;

    tabix_vcf("$file_output_is_denovo") ;

    $log->write("Paired sample denovo variants detect finish") ;
    push @$col_list_for_wf_union_set, "is_denovo_$jobid" ;
    
    my $file_gz_inanalysis = "$dir_results_output/inanalysis_$jobid.txt.gz" ;
    my $file_vcfgz_inanalysis = "$dir_results_output/inanalysis_$jobid.vcf.gz" ;

    $vcf_input_for_next_step = $file_vcfgz_is_denovo ;

    my $in_analysis_file_list = {
	$file_output_isLOH => 5,
	$file_output_is_denovo => 5,
    } ;

    my $file_gz_inanalysis = "$dir_results_output/inanalysis_$jobid.txt.gz" ;
    my $file_vcfgz_inanalysis = "$dir_results_output/inanalysis_$jobid.vcf.gz" ;

    get_union($in_analysis_file_list , $file_gz_inanalysis ) ;
    tabix_vcf($file_gz_inanalysis) ;

    my $cmd_vannot_inanalysis = "zcat $vcf_input_for_next_step | vcf-annotate " ;
    $cmd_vannot_inanalysis .= " -a $file_gz_inanalysis " ;
    $cmd_vannot_inanalysis .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/in_analysis_$jobid")) ;
    $cmd_vannot_inanalysis .= " -d key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='in_analysis_$jobid' " ;
    $cmd_vannot_inanalysis .= "2> $dir_log/stderr_vannot_inanalysis_$jobid.log " ;
    $cmd_vannot_inanalysis .= " |bgzip -@ $threads_enable -c > $file_vcfgz_inanalysis " ;

    $log->andRun($cmd_vannot_inanalysis) ;
    
    $vcf_input_for_next_step = $file_vcfgz_inanalysis ;

    my $collist_pairworkflow_extract = [qw/is_denovo denovo_samples is_LOH LOH_samples/ , "in_analysis_$jobid"] ;
    my $collist_pairworkflow_coladd = [qw/is_denovo denovo_samples is_LOH LOH_samples in_analysis/] ;

    my $cmd_vslannot_paired = "$script_vslannot" ;
    $cmd_vslannot_paired .= " -a extract " ;
    $cmd_vslannot_paired .= " -f $file_vcfgz_inanalysis " ;
    $cmd_vslannot_paired .= " -e " . join("," , @$collist_pairworkflow_extract)  ;
    $cmd_vslannot_paired .= " -t " . join("," , qw/boolean text boolean text boolean/) ;
    $cmd_vslannot_paired .= " -c " . join("," , (map {"$_\_$jobid"} @$collist_pairworkflow_coladd, )) ;
    $cmd_vslannot_paired .= " -b varselect_variants " ;
    $cmd_vslannot_paired .= " -d $file_db" ;

#    $log->andRun($cmd_vslannot_paired) ;
    $log->write("enqueue $cmd_vslannot_paired", 1) ;
    $vslannotQueue->enqueue($cmd_vslannot_paired) ;
    $log->write("Paired case/control workflow finish" , 1) ;

} else {
    # none
    $log->write("None workflow start",1) ;
    my $file_gz_inanalysis = "$dir_results_output/inanalysis_$jobid.txt.gz" ;
    my $file_vcfgz_inanalysis = "$dir_results_output/inanalysis_$jobid.vcf.gz" ;

    my $cmd_get_annotate_file_for_none = "zcat $vcf_input_for_next_step | cut -f 1,2,4,5 | sed 's/\$/\\t1/g' | bgzip -c > $file_gz_inanalysis" ;
    $log->andRun($cmd_get_annotate_file_for_none) ;

    tabix_vcf($file_gz_inanalysis) ;

    my $cmd_vannot_inanalysis_none = "zcat $vcf_input_for_next_step | vcf-annotate " ;
    $cmd_vannot_inanalysis_none .= " -a $file_gz_inanalysis " ;
    $cmd_vannot_inanalysis_none .= " -c " . join("," , ('CHROM', 'POS','REF','ALT' , "INFO/in_analysis_$jobid") ) ;
    $cmd_vannot_inanalysis_none .= " -d key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='in_analysis_$jobid' " ;
    $cmd_vannot_inanalysis_none .= "2> $dir_log/stderr_vannot_inanalysis_$jobid.log " ;
    $cmd_vannot_inanalysis_none .= " |bgzip -@ $threads_enable -c > $file_vcfgz_inanalysis " ;

    $log->andRun($cmd_vannot_inanalysis_none)  ;
    tabix_vcf($file_vcfgz_inanalysis) ;

    $vcf_input_for_next_step = $file_vcfgz_inanalysis ;
    
    my $cmd_vslannot_nonemodel = "$script_vslannot" ;
    $cmd_vslannot_nonemodel .= " -a extract " ;
    $cmd_vslannot_nonemodel .= " -f $file_vcfgz_inanalysis " ;
    $cmd_vslannot_nonemodel .= " -e in_analysis_$jobid " ;
    $cmd_vslannot_nonemodel .= " -t boolean " ;
    $cmd_vslannot_nonemodel .= " -c in_analysis_$jobid " ;
    $cmd_vslannot_nonemodel .= " -b varselect_variants " ;
    $cmd_vslannot_nonemodel .= " -d $file_db " ;

#    $log->andRun($cmd_vslannot_nonemodel) ;
    $log->write("enqueue $cmd_vslannot_nonemodel" , 1) ;
    $vslannotQueue->enqueue($cmd_vslannot_nonemodel) ;
    $log->write("None workflow finish",1) ;
}

#close $CONFIG ;

$log->write("Thread finish" , 1) ;
$vslannotQueue->enqueue(undef) ;
$vslannot_thread->join ;

# CREATE VIEW
#=======================================
my $cmd_createview = "$dir_script/create_view_for_varselect.pl -j $jobid -d $file_db" ;
$log->andRun($cmd_createview, 1) ;

$log->write("VSanlz finish") ;

sub wf_vannot {
    my $file_annotation = shift ;
    my $tag = shift ;
    my $file_vcf_input = shift ;
    my $file_vcf_output = shift ;

    my $cmd_vannot = "vcf-annotate " ;
    $cmd_vannot .= " -a $file_annotation " ;
    $cmd_vannot .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/$tag")) ;
    $cmd_vannot .= " -d key=INFO,ID=$tag,Number=1,Type=Integer,Description='genetic model $tag' " ;
    $cmd_vannot .= "$file_vcf_input 2> $dir_log/stderr_vannot_workflow_$tag.log " ;
    $cmd_vannot .= " |bgzip -c > $file_vcf_output " ;

    $log->andRun($cmd_vannot);

    tabix_vcf($file_vcf_output) ;
}

sub running_workflow {
    my $script = shift ;
    my $file_vcf = shift ;
    my $file_ped = shift ;
    my $file_stderr = shift ;
    my $file_gz_output_wf = shift ;

    my $script_compare = shift ;
    my $db = shift || 0 ;

    my $output_tmp = "$file_gz_output_wf.tmp" ;

    # workflow step 1 of genetic model
    my $cmd_runworkflow = "$script " ;
    $cmd_runworkflow .= " -v $file_vcf " ;
    $cmd_runworkflow .= " -p $file_ped " ;
    $cmd_runworkflow .= " -d $db " if ($db) ;
    $cmd_runworkflow .= " 2> $file_stderr " ;
    $cmd_runworkflow .= " > $output_tmp" ;

    $log->andRun($cmd_runworkflow) ;

    if (-e $script_compare) {
	# workflow step 2 of genetic model
        my $cmd_step2 = "$script_compare " ;
        $cmd_step2 .= " -v $file_vcf " ;
        $cmd_step2 .= " -p $file_ped " ;
        $cmd_step2 .= " -c $output_tmp " ;
        $cmd_step2 .= " | bgzip -c > $file_gz_output_wf" ;
        $log->andRun($cmd_step2) ;

    } else { # without Compare.py, run step 1 ONLY
        my $cmd_bgzip_output_tmp = "bgzip -@ $threads_enable -c $output_tmp > $file_gz_output_wf" ;
        $log->andRun($cmd_bgzip_output_tmp) ;
    }

    tabix_vcf($file_gz_output_wf) ;
}


sub get_union {
    my $files = shift ;
    my $file_union = shift ;

    my $union_set = {} ;
    foreach my $file (keys %$files) {

	my $col1 = $files->{$file} ;

	my $SRC ;

	if ($file =~ /gz$/) {
	    open ($SRC,"zcat $file|") || die "can't open $file .\n$!" ;
	} else {
	    open ($SRC,"$file") || die "can't open $file .\n$!" ;
	}

	while(<$SRC>) {
	    chomp ;
	    my @data = split /\t/ ;
	    my $col = $col1 - 1 ;

	    if ($data[$col]) {
		my $tag = join("\t" , ($data[0],$data[1],$data[2],$data[3])) ;

		push @{$union_set->{$tag}} , $file ;

	    }
	}
	close $SRC ;
    }

    open (my $TGT,">$file_union.tmp") ;

    foreach my $tag (keys %$union_set) {

	my $data = join("|" , @{$union_set->{$tag}}) ;
	print $TGT "$tag\t1\t$data\n" ;
    }

    close $TGT ;


#    `vcf-sort $file_union.tmp |bgzip -c > $file_union` ;
    my $cmd_sort_output = "grep '^#' $file_union.tmp > $file_union.tmp.header ; grep -v '^#' $file_union.tmp | sort -V | cat $file_union.tmp.header - | bgzip -c > $file_union" ;
    $log->andRun($cmd_sort_output) ;
}
