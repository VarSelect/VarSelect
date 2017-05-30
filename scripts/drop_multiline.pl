#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;
use Data::Dumper ;

# default global setting
my $debug = 0 ;
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
my $opts = {} ;
getopts("j:i:", $opts) ;
my $jobid      = $opts->{j} ;
my $file_input = $opts->{i} || die "-i VCF file is required" ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;


my ($fname_input, $fpath_input , $fileext_input) = fileparse($file_input, qw/.vcf .vcf.gz/) ;

$dir_log = "./" if ($debug) ;
my $file_log = "$dir_log/drop_multiline_vcf_$fname_input\_$jobid.log" ;

open (my $LOG, ">$file_log") ;
print $LOG "FILE\t$file_input\n" ;
print $LOG "LOG\t$file_log\n" ;
print $LOG "==========\n" ;

my $vcf_parser = Vcf_parser->new (file => $file_input) ;
print $vcf_parser->header() ;
my $samples = $vcf_parser->{samples} ;

my $prevar ;
my $count = {} ;
my $flag_new =0 ;
while (my $var = $vcf_parser->next_var) {
    $count->{all} ++ ;

    if ($prevar && $var->{POS} == $prevar->{POS} && $var->{CHROM} eq $prevar->{CHROM}) {
	$flag_new = 0 ;

	print $LOG "DROP\t$prevar->{line}\n" ;
	print $LOG "DROP\t$var->{line}\n\n" ;

	$count->{multiline} ++ ;

    } elsif ($prevar) {

	if (validate_format_num($prevar) ) {
	    if ($flag_new) {
		my $output_line = $prevar->{line} ;
		$output_line =~ s/^chrM/MT/ ;
		$output_line =~ s/^chr// ;
		print "$output_line\n" ;
	    }

	} else {
	    if ($flag_new) {
		print $LOG "DROP\t$prevar->{line}\n\n"  ;
		$count->{format_num_invalid} ++ ;
	    }
	}

	$flag_new = 1 ;
    }

    $prevar = $var ;
}

# final var
if(validate_format_num($prevar) )  {
    if ($flag_new) {
	my $output_line = $prevar->{line} ;
	$output_line =~ s/^chrM/MT/ ;
	$output_line =~ s/^chr// ;
	print "$output_line\n" ;
    }
} else {
    print $LOG "DROP\t$prevar->{line}\n\n" ;
    $count->{format_num_invalid} ++ ;
}

foreach my $idx (sort keys %$count) {
    my $perc = sprintf ("%5.2f" , $count->{$idx} / $count->{all} * 100) if ($count->{all}) ;
    print $LOG "STAT\t$idx\t$count->{$idx}\t$perc %\n" ;
}


sub validate_format_num {
    my $var = shift ;

    my $output = 1 ;

    my $alt_num = $var->get_alt_allele_num ;
    my $gt_combine_num = $var->get_gt_combination_num ;

    if ($var->{CHROM} eq 'Y' || $var->{CHROM} eq 'chrY' || $var->{CHROM} eq 'MT' || $var->{CHROM} eq 'chrM') {
	$gt_combine_num = $var->get_gt_combination_num(1) ;
    }

    foreach my $sample (@{$var->{samples}}) {
	my $format_value_pair_of_sample = $var->{sample_val}->{$sample} ;

	foreach my $tag (keys %$format_value_pair_of_sample ) {
	    my $value_num = scalar split("\,",$format_value_pair_of_sample->{$tag}) ;
	    
	    my $header_setting = $vcf_parser->{header}->{FORMAT}->{$tag}->{Number} ;

	    my $expect_val_number  ;
	    if ($header_setting eq 'G') {
		$expect_val_number = $gt_combine_num ;
	    
	    } elsif ($header_setting eq 'A') {
		$expect_val_number = $alt_num ;
	    
	    } elsif ($header_setting eq 'R') {
		$expect_val_number = $alt_num + 1 ;

	    } elsif ($header_setting eq '.') {
		$expect_val_number = '.' ;

	    } else {
		$expect_val_number = $header_setting ;
	    }
	    
	    $output = 0 if ($header_setting ne '.' && $value_num != $expect_val_number) ;

	}
    }

    return $output ;
}
