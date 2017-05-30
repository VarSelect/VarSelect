#!/usr/bin/perl
#
# Perl partial port of gemini_annotate.py
# enable annotate to non-variants table
#
# gmn_annotate.pl -f other.vcf.gz 
#                 -c other_allele_freq 
#                 -t float 
#                 -e AF 
#                 -b table
#                 -d db
#                 -k Chunk size

use strict ;
use DBI ;

use Getopt::Std ;
use File::Basename ;
use DBI;
use Data::Dumper ;
use Cwd qw(abs_path) ;
use lib dirname(dirname( abs_path($0) )) . '/lib' ;
use VarSelect ;
use VarSelect::Vcf_parser ;

my $DEFAULT_CHUNK_SIZE = 100000 ;

my $opts = {} ;
getopts("j:f:a:e:t:c:o:d:b:k:",$opts) ;

my $jobid		= $opts->{j} ;
my $file_vcfgz_input	= $opts->{f} ;
#my $type_ganno		= $opts->{a} || 'extract' ;  # should support extract, boolean , count but now extract only
my $list_cols_extract	= [split(/\,/, $opts->{e})] ;
my $list_coltypes	= [split(/\,/, $opts->{t})] ;
my $list_cols		= [split(/\,/, $opts->{c})] ;
my $list_col_operations = [split(/\,/, $opts->{o})] ; 
my $file_db		= $opts->{d} || die "-d file_db is required!" ;
my $target_table	= $opts->{b} || die "-b target_table is required!" ;
my $CHUNK_SIZE		= $opts->{k} || $DEFAULT_CHUNK_SIZE ;

my $cmd_tabix_input_anyway = "tabix -f -s 1 -b 2 -e 2 -c # $file_vcfgz_input " ;
`$cmd_tabix_input_anyway` ;


my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","");

add_requested_columns($list_cols,$list_coltypes,$target_table) ;

my $sql_update_table = "update $target_table SET " ;
$sql_update_table .= join (" , " , map { "'$_'=?" } @$list_cols) ;
$sql_update_table .= " where variant_id = ? " ;

print "$sql_update_table\n" ;
my $sth_updatetable = $dbh->prepare("$sql_update_table") ;

get_annotate_variant($dbh, $target_table , $file_vcfgz_input) ;


sub add_requested_columns {
    my $col_names = shift ;
    my $col_types = shift ;
    my $table = shift ;

    my $i = 0 ;

    foreach my $col_name (@$col_names) {
	my $col_type = $col_types->[$i] ;

	if ( col_is_in_table($col_name,$target_table,$file_db)) { 
	    print STDERR "Column \"$col_name\"  already exists in $target_table table. Overwriting values.\n" ;

	} else {
	    my $sql_add_col = "ALTER table $table ADD COLUMN '$col_name' $col_type DEFAULT NULL" ;
	    print STDERR "$sql_add_col \n" ;

	    $dbh->do("$sql_add_col") ;
	}

	$i ++ ;
    }
}


sub col_is_in_table {
    my $col  = shift ;
    my $table  = shift ;
    my $db = shift ;

    my $sql_PRAGMA = "PRAGMA table_info($table) " ;

    open (my $SRC_tableinfo,"sqlite3 $db '$sql_PRAGMA' |") ;
    while (<$SRC_tableinfo>) {
	chomp ;
	my ($sn,$col_name,$col_type,$othercols) = split(/\|/ , $_) ;

	return 1 if ($col_name eq $col) ;
    }
    close $SRC_tableinfo ;

    return 0 ;
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

sub get_annotate_variant {

    my $dbh = shift ;
    my $table = shift ;
    my $file_vcfgz_input = shift ;
    my $list_to_extract = shift ;

    my $i = 0 ;
    my $sql_select_table = "select chrom, start, end, ref, alt, variant_id from $table " ;
    my $sth_select_table = $dbh->prepare($sql_select_table) ;
    $sth_select_table->execute() ;

    my $data_tobe_processed = {} ;
    while (my $row = $sth_select_table->fetchrow_hashref) {
	$i ++ ;

	$row->{converted_chrom} = CHROM_nochr($row->{chrom}) ;
	$row->{start1} = $row->{start} + 1 ;  # change gemini 0-based to vcf 1-based
	$row->{end1} = $row->{end} + 1 ;

	push @{ $data_tobe_processed->{$row->{converted_chrom}} } , $row ;

	if ($i == $CHUNK_SIZE) {
	    finding_match($data_tobe_processed , $file_vcfgz_input) ; # and update

	    # initialize for next chunk
	    $i = 0 ;
	    $data_tobe_processed = {} ;
	}
    }

    finding_match($data_tobe_processed, $file_vcfgz_input) ; # for data after last chumk
}

sub finding_match {
    my $data_tobe_processed = shift ;
    my $file_vcfgz_input = shift ;

    $dbh->begin_work ;
    my $match = 0 ;

    my $vcf_region_string ;
    my $flag_first_tabix_with_header = 1 ;
    foreach my $chr_tobe_processed (sort keys %$data_tobe_processed) {
	my $frow = $data_tobe_processed->{$chr_tobe_processed}->[0] ; # first row
	my $lrow = $data_tobe_processed->{$chr_tobe_processed}->[-1] ; # last row
	print "$frow->{chrom}, $frow->{start}, $frow->{ref}, $frow->{alt}, $frow->{variant_id}\n" ;
	print "$lrow->{chrom}, $lrow->{start}, $lrow->{ref}, $lrow->{alt}, $lrow->{variant_id}\n" ;

	# Load VCF in the region
	my $cmd_tabix = "tabix " ;
	$cmd_tabix .= " -h " if ($flag_first_tabix_with_header) ;
	$cmd_tabix .= " $file_vcfgz_input $frow->{converted_chrom}:$frow->{start}-$lrow->{end} " ;

	$vcf_region_string .= `$cmd_tabix` ;

    }

    my $vcf_records = {} ;
    my $vcf_region_parser = Vcf_parser->new (string => $vcf_region_string) ;
    while (my $var = $vcf_region_parser->next_var) {
	$var->{converted_CHROM} = CHROM_nochr($var->{CHROM}) ;
	my $var_tag = join("_" , map { $var->{$_} } qw/converted_CHROM POS REF ALT/ ) ;
	$vcf_records->{$var_tag} = $var ;
    }

    my $count_realupdate = 0 ;
    my $count_skip = 0 ;
    foreach my $chr (sort keys %$data_tobe_processed) {
	foreach my $data (@{$data_tobe_processed->{$chr}}) {
	    my $data_tag = join("_", map {$data->{$_}} qw/converted_chrom start1 ref alt/ ) ;

	    my $vid = $data->{variant_id} ;
	    if (exists $vcf_records->{$data_tag}) { # match
		$match++ ;
		my $matched_var = $vcf_records->{$data_tag} ;

		if (update_variant($sth_updatetable, $list_cols_extract, $list_cols , $matched_var,$vid) ) {
		    # real update
		    $count_realupdate ++ ;
		    
		} else {
		    # empty in vcf-INFO-value, skip
		    $count_skip ++ ;
		}

		my $window_verbose = $CHUNK_SIZE / 10 ;
		if ($match % $window_verbose == 0) {
		    my $ts = localtime(time) ;
		    print "$ts\tupdate $match variants\t$count_realupdate : $count_skip\n" ;
		    $count_realupdate = 0 ;
		    $count_skip = 0 ;
		}
	    }
	}
    }

    print "match $match records\n" ;
    $dbh->commit ;
}

sub update_variant {
    my $sth_updatetable = shift ;
    my $list_cols_extract = shift ;
    my $list_cols = shift ;
    my $var = shift ;
    my $vid = shift ;

    my $i = 0 ;
    my $var_tobe_execute = [] ;
    my $flag_needto_update = 0 ;
    foreach my $col (@$list_cols_extract) {
	push @$var_tobe_execute , $var->{INFO}->{$col} ;
	$flag_needto_update = 1 if ($var->{INFO}->{$col} && $var->{INFO}->{$col} ne '.') ;
	$i ++ ;
    }
    push @$var_tobe_execute , $vid ;
    
#    print join("\t",@$var_tobe_execute) . "\n";

    if ($flag_needto_update) {
	$sth_updatetable->execute(@$var_tobe_execute) ;
	return 1 ;
    } else {
	return 0 ;
    }
}
