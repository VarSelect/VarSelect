#!/usr/bin/perl
#
use strict ;
use Getopt::Std ;
use File::Path qw(make_path remove_tree);
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


my $dir_vep = $Setting->{dir_vep} ;
my $exe_vep = $Setting->{vep} ;
my $file_ref = $Setting->{file_ref} ;
my $dir_vep_cache = $Setting->{dir_vep_cache} ;
my $file_dbNSFP = $Setting->{file_dbNSFP} ;



# Options handling
my $opts = {} ;
getopts("i:o:l:n:t:", $opts) ;

my $file_input		= $opts->{i} ;
my $file_output		= $opts->{o} ;	# should be vcf
my $dir_log		= $opts->{l} ;
my $threads_available	= $opts->{n} ;
my $ts			= $opts->{t} || time ;

my $file_log = "$dir_log/running_vep.$ts.log" ;
my $file_stderr = "$dir_log/stderr.vep.$ts.log" ;

open (my $LOG,">$file_log") ;


    my $cmd_vep = "  export PATH=$dir_vep:\$PATH && $exe_vep " ;
    $cmd_vep .= " -i $file_input " ;
    $cmd_vep .= " --dir $dir_vep_cache " ;
    $cmd_vep .= " --cache " ;
    $cmd_vep .= " --offline " ;
    $cmd_vep .= " --fasta $file_ref " ;
    $cmd_vep .= " --port 3337 " ;
    $cmd_vep .= " --fork $threads_available " ;
    $cmd_vep .= " --everything " ;
    $cmd_vep .= " --total_length " ;
    $cmd_vep .= " --plugin dbNSFP,$file_dbNSFP" ;
        $cmd_vep .= ",Ancestral_allele,AltaiNeandertal,Denisova" ;  #column of dbNSFP 
        $cmd_vep .= ",Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid" ;
        $cmd_vep .= ",SIFT_score,SIFT_converted_rankscore,SIFT_pred" ;
        $cmd_vep .= ",LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega" ;
        $cmd_vep .= ",MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE" ;
        $cmd_vep .= ",MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_pred" ;
        $cmd_vep .= ",FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred" ;
        $cmd_vep .= ",PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred" ;
        $cmd_vep .= ",fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group" ;
        $cmd_vep .= ",MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred" ;
        $cmd_vep .= ",MetaLR_score,MetaLR_rankscore,MetaLR_pred" ;
        $cmd_vep .= ",Reliability_index" ;
        $cmd_vep .= ",integrated_fitCons_score,integrated_fitCons_score_rankscore,integrated_confidence_value" ;
        $cmd_vep .= ",GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value" ;
        $cmd_vep .= ",H1-hESC_fitCons_score,H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value" ;
        $cmd_vep .= ",HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value" ;
        $cmd_vep .= ",GERP++_NR,GERP++_RS,GERP++_RS_rankscore" ;
	$cmd_vep .= ",phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore" ;
	$cmd_vep .= ",phyloP20way_mammalian,phyloP20way_mammalian_rankscore" ;
	$cmd_vep .= ",phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore" ;
	$cmd_vep .= ",phastCons20way_mammalian,phastCons20way_mammalian_rankscore" ;
        $cmd_vep .= ",SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore" ;
        $cmd_vep .= ",1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF" ;
        $cmd_vep .= ",TWINSUK_AC,TWINSUK_AF" ;
        $cmd_vep .= ",ALSPAC_AC,ALSPAC_AF" ;
        $cmd_vep .= ",ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF" ;
        $cmd_vep .= ",ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF" ;
        $cmd_vep .= ",clinvar_rs,clinvar_clnsig,clinvar_trait" ;
        $cmd_vep .= ",Interpro_domain" ;
        $cmd_vep .= ",GTEx_V6_gene,GTEx_V6_tissue" ;
	$cmd_vep .= ",cds_strand" ;
#	$cmd_vep .= ",'hg18_pos(1-based)'" ;
	$cmd_vep .= ",Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2" ;	# academic version only
	$cmd_vep .= ",Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred" ;	# academic version only
	$cmd_vep .= ",Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred" ;	# academic version only
	$cmd_vep .= ",Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore" ;	# academic version only 
        $cmd_vep .= ",Eigen-raw,Eigen-phred,Eigen-raw_rankscore,Eigen-PC-raw,Eigen-PC-raw_rankscore" ; # academic version only 
	$cmd_vep .= ",GenoCanyon_score,GenoCanyon_score_rankscore" ;				# academic version only


$cmd_vep .= " -o $file_output " ;
$cmd_vep .= " --vcf " ;
$cmd_vep .= " --force_overwrite " ;
$cmd_vep .= " &>$dir_log/stderr.vep.$ts " ;

tlog($cmd_vep) ;
`$cmd_vep` ;
close $LOG ;

sub timelog {
    my $text = shift ;
    my $time = localtime(time) ;

    return "$time\t$text" ;
}

sub logit {
    my $msg = shift ;
    my $fh = shift ;

    print $fh "$msg\n" ;
}

sub tlog {  # timelog + logit
    my $msg = shift ;
    logit(timelog("$msg"),$LOG) ;
}

