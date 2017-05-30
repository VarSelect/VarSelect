#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use File::Path qw(make_path remove_tree);
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;
use DBI qw(:sql_types) ;

use threads ;
use threads::shared ;
use Thread::Queue ;


use IO::Handle ;
autoflush STDERR 1 ;
autoflush STDOUT 1 ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db          = "$dir_vsHOME/db" ;
my $dir_lib         = "$dir_vsHOME/lib" ;
my $dir_workflows   = "$dir_vsHOME/workflows" ;

my $script_vslannot = "$dir_script/vsl_annot.pl" ;

# options handle
my $opts = {} ;
getopts("v:p:j:n:kd:",$opts) ;

my $file_vcfgz_input = $opts->{v} ;
my $jobid = $opts->{j} ;
my $file_ped = $opts->{p} ;
my $threads_enable = $opts->{n} || 16 ;
my $threads_active = $threads_enable ;
my $file_db       = $opts->{d} ;

# jobs specific 
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log  = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $dir_temp_for_gemini_vep = "$dir_work/temp_for_gemini_vep" ;
make_path($dir_temp_for_gemini_vep , {chmod => 0755 ,}) unless (-e $dir_temp_for_gemini_vep) ;

my $file_config = "$dir_results_output/varselect.config" ;
my $file_log = "$dir_log/log_vsann2_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

# Threads Starting
my @thr ;

# InitAnnot VEP
#=======================================
my $file_prefix = "$dir_results_output/VsResult_$jobid" ;
my $file_vcf_vep = $file_prefix . "_vep.vcf" ;
my $file_vcfgz_vep = $file_prefix . "_vep.vcf.gz" ;
my $script_runvep = "$dir_script/run_vep.pl" ;
my $file_log_vep = "$dir_log/running_vep.$jobid.log" ;
my $file_vcfgz_pathway = "$file_prefix\_pathway.vcf.gz" ;

my $thread_vep = async {
    $log->write ("Running VEP start",1) ;
    my $cmd_run_vep = "$script_runvep -i $file_vcfgz_input -o $file_vcf_vep -t $jobid -n $threads_enable -l $dir_log " ;
    $log->andRun($cmd_run_vep) ;

    my $cmd_bgzip_vep = "bgzip -c -@ $threads_enable $file_vcf_vep > $file_vcfgz_vep" ;
    $log->andRun($cmd_bgzip_vep) ;

    tabix_vcf($file_vcfgz_vep) ;
    $log->write ("Running VEP finish",1) ;

    # InitAnnot VEP: Gemini load VEP
    #=======================================
    $log->write ("Gemini load vcf to db start",1) ;

    my $cmd_vep_geminiload = "gemini load " ;
    $cmd_vep_geminiload .= " -t VEP " ;
    $cmd_vep_geminiload .= " -v $file_vcfgz_vep " ;
    $cmd_vep_geminiload .= " -p $file_ped " ;
    $cmd_vep_geminiload .= " --cores $threads_enable " ;
    $cmd_vep_geminiload .= " --tempdir $dir_temp_for_gemini_vep " ;
    $cmd_vep_geminiload .= " $file_db 2> $dir_log/stderr_geminiload_$jobid.log" ;

    $log->andRun($cmd_vep_geminiload) ;
    $log->write("Gemini load vcf to db finish",1) ;

    # InitAnnot get pathways from variants
    #=======================================
    $log->write("Gemini pathways start" , 1) ;

    my $cmd_running_pathways = "$dir_script/pathway_parse.pl " ;
    $cmd_running_pathways   .= " -j $jobid " ;
    $cmd_running_pathways   .= " -d $file_db " ;
    $cmd_running_pathways   .= " -v $file_vcfgz_input " ;
    $cmd_running_pathways   .= " -l $dir_log " ;
    $cmd_running_pathways   .= " -o $file_vcfgz_pathway " ;
    $cmd_running_pathways	.= " 2> $dir_log/stderr_gemini_pathways_$jobid.log " ;
    $log->andRun($cmd_running_pathways) ;

    $log->write("Gemini pathways finish" , 1) ;

    # InitAnnot Gemini annotate varapp
    #=======================================
    $log->write("Extract 6 info fields from vcf for Varapp (AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR) start") ;
    my $cmd_gemini_annotate = "$script_vslannot -f $file_vcfgz_vep -a extract " ;
    $cmd_gemini_annotate .= " -e AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR " ; # extract from INFO, for Varapp compatiable
    $cmd_gemini_annotate .= " -t float,float,float,float,float,float " ;
    $cmd_gemini_annotate .= " -c AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR " ;
    $cmd_gemini_annotate .= " -o mean,mean,mean,mean,mean,mean " ;
    $cmd_gemini_annotate .= " -d $file_db " ;
    $cmd_gemini_annotate .= " -b variants " ;
    $cmd_gemini_annotate .= " 2> $dir_log/stderr.vslannot_forvarapp_$jobid.log" ;

    $log->andRun($cmd_gemini_annotate) ;
    $log->write("Extract 6 info fields from vcf for Varapp finish") ;

    # InitAnnot VEP: Rename to vep_variants
    #=======================================
    my $table_vep = "vep_variants" ;
    my $cmd_variants2_vepvar = "sqlite3 $file_db 'alter table variants rename to $table_vep' " ;
    $cmd_variants2_vepvar .= "2> $dir_log/stderr.tbl_rename_$table_vep\_$jobid.log " ;
    $log->andRun($cmd_variants2_vepvar , 1) ;
} ;

push @thr , $thread_vep ;

# InitAnnot running snpEff
#=======================================
my $file_vcfgz_snpeff = $file_prefix . "_snpeff.vcf.gz" ;
my $script_run_snpeff = "$dir_script/run_snpeff.pl" ;
my $file_log_snpeff = "$dir_log/running_snpeff.$jobid.log" ;
my $file_snpeffdb = $file_prefix ."_snpeff.db" ;

my $thread_snpeff = async {
    $log->write("Running snpEff start" ,1) ;
    my $cmd_run_snpeff = "$script_run_snpeff -i $file_vcfgz_input -o $file_vcfgz_snpeff -j $jobid -n $threads_active -d $file_snpeffdb" ;
    $log->andRun($cmd_run_snpeff) ;
    $log->write("Running snpEff finish",1) ;
} ;

push @thr, $thread_snpeff ;

# InitAnnot Running ANNOVAR
#=======================================
my $file_vcfgz_annovar = "$file_prefix\_annovar.vcf.gz" ;
my $file_log_running_annovar = "$dir_log/running_annovar.log" ;

my $thread_annovar = async {
    $log->write("Running ANNOVAR start" , 1) ;
    my $cmd_runannovar = "$dir_script/run_annovar.pl" ;
    $cmd_runannovar .= " -i $file_vcfgz_input " ;
    $cmd_runannovar .= " -o $file_vcfgz_annovar " ;
    $cmd_runannovar .= " -j $jobid " ;

    $log->andRun($cmd_runannovar) ;

    tabix_vcf($file_vcfgz_annovar) ;
    $log->write("Running ANNOVAR finish" , 1) ;
} ;

push @thr , $thread_annovar ;


foreach (@thr) {
    $_ -> join ;
}

# InitAnnot VarSelect: CREATE TABLE varselect_info
#=======================================

my $table_varselectinfo = 'varselect_info' ;
$log->write("CREATE & INSERT TABLE $table_varselectinfo start",1) ;

my $dbh = DBI->connect("dbi:SQLite:dbname=$file_db","","");
$dbh->{RaiseError} = 1 ;

# CREATE TABLE varselect_info
my $sql_varsel_info = "CREATE TABLE IF NOT EXISTS $table_varselectinfo ( " ;
$sql_varsel_info .= " sn INTEGER " ;
$sql_varsel_info .= " , vcfgz_for_analysis BLOB " ; #vcf_input
$sql_varsel_info .= " , sample_list TEXT " ;
$sql_varsel_info .= " , caller_list TEXT " ;
$sql_varsel_info .= " , multicaller BOOLEAN " ;
$sql_varsel_info .= ")" ; 

$log->write($sql_varsel_info,1) ;
$dbh->do($sql_varsel_info) ;

# Load config
my $config = {} ;
open (my $SRC_config , "$file_config") ;
while (<$SRC_config>) {
    chomp ;
    my ($idx,$val) = split /\t/ ;
    $config->{$idx} = $val ;
}
close $SRC_config ;

# Load bgziped VCF blob from file
my $buff ;
my $vcf_blob = [] ;
#my $limit_datasize = 1000000000 ; # 1G
my $limit_datasize = 100000000 ; # 100M
#my $limit_datasize = 1000000 ; # 1M
my $reading_chunk_size = 1024 ;
my $pivot_insert_block = 0 ;
my $insert_sum = 0 ;

open (my $SRC_vcfinput,"$file_vcfgz_input") ;
while (read $SRC_vcfinput, $buff , $reading_chunk_size) {
    $insert_sum += $reading_chunk_size ;
    if ($insert_sum >= $limit_datasize) {
	$insert_sum = 0 ;
	$pivot_insert_block ++ ;
    }

    $vcf_blob->[$pivot_insert_block] .=  $buff ;
}
close $SRC_vcfinput ;

# Save to varselect_info Table
my $multicaller_mode = ($config->{multicaller})?1:0 ;
my $sql_insert_varselinfo = "INSERT INTO $table_varselectinfo " ;
$sql_insert_varselinfo .= " (sn, vcfgz_for_analysis , sample_list, caller_list, multicaller) values " ;
$sql_insert_varselinfo .= " (?, ?, ?, ?, ?)" ;

my $sth = $dbh->prepare($sql_insert_varselinfo) ;

for (my $i = 0 ; $i <= $pivot_insert_block ; $i ++) {
    $sth->bind_param(1, $i) ;
    $sth->bind_param(2, $vcf_blob->[$i] , SQL_BLOB) ;
    $sth->bind_param(3, "$config->{sample_list}" ) ;
    $sth->bind_param(4, "$config->{caller_list}" ) ;
    $sth->bind_param(5, $multicaller_mode) ;
    $sth->execute() ;
}

$log->write("Store config and vcf file into Table $table_varselectinfo") ;

$log->write("CREATE & INSERT TABLE $table_varselectinfo finish",1) ;


# InitAnnot snpEff processing
#=======================================

my $table_snpeff = "snpeff_variants" ;
my $file_sql_snpeff_create = "$dir_work/create_s_vars.sql" ;
my $file_sql_attach = "$dir_work/snpeff_atach.sql" ;

my $cmd_get_snpeff_schema = "sqlite3 $file_snpeffdb '.schema variants' > $file_sql_snpeff_create " ;
$cmd_get_snpeff_schema .= " 2> $dir_log/stderr_schema_variants_$jobid.log" ;
$log->andRun($cmd_get_snpeff_schema) ;

my $cmd_sed_replace  = "sed -i 's/variants/snpeff_variants/' $file_sql_snpeff_create" ;
$cmd_sed_replace .= " 2> $dir_log/stderr_replace_sql_variants_$jobid.log" ;
$log->andRun($cmd_sed_replace) ;

my $cmd_sed_replace2 = "sed -i 's/idx/sidx/' $file_sql_snpeff_create" ;
$cmd_sed_replace2 .= " 2> $dir_log/stderr_replqce_sql_idx_$jobid.log " ;
$log->andRun($cmd_sed_replace2) ;

my $cmd_create_table = "sqlite3 $file_db < $file_sql_snpeff_create" ;
$cmd_create_table .= " 2> $dir_log/stderr_tbl_create_snpeffvariants_$jobid.log"  ;
$log->andRun($cmd_create_table) ;

open (my $TGT_attachsql , ">$file_sql_attach") ;
print $TGT_attachsql "attach database '$file_snpeffdb' as dbsnpeff ;" ;
print $TGT_attachsql "insert into snpeff_variants select * from dbsnpeff.variants ;" ;
close $TGT_attachsql ;

my $cmd_attach_insert = "sqlite3 $file_db < $file_sql_attach" ;
$cmd_attach_insert .= " 2> $dir_log/stderr_attachinsert_$jobid.log" ;
$log->andRun($cmd_attach_insert) ;


# InitAnnot ANNOVAR processing
#=======================================
# Create TABLE annovar_variants
my $table_annovar = "annovar_variants" ;
$log->write("Create table annovar_variants start" , 1) ;

my $sql_create_table_anvvars = "CREATE TABLE $table_annovar as select variant_id,chrom,start,\"end\",ref,alt from vep_variants " ;
my $cmd_addtbl_anvvars = "sqlite3 $file_db '$sql_create_table_anvvars'" ;
$cmd_addtbl_anvvars .= " 2> $dir_log/stderr_tbl_create_$table_annovar\_$jobid.log" ;
$log->andRun($cmd_addtbl_anvvars) ;

my $sql_create_index_anvvars = "CREATE INDEX anvvars_vid_idx ON $table_annovar(variant_id) " ;
my $cmd_addidx_anvvars = "sqlite3 $file_db '$sql_create_index_anvvars' " ;
$cmd_addidx_anvvars .= " 2> $dir_log/stderr_idx_create_$table_annovar\_$jobid.log " ;
$log->andRun($cmd_addidx_anvvars) ;
$log->write("Create table annovar_variants finish" , 1) ;


# Exetract INFO from ANNOVAR
$log->write("Extract info fields from ANNOVAR finish",1) ;
my $annovar_header_parser = Vcf_parser->new (file => $file_vcfgz_annovar) ;

my $list_extract = [] ;
my $list_type = [] ;
my $list_column = [] ;
my $list_operation = [] ;

foreach my $info_tag (sort keys %{$annovar_header_parser->{header}->{INFO}} ) {
    next if ($info_tag eq 'lines') ;  # dirty hack for Vcf_parser bug

    my $desc = $annovar_header_parser->{header}->{INFO}->{$info_tag}->{Description} ;

    if ($desc =~ /provided by ANNOVAR/) {
	my $col_tag = $info_tag ;
	$col_tag =~ s/\-/\_/g ;
	$col_tag =~ s/\+/x/g ;

	next if ($info_tag eq 'esp6500siv2_all') ; # dirty hack for duplicate INFO tag, there is ESP6500siv2_ALL already

	push @$list_extract , $info_tag ;
	push @$list_type , 'text' ;
	push @$list_column , "$col_tag" ;
	push @$list_operation , 'list' ;
   }
}

my $cmd_gannot_annovar = "$script_vslannot " ;
$cmd_gannot_annovar .= "-f $file_vcfgz_annovar -a extract " ;
$cmd_gannot_annovar .= " -e " . join("," , @$list_extract) ;
$cmd_gannot_annovar .= " -t " . join("," , @$list_type) ;
$cmd_gannot_annovar .= " -c " . join("," , map {"'" . $_ . "'"} @$list_column) ;
$cmd_gannot_annovar .= " -d $file_db " ;
$cmd_gannot_annovar .= " -b $table_annovar " ;
$cmd_gannot_annovar .= " 2> $dir_log/stderr_gemini_annovar_$jobid.log " ;
#$cmd_gannot_annovar .= " -o " . join("," , @$list_operation) ;
#$cmd_gannot_annovar .= " $file_db 2> $dir_log/stderr.gannot_annovar.$jobid.log" ;
$log->andRun($cmd_gannot_annovar) ;

$log->write("Extract info fields from ANNOVAR finish",1) ;

#
# Create TABLE varselect_variants
#=======================================
$log->write("Create table varselect_variants start" , 1) ;

my $table_varselect = "varselect_variants" ;
my $sql_create_table_vslvars = "CREATE TABLE $table_varselect as select variant_id,chrom,start,\"end\",ref,alt from vep_variants " ;
my $cmd_addtbl_vslvars = "sqlite3 $file_db '$sql_create_table_vslvars'" ;
#$log->write($cmd_addtbl_vslvars , 1) ;
$log->andRun($cmd_addtbl_vslvars , 1) ;

my $sql_index_table_vslvars = "CREATE INDEX vslvars_vid_idx ON $table_varselect(variant_id)" ;
my $cmd_addidx_vslvars = "sqlite3 $file_db '$sql_index_table_vslvars' " ;
$cmd_addidx_vslvars .= " 2> $dir_log/stderr_idx_create_vslvars_$table_varselect\_$jobid.log " ;
$log->andRun($cmd_addidx_vslvars,1) ;

$log->write("Create table varselect_variants finish" , 1) ;

# vsl_annotate pathways
#=======================================
#my $file_vcfgz_pathway = "$file_prefix\_pathway.vcf.gz" ;

my $cmd_vsl_annot_pathway = "$script_vslannot " ;
$cmd_vsl_annot_pathway .= " -j $jobid " ;
$cmd_vsl_annot_pathway .= " -f $file_vcfgz_pathway " ;
$cmd_vsl_annot_pathway .= " -a extract " ;
$cmd_vsl_annot_pathway .= " -e pathway " ;
$cmd_vsl_annot_pathway .= " -t text " ;
$cmd_vsl_annot_pathway .= " -c pathway " ;
$cmd_vsl_annot_pathway .= " -d $file_db " ;
$cmd_vsl_annot_pathway .= " -b $table_varselect " ;
$cmd_vsl_annot_pathway .= " 2> $dir_log/stderr_vsl_pathway_$jobid.log " ;
#$cmd_vsl_annot_pathway .= " -o list " ;

$log->andRun($cmd_vsl_annot_pathway,1) ;

# InitAnnot GO parsing 
#=======================================
my $file_vcfgz_GO = "$file_prefix\_go.vcf.gz " ;
$log->write("GO parsing start",1) ;
my $cmd_parse_GO = "$dir_script/GO_parse.pl " ;
$cmd_parse_GO   .= " -j $jobid " ;
$cmd_parse_GO   .= " -d $file_db " ;
$cmd_parse_GO   .= " -v $file_vcfgz_input " ;
$cmd_parse_GO   .= " -l $dir_log " ;
$cmd_parse_GO   .= " -o $file_vcfgz_GO " ;
$log->andRun($cmd_parse_GO) ;

# InitAnnot GO processing 
#=======================================
my $cmd_vsl_annot_GO = "$script_vslannot " ; 
$cmd_vsl_annot_GO .= " -j $jobid " ;
$cmd_vsl_annot_GO .= " -f $file_vcfgz_GO " ;
$cmd_vsl_annot_GO .= " -a extract " ;
$cmd_vsl_annot_GO .= " -e GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function " ;
$cmd_vsl_annot_GO .= " -t text,text,text " ;
$cmd_vsl_annot_GO .= " -c GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function " ;
$cmd_vsl_annot_GO .= " -d $file_db " ;
$cmd_vsl_annot_GO .= " -b $table_varselect " ;
$cmd_vsl_annot_GO .= " 2> $dir_log/stderr_vsl_GO_$jobid.log" ;

$log->andRun($cmd_vsl_annot_GO,1) ;

$log->write("GO parse finish",1) ;


# CREATE VIEW
#=======================================
my $cmd_createview = "$dir_script/create_view_for_varselect.pl -j $jobid -d $file_db 2> $dir_log/stderr_view_create_$jobid.log" ;
$log->andRun($cmd_createview, 1) ;

