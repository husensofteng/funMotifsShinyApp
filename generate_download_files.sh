downloadhtmlfile="../../download.html.funmotifs"
db_name="funmotifsdb"
tables_db_filename="table_names_funmotifsdb.txt"
output_dir="www/"
output_dir2="bed/"
motifs_table="motifs"
all_tissues_table="all_tissues"
col_names_motifs="chr,motifstart,motifend,name,score,strand,pval"
col_names_all_tissues="blood,brain,breast,cervix,colon,esophagus,kidney,liver,lung,myeloid,pancreas,prostate,skin,stomach,uterus"
col_names_tissue_table="fscore, chromhmm, contactingdomain, dnase__seq, fantom, loopdomain, numothertfbinding, othertfbinding, replidomain, tfbinding, tfexpr"
limit_stmt=""

cd $output_dir$output_dir2
#get tissue names from the database
psql -d $db_name -c "SELECT distinct(table_name) as table_name FROM information_schema.tables where table_schema='public' and table_name not like '%motif%' order by table_name;" | awk '$0!~"table" && $0!~"---" && $0!~"row" && $0!="" && $0!~"all_tissues"'> $tables_db_filename

echo "" > $downloadhtmlfile
#GENERATE filtered output (potential functional motifs)
echo "<h4><p><b>Potential functional motifs per tissue type</b></p></h4><hr/><ul>" >> $downloadhtmlfile 
#export a file for each tissue table
while read -r name
do
        echo $output_dir$name.funmotifs
        psql -d $db_name -A -F $'\t' --pset footer -c "select $col_names_motifs,$col_names_tissue_table from $motifs_table,$name where $motifs_table.mid=$name.mid and dnase__seq>0 and (tfexpr>0 or tfexpr='nan') and (fscore>2.55 or (tfbinding>0 and tfbinding!='NaN'))$limit_stmt order by chr,motifstart,motifend;" -o $name.funmotifs_sorted.bed
bgzip $name.funmotifs_sorted.bed
tabix -p bed --skip-lines 1 $name.funmotifs_sorted.bed.gz

#tar -czvf $name.funmotifs.bed.tar.gz $name.funmotifs.bed.gz
echo "<li><a href='$output_dir2$name.funmotifs_sorted.bed.gz' download>$name.funmotifs</a> (<a href='$output_dir2$name.funmotifs_sorted.bed.gz.tbi'>.tbi</a>)<br/></li>" >> $downloadhtmlfile
done < $tables_db_filename
echo "</ul>" >>  $downloadhtmlfile



#get tissue names
tissues=$(
while read line
do
  echo -n ",${line}"
done < $tables_db_filename)

########variants######
echo "<hr/><h4><p><b>Annotated variatns using funMotifsDB</b></p></h4><hr/>
<ul>
<li>Putative regulatory variants from the 1000 Genomes project: <a href='../1000GenomeProject/functional/candidates1000G_tfexpr_entropyabs0.3_dnase_tfbindOrFscore2.55_grouped.tar.gz' download>candidate_variants_1000G.tar.gz</a>.</li>
<li>Putative regulatory variants from the GTEx project (eQTLs): <a href='../GTEx_eQTLs/GTEx_eQTLs_candidates.tar.gz' download>GTEx_eQTLs_candidates.tar.gz</a>.</li>
<li>Putative regulatory variants from GWAS : <a href='../GWAS_data/GWAS_SNPs_LD_candiadtes.tar.gz' download>GWAS_SNPs_LD_candiadtes.tar.gz</a>.</li>

</ul>" >> $downloadhtmlfile

##########get motif tables
tables_motifs_db_filename="motifs_tables.txt"
psql -d $db_name --pset footer -c "SELECT distinct(table_name) as table_name FROM information_schema.tables where table_schema='public' and table_name like '%motif%' and table_name like '%chr%' order by table_name;" | awk '$0~"chr" && $0~"motifs"' > $tables_motifs_db_filename

echo "<br/><hr/><h4><p><b>All scored motifs per tissue type</b></p></h4><hr/>" >> $downloadhtmlfile

#write the column headers to an outputfile
psql -d $db_name -A -F $'\t' --pset footer -c "select $col_names_motifs$tissues from $motifs_table,$all_tissues_table where $motifs_table.mid=$all_tissues_table.mid limit 0;" > $all_tissues_table.funmotifs_sorted.bed

while read -r name
do
        echo $output_dir$name
psql -d $db_name -A -F $'\t' -t -c "select $col_names_motifs$tissues from $name,$all_tissues_table where $name.mid=$all_tissues_table.mid$limit_stmt order by chr, motifstart, motifend;" >> $all_tissues_table.bed
done < $tables_motifs_db_filename
#tar -czvf $all_tissues_table.bed.tar.gz $all_tissues_table.bed
sort -k1,1n -k2,2n -k3,3n $all_tissues_table.bed >> $all_tissues_table.funmotifs_sorted.bed
bgzip $all_tissues_table.funmotifs_sorted.bed
tabix -p bed --skip-lines 1 $all_tissues_table.funmotifs_sorted.bed.gz

echo "<ul><li><a href='$output_dir2$all_tissues_table.funmotifs_sorted.bed.gz' download>All motifs: fscore per tissue type</a> (<a href='$output_dir2$all_tissues_table.funmotifs_sorted.bed.gz.tbi'>.tbi</a>)</li></ul><br/>" >> $downloadhtmlfile




#####data_files.tar.gz
echo "<hr/><h4><p><b>Data files used to build funMotifsDB</b></p></h4><hr/><ul><li><a href='../datafiles.tar.gz' download>datafiles.tar.gz</a>: information about the contents of this archive is listed <a href='https://github.com/husensofteng/funMotifs/tree/master/ReadMe'>here</a>.</li></ul>" >> $downloadhtmlfile



