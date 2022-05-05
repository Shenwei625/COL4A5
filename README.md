# COL4A5
### 数据下载
+ GnomAD
```bash
# 利用gsutil下载数据
pip3 install gsutil

# 检索gnomAD数据库中的数据
gsutil ls gs://gcp-public-data--gnomad/release/
gsutil ls gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/
gsutil du gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz # 查看文件大小
gsutil cp gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz . # 下载

# 利用bcftools筛选下载的数据

```

+ Clinvar
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20220430.vcf.gz

# md5
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20220430.vcf.gz.md5
cat clinvar_20220430.vcf.gz.md5
# 8aeb3aa191632bf06556948c853b5f06  /panfs/pan1/clintest/ftp_test_MGG_REPORT/vcf/vcf_GRCh38/clinvar_20220430.vcf.gz 待测试文件路径有点问题，修改一下
sed -i 's/\/.\+\///g' clinvar_20220430.vcf.gz.md5
md5sum --check clinvar_20220430.vcf.gz.md5 # md5检验
```
+ 筛选出col4a5中的单碱基替换位点
```bash
gzip -d clinvar_20220430.vcf.gz # 解压
cat clinvar_20220430.vcf.gz | grep "##" > head.tsv # 保存注释区
cat clinvar_20220430.vcf_ | grep -v "##" > clinvar.tsv # 删除注释

# 利用bcftools筛选下载的数据




wc -l clinvar.col4a5.tsv 
# 1252
# X       108580280       1678410 A       G
# X       108598691       1678532 T       A
# X       108694845       1678593 T       C
```

### 数据集合并
+ GnomAD处理
```bash
cat gnomAD.csv | tr "," "\t" > gnomeAD.tsv

```

+ Clinvar处理
```bash
cat clinvar.col4a5.tsv | perl -e' while(<>){
    if (/\s(\d*)\s\d*\s([A-Z])\s([A-Z])\s/) {
        print "X:$1:$2:$3"
    }
    
    if (/(CLNSIG=.*?);/) {
        print "\t$1\n"
    }
}' > clinvar.filter.tsv

# 添加致病与否的信息（致病T，不致病F）
cat clinvar.filter.tsv | perl -e' while(<>){
    chomp($_);
    if (/Pathogenic/) {
        print "$_\tT\n";
    }else {
        print "$_\tF\n";
    }
}' > tem&&
    mv tem clinvar.filter.tsv

# 统计
tsv-summarize --count -g 2 clinvar.filter.tsv
CLNSIG=Pathogenic       252
CLNSIG=Likely_benign    390
CLNSIG=Uncertain_significance   141
CLNSIG=Benign   123
CLNSIG=Conflicting_interpretations_of_pathogenicity     31
CLNSIG=Likely_pathogenic        260
CLNSIG=Benign/Likely_benign     25
CLNSIG=Pathogenic/Likely_pathogenic     28

tsv-summarize --count -g 3 clinvar.filter.tsv
T       280
F       970
```

+ 合并
```bash
cat clinvar.tsv | grep -f gnomAD.tsv | wc -l
303
cat clinvar.tsv | grep -f gnomAD.tsv > both.tsv
cat both.tsv | cut -f 1 > both.marker.tsv

cat clinvar.tsv | cut -f 1 > clinvar.marker.tsv
cat gnomAD.tsv | cut -f 1 > gnomAD.marker.tsv
cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv | wc -l
304

cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv | grep -v -f both.marker.tsv
X:108559155:T:C

for i in clinvar.tsv gnomAD.tsv;do
    cat $i | grep -v -f both.tsv | grep -v "X:108559155:T:C" >> merge.tsv
done
cat both.tsv >> merge.tsv

# 统计
tsv-summarize --count -g 2 merge.tsv
T       279
F       2137
```

### VEP的使用
```bash
# 安装
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl

# 注释

```





