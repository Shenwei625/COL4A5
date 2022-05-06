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

mv gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz sites.chrX.vcf.bgz #重命名
# 利用bcftools筛选下载的数据
bcftools index --threads 4 sites.chrX.vcf.bgz
bcftools filter sites.chrX.vcf.bgz --regions chrX:108439838-108697545 > gnomAD.x.vcf
bcftools view -v snps gnomAD.x.vcf > gnomAD.col4a5.vcf
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
mv clinvar_20220430.vcf.gz clinvar.vcf.gz
# 利用bcftools筛选下载的数据
bcftools index --threads 4 clinvar.vcf.gz
bcftools filter clinvar.vcf.gz --regions X:108439838-108697545 > clinvar.x.vcf
bcftools view -v snps clinvar.x.vcf > clinvar.col4a5.vcf

wc -l clinvar.col4a5.vcf 
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
cat clinvar.col4a5.vcf | grep -v "##" > tem&&
    mv tem clinvar.col4a5.vcf

cat clinvar.col4a5.vcf | perl -e' while(<>){   
    if (/CLNSIG=/) {
        if (/\s(\d*)\s\d*\s([A-Z])\s([A-Z])\s/) {
            print "X:$1:$2:$3";
        }
    
        if (/(CLNSIG=.*?);/) {
            print "\t$1\n";
        }
    }else {
        next;
    }    
}' > clinvar.filter.tsv

# 添加致病与否的信息（致病T，不致病F）
cat clinvar.filter.tsv | perl -e' while(<>){
    chomp($_);
    if (/Pathogenic|Likely_pathogenic/) {
        print "$_\tT\n";
    }else {
        print "$_\tF\n";
    }
}' >tem&&
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
T       540
F       710
```

+ 合并
```bash
mkdir merge
cd merge
cat ../clinvar.filter.tsv | cut -f 1.3 > clinvar.tsv



cat clinvar.tsv | grep -f gnomAD.tsv | wc -l
300
cat clinvar.tsv | grep -f gnomAD.tsv > both.tsv
cat both.tsv | cut -f 1 > both.marker.tsv

cat clinvar.tsv | cut -f 1 > clinvar.marker.tsv
cat gnomAD.tsv | cut -f 1 > gnomAD.marker.tsv
cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv | wc -l
304

cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv | grep -v -f both.marker.tsv
X:108559155:T:C
X:108601959:A:G
X:108624316:G:A
X:108677590:G:C
cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv | grep -v -f both.marker.tsv > discard.tsv

for i in clinvar.tsv gnomAD.tsv;do
    cat $i | grep -v -f both.tsv | grep -v -f discard.tsv >> merge.tsv
done
cat both.tsv >> merge.tsv

# 统计
tsv-summarize --count -g 2 merge.tsv
T       536
F       1877

# 转换格式
cat merge.tsv | cut -f 1 | perl -e ' while(<>){
    chomp($_);
    if (/X:(\d*):([A-Z]):([A-Z])/) {
        print "X\t$1\t$1\t$2/$3\n"
    }
}'   > format.tsv
```

### VEP的使用
```bash
# 安装
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl

# 注释



# 在VEP结果后面添加致病性
# 将vep文件的第一列改为唯一标识符
HEAD=$(head -n 1 vep.tsv) 
sed -i '1d' vep.tsv
sed -i 's/_/:/1' vep.tsv
sed -i 's/_/:/1' vep.tsv
sed -i 's/\//:/1' vep.tsv
(echo $HEAD | tr " " "\t" && cat vep.tsv) > tem&&
    mv tem vep.tsv

(echo -e "#Uploaded_variation\tPathogenicity" && cat merge.tsv) > tem&&
    mv tem merge.tsv # 添加表头  
    
tsv-join --filter-file merge.tsv --H --key-fields 1 --append-fields Pathogenicity vep.tsv > tem&&
    mv tem vep.tsv
```





