# COL4A5
### 数据下载
+ GnomAD
```bash
# 利用gsutil下载数据
sudo apt upgrade
sudo apt install gsutil

# 检索gnomAD数据库中的数据
gsutil ls gs://gcp-public-data--gnomad/release/
gsutil ls gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/
gsutil du gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz # 查看文件大小
gsutil cp gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz . # 下载
# 文件太大，是否可以下载特定转录本、基因的突变信息？？
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
cat clinvar_20220430.vcf_ | grep -v -f head.tsv > clinvar.tsv # 删除注释

# X染色体
tsv-filter -H --str-eq 1:X clinvar.tsv > clinvar.x.tsv
# 位置：108439838-108697545
tsv-filter -H --ge 2:108439838 clinvar.x.tsv | 
  tsv-filter -H --le 2:108697545 > clinvar.x.loc.tsv
# 单碱基替换
tsv-filter -H --char-len-eq 4:1 clinvar.x.loc.tsv | 
  tsv-filter -H --char-len-eq 5:1 > clinvar.col4a5.tsv 

wc -l clinvar.col4a5.tsv 
# 1252
# X       108580280       1678410 A       G
# X       108598691       1678532 T       A
# X       108694845       1678593 T       C
```

### 数据集合并
+ GnomAD处理
```bash

```

+ Clinvar处理
```bash
# 构建唯一标识符，vcf文件为例
cat clinvar.vcf | grep -v "##" > tem&&
  mv tem clinvar.vcf
LINE=$(wc -l clinvar.vcf)
head -n 1 clinvar.vcf
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
for i in $(seq $LINE);do
  CHROM=$(sed -n "${i}p" clinvar.vcf | cut -f 1)
  POS=$(sed -n "${i}p" clinvar.vcf | cut -f 2)
  REF=$(sed -n "${i}p" clinvar.vcf | cut -f 4)
  ALT=$(sed -n "${i}p" clinvar.vcf | cut -f 5)
  echo -e "$CHROM:$POS:$REF:$ALT" >> marker.tsv
done  
sed -i '1d' clinical.marker.tsv
```



### VEP的使用
```bash
# 安装
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl

# 注释



```





