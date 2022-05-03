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
+ 删选出col4a5中的snv位点
```bash

```



### 数据集合并
+ GnomAD处理
```bash

```

+ Clinvar处理
```bash

```



### VEP的使用






