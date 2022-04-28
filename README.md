# COL4A5
### 数据下载

### 数据集合并
+ GnomAD处理
```bash

```

+ Clinvar处理
```bash
cat clinvar_result.txt | cut -f 11,15,5 > clinvar.filter.tsv    #选择对我们有用的列
# 构建唯一标识符
cat clinvar.filter.tsv | cut -f 2 > clinvar.location.tsv #找到位置
cat clinvar_result.txt | cut -f 10 > clinvar.x.tsv #找到所在的染色体
cat clinvar.filter.tsv | cut -f 3 | cut -d ":" -f  3,4 > clinvar.allele.tsv #找到对应snv
paste -d ":" clinvar.x.tsv clinvar.location.tsv clinvar.allele.tsv > clinvar.marker.tsv #以冒号为连接符连接
sed -i '1d' clinvar.marker.tsv #删去表头
rm clinvar.x.tsv clinvar.allele.tsv clinvar.location.tsv #删去过程文件

# 标注是否致病（Clinical significance列含有Pathogenic即为致病）
cat clinvar.filter.tsv | grep -i "pathogenic" | grep -v "Conflicting interpretations of pathogenicity" | cut -f 2,3 > clinvar.T.tsv | wc -l  #找到致病位点的相关信息
cat clinvar.T.tsv | wc -l #541
# 构建致病位点的唯一标识符
cat clinvar.T.tsv | cut -f 1 > clinvar.T.location.tsv
cat clinvar.T.tsv | cut -f 2 | cut -d ":" -f 3,4 > clinvar.T.allele.tsv
paste -d ":" clinvar.T.location.tsv clinvar.T.allele.tsv > clinvar.T.marker.tsv
rm clinvar.T.tsv clinvar.T.location.tsv clinvar.T.allele.tsv
# 标注
awk '{print $0,"\tT"}' clinvar.T.tsv  > tem&&
  mv tem clinvar.T.tsv
awk '{print $0,"\tF"}' clinvar.F.tsv  > tem&& 
  mv tem clinvar.F.tsv
cat clinvar.F.tsv >> clinvar.T.tsv
cat clinvar.T.tsv | wc -l #1248
mv clinvar.T.tsv clinvar.tsv
rm clinvar.F.tsv clinvar.T.marker.tsv clinvar.filter.tsv clinvar.marker.tsv
```
