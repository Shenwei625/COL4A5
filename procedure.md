# COL4A5
### 数据下载

```bash

```




### 数据集合并
+ GnomAD处理
```bash
### 表头，利用tsv-select提取列，利用perl单行生成唯一标识符，将awk换成sed
cat gnomAD.csv | tr "," "\t" | cut -f 1,2,4,5,14 > gnomAD.filter.tsv # 提取有用的列
cat gnomAD.csv | tr "," ":" | cut -d ":" -f 1,2,4,5 > gnomAD.marker.tsv # 整理出唯一标识符
sed -i '1d' gnomAD.marker.tsv # 删除表头
cat gnomAD.marker.tsv | wc -l # 1471

# 标注致病信息(Pathogenic和Pathogenic/Likely pathogenic)
cat gnomAD.filter.tsv | grep "Pathogenic" | cut -f 1,2,3,4 | tr "\t" ":" > gnomAD.T.tsv # 五个致病点的唯一标识符
cat gnomAD.marker.tsv | grep -v -f gnomAD.T.tsv > gnomAD.F.tsv
cat gnomAD.F.tsv | wc -l # 1466
awk '{print $0,"\tT"}' gnomAD.T.tsv  > tem&& # 标注致病信息，致病为T，不致病为F
  mv tem gnomAD.T.tsv
awk '{print $0,"\tF"}' gnomAD.F.tsv  > tem&&
  mv tem gnomAD.F.tsv  
cat gnomAD.F.tsv >> gnomAD.T.tsv # 合并
mv gnomAD.T.tsv gnomAD.tsv # 重命名
rm gnomAD.F.tsv gnomAD.filter.tsv gnomAD.marker.tsv  # 删除多余的过程文件
```

+ Clinvar处理
```bash
cat clinvar_result.txt | cut -f 11,15,5 > clinvar.filter.tsv    # 选择对我们有用的列
# 构建唯一标识符
cat clinvar.filter.tsv | cut -f 2 > clinvar.location.tsv # 找到位置
cat clinvar_result.txt | cut -f 10 > clinvar.x.tsv # 找到所在的染色体
cat clinvar.filter.tsv | cut -f 3 | cut -d ":" -f  3,4 > clinvar.allele.tsv # 找到对应snv
paste -d ":" clinvar.x.tsv clinvar.location.tsv clinvar.allele.tsv > clinvar.marker.tsv # 以冒号为连接符连接
sed -i '1d' clinvar.marker.tsv # 删去表头
rm clinvar.x.tsv clinvar.allele.tsv clinvar.location.tsv # 删去过程文件

# 标注是否致病（Clinical significance列含有Pathogenic即为致病）
### 利用tsv-summarize统计一下不同致病信息的数量（gnomAD和clinvar分开统计）
cat clinvar.filter.tsv | grep -i "pathogenic" | grep -v "Conflicting interpretations of pathogenicity" | cut -f 2,3 > clinvar.T.tsv   # 找到致病位点的相关信息
cat clinvar.T.tsv | wc -l #541
# 构建致病位点的唯一标识符
cat clinvar.T.tsv | cut -f 1 > clinvar.T.location.tsv
cat clinvar.T.tsv | cut -f 2 | cut -d ":" -f 3,4 > clinvar.T.allele.tsv
paste -d ":" clinvar.T.location.tsv clinvar.T.allele.tsv > clinvar.T.marker.tsv
grep -f clinvar.T.marker.tsv clinvar.marker.tsv > tem&&
  mv tem clinvar.T.tsv
cat clinvar.marker.tsv | grep -v -f clinvar.T.marker.tsv > clinvar.F.tsv
rm clinvar.T.location.tsv clinvar.T.allele.tsv clinvar.T.marker.tsv
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

+ 合并两个数据集
```bash
cat gnomAD.tsv | grep -f clinvar.tsv | wc -l #300
cat gnomAD.tsv | grep -f clinvar.tsv > both.tsv # 两个数据集的交集

# 找到两个数据集中突变位置相同，但是致病性注释不一样的点（舍去）
cat clinvar.tsv | cut -f 1 > clinvar.marker.tsv
cat gnomAD.tsv | cut -f 1 > gnomAD.marker.tsv
cat both.tsv | cut -f 1 > both.marker.tsv
cat clinvar.marker.tsv | grep -f gnomAD.marker.tsv > samemarker.tsv
cat samemarker.tsv | grep -v -f both.marker.tsv > discard.tsv # 两个数据集中注释不同的点，需要删除

# 取并集
cat clinvar.tsv | grep -v -f both.tsv | grep -v -f discard.tsv > clinvar.rest.tsv 
cat gnomAD.tsv | grep -v -f both.tsv | grep -v -f discard.tsv > gnomAD.rest.tsv
for i in clinvar.rest.tsv gnomAD.rest.tsv both.tsv;do
  cat $i >> terminal.tsv
done
cat terminal.tsv | wc -l #2411

# 删除多余的过程文件
rm both.marker.tsv clinvar.marker.tsv clinvar.rest.tsv gnomAD.marker.tsv samemarker.tsv gnomAD.rest.tsv
```

### VEP注释
+ 更改格式
```bash
cat terminal.tsv | cut -f 1 | cut -d ":" -f 1 > terminal.chromosome.tsv # 提取染色体信息
cat terminal.tsv | cut -f 1 | cut -d ":" -f 2 > terminal.location.tsv # 提取位置信息
cat terminal.tsv | cut -f 1 | cut -d ":" -f 3 > terminal.ref.tsv # 提取等位基因信息
cat terminal.tsv | cut -f 1 | cut -d ":" -f 4 > terminal.alt.tsv
# 合并
paste -d "\t" terminal.chromosome.tsv terminal.location.tsv terminal.location.tsv terminal.ref.tsv > format.tsv 
paste -d "/" format.tsv terminal.alt.tsv > tem&&
  > mv tem format.tsv
rm terminal.chromosome.tsv terminal.location.tsv terminal.ref.tsv terminal.alt.tsv #删除过程文件  
```

```bash
### 命令行使用vep
```





### 氨基酸保守位点(改成独立的脚本)
+ muscle软件的安装
```bash
brew install muscle
muscle --help
# 用法：muscle -in seqs.fa -out seqs.afa
```

+ 数据前处理
```bash
muscle -in lingzhang.fa -out lingzhang.afa
# faops size lingzhang.fa
# faops size lingzhang.afa    1721

cat lingzhang.afa | grep ">" | cut -d ">" -f 2 > lingzhang.lst # 保存序列名称，方便后续查找
faops some -l 1721 lingzhang.afa lingzhang.lst lingzhang.align.tsv
# 删去序列名称，方便后续统计
cat lingzhang.align.tsv | grep -v ">" > tem&&
  mv tem lingzhang.tsv
head lingzhang.align.tsv # 查看比对后序列

# 在每个字符间用tab分隔
sed 's/./&\t/g' lingzhang.tsv > lingzhang.tab.tsv
cat lingzhang.tab.tsv | wc -l #20

#转置
for i in `seq $(head -n 1 lingzhang.tab.tsv | awk '{print NF}')`; do
  cut -f $i lingzhang.tab.tsv | tr "\n" "\t"| sed '$ s/$/\n/' >> lingzhang_T.tsv 
done
cat lingzhang_T.tsv | wc -l # 1721
(echo -e "1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20" && cat lingzhang_T.tsv) > tem&&
  mv tem lingzhang_T.tsv # 添加表头
rm lingzhang.afa lingzhang.tsv # 删除多余文件
```

+ 统计
```bash
# 统计每一行除了-之外的字符数量
echo "total" >> number.tsv
for i in `seq $(head -n 1 lingzhang.tab.tsv | awk '{print NF}')`; do
  NUMBER=$(cut -f $i lingzhang.tab.tsv | grep -v "-" | wc -l)
  echo -e "$NUMBER"  >> number.tsv
done 

paste lingzhang_T.tsv number.tsv > tem&&
  mv tem lingzhang_T.tsv
rm number.tsv

# 找到人类col4a5每个位点氨基酸类型
cat lingzhang_T.tsv | cut -f 12 > human.tsv
sed -i '1d' human.tsv
(echo "human" && cat human.tsv) > tem&&
  mv tem human.tsv
paste lingzhang_T.tsv human.tsv > tem&&
  mv tem lingzhang_T.tsv
rm human.tsv

# 统计每个位置人类氨基酸出现的次数
for i in `seq $(head -n 1 lingzhang.tab.tsv | awk '{print NF}')`; do
  cat lingzhang.tab.tsv | cut -f $i > $i.tsv
  HUMAN=$(sed -n '12p' $i.tsv)
  STATISTIC=$(cat $i.tsv | grep $HUMAN | wc -l)
  echo -e "$STATISTIC" >> statistic.tsv
  rm $i.tsv
done

(echo "statistic" && cat statistic.tsv) >tem&&
  mv tem statistic.tsv #添加表头
paste lingzhang_T.tsv statistic.tsv > tem&&
  mv tem lingzhang_T.tsv #添加新的列

#过滤人类中的-
cat lingzhang_T.tsv | tsv--filter  --str-ne 12:- | wc -l # 1692
faops size lingzhang.fa # 人类col4a5蛋白序列长度为1691
cat lingzhang_T.tsv | tsv-filter  --str-ne 12:- > lingzhang.human.tsv # 因为要关注人类COL4A5氨基酸保守位点，所以要删去人类中-的行
```

+ 计算人类氨基酸在灵长目中的频率,保留小数点后三位
```bash
sed -n '1p' lingzhang.human.tsv > head.tsv # 保存表头
sed -i '1d' lingzhang.human.tsv # 删除统计表中的表头

cat lingzhang.human.tsv | cut -f 22 > total.tsv
cat lingzhang.human.tsv | cut -f 24 > statistic.tsv
for i in `seq $(cat lingzhang.human.tsv | wc -l)`; do
  TOTAL=$(eval sed -n '${i}p' total.tsv)
  STATISTIC=$(eval sed -n '${i}p' statistic.tsv)
  FREQUENCY=$(printf "%.3f" `echo "scale=3;$STATISTIC/$TOTAL" | bc`)
  echo -e "$FREQUENCY" >> frequency.tsv
done # 计算

(echo "frequency" && cat frequency.tsv) >tem&&
  mv tem frequency.tsv # 添加表头
(cat head.tsv && cat lingzhang.human.tsv) >tem&&
  mv tem lingzhang.human.tsv
paste lingzhang.human.tsv frequency.tsv >tem&&
  mv tem lingzhang.human.tsv #合并
rm head.tsv frequency.tsv total.tsv statistic.tsv
```

+ 筛选出保守位点
```bash
# 标记出位点
echo "location" > location.tsv
seq $(($(cat lingzhang.human.tsv | wc -l) - 1 )) >> location.tsv
paste lingzhang.human.tsv location.tsv > tem&&
  mv tem lingzhang.human.tsv
rm location.tsv

# 筛选 
cat lingzhang.human.tsv | tsv-filter -H --ge 25:0.950 | cut -f 26 > conserve.tsv
sed -i '1d' conserve.tsv
```













