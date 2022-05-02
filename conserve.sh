### 氨基酸保守位点的确定
#! usr/bin/bash
file=$1 #从命令行中传递文件名称
file_name=$(echo $file | cut -d "." -f 1) #获取不带扩展名的文件名称

muscle -in $file -out $file_name.afa # 利用musucle比对文件


## 数据前处理
# 将比对后的序列一一对应，排布到tsv文件中
cat $file_name.afa | grep ">" | cut -d ">" -f 2 > $file_name.lst # 保存序列名称，方便后续查找
human_line=$(cat $file_name.lst | grep -n "Homo_sapiens" | cut -d ":" -f 1) # 保存人类序列所在的列，便于后续统计
align_length=$(faops size $file_name.afa | cut -f 2 | uniq) # 保存序列比对后长度
faops some -l $align_length $file_name.afa $file_name.lst $file_name.align.tsv # 将所有的的序列放在一列中，方便后续统计比较
cat $file_name.align.tsv | grep -v ">" > tem&&
  mv tem $file_name.tsv # 删去序列名称，方便后续统计

# 为了方便后续统计，将序列以列的方式排布，相邻两个氨基酸位点用tab分隔，并添加表头
sed 's/./&\t/g' $file_name.tsv > $file_name.tab.tsv # 在每个字符间用tab分隔
for i in `seq $align_length`; do
  cut -f $i $file_name.tab.tsv | tr "\n" "\t"| sed '$ s/$/\n/' >> $file_name.T.tsv 
done # 转置
sequence_number=$(cat $file_name.lst | wc -l) # 保存序列的数量
seq $sequence_number | tr "\n" "\t" | sed '$ s/$/\n/' > head.tsv # 新建表头文件
(cat head.tsv && cat $file_name.T.tsv) > tem&&
    mv tem $file_name.T.tsv # 在统计文件中添加表头


## 统计
# 统计每一行除了-之外的字符数量
for i in `seq $align_length`; do
  NUMBER=$(cut -f $i $file_name.tab.tsv | grep -v "-" | wc -l)
  echo -e "$NUMBER"  >> number.tsv
done  # 计算
(echo "total" && cat number.tsv) >tem&&
  mv tem number.tsv # 添加表头
paste $file_name.T.tsv number.tsv > tem&&
  mv tem $file_name.T.tsv # 将统计结果追加到文件末尾（新的一列）

# 统计每个位置人类氨基酸出现的次数
for i in `seq $align_length`; do
  cat $file_name.tab.tsv | cut -f $i > $i.tsv
  HUMAN=$(eval sed -n '${human_line}p' $i.tsv)
  STATISTIC=$(cat $i.tsv | grep "$HUMAN" | wc -l)
  echo -e "$STATISTIC" >> statistic.tsv
  rm $i.tsv
done # 计算
(echo "statistic" && cat statistic.tsv) >tem&&
  mv tem statistic.tsv # 添加表头
paste $file_name.T.tsv statistic.tsv > tem&&
  mv tem $file_name.T.tsv # 将统计内容追加到文件末尾（新的一列）

# 计算每个人类氨基酸的频率,保留小数点后三位
cat $file_name.T.tsv | tsv-filter --str-ne 12:- > $file_name.human.tsv # 因为要关注人类COL4A5氨基酸保守位点，所以要删去人类中-的行
sed -n '1p' $file_name.human.tsv > head.tsv # 保存表头
sed -i '1d' $file_name.human.tsv # 删除统计表中的表头
total_col=$(( $sequence_number + 2 )) # 保存total列数
statistic_col=$(( $sequence_number + 3 )) # 保存statistic列数
sed -i '1d' number.tsv 
sed -i '1d' statistic.tsv
for i in `seq $(cat $file_name.human.tsv | wc -l)`; do
  TOTAL=$(eval sed -n '${i}p' number.tsv)
  STATISTIC=$(eval sed -n '${i}p' statistic.tsv)
  FREQUENCY=$(printf "%.3f" `echo "scale=3;$STATISTIC/$TOTAL" | bc`)
  echo -e "$FREQUENCY" >> frequency.tsv
done # 计算
(echo "frequency" && cat frequency.tsv) >tem&&
  mv tem frequency.tsv # 添加表头
(cat head.tsv && cat $file_name.human.tsv) >tem&&
  mv tem $file_name.human.tsv
paste $file_name.human.tsv frequency.tsv >tem&&
  mv tem $file_name.human.tsv # 合并


## 筛选出氨基酸保守位点（频率大于0.950）
# 标记出位点
echo "location" > location.tsv
seq $(($(cat $file_name.human.tsv | wc -l) - 1 )) >> location.tsv
paste $file_name.human.tsv location.tsv > tem&&
  mv tem $file_name.human.tsv
rm location.tsv

# 筛选
frequency_col=$(( $sequence_number + 4 )) # 保存频率列
location_col=$(( $sequence_number + 5 )) # 保存位置列
cat $file_name.human.tsv | tsv-filter -H --ge $frequency_col:0.950 | cut -f $location_col > conserve.tsv
sed -i '1d' conserve.tsv

echo "==> Complete! Conservative sites are saved in conserve.tsv"





