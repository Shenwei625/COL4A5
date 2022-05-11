#! usr/bin/bash
file1=$1
file2=$2

# 构建结果文件
echo -e "参数\ncDNA_position\nCDS_position\nProtein_position\nSIFT\nPolyPhen" >> result.tsv

for F in $1 $2; do 
    file_name=$(echo $F | cut -d "." -f 1)
    # 对输入的文件中的SIFT和PolyPhen两列进行预处理（需要改进）
    tsv-select -H --fields SIFT,PolyPhen,cDNA_position,CDS_position,Protein_position $F > $file_name.filter.tsv
    HEAD=$(head -n 1 $file_name.filter.tsv)
    sed -i '1d' $file_name.filter.tsv
    sed -i 's/_//g' $file_name.filter.tsv # 删除下划线
    sed -i 's/[a-z]\+//g' $file_name.filter.tsv # 删除括号前面的字母
    # 删除括号
    sed -i 's/(//g' $file_name.filter.tsv
    sed -i 's/)//g' $file_name.filter.tsv
    (echo $HEAD | tr " " "\t" && cat $file_name.filter.tsv) > tem&&
        mv tem $file_name.filter.tsv


    # 计算结果
    for i in cDNA_position CDS_position Protein_position SIFT PolyPhen;do
        tsv-filter -H --str-ne $i:- $file_name.filter.tsv | 
        tsv-select -H --fields $i |
        keep-header -- sort -n > $file_name.$i.tsv # 取出对应的列并排序

        # 获取中位数
        LINE=$(sed '1d' $file_name.$i.tsv | wc -l)
        if [ $(($LINE%2)) == 0 ] ; then
            MEDIAN_LINE=$(( $LINE / 2 ))
            MEDIAN_LINE_AFTER=$(( $MEDIAN_LINE + 1 ))
            export AFTER=$(sed '1d' $file_name.$i.tsv | sort -n | eval sed -n '${MEDIAN_LINE_AFTER}p')
            export FRONT=$(sed '1d' $file_name.$i.tsv | sort -n | eval sed -n '${MEDIAN_LINE}p')
            MEDIAN=$(perl -e '
                my $F = $ENV{"FRONT"};
                my $A = $ENV{"AFTER"};
                my $M = ( $F + $A ) / 2;
                print "$M\n"
                ')
#            echo -e "Number of SNV mutation sites in $file_name.$i is 偶数"
#            echo -e "The median of $file_name.$i is $MEDIAN\n"
            echo $MEDIAN >> $file_name.median.tsv 
        else
            MEDIAN_LINE=$(( $(( $LINE + 1 )) / 2 ))
            MEDIAN=$(sed '1d' $file_name.$i.tsv | sort -n | eval sed -n '${MEDIAN_LINE}p')
#            echo -e "Number of SNV mutation sites in $file_name.$i is 奇数"
#            echo -e "The median of $file_name.$i is $MEDIAN\n"
            echo $MEDIAN >> $file_name.median.tsv     
        fi

        # 获取最大值和最小值
#        echo -e "The number of SNV mutation sites in $file_name.$i is $LINE\n"
        MAX=$(sed '1d' $file_name.$i.tsv | sort -n | eval sed -n '${LINE}p')
#        echo -e "The max of $file_name.$i is $MAX\n"
        echo $MAX >> $file_name.max.tsv
        
        MIN=$(sed '1d' $file_name.$i.tsv | sort -n | eval sed -n '1p')
#        echo -e "The min of $file_name.$i is $MIN\n"
        echo $MIN >> $file_name.min.tsv    
        done        
done

# 合并所有的结果信息
for M in median max min; do
    for F in $1 $2; do
        file_name=$(echo $F | cut -d "." -f 1)
        (echo -e "$file_name.$M" && cat $file_name.$M.tsv) > tem&&
            mv tem $file_name.$M.tsv
        paste -d "\t" result.tsv $file_name.$M.tsv >> tem&&
            mv tem result.tsv  
        rm $file_name.$M.tsv
    done
done    

#美化 
tsv-pretty -f -u result.tsv > tem&&
    mv tem result.tsv  

echo -e "==> Done, the results are in the result.tsv"