#! usr/bin/bash

file1=$1

# 对vep文件进行预处理
sed -i 's/#//g' $file1
# Consequence列中存在逗号,进行更换
sed -i 's/,/;/g' $file1
# SIFT和Codons列不是字符串形式,进行更改,保留括号中的数字
tsv-select -H --fields Uploaded_variation,SIFT,PolyPhen $file1 > merge.tsv
HEAD=$(head -n 1 merge.tsv)
sed -i '1d' merge.tsv
sed -i 's/_//g' merge.tsv
sed -i 's/[a-z]*(//g' merge.tsv
sed -i 's/)//g' merge.tsv
(echo $HEAD | tr " "  "\t" && cat merge.tsv) > tem&&
    mv tem merge.tsv
tsv-select -H -e SIFT,PolyPhen $file1 > tem&&
    mv tem $file1
tsv-join --filter-file merge.tsv --H --key-fields 1 --append-fields SIFT,PolyPhen $file1 > tem&&
    mv tem $file1
rm merge.tsv

for COL in MAF renke_conserve xiabixiamu_conserve jianbiyamu_conserve lingzhang_conserve lingzhangzongmu_conserve beifangshoulei_conserve Consequence IMPACT EXON INTRON cDNA_position CDS_position topmed_AF Amino_acids Codons SIFT PolyPhen Protein_position;do
echo -e "==> 正在对$COL数据进行logistic分析"
# 提取出我们需要分析的数据
tsv-filter -H --str-ne $COL:- $file1 |
    tsv-select -H --fields Pathogenicity,$COL > $COL.logi.lst
sed -i 's/^T/1/g' $COL.logi.lst
sed -i 's/^F/0/g' $COL.logi.lst

# logistic回归分析
Rscript -e '
    FILE1 <- list.files(path = ".", pattern = "*.lst")
    DATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
    C <- colnames(DATA1[2])
    mod <- glm(Pathogenicity~get(C),data = DATA1,family=binomial)
    #取P值
    p<-summary(mod)$coefficients[,4]
    #wald值
    wald<-summary(mod)$coefficients[,3]^2
    #B值
    valueB<-coef(mod)
    #OR值
    valueOR<-exp(coef(mod))
    #OR值得95%CI
    confitOR<-exp(confint(mod))
    #构建表格
    FRAME <- data.frame(
        B=round(valueB,3),
        Wald=round(wald,3),
        OR_with_CI=paste(round(valueOR,3),"(",
                round(confitOR[,1],3),"~",round(confitOR[,2],3),")",sep=""),
        P=format.pval(p,digits = 3,eps=0.001)
    )

    NAME <- paste(sep = "", C, ".csv")
    write.csv(FRAME, file = NAME)
'
rm $COL.logi.lst
cat $COL.csv | tr "," "\t" | sed 's/"//g' > $COL.logi.tsv
sed -i 's/get(C)/'$COL'/g' $COL.logi.tsv
cat $COL.logi.tsv | grep -v "Intercept" > tem&&
    mv tem $COL.logi.tsv
# 添加显著标记(0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1)
sed '1d' $COL.logi.tsv | perl -na -F"\t" -e'
    if ( 0.1 < @F[4] <= 1 ) {
        print "\n";
    }elsif ( 0.05 < @F[4] <= 0.1 ) {
        print ".\n";
    }elsif ( 0.01 < @F[4] <= 0.05 ) {
        print "*\n";
    }elsif ( 0.001 < @F[4] <= 0.01 ) {
        print "**\n";
    }else {
        print "***\n";
    }
' > $COL.merge.lst
(echo -e "\n" && cat $COL.merge.lst) > tem&&
    mv tem $COL.merge.lst
cat $COL.logi.tsv | mlr --itsv --omd cat > tem&&
    mv tem $COL.logi.tsv
sed -i 's/|$//g' $COL.logi.tsv
paste -d "" $COL.logi.tsv $COL.merge.lst > tem&&
    mv tem $COL.logi.tsv
cat $COL.logi.tsv | perl -ne '
    chomp($_);
    print "$_|\n";
' > tem&&
    mv tem $COL.logi.tsv
rm $COL.merge.lst

echo -e "\n显著因素:\n" >> $COL.logi.tsv
cat $COL.csv | tr "," "\t" | 
    sed 's/"//g' | sed 's/get(C)/'$COL'/g' |
    grep -v "<" | 
    tsv-filter -H --le 5:0.1 >> $COL.sig.tsv
cat $COL.csv | tr "," "\t" | 
    sed 's/"//g' | sed 's/get(C)/'$COL'/g' |
    grep "<" >> $COL.sig.tsv
cat $COL.sig.tsv | keep-header -- sort -n -k5,5 > tem&&
    mv tem $COL.sig.tsv
cat $COL.sig.tsv | grep -v "Intercept" > tem&&
    mv tem $COL.sig.tsv
sed '1d' $COL.sig.tsv | perl -na -F"\t" -e'
    if ( 0.1 < @F[4] <= 1 ) {
        print "\n";
    }elsif ( 0.05 < @F[4] <= 0.1 ) {
        print ".\n";
    }elsif ( 0.01 < @F[4] <= 0.05 ) {
        print "*\n";
    }elsif ( 0.001 < @F[4] <= 0.01 ) {
        print "**\n";
    }else {
        print "***\n";
    }
' > $COL.merge.lst
(echo -e "\n" && cat $COL.merge.lst) > tem&&
    mv tem $COL.merge.lst
cat $COL.sig.tsv | mlr --itsv --omd cat > tem&&
    mv tem $COL.sig.tsv
sed -i 's/|$//g' $COL.sig.tsv
paste -d "" $COL.sig.tsv $COL.merge.lst > tem&&
    mv tem $COL.sig.tsv
cat $COL.sig.tsv | perl -ne '
    chomp($_);
    print "$_|\n";
' > tem&&
    mv tem $COL.sig.tsv
rm $COL.merge.lst
cat $COL.sig.tsv >> $COL.logi.tsv
rm $COL.csv $COL.sig.tsv
done


# 构建MAF和topmed_AF的最大值、最小值和中位数
echo -e "开始计算MAF和topmed_AF的最大值、中位数和最小值\n"
for i in topmed_AF MAF;do
    # 提取出我们需要分析的数据
    tsv-filter -H --str-ne $i:- $file1 |
        tsv-select -H --fields Pathogenicity,$i | 
        sort -n -k2,2 > $i.logi.lst
    sed -i 's/^T/1/g' $i.logi.lst
    sed -i 's/^F/0/g' $i.logi.lst

    # 计算最大值和最小值
    MIN=$(cut -f 2 $i.logi.lst | sed -n '2p')
    MAX_LINE=$(cat $i.logi.lst | wc -l) 
    MAX=$(eval sed -n '${MAX_LINE}p' $i.logi.lst | cut -f 2)

    # 计算中位数
    NUM_LINE=$(( $MAX_LINE - 1 ))
    if [ $(($NUM_LINE%2)) == 0 ] ; then
        MEDIAN_LINE=$(( $NUM_LINE / 2 ))
        MEDIAN_LINE_AFTER=$(( $MEDIAN_LINE + 1 ))
        export AFTER=$(sed '1d' $i.logi.lst | eval sed -n '${MEDIAN_LINE_AFTER}p' | cut -f 2)
        export FRONT=$(sed '1d' $i.logi.lst | eval sed -n '${MEDIAN_LINE}p' | cut -f 2)
        MEDIAN_RAW=$(perl -e '
                my $F = $ENV{"FRONT"};
                my $A = $ENV{"AFTER"};
                my $M = ( $F + $A ) / 2;
                print "$M\n"
                ')
        MEDIAN=$(printf "%*.*f\n" 16 15 $MEDIAN_RAW)
    else
        MEDIAN_LINE=$(( $(( $NUM_LINE + 1 )) / 2 ))
        MEDIAN_RAW=$(sed '1d' $i.logi.lst | eval sed -n '${MEDIAN_LINE}p' | cut -f 2)
        MEDIAN=$(printf "%*.*f\n" 16 15 $MEDIAN_RAW)
    fi
    (echo -e "$i的最大值是:$MAX\n\n$i的中位数是:$MEDIAN\n\n$i的最小值是:$MIN\n" && cat $i.logi.tsv) > tem&&
        mv tem $i.logi.tsv
        rm $i.logi.lst
done