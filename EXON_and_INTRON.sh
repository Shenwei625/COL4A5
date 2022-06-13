#! usr/bin/bash
file1=$1
file2=$2
file3=$3

echo -e "==> Construction of mutation site distribution map."
for F in $1 $2 $3; do 
    file_name=$(echo $F | cut -d "." -f 1)
    
    for i in EXON INTRON;do
        tsv-filter -H --str-ne $i:- $F | 
        tsv-select -H --fields $i |
        keep-header -- sort -n -k1,1 > $file_name.$i.tsv
        
        # 判断输入的文件是致病还是不致病
        if [ $file_name == T ] ; then
            export PCONTENT="Pathogenic"
        elif [ $file_name == F ] ; then
            export PCONTENT="Non-pathogenic"
        else
            export PCONTENT="Unknown"    
        fi
        
        # 构建作图文件
        cat $file_name.$i.tsv | perl -e' while (<>) {
            chomp($_);
            if (/^\d/) {
                print "$_\t" . $ENV{"PCONTENT"} . "\n"
            }else {
                print "$_\tgroup\n"
            }
        }' > tem&&
            mv tem $file_name.$i.tsv   
        sed -i 's/\/\S\S//g' $file_name.$i.tsv

    done     
done

# 作图(需要改进)
Rscript -e '
    library(ggplot2)
    JOB <- c("EXON","INTRON")
    Pathogenicity <- c("T.","F.","unknown.")
    for(J in JOB) {
        for(P in Pathogenicity) {
            FILE1 <- paste(sep="", P, J, ".tsv")
            HISTABLE1 <- read.table(FILE1,header=TRUE,sep = "\t")
            LABEL1 <- paste(sep="", P, J)
            ggplot() +
            geom_histogram(aes(x=get(J),fill=group,y=..density..), HISTABLE1, binwidth = 1,colour="black", position="identity") +
            scale_x_continuous(breaks = seq(55)) + xlab(LABEL1) + ylab("proportion") + theme(panel.grid = element_blank()) +
            geom_vline(aes(xintercept=median(HISTABLE1[,J], na.rm=T)),color="black", linetype="dashed", size=.5) +
            theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold")) 
            NAME <- paste(sep="", P, J, ".png" )
            ggsave(file=NAME, width=10, height=3, dpi=300)
        }
    }
'

echo -e "==> Construct the distribution density map of mutation sites."
# 计算内含子与外显子上的突变密度
# 构建内含子和外显子长度文件
sed '2~2d' length.tsv | cut -f 1,3 > EXON.LENGTH.tsv
sed '1~2d' length.tsv | cut -f 3 > INTRON.LENGTH.tsv
(echo -e "EXON\tlength" && cat EXON.LENGTH.tsv) > tem&& 
    mv tem EXON.LENGTH.tsv
LINE=$(cat INTRON.LENGTH.tsv | wc -l)
seq $LINE > INTRON.tem.tsv
paste -d "\t" INTRON.tem.tsv INTRON.LENGTH.tsv > tem&&
    mv tem INTRON.LENGTH.tsv
rm INTRON.tem.tsv
(echo -e "INTRON\tlength" && cat INTRON.LENGTH.tsv) > tem&& 
    mv tem INTRON.LENGTH.tsv

# 计算密度（每1kb发生的单核苷酸变异数量）
for P in T F unknown;do
    for i in EXON INTRON;do
        tsv-summarize -H --count -g 1 $P.$i.tsv > $P.$i.count.tsv
        tsv-join -H --filter-file $i.length.tsv --key-fields $i --append-fields length $P.$i.count.tsv > tem&&
            mv tem $P.$i.count.tsv
        cat $P.$i.count.tsv | perl -e 'while(<>) {
            chomp($_);
            if (/^(\d+)\t(\d+)\t(\d+)/) {
                my $result = ( $2 / $3 ) * 1000 ;
                printf "$_\t%.02f\n", $result;
            }else {
                print "$_\tSNP\n";
            }  
        }
        ' > tem&&
        mv tem $P.$i.count.tsv
        rm $P.$i.tsv
    done
done            
rm *.LENGTH.tsv

Rscript -e '
    library(ggplot2)
    JOB <- c("EXON","INTRON")
    Pathogenicity <- c("T.","F.","unknown.")
    for(J in JOB) {
        for(P in Pathogenicity) {
            FILE1 <- paste(sep="", P, J, ".count.tsv")
            COLDATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
            LABEL1 <- paste(sep="", P, J)
            ggplot(COLDATA1,aes(x=COLDATA1[,J], y=SNP)) +
            geom_bar(stat = "identity", width = 1,colour="black",fill="34170212") +
            scale_x_continuous(breaks = seq(55)) + ylab("SNP( per kilo base)") + xlab(LABEL1) + theme(panel.grid = element_blank()) +
            theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold")) 
            NAME <- paste(sep="", P, J, ".SNP.png" )
            ggsave(file=NAME, width=10, height=3, dpi=300)
        }
    }
'
rm *.count.tsv




