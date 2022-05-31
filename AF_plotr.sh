#! usr/bin/bash
file1=$1
file2=$2

# 获得绘图文件
for F in $1 $2; do 
    file_name=$(echo $F | cut -d "." -f 1)
    tsv-filter -H --str-ne topmed_AF:- $F | 
    tsv-select -H --fields topmed_AF |
    keep-header -- sort -n -k2,2 > $file_name.topmed_AF.tsv
done

# 构建作图文件
for F in $1 $2; do
    file_name=$(echo $F | cut -d "." -f 1)
    # 判断输入的文件是致病还是不致病
    if [ $file_name == T ] ; then
        export PCONTENT="Pathogenic"
    else
        export PCONTENT="Non-pathogenic"    
    fi
    
    cat $file_name.topmed_AF.tsv | perl -e' while (<>) {
        chomp($_);
        if (/^\d/) {
            print "$_\t" . $ENV{"PCONTENT"} . "\n"
        }else {
            print "$_\tgroup\n"
            }
        }' > tem&&
        mv tem $file_name.topmed_AF.tsv 

#    NAME=$file_name.topmed_AF
#    FILE=$NAME.pdf
#    plotr hist $file_name.topmed_AF.tsv \
#        -g 2 \
#        -c 1 \
#        -p \
#        --xl "$NAME" \
#        --yl "propertion" \
#        --device pdf -o $FILE 
done


# 作图
Rscript -e '
    library(ggplot2)
    Pathogenicity <- c("T.","F.")
    for(P in Pathogenicity) {
            FILE1 <- paste(sep="", P, "topmed_AF", ".tsv")
            HISTABLE1 <- read.table(FILE1,header=TRUE,sep = "\t")
            LABEL1 <- paste(sep="", P, "topmed_AF")
            ggplot() +
            geom_histogram(aes(x=topmed_AF,fill=group,y=..density..), HISTABLE1,colour="black", position="identity") +
            xlab(LABEL1) + ylab("proportion") + 
            theme(panel.grid = element_blank()) +
            geom_vline(aes(xintercept=median(HISTABLE1[,1], na.rm=T)),color="black", linetype="dashed", size=.5) +
            scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
            theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold"))
            NAME <- paste(sep="", P, "topmed_AF", ".png" )
            ggsave(file=NAME, width=8, height=3, dpi=300)
        }
'







