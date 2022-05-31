#! usr/bin/bash
file1=$1
file2=$2

# 构建作图文件
for F in $1 $2; do 
    file_name=$(echo $F | cut -d "." -f 1)
    tsv-filter -H --str-ne MAF:- $F | 
    tsv-select -H --fields MAF |
    keep-header -- sort -n -k2,2 > $file_name.MAF.tsv
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

    cat $file_name.MAF.tsv | perl -e' while (<>) {
        chomp($_);
        if (/^\d/) {
            print "$_\t" . $ENV{"PCONTENT"} . "\n"
        }else {
            print "$_\tgroup\n"
            }
        }' > tem&&
        mv tem $file_name.MAF.tsv 
done

# 作图
Rscript -e '
    library(ggplot2)
    Pathogenicity <- c("T.","F.")
    for(P in Pathogenicity) {
            FILE1 <- paste(sep="", P, "MAF", ".tsv")
            HISTABLE1 <- read.table(FILE1,header=TRUE,sep = "\t")
            LABEL1 <- paste(sep="", P, "MAF")
            ggplot() +
            geom_histogram(aes(x=MAF,fill=group,y=..density..), HISTABLE1,colour="black", position="identity") +
            xlab(LABEL1) + ylab("proportion") + 
            scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
            theme(panel.grid = element_blank()) +
            geom_vline(aes(xintercept=median(HISTABLE1[,1], na.rm=T)),color="black", linetype="dashed", size=.5)
            NAME <- paste(sep="", P, "MAF", ".png" )
            ggsave(file=NAME, width=8, height=3, dpi=300)
        }
'