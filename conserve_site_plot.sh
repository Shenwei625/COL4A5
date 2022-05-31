#! usr/bin/bash
file1=$1

tsv-select -H --fields '*_conserve' $1 > raw_conserve.tsv
JOB=$(head -n 1 raw_conserve.tsv)
rm raw_conserve.tsv
for J in $JOB;do
    echo -e "==>Draw a pie chart of conserved amino acid sites of COL4A5 in $J\n"
    
    # 构建作图文件
    tsv-filter -H --str-ne $J:- $1 |
    tsv-select -H --fields $J |
    tsv-summarize -H --count -g $J > $J.lst

    # 作图
    Rscript -e '
    library(ggplot2)
    
    FILE1 <- list.files(path = ".", pattern = "*.lst")
    PIEDATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
    LABEL1 <- colnames(PIEDATA1[1])
    myLabel1 <- paste(round(PIEDATA1$count / sum(PIEDATA1$count) * 100, 2), "%")

    ggplot(PIEDATA1,aes(x="", y=count, fill=get(LABEL1))) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") + 
    labs(x = "", y = "", title = LABEL1) +
    theme(axis.ticks = element_blank(), plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x = element_blank()) + 
    theme(panel.grid = element_blank()) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
    theme(legend.title = element_blank(), legend.position = "bottom") +
    geom_text(aes(1, label = myLabel1),position = position_stack(vjust = 0.5),size = 4) 
    
    NAME <- paste(sep="", LABEL1, ".png" )
    ggsave(file=NAME, width=6, height=7.5, dpi=300)
    '
    rm $J.lst
done
rm Rplots.pdf





