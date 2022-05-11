#! usr/bin/bash
file1=$1
file2=$2

# 获得绘图文件
for F in $1 $2; do 
    file_name=$(echo $F | cut -d "." -f 1)
    for i in Amino_acids Codons; do
        tsv-filter -H --str-ne $i:- $F | 
        tsv-select -H --fields $i |
        tsv-summarize -H --count -g $i |
        keep-header -- sort -nr -k2,2 > $file_name.$i.tsv # 取出对应的列并排序
    done
done


# 作图
Rscript -e '
    library(treemapify)
    library(ggplot2)
    library(patchwork)

    JOB <- c("Amino_acids","Codons")
    for(J in JOB) {
        FILE1 <- paste(sep="","T.", J, ".tsv")
        TREEDATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
        TITLE1 <- paste(sep="","T.", J)
        p1 <- ggplot(TREEDATA1, aes(area = count, fill = get(J), label = paste(get(J), count,sep = "," ),subgroup = get(J))) +
        geom_treemap() +
        geom_treemap_text(colour = "white",place = "bottomleft",size = 10) +
        geom_treemap_subgroup_border(colour = "white", size = 2) +
        labs(title = TITLE1) +
        theme(legend.position = "none", plot.title=element_text(hjust=0.5))

        FILE2 <- paste(sep="","F.", J, ".tsv")
        TREEDATA2 <- read.table(FILE2,header=TRUE,sep = "\t")
        TITLE2 <- paste(sep="","F.", J)
        p2 <- ggplot(TREEDATA2, aes(area = count, fill = get(J), label = paste(get(J), count,sep = "," ),subgroup = get(J))) +
        geom_treemap() +
        geom_treemap_text(colour = "white",place = "bottomleft",size = 10) +
        geom_treemap_subgroup_border(colour = "white", size = 2) +
        labs(title = TITLE2) +
        theme(legend.position = "none", plot.title=element_text(hjust=0.5))
        p1 / p2
        NAME <- paste(sep="", J, ".pdf" )
        ggsave(file=NAME, width=20, height=12, dpi=300)
    }
'





