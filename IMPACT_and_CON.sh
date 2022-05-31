#! usr/bin/bash
file1=$1
file2=$2

# 获得绘图文件
for F in $1 $2; do 
    file_name=$(echo $F | cut -d "." -f 1)
    for i in Consequence IMPACT; do
        tsv-filter -H --str-ne $i:- $F | 
        tsv-select -H --fields $i |
        tsv-summarize -H --count -g $i |
        keep-header -- sort -nr -k2,2 > $file_name.$i.tsv # 取出对应的列并排序
    done
done

# 绘图
Rscript -e '
    library(ggplot2)
    library(patchwork)
    JOB <- c("IMPACT","Consequence")

    for(J in JOB) {
        FILE1 <- paste(sep="","T.", J, ".tsv")
        PIEDATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
        TITLE1 <- paste(sep="","T.", J)
        FENZI <- PIEDATA1$count
        FENMU <- sum(PIEDATA1$count)
        myLabel1 = as.vector(PIEDATA1[,J])
        myLabel1 = paste(myLabel1, "(", FENZI, "/", FENMU, ",", round(PIEDATA1$count / sum(PIEDATA1$count) * 100, 2), "%" , ")")
        PLABEL1 = as.vector(PIEDATA1[,J])

        p1 <- ggplot(PIEDATA1,aes(x="", y=count, fill=get(J))) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar(theta = "y") + 
            labs(x = "", y = "", title = TITLE1) +
            theme(axis.ticks = element_blank(), plot.title=element_text(hjust=0.5)) +
            theme(axis.text.x = element_blank()) + 
            theme(panel.grid = element_blank()) + 
            theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
            theme(legend.title = element_blank(), legend.position = "bottom") +
            theme(legend.direction = "vertical") +
            scale_fill_discrete( breaks = PIEDATA1[,J], labels = myLabel1 ) +
            theme(legend.text=element_text(size=10)) +
            theme(text=element_text(size = 20, face = "bold")) +
            theme(legend.text=element_text(size=12))
    
    
        FILE2 <- paste(sep="","F.", J, ".tsv")
        PIEDATA2 <- read.table(FILE2,header=TRUE,sep = "\t")
        TITLE2 <- paste(sep="","F.", J)
        FENZI <- PIEDATA2$count
        FENMU <- sum(PIEDATA2$count)
        myLabel2 = as.vector(PIEDATA2[,J])
        myLabel2 = paste(myLabel2, "(", FENZI, "/", FENMU, ",", round(PIEDATA2$count / sum(PIEDATA2$count) * 100, 2), "%" , ")")
        PLABEL2 = as.vector(PIEDATA2[,J])

        p2 <- ggplot(PIEDATA2,aes(x="", y=count, fill=get(J))) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar(theta = "y") + 
            labs(x = "", y = "", title = TITLE2) +
            theme(axis.ticks = element_blank(), plot.title=element_text(hjust=0.5)) +
            theme(axis.text.x = element_blank()) + 
            theme(panel.grid = element_blank()) + 
            theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
            theme(legend.title = element_blank(), legend.position = "bottom") +
            theme(legend.direction = "vertical") +
            scale_fill_discrete( breaks = PIEDATA2[,J], labels = myLabel2 ) +
            theme(legend.text=element_text(size=10)) +
            theme(text=element_text(size = 20, face = "bold")) +
            theme(legend.text=element_text(size=12))
    
        p1 + p2
        NAME <- paste(sep="", J, ".png" )
        ggsave(file=NAME, width=17, height=12, dpi=300)
    }
'

