#! usr/bin/bash
file1=$1
file2=$2
file3=$3

# 对SIFT PolyPhen进行处理（保留括号里面的数字）
for F in $1 $2 $3; do 
    sed -i 's/#//g' $F
    tsv-select -H --fields Uploaded_variation,SIFT,PolyPhen $F > merge.tsv
    HEAD=$(head -n 1 merge.tsv)
    sed -i '1d' merge.tsv
    sed -i 's/_//g' merge.tsv
    sed -i 's/[a-z]*(//g' merge.tsv
    sed -i 's/)//g' merge.tsv
    (echo $HEAD | tr " "  "\t" && cat merge.tsv) > tem&&
        mv tem merge.tsv
    tsv-select -H -e SIFT,PolyPhen $F > tem&&
        mv tem $F
    tsv-join --filter-file merge.tsv --H --key-fields 1 --append-fields SIFT,PolyPhen $F > tem&&
        mv tem $F
done
rm merge.tsv

# 每个参数做一个箱线图,致病与非致病作为对比放在一张图里
for i in cDNA_position CDS_position Protein_position SIFT PolyPhen;do 
    echo -e "Making plot for $i\n"
    for F in $1 $2 $3; do 
        file_name=$(echo $F | cut -d "." -f 1)
        tsv-filter -H --str-ne $i:- $F |
        tsv-select -H --fields Pathogenicity,$i > $file_name.$i.lst 
    done 
    sed '1d' F.$i.lst >> T.$i.lst
    sed '1d' unknown.$i.lst >> T.$i.lst
    rm F.$i.lst unknown.$i.lst
    mv T.$i.lst $i.lst

    #作图
    Rscript -e '
    library(ggplot2)
    FILE1 <- list.files(path = ".", pattern = "*.lst")
    DATA1 <- read.table(FILE1,header=TRUE,sep = "\t")
    C <- colnames(DATA1[2])
    ggplot(DATA1, aes(x=Pathogenicity, y=get(C), fill=Pathogenicity)) +
    geom_boxplot() +
    theme(text=element_text(size=15, face = "bold"),axis.text=element_text(size=15, face = "bold"))+
    ylab(C)
    NAME <- paste(sep="", C, ".png" )
    ggsave(file=NAME, width=5, height=5, dpi=300)
    '
    rm $i.lst
done
rm Rplots.pdf





