#! usr/bin/bash
file1=$1

tsv-filter -H --str-ne topmed_AF:- $file1 |
    tsv-select -H --fields '#Uploaded_variation',topmed_AF > raw.tsv
sed -i '1d' raw.tsv
sed -i 's/:/\t/3' raw.tsv
(echo -e "REF\tALT\tAF" && cat raw.tsv) > tem&&
    mv tem raw.tsv

# 对科学计数法进行转化
cat raw.tsv | perl -na -F"\t" -e '
    chomp($_);
    if ( @F[2] =~ /E/) {
        print "@F[0]\t@F[1]\t";
        printf "%*.*f\n", 11, 10, @F[2];
    }else {
        print "$_\n";
    }
' > tem&&
    mv tem raw.tsv

tsv-summarize -H --count -g REF raw.tsv > marker.tsv
ALT=$(tsv-summarize -H --count -g REF raw.tsv | tsv-summarize -H --count -g count | sed '1d' | cut -f 1)

for A in $ALT;do
    tsv-filter -H --str-eq count:$A marker.tsv | cut -f 1 > $A.ALT.marker.tsv
done
rm marker.tsv

for A in $ALT;do
    cat raw.tsv | grep -f $A.ALT.marker.tsv > $A.filter.tsv
done
rm *.marker.tsv

# 对一个位点有两个碱基的ALT的MAF进行计算（小于0.5）
cat 1.filter.tsv | perl -ne '
    chomp($_);
    if (/^REF/) {
        print "$_\tMAF\n";
    }elsif (/(\S+)\s(\S)\s(0.\d+)/) {
        if ( $3 <= 0.5 ) {
            print "$_\t$3\n";
        }else {
            my $MAF = ( 1 - $3 );
            print "$_\t$MAF\n";
        }
    }
' > merge.tsv
tsv-select -H --fields REF,ALT,MAF merge.tsv > tem&&
    mv tem merge.tsv
rm 1.filter.tsv

# 对一个位点有三个碱基的情况的MAF计算（取第二小）
if [ -e 2.filter.tsv ]; then
    JOB=$(cat 2.filter.tsv | sed '1d' | cut -f 1 | uniq )
    for J in $JOB;do
        NAME=$(echo $J | tr ":" "_")
        cat 2.filter.tsv | grep "$J" > $NAME.tsv
        cat $NAME.tsv | 
            tr "\n" "\t" | perl -a -F"\t" -e '
            my $REST = ( ( 1 - @F[2] ) - @F[5] );
            ( $REF = @F[0] ) =~ s/^X:\d+://;
            print "@F[o]\t@F[1]\t@F[2]\n@F[3]\t@F[4]\t@F[5]\n@F[0]\t$REF\t$REST\n";
        ' | sort -n -k3,3 > tem&&
            mv tem $NAME.tsv
        MAF=$(sed -n '2p' $NAME.tsv | cut -f 3 )
        cat 2.filter.tsv | grep "$J" > $NAME.tsv
        for i in 1 2;do
            echo -e "$MAF" >> paste2.tsv
        done
        cut -f 1,2 $NAME.tsv > paste1.tsv
        paste -d "\t" paste1.tsv paste2.tsv >> merge.tsv
        rm paste1.tsv paste2.tsv $NAME.tsv
    done    
    rm 2.filter.tsv
else 
    echo -e "没有一个位点存在三个碱基的情况\n"
fi

# 对一个位点有四个碱基的情况的MAF的计算（取第二小）
if [ -e 3.filter.tsv ]; then
    JOB=$(cat 3.filter.tsv | sed '1d' | cut -f 1 | uniq )
    for J in $JOB;do
        NAME=$(echo $J | tr ":" "_")
        cat 3.filter.tsv | grep "$J" > $NAME.tsv
        cat $NAME.tsv | 
            tr "\n" "\t" | perl -a -F"\t" -e '
            my $REST = ( ( ( 1 - @F[2] ) - @F[5] ) - @F[8] );
            ( $REF = @F[0] ) =~ s/^X:\d+://;
            print "@F[o]\t@F[1]\t@F[2]\n@F[3]\t@F[4]\t@F[5]\n@F[6]\t@F[7]\t@F[8]\n@F[0]\t$REF\t$REST\n";
        ' | sort -nr -k3,3 > tem&&
            mv tem $NAME.tsv
        MAF=$(sed -n '2p' $NAME.tsv | cut -f 3 )
        cat 3.filter.tsv | grep "$J" > $NAME.tsv
        cut -f 1,2 $NAME.tsv > paste1.tsv
        for i in 1 2 3;do
            echo -e "$MAF" >> paste2.tsv
        done
        paste -d "\t" paste1.tsv paste2.tsv >> merge.tsv
        rm paste1.tsv paste2.tsv $NAME.tsv
    done
    rm 3.filter.tsv
else 
    echo -e "没有一个位点存在四个碱基的情况\n"
fi



# 添加没有AF的位点
tsv-filter -H --str-eq topmed_AF:- $file1 |
    tsv-select -H --fields '#Uploaded_variation',topmed_AF > blank.tsv
sed -i 's/topmed_AF/MAF/' blank.tsv

# 合并到merge.tsv
sed -i "s/\t/:/1" merge.tsv
sed '1d' merge.tsv >> blank.tsv
rm merge.tsv
mv blank.tsv merge.tsv

tsv-join --filter-file merge.tsv --H --key-fields 1 --append-fields MAF  $file1 > tem&&\
    mv tem $file1
rm merge.tsv raw.tsv