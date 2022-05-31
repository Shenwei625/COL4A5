#! usr/bin/bash

# 将模板中需要替换的部分替换
cat template.md |
    perl -ne '
        chomp($_);
        if (/^sub/) {
            s/sub//;
            print "$_\n";
        }
    ' > JOB.lst

JOB=$(cat JOB.lst)
rm JOB.lst

for J in $JOB;do
    echo -e "开始对$J进行替换\n"
    export PA=$(find -maxdepth 2 -type f -name "$J")
    export SUB=sub$J
    cat template.md |
        perl -n -MPath::Tiny -e '
            my $SUBS = $ENV{"SUB"};
            $data = path($ENV{"PA"})->slurp;
            s/$SUBS/$data/;
            print "$_";
        ' > tem&&
        mv tem template.md
done

pandoc -i template.md -o report.html

