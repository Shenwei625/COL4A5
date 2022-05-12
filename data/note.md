### 笔记
# 5月9日
# sed -i 's/指定的字符/要插入的字符&/' 文件
```bash
sed -i 's/[A-Z]:[A-Z]/:&/' LOVD.tsv #在指定的字符前面插入字符
```

# R绘制分布图
```R
library(ggplot2)
HISTABLE <- "F.EXON.tsv"
data <- read.table(HISTABLE,header=TRUE,sep = "\t")
# 单独做图
p1 <- ggplot(data, aes(x=EXON,fill=group)) +
+ geom_histogram(binwidth = 1,colour="black", position="identity")
# 或(更改y后面的参数count或者density可以选择作图的形式)
ggplot() +
+ geom_histogram(aes(x=EXON,fill=group,y=..density..), data, binwidth = 1,colour="black", position="identity") +

# 合并作图
p <- ggplot() +
+ geom_histogram(aes(x=EXON,fill=group,y=..density..), data, binwidth = 1,colour="black", position="identity", fill="#00F5FF") + # 设置颜色
+ geom_histogram(aes(x=EXON,fill=group,y=-..density..), data2, binwidth = 1,colour="black", position="identity")

# 修饰
p + scale_x_continuous(breaks = seq(55)) + ylab("density") + theme(panel.grid = element_blank())# 坐标轴刻度+y标签+删去网格线

EXON <- data[,1]
ggplot() +
+ geom_histogram(aes(x=EXON,fill=group,y=..density..), data, binwidth = 1,colour="black", position="identity") +
+ geom_vline(aes(xintercept=median(EXON, na.rm=T)),color="black", linetype="dashed", size=.5) # 添加中位数线
```





# 5月11日
# 绘制树形图
+ 可以利用paste()函数来构建一些变量
```R
install.packages("treemapify")
library(treemapify)
library(ggplot2)
library(patchwork) # 可以将多张图组合在一起(p1 + p2左右组合，或者p1 / p2上下组合)

COUNTF <- "F.Amino_acids.tsv" 
TREEDATA <- read.table(COUNTF,header=TRUE,sep = "\t")
ggplot(TREEDATA, aes(area = count, fill = Amino_acids, label = paste(Amino_acids, count,sep = "," ),subgroup = Amino_acids)) +  
# area后面添加计数列，fill后面添加被统计的列
    geom_treemap() +
    geom_treemap_text(colour = "white",place = "bottomleft",size = 10) +   
# 添加标签,place设置位置（top,bottom.center,left,rght）,前面的label选项可以设置标签内容
    geom_treemap_subgroup_border(colour = "white", size = 3) 
# 设置每个区块之间的间隔，对应于前面的subgroup

ggsave(filename = <variants>,width=20, height=10, dpi=300) # 保存图片，filename后面可以直接使用变量

labs(title = "content") # 设置图片标题，filename后面可以直接使用变量

theme(legend.position = "none",plot.title=element_text(hjust=0.5)) # 删除图例,标题居中
```

# 绘制饼状图
+ 通过ggplot2利用极坐标将柱状图折叠成饼状图
```R
library(ggplot2)
FILE1 <- "T.Consequence.tsv"
PIEDATA <- read.table(FILE1,header=TRUE,sep = "\t")
TITLE1 <- "T.Consequence"
J <-"Consequence"
PLABEL = as.vector(PIEDATA[,J])
myLabel = as.vector(PIEDATA[,J])
myLabel = paste(myLabel, "(", round(PIEDATA$count / sum(PIEDATA$count) * 100, 2), "%" , ")") # 重新定义标签

ggplot(PIEDATA,aes(x="", y=count, fill=get(J))) +  #只需要一个柱子所以x="",y为计数列，fill为分类
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +   # 把柱状图折叠成饼状图
    labs(x = "", y = "", title = TITLE1) + # 删去饼图旁边的标签，原柱状图x,y轴坐标
    theme(axis.ticks = element_blank(), plot.title=element_text(hjust=0.5)) + # 删去坐标轴刻度,标题居中
    theme(axis.text.x = element_blank()) + # 删去x轴标签
    theme(panel.grid = element_blank()) + # 删除网格线
    theme(panel.background = element_rect(fill = "transparent",colour = NA)) + # 删除灰色的背景
    theme(legend.title = element_blank(), legend.position = "bottom") + # 删除图例标题
    scale_fill_discrete( breaks = PIEDATA[,J], labels = myLabel ) + # 重新更换图例中的标签
    theme(legend.text=element_text(size=10)) + # 图例中字体大小
    geom_text(aes(2, label = myLabel),position = position_stack(vjust = 0.5),size = 4) # 在饼图中添加标签

```

