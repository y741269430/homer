# homer

---

具体流程
- http://homer.ucsd.edu/homer/ngs/index.html
- [利用HOMER预测目标序列的motif（从运行程序到结果解读，以及注意事项）](https://www.jianshu.com/p/467d970ec097)    

安装过程
- http://homer.ucsd.edu/homer/introduction/install.html

方法参考
- [Nat Struct Mol Biol 30, 948–957 (2023)](https://www.nature.com/articles/s41594-023-01021-8#Sec15)

## 0. 配置环境
```bash
mamba create -n homer
mamba activate homer 
mamba install homer
mamba install wget samtools r-essentials bioconductor-deseq2 bioconductor-edger 
```

## 1. 构建或下载基因组
- http://homer.ucsd.edu/homer/introduction/configure.html

## 1.1 直接下载构建好的基因组
```bash

# 使用该命令查看可供下载的基因组（如果是自己构建的，也可以看的到）
perl ~/.conda/envs/homer/share/homer/configureHomer.pl -list

# 下载mm39基因组（速度很慢建议自己构建）
nohup perl ~/.conda/envs/homer/share/homer/configureHomer.pl -install mm39 &
```

## 1.2 预先下载mm39，然后自己构建

```bash
loadPromoters.pl -name myMM39-p -org mouse -id refseq -fasta ~/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa -offset 2000 &
loadGenome.pl -name myMM39 -org mouse -fasta ~/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa -gtf ~/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf -promoters &

# 构建完成后通过以下命令查看是否成功
cat ~/.conda/envs/homer/share/homer/config.txt
perl ~/.conda/envs/homer/share/homer/configureHomer.pl -list
```

创建homer文件夹，并将数据存储在里面

## 2.转换macs3的peak文件为homer的bed文件

这一步之前，应该要先去除blacklist [remove-blacklist](https://github.com/y741269430/ATAC-seq#8remove-blacklist)

```bash
vim nar2homer.sh

#!/bin/bash
## narrowPeak to homer ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

# 读取文件名列表
cat filenames | while read i; do
    # 构造输入和输出文件路径
    # input_file="$input_dir/${i}_peaks.narrowPeak"
    input_file="$input_dir/${i}_rmBL.narrowPeak"
    output_file="$output_dir/${i}_homer.bed"

    # 运行 
    awk -F'\t' '{print $1, $2, $3, $4, $9, $6}' OFS='\t'  "$input_file" > "$output_file" &
done
```

## 3.homer 预测motif

```bash
vim h1_homer.sh

#!/bin/bash
## homer find motif mm39 ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

cat filenames | while read i; 
do
input_file="$input_dir/${i}_homer.bed"
output_file="$output_dir/MotifOutput_${i}/"

nohup findMotifsGenome.pl "$input_file" mm39 "$output_file" -size 200 -mask &

done
```

## 4.homer 预测motif (RNAseq)
使用一个csv列表，以基因为列，去除行名列名，放到mouse数据库中，进行motif查找     
```bash
findMotifs.pl h.csv mouse MotifOutput/ -rna -len 8 &
```

## 5.R语言统计motif结果

以下是使用的函数
```r
subString <- function(strings, idx, sep = NA){
  
  strings = as.character(strings)
  if(is.na(sep)){
    res = as.character(lapply(strings, function(x) paste(strsplit(x, "")[[1]][idx], collapse = "")))
  } else{
    res = sapply(strsplit(strings, sep), function(x) x[idx])
  }
  return(res)
}
summaryHomer <- function(outFolder){
  
  homerFolder = paste0(outFolder, "/homerResults")
  xFiles = list.files(homerFolder, ".motif$")
  xFiles = xFiles[-grep("similar", xFiles)]
  xFiles = xFiles[-grep("RV", xFiles)]
  xFiles = xFiles[order(as.numeric(gsub("\\.", "", gsub("motif", "", xFiles))))]
  texts  = sapply(paste0(homerFolder, "/", xFiles), readLines)
  chunks = sapply(texts, function(x) strsplit(x[1], "[\t]"))
  
  motif = sapply(chunks, function(x) subString(x[1], 2, ">"))
  match = sapply(chunks, function(x) subString(subString(x[2], 2, "BestGuess:"),  1, "/"))
  score = sapply(chunks, function(x) rev(strsplit(x[2], "[()]")[[1]])[1])
  count = sapply(chunks, function(x) subString(x[6], 3, "[T:()]"))
  ratio = sapply(chunks, function(x) subString(x[6], 2, "[()]"))
  p_value = sapply(chunks, function(x) subString(x[6], 2, "P:"))
  
  xresT = data.frame(motif, 
                     match, 
                     score = as.numeric(score), 
                     count = as.numeric(count),
                     ratio_perc = as.numeric(gsub("%", "", ratio)), 
                     p_value = as.numeric(p_value)
  )
  rownames(xresT) = gsub(".motif", "", basename(rownames(xresT)))
  return(xresT)
}
summaryHomerKnown <- function(outFolder){
  
  knownFolder = paste0(outFolder, "/knownResults")
  xFiles = list.files(knownFolder, ".motif$")
  xFiles = xFiles[order(as.numeric(gsub("\\.motif", "", gsub("known", "", xFiles))))]
  texts  = sapply(paste0(knownFolder, "/", xFiles), readLines)
  chunks = sapply(texts, function(x) strsplit(x[1], "[\t]"))
  
  motif = sapply(chunks, function(x) subString(x[1], 2, ">"))
  TF    = sapply(chunks, function(x) subString(x[2], 1, "/"))
  count = sapply(chunks, function(x) subString(x[6], 3, "[T:()]"))
  ratio = sapply(chunks, function(x) subString(x[6], 2, "[()]"))
  p_value = sapply(chunks, function(x) subString(x[6], 2, "P:"))
  
  xresT = data.frame(motif, 
                     TF, 
                     count = as.numeric(count),
                     ratio_perc = as.numeric(gsub("%", "", ratio)), 
                     p_value = as.numeric(p_value)
  )
  rownames(xresT) = gsub("\\.motif", "", basename(rownames(xresT)))
  return(xresT)
}
get_TF_Vennlist <- function(x){
  
  inter <- VennDiagram::get.venn.partitions(x)
  inter$values <- sapply(inter$..values.., paste, collapse = "~")
  output <- lapply(inter$values, function(x){ x <- unlist(strsplit(x, "~"))})
  
  return(output)
}

```
#### 寻找已知的转录因子结合区域 ####
```r
known_find <- list(CON = summaryHomerKnown('MotifOutput_merge_CON/'),
                   Tre = summaryHomerKnown('MotifOutput_merge_Tre/') )

known_find2 <- lapply(known_find, function(x){x <- x$TF })
ggvenn(known_find2)

# 获取venn图TF集合
a1 <- get_TF_Vennlist(known_find2[c(1,2)])
TF_tre_vs_con <- list(only_con = known_find[[1]][known_find[[1]]$TF %in% a1[[3]], ],
                      overlap_in_con = known_find[[1]][known_find[[1]]$TF %in% a1[[1]], ],
                      only_tre = known_find[[2]][known_find[[2]]$TF %in% a1[[2]], ],
                      overlap_in_tre = known_find[[2]][known_find[[2]]$TF %in% a1[[1]], ])
```

#### 寻找潜在的转录因子结合区域 ####
```r
homer_find <- list(CON = summaryHomer('MotifOutput_merge_CON/'),
                   Tre = summaryHomer('MotifOutput_merge_Tre/'),)

homer_find2 <- lapply(homer_find, function(x){x <- x$TF })

```

---
## 不太记得下面的代码是用来做什么的了，之前好像是用homer直接call peak

## 2.Creation Tag directories, quality control, and normalization. (makeTagDirectory)
- [Creating a "Tag Directory" with makeTagDirectory](http://homer.ucsd.edu/homer/ngs/tagDir.html)
- https://mp.weixin.qq.com/s?src=11&timestamp=1727681189&ver=5537&signature=HcXSefB-Uf3WAeJ3qYE-jtoh59ts2mjCGaob58KwIyjbyxt2uhKCbU3pmbZS8KFiviCvWO2iM*GZ83QcTohGInsEuUwc1vgdwlnTnypBOYeN24QTD2Fn-lqL0oK8BfGz&new=1

```bash
nohup cat filenames | while read i; do makeTagDirectory ./homer/${i} -genome myMM39 -checkGC ./bam/${i}.last.bam; done
```

## 3.Peak finding / Transcript detection / Feature identification (findPeaks, getDifferentialPeaksReplicates.pl)

```bash
nohup cat filenames | while read i; do findPeaks ./homer/${i} -style factor -o auto; done &
```

## 转换macs3的peak文件为homer文件
- http://homer.ucsd.edu/homer/ngs/peaks.html

```bash
nohup cat filenames | while read i; do bed2pos.pl ../macs3/${i}_peaks.narrowPeak > ./${i}/peaks.txt; done &
```


## 9.Peak finding / Differential Peak calling with Replicates (getDifferentialPeaksReplicates.pl)
```bash
nohup getDifferentialPeaksReplicates.pl \
-t SUS_1/ SUS_2/ \
-i CON_1/ CON_2/ -genome myMM39 > deg_SUSvsCON.csv -all -f 1.5 &

nohup getDifferentialPeaksReplicates.pl \
-t CONP_1/ CONP_2/ \
-i CON_1/ CON_2/ -genome myMM39 > deg_CONPvsCON.csv -all -f 1.5 &

nohup getDifferentialPeaksReplicates.pl \
-t SUSP_1/ SUSP_2/ \
-i SUS_1/ SUS_2/ -genome myMM39 > deg_SUSPvsSUS.csv -all -f 1.5 &

nohup getDifferentialPeaksReplicates.pl \
-t SUSP_1/ SUSP_2/ \
-i CONP_1/ CONP_2/ -genome myMM39 > deg_SUSPvsCONP.csv -all -f 1.5 &
```





