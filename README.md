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
这里重点展示两类最常用的结果：knownResults.html 和 homerResults.html。

一、homerResults.html — 从头发现的 Motif（De Novo）    

这个 HTML 展示的是 HOMER 不依赖先验数据库、直接从 Peak 序列里挖掘出来的重复模式。和 knownResults 的区别是：knownResults 是“拿已知 motif 去数据里找”，而 homerResults 是“让数据自己告诉你可能存在什么 motif”。    

| 名称 | 解释 |
|:-----:|:-----:|
| Motif | 从头发现的 motif 序列 logo |
| P-value / Log P-value | 该 motif 的富集显著性 |
| % Target | Peak 序列中包含该 motif 的比例 |
| % Background | 背景序列中包含该 motif 的比例 |
| Best Match | 和已知 motif 数据库最相似的匹配结果 |
| Match Score | 匹配得分，通常越接近 1 说明越相似 |

怎么判断：
先看 `Best Match`，判断 de novo motif 最像哪个已知转录因子。
`Match Score > 0.9` 时，通常说明这个 de novo motif 和对应已知 motif 非常接近。
`Match Score` 中等时，可能代表同一家族 motif 的变体。
如果 de novo 结果和 knownResults 中排名靠前的 motif 彼此一致，通常说明结果更可靠。
motif 富集结果本身不能单独证明直接结合关系，仍建议结合 peak 注释、功能富集和 IGV 信号综合解释。

二、knownResults.html — 已知 Motif 富集结果    

用浏览器打开这个 HTML 文件，可以看到一个表格，展示的是：在你的 Peak 区域中，哪些已知转录因子结合 motif 被显著富集了。HOMER 会把目标 Peak 序列和背景序列进行比较，如果某个已知 motif 在目标序列中出现得明显更多，就说明它被富集了。

| 名称 | 解释 |
|:-----:|:-----:|
| Motif | motif 的序列 logo 图，用来可视化展示碱基偏好 |
| Name | motif 名称或对应的转录因子 |
| P-value | 富集显著性，越小通常越显著 |
| Log P-value | P 值的对数形式，用于排序和展示 |
| q-value | 多重检验校正后的显著性水平 |
| # / % Target with Motif | 目标 Peak 序列中包含该 motif 的数量和比例 |
| # / % Background with Motif | 背景序列中包含该 motif 的数量和比例 |

怎么判断：
`P-value` 越小，说明该已知 motif 在 Peak 中富集越显著。
`% Target` 明显高于 `% Background`，说明该 motif 更可能在目标 Peak 中特异出现。
排在最前面的 motif 通常最值得优先关注。
富集到的转录因子是否合理，最好结合你的细胞类型、靶蛋白和实验背景一起解释。


## 4.homer 预测motif (RNAseq)
使用一个csv列表，以基因为列，去除行名列名，放到mouse数据库中，进行motif查找     
```bash
findMotifs.pl h.csv mouse MotifOutput/ -rna -len 8 &
```
---
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

