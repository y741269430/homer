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