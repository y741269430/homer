# homer

- 0.Build source
- 1.homer peak files (R)
- 2.Perform the bash

## 0.Build source  

Firstly, we create conda source to perform homer analysis.  

    conda create -n homer
    conda activate homer
    conda install -c bioconda homer
    ...  

And then we download the genome files. (It depends on your internet speed.)  

    perl /home/yangjiajun/miniconda3/envs/homer/share/homer/.//configureHomer.pl -list
    perl /home/yangjiajun/miniconda3/envs/homer/share/homer/.//configureHomer.pl -install mm10

## 1.homer peak files (R)  

    projPath = "/home/yangjiajun/old/cutcfa/"
    load(paste0(projPath, "featurecounts/cut_Anno_df.RData"))
        
And then we divided the chromatin into promoter regions and gene body regions.  

![cut%26tag_peak_anno.png](https://github.com/y741269430/homer/blob/main/cut%26tag_peak_anno.png)  

    rmdis_bed <- lapply(peakAnno_df[c(1:4,10,11)], function(x){
      x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
      x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
      x <- x[, c(8,1:3)]
      x$Strand <- 0
      colnames(x)[1] <- 'Peak ID'
      return(x)
    })
    
    for (i in 1:length(rmdis_bed)) {
      write.table(rmdis_bed[i],
                  paste0(projPath, 'homer/', names(rmdis_bed[i]), "_homer.peak"),
                  sep = "\t", row.names = F, col.names = F, quote = F)
    }

![homer_peak_anno.png](https://github.com/y741269430/homer/blob/main/homer_peak_anno.png)  

## 2.Perform the bash

    vim h1_homer.sh
  
    #!/bin/bash
    ## Alignment to mm10 ##
  
    cat filenames | while read i; 
    do
    nohup findMotifsGenome.pl ${i}_homer.peak mm10 MotifOutput_${i}/ -size 200 -mask &
    done

