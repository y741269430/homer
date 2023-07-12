# homer

## 0.Build source  

Firstly, we create conda source to perform homer analysis.  

    conda create -n homer
    conda activate homer
    conda install -c bioconda homer
    ...  

And then we download the genome files. (It depends on your internet speed.)  

    perl /home/yangjiajun/miniconda3/envs/homer/share/homer/.//configureHomer.pl -list
    perl /home/yangjiajun/miniconda3/envs/homer/share/homer/.//configureHomer.pl -install mm10

Perform the bash.

    vim h1_homer.sh
  
    #!/bin/bash
    ## Alignment to mm10 ##
  
    cat filenames | while read i; 
    do
    nohup findMotifsGenome.pl ${i}_homer.peak mm10 MotifOutput_${i}/ -size 200 -mask &
    done

