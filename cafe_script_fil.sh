#!shell
date

#specify data file, p-value threshold, # of threads to use, and log file
load -i filtered_cafe_input_prueba.txt -filter -t 1 -l log_run_scriptCAFE_filter.txt

#the phylogenetic tree structure with branch lengths
tree (Pi:35,(((Date:0.1,Dati:0.1):30,((Cag:1.3,Cam:1.3):17.7,(Stu:7.9,(Spe:3.6,(Sly:1.5,Spi:1.5):2.1):4.3):11.1):11):1,(Nto:10,(Nat:7,(Nsy:4.2,Ntab:4.2):2.7):3):21):4)

#search for 1 parameter model
lambda -s -t (1,(((1,1)1,((1,1)1,(1,(1,(1,1)1)1)1)1)1,(1,(1,(1,1)1)1)1)1)

report resultfileCAFE_run
