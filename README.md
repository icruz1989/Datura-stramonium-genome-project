# *Datura stramonium* genome project (the devil weed)

* This is the workflow that was followed to assembly and annotate the genome of *Datura stramonium* and for comparative genomics with 13 Solanaceae species. Many programs have been used and links of some of them have been stored here. Also, other pipelines have been followed and a link of each one will be found in this workflow. 


### Firts step is trimming poor quality sequences

    java -jar trimmomatic-0.33.jar PE -phred30 Tic23_S155_L006_R1_001.fastq Tic23_S155_L006_R2_001.fastq output_Tic23_S155_L006_R1_001_paired.fastq output_Tic23_S155_L006_R1_001_unpaired.fastq output_Tic23_S155_L006_R2_001_paired.fastq output_Tic23_S156_L006_R2_001_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36

### Fastqc program was used to visualize the quality of the sequences

    fastqc output_Tic23_S155_L006_R1_001_paired.fastq output_Tic23_S155_L006_R2_001_paired.fastq

### Genome size estimation with Kmergenie

   reads_file contains the name of the files of PE sequences (trimmed) from Illumina (per line)
   
    kmergenie-1.7048/kmergenie reads_file.txt 

   One this step is accomplished, now work on PacBio Sequences
   
### Transform PacBio subreads from .bam format to fastq or fasta, fasta files are much bigger in size

    bamtools convert -format fastq -in output_tic_pacbio_merged.bam -out tic23.subreads.fastq

### Firts, we produced an assembly using only the longreads, Canu program is the option. Canu detects the number of cores 

    canu -d ensamble_tic_pacbio -p datura_tic23 genomeSize=1.5g gridOptionscormhap="--mem=40g" merylMemory=62 batMemory=62         corMhapSensitivity=high correctedErrorRate=0.105 corOutCoverage=100 corMinCoverage=0 gridOptions="--time=168:00:00 --         partition=FAST" gnuplotTested=true -pacbio-raw tic23.subreads.fastq 1>run2.log

### At the same time, generate contigs only with the Illumina PE sequences using the SparseAssembler program. This will generate an assembly only with Illumina sequences. Kmergenie program used in the step above gives the best kmer to produce an assembly. Use this one to feed SparseAssembler

    SparseAssembler LD 0 k 73 g 15 NodeCovTh 1 EdgeCovTh 0 GS 15000000 i1 output_Tic23_S155_L006_R1_001_paired.fastq i2 output_Tic23_S155_L006_R2_001_paired.fastq

### Once both assemblies are done (Canu and SparseAssembler), now do an Hybrid assembly but using raw reads from PacBio and contigs.txt file from SparseAssembler. This program uses De bruin graph and overlap layout consensus algorithms

    DBG2OLC LD 1 k 17 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 10 RemoveChimera 1 Contigs Contigs.txt f tic23.subreads.fastq
  
 A file called scaffolds.fasta is generated, this is the hybrid assembly

### In this step we want to align Canu assembly (Canu) to the Hybrid assembly (DBG2OLC)

    nucmer --mumreference -l 100 self.fasta hybrid.fasta
    
 self.fasta corresponds to the Canu assembly file and hybrid.fasta corresponds to the hybrid assembly from DBG2OLC

### Once alignment from above step is done, now use Quickmerge, Ncmer gives you an out.delta file and use again the hybrid assembly and the Canu assembly

    quickmerge -d out.delta -q hybrid.fasta -r self.fasta -hco 3.0 -c 1.1 -l 70000 -ml 7000 -o genome_tic.fasta
   
   A final merged assembly is obtained. This is the BackBone

## Polishing and scaffolding

### Bowtie2 was used to index the genome and we aligned the raw Illumina reads to the merged genome
   
   Firts, index the genome
   
    bowtie2-build --threads 30 genome_tic.fasta bt2_index_genomeTic
   
   Align the Illumina PED sequences to the BackBone
    
    bowtie2 -x bt2_index_genomeTic -1 output_Tic23_S155_L006_R1_001_paired.fastq -2 output_Tic23_S155_L006_R2_001_paired.fastq     -S tic_gen.sam
   
   Convert the aligned file; sam to bam.
    
    samtools view -Sb tic_gen.sam > tic_gen.bam
   
   Index the aligned file and sort
   
    samtools index tic_gen.bam | sort > tic_gen.sorted.bam
    
### Now, its time to use Pilon for polishing the genome. Three rounds of polishing are recommended. Use the aligned file and Backbone to feed Pilon
    
    java -jar pilon-1.22.jar --genome genome_tic.fasta —-b tic_gen.sorted.bam —-fix bases > genome_tic_pilon1.fasta
   
   In this step, align the raw PacBio sequences to BackBone already polished with Pilon
     
    pbalign --minAnchorSize 15 --maxMatch 20 --nproc 15 output_tic_pacbio_merged.bam genome_tic_pilon1.fasta pbalign_tic.bam
    
   Now use arrow to obtain a consensus sequence and a polished genome
    
    arrow pbalign_tic.bam -r genome_tic_pilon1.fasta -o consensus_tic_pilon1_arrow.fasta

### All is ready for scaffolding, OPERA-LG uses the raw Illumina PE sequences, the genome and long reads

    perl OPERA-long-read.pl --short-read-maptool bowtie2 --opera /LUSTRE/Genetica/ivan/bin_app/OPERA-LG_v2.0.6/bin —num-of-       processors 5 --kmer 17 --contig-file consensus_tic_pilon1_arrow.fasta --illumina-read1 reads_1.fasta --illumina-read2         reads_2.fasta --long-read-file raw_reads_pacbio_tic.fasta --output-prefix opera_lr --output-directory RESULTS

### Two polishing steps with PILON were carried out but now to the last version of the genome (output from OPERA-LG). Seesteps above

### Now its time to check the quality and completness of the final assembly
   
   Quast gives you summary statistics
   
    quast.py genome.polished.draft.fasta
   
   BUSCO looks in a specific database set by the user the number of Single copy Orthologs ans gives you the percentage of genome completness
   
    python run_BUSCO.py -r -i final_genome_teotihuacan.fasta -o busco_finaldraft_genome_teotihuacan -l solanaceae_odb10 -m geno

# Repetitive elements analysis

Some script used in this section are already stored in RepeatModeler or RepeatMasker programs. See manual of these to more information

see http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction--Basic

### The firts step for repetitive elements identification is build a databse

    BuildDatabase -engine ncbi -name datura GenomeTic.fa

### Then use RepeatModeler to obtain specific repetitive elements ofe Datura genome
    
    RepeatModeler -database datura -engine ncbi -pa 32
    
### Repeatmodeler produces a consensi.fa.classified file with knows and unknows repetitive elements. Thus, to separete both of them we used

    perl repeatmodeler_parse.pl --fastafile consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta  \
    --identities repeatmodeler_identities.fasta 

### Thus we obtained a file with only unknows sequences. These are searched against a transposase database and sequences matching transposase are considered as transposons belonging to the relevant superfamily and are incorporated into repeatmodeler_identities and excluded from repeatmodeler_unknowns. Once filtering is complete using the script transposon_blast_parse.pl the libraries ModelerUnknown.lib and ModelerID.lib are created. The trans database can be obtained downloaded from http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced#3._Collecting_repetitive_sequences_by_RepeatModeler

    makeblastdb -in Tpases020812  -dbtype prot
    
    blastx -query repeatmodeler_unknowns.fasta -db Tpases020812  -evalue 1e-10 -num_descriptions 10 -out                           modelerunknown_blast_results.txt
    
    transposon_blast_parse.pl --blastx <blastx output file> --modelerunknown <full path of the modeler unknown file>
    
   Outputs:

   identified_elements.txt

   unknown_elements.txt

   Rename/combine files from RepeatModeler using following commands:
   
    mv  unknown_elements.txt  ModelerUnknown.lib
    cat identified_elements.txt  repeatmodeler_identities.fasta  > ModelerID.lib
   
### All repeats collected so far are used to search against a plant protein database where proteins from transposons are excluded. Elements with significant hits to genes are removed, along with 50 bp upstream and downstream of the blast hit. Remaining sequence that is less than 50 bp is removed completely. Outputs from this script are elements with no significant blast hits to the protein database and the remaining sequence from elements with blast hits that is greater than 50 bp. Use ProExcluder prograam for this purpose
   
   Database alluniRefprexp070416 can be downloaded here 
   
  http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-           Advanced#3._Collecting_repetitive_sequences_by_RepeatModeler
    
    blastx -query ModelerUnknown.lib -db alluniRefprexp070416  -evalue 1e-10 -num_descriptions 10 \
    -out ModelerUnknown.lib_blast_results.txt
    
    ProtExcluder.pl -option ModelerUnknown.lib_blast_results.txt  ModelerUnknown.lib
    
 ### Now obtain a fasta file from RepeatMasker. This program uses the database Dfam and RepeatDabase. Use all to search in the entire database
 
    queryRepeatDatabase.pl -species all > repeatmasker.all.fa
 
   Concatenate now the RepeatMasker database with the one generated manually above
  
    tail -n+2 repeatmasker.taxon.fa | cat - taxon.consensi.fa.classified > taxon.repeatlib.fa
 
 ### Finally run RepeatMasker with this databse
 
    RepeatMasker -pa 20 -q -lib taxon.repeatlib.fa -gff final_genome_tic23.fasta
   
# Gene annotation 

### Firts we obtained and trained gene models using Augustus inside BUSCO program. See BUSCO documentation for more information. Sugustus produces gene models and this models will be used in MAKER

    python run_BUSCO.py -r -i final_genome_tic23.fasta -o busco_finaldraft_genome_teotihuacan -l sonaceae_odb10 -m geno -long      
### MAKER was run following this pipeline 

http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018

### Running MAKER 
   
   MAKER uses several other programs that run inside MAKER and also uses three main configuration files
   
   - maker_exe.ctl - contains the path information for the underlying executables    
   
   - maker_bopt.ctl - contains filtering statistics for BLAST and Exonerate
   
   - maker_opt.ctl - contains all other information for MAKER, including the location of the input genome filem input proteins      and transcripts and gene models

### Once these configuration files are filled, running MAKER is pretty easy. From the folder where you have stored the configuration files just type

    maker -base maker_gff_final -g draffinal_genome_ticuman.fasta -fix_nucleotides

The outputs from MAKER are:

   - The maker_opts.log, maker_exe.log, and maker_bopts.log files are logs of the control files used for this run of MAKER
   
   - The mpi_blastdb directory contains FASTA indexes and BLAST database files created from the input EST, protein, and repeat      databases
   
   - The dpp_contig_master_datastore_index.log contains information on both the run status of individual contigs and information    on where individual contig data is stored
   
   - The dpp_contig_datastore directory contains a set of subfolders, each containing the final MAKER output for individual        contigs from the genomic fasta file

MAKER was run four times, each time was changed the gene models predicted from the prior run 
      
    Now, to obtain the proteins, transcript and a gff file with coordinates of these just use this script and the path of the directory created by MARKER
     
    fasta_merge -d dpp_contig_master_datastore_index.log

    gff3_merge -d dpp_contig_master_datastore_index.log

### AED quality filtering was done using the script from MAKER: quality_filter.pl 

available here: https://groups.google.com/forum/#!searchin/maker-devel/quality_filter.pl%7Csort:relevance/maker-devel/LC4STWWlwgo/XV4nhGiHsfIJ

### Functional annotation was done using the MAKER proteins and transcripts. Names of the genes were edited using the program AHRD, alternative annotation was done using Mercartor with the database MapMan4

A comprenhensive tutorial for funtional annotations is found in below links

Link for MAKER pipeline: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018 

Link for AHRD program: https://github.com/asishallab/AHRD https://github.com/asishallab/AHRD

### Alternative, these scripts are very useful to obtain proteins and CDS from gff3 format files 

scripts are part of https://github.com/NBISweden/GAAS

    gff3_sp_extract_sequences.pl -gff Final_just_genemodels_tic23_AED_0.5 -f draffinal_genome_ticuman.fasta -p -o mkr_snap3_final.all.maker.proteins_extractions.fasta 

    gff3_sp_extract_sequences.pl -t cds --cfs -g s.pimp.mgm.FINAL.gff -f tomato-scaffolds.abyss.77.fasta -o cds_solpi.fasta &

### InterproScan command to annotate domains of proteins 

    interproscan.sh -appl TIGRFAM,SUPERFAMILY,Coils,ProSiteProfiles,SMART,PRINTS,ProSitePatterns,Pfam,ProDom -dp -f TSV -goterms -iprlookup -pa -t p -i Capsicum_annuum_cvCM334.fasta Capsicum_annuum_cvCM334.iprscan >& run1.out

# Comparative analyses

### OrthoFinder program to construct Orthologs and paralogs, you can have proteins families from this program

    orthofinder -S diamond -T fasttree -M msa -t 8 -a 8 -f /home/icruz/data/orthofinder_all/all_solanaceas__ortho_correction_CAFE/cafe_headers >& run.out &

### Constructing an ultra metric three 

Make your three ultrametric with the program of your preference, could be with the R package Ape and its Chronos function

### Protein families expansions and contractions 

    cafetutorial_clade_and_size_filter.py -i cafe_all_sol.tab -o filtered_cafe_input_prueba.tab -s

    cafe cafe_script_fil.sh

    cafetutorial_report_analysis.py -i resultfile_run2.cafe -o summary_run_2_cafe

### Families with signal of expansion and contraction can be retrieve from CAFE outputs, then this script will carry out enrichment test in protein families with significant expansion and contraction 

For InterPro function annotation (domains) 

https://github.com/asishallab/SlydGeneFamsAnalyses/blob/icruz/exec/enrichedAnnosInExpContrFams.R

For Mapman4 enrichment funtion annotation

https://github.com/asishallab/SlydGeneFamsAnalyses/blob/icruz/exec/enrichedMapManInMappFams.R 

### Running MAPP for physicochemical divergence 

##### first we have to retrieve all of the multiple sequence alignments and gene therefore each protein family from Orthofinder results. OrthoFinder create two folders with all of the MSA and genes three. This script will do that and create a folder for each protein family storing both files (MSA and the gene three)

https://github.com/asishallab/SlydGeneFamsAnalyses/blob/icruz/exec/readAndParseOrthogroupsTxt.R

### Running MAPP en a loop mode for all gene families

    for d in *; do
    if [ -d "$d" ]; then         # or:  if test -d "$d"; then
    ( cd "$d" && java -jar /home/icruz/data/MAPPrun/MAPPprogram/MAPP.jar -f ${d}_edited_AA_MSA.fa -t ${d}_phyl_tree.newick -o ${d}_mapp.out )
    fi;
    done

### Once MAPP is done use this script to read the results and also obtain a list with enriched proteins and its function 

    readMappResults.R --workDir /working_dir --outDir resulstMAPP --allInterProTable "here load the Interproannotations of the proteins"/"Also load MapMan4 anotation"

### Positive selection analyses with FUBAR, this pipeline require CDS and the multiple sequences aligments from OrthoFinder. The CDS sequences need to have the same ID name of the alignments, protein alignments are already stored in each folder from previous analyses. Each protein family needs its CDS to construct a CDS alignment in each family. 

### Remove stop codons from proteins alignments in each protein family

    for f in *.fa; do sed 's/*//g' $f > ${f%\.*}_cleanedAs.fa; done

# Concatenate all the CDS

    cat cds.sl.fasta cds.ds.fasta. cds.spe.fasta alltheCDS.fasta

### Retrieve all the ID names from each protein aligment in all the folders

    for d in *; do if [ -d "$d" ]; then ( cd "$d" && grep "^>" *_AA_MSA_orig_gene_ids.fa > ${d}_listIDS.txt ); fi; done &> runningheaderextrac.txt &

### Extract the CDS specific for each protein family in all the folders

    for d in *; do if [ -d "$d" ]; then  ( cd "$d" && perl seqextr.pl $f Solanumtuberosum.fasta > ${f%\.*}_CDS.fa; done) fi; done

### Running PAL2NAL to generate a CDS alignment for all folders

    for d in *; do if [ -d "$d" ]; then ( cd "$d" && perl /home/icruz/apps/pal2nal.v14/pal2nal.pl *_AA_MSA_orig_gene_ids.fa *_CDS.fa -codontable 1 -output fasta > ${d}_pal2nal_CDS.fa ); fi; done

### Running FUBAR for all folders 

    for d in *; do if [ -d "$d" ]; then ( cd "$d" && hyphy fubar --alignment ${d}_pal2nal_CDS.fa --tree ${d}_phyl_tree_orig_gene_ids.newick ); fi; done &> running_FUBAR.txt &

### Parsing FUBAR outputs

    myext = phyphy.Extractor("/home/icruz/data/PRUEBA_MEME_fams/PRUEBA_MEME_withnewpipeline/working_dir/MEME_22001_23000/OG0022001/OG0022001_pal2nal_CDS.fa.MEME.json")
    myext.extract_csv("fest.tsv", delim = "\t")

### Parsing .json format for all directories

    for d in *; do if [ -d "$d" ]; then ( cd "$d" && python /LUSTRE/Genetica/ivan/prueba_parsing_fubar/parsejason_v2.py -f *_pal2nal_CDS.fa.FUBAR.json ); fi; done &> running_parseJSON.txt &

### Now you will have a txt file in each directory or folder with the sites under positive selection, the next step is to parse that file to just keep the sites with a byes factor > 100, significant posterior probabilities ≥ 0.98. Use this script for the purpose. 

https://github.com/asishallab/GeneFamilies/blob/master/exec/loadFubarResults.R

### Enrichment tests can be done using InterProscan and Mapman4 annotations with this script 

Read also these scripts enrichedMapManInMappFams.R and enrichedAnnosInExpContrFams.R

https://github.com/asishallab/SlydGeneFamsAnalyses/blob/icruz/exec/identifyDomainsAtSelectedSites.R
