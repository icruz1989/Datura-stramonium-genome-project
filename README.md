# Datura-stramonium-genome-project
This repository contains the commands used in the assembly of Datura stramonium genome and its comparative genomics analyses with other Solanaceae species

# Datura stramonium genome project (the devil weed)

* This is the workflow that was followed to assembly and annotate the genome of Datura stramonium and for comparative genomics with 13 Solanaceae species. Many programs have been used and links of some of them have been stored here. Also, other pipelines have been followed and a link of each one will be found in this workflow.


### Firts step is trimming poor quality sequences

$ java -jar ~/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred30 Tic23_S155_L006_R1_001.fastq Tic23_S155_L006_R2_001.fastq output_Tic23_S155_L006_R1_001_paired.fastq output_Tic23_S155_L006_R1_001_unpaired.fastq output_Tic23_S155_L006_R2_001_paired.fastq output_Tic23_S156_L006_R2_001_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36

### Fastqc program was used to visualize the quality of the sequences

$ fastqc output_Tic23_S155_L006_R1_001_paired.fastq output_Tic23_S155_L006_R2_001_paired.fastq

### Genome size estimation with Kmergenie

$ /LUSTRE/Genetica/ivan/bin_app/kmergenie-1.7048/kmergenie reads_file.txt #reads_files contain PE reads from Illumina

### Transform subreads in .bam format of PacBio Sequel to fastq

$  bamtools convert -format fastq -in output_tic_pacbio_merged.bam -out tic23.subreads.fastq

### Canu assembly

$ canu -d ensamble_teo1_pacbio -p datura_tic23 genomeSize=1.5g gridOptionscormhap="--mem=40g" merylMemory=62 batMemory=62 corMhapSensitivity=high correctedErrorRate=0.105 corOutCoverage=100 corMinCoverage=0 gridOptions="--time=168:00:00 --partition=FAST" gnuplotTested=true -pacbio-raw /home/icruz/data/secuencias_pacbio/bamfiles_tic_g/tic23.subreads.fastq 1>run2.log &

### Generate contigs only with the Illumina PE reads using SparseAssembler

$ SparseAssembler LD 0 k 73 g 15 NodeCovTh 1 EdgeCovTh 0 GS 15000000 i1 /home/icruz/data/sec_illumina_tic23/output_Tic23_S155_L006_R1_001_paired.fastq i2 /home/icruz/data/illumina_tic23/output_Tic23_S155_L006_R2_001_paired.fastq

### Hybrid assembly using raw reads from PacBio and contigs.txt file from SparseAssembler 

$ DBG2OLC LD 1 k 17 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 10 RemoveChimera 1 Contigs /home/icruz/data/sec_illumina_tic23/Contigs.txt f /home/icruz/data/sec_illumina_tic23/tic23.subreads.fastq 

### alignment between Canu (query) assembly and DBG2OLC (reference) assembly

$ /home/icruz/MUMmer3.23/nucmer --mumreference -l 100 self.fasta hybrid.fasta 

### Quickmerge between both assemblies using the output out.delta from Nucmer

$ quickmerge -d out.delta -q hybrid.fasta -r self.fasta -hco 3.0 -c 1.1 -l 70000 -ml 7000

### Polishing and scaffolding

### Bowtie2 was used to index the genome and we aligned the raw Illumina reads to the merged genome

$ bowtie2-build --threads 30 genome_tic.fasta bt2_index_genomeTic

$ bowtie2 -x bt2_index_genomeTic -1 ../../../../../../sec_illumina_tic23/output_Tic23_S155_L006_R1_001_paired.fastq -2 ../../../../../../sec_illumina_tic23/output_Tic23_S155_L006_R2_001_paired.fastq -S tic_gen.sam

$ samtools view -Sb tic_gen.sam > tic_gen.bam

$ samtools index tic_gen.bam | sort > tic_gen.sorted.bam

$ java -jar pilon-1.22.jar --genome genome_tic.fasta —-b tic_gen.sorted.bam —-fix bases > genome_tic_pilon1.fasta

$ pbalign --minAnchorSize 15 --maxMatch 20 --nproc 15 output_tic_pacbio_merged.bam genome_tic_pilon1.fasta pbalign_tic.bam

$ arrow pbalign_tic.bam -r genome_tic_pilon1.fasta -o consensus_tic_pilon1_arrow.fasta

$ /LUSTRE/Genetica/ivan/bin_app/OPERA-LG_v2.0.6/bin/OPERA-long-read.pl --short-read-maptool bowtie2 --opera /LUSTRE/Genetica/ivan/bin_app/OPERA-LG_v2.0.6/bin —num-of-processors 5 --kmer 17 --contig-file /LUSTRE/Genetica/ivan/pacbio_sequences/bamfiles_tic/pbalign/OPERA/consensus_tic_pilon1_arrow.fasta --illumina-read1 /LUSTRE/Genetica/ivan/pacbio_sequences/bamfiles_tic/pbalign/OPERA/reads_1.fasta --illumina-read2 /LUSTRE/Genetica/ivan/pacbio_sequences/bamfiles_tic/pbalign/OPERA/reads_2.fasta --long-read-file /LUSTRE/Genetica/ivan/pacbio_sequences/bamfiles_tic/pbalign/OPERA/raw_reads_pacbio_tic.fasta --output-prefix opera_lr --output-directory RESULTS

### A second polishing step with PILON but now to the last version of the genome

### Nuclear genome validation

$ quast.py genome.polished.draft.fasta

$ python run_BUSCO.py -r -i /home/icruz/data/ensamble_hibrido_teo1/merge/maker/final_genome_teotihuacan.fasta -o busco_finaldraft_genome_teotihuacan -l /home/icruz/data/sec_illumina_tic23/merge/quickmerge2/quickmerge3/quickmerge4/finish/MARKER/busco/solanaceae_odb10 -m geno

### Also a last alignment of the raw Illumina reads to the last version of the genome was carried out to obtain the alignment rates

### Repetitive elements analysis

see http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction--Basic

### Then repeat masker was run 

RepeatMasker -pa 20 -dir finalrepeatmasker_genometicuman -q -species eukaryota -html -gff draffinal_genome_ticuman.fasta

### plotting repeat landscape

$ RepeatMasker/util/calcDivergenceFromAlign.pl -s teo1.divsum -a teo1.align final_genome_teotihuacan.fasta.cat.gz &

$ createRepeatLandscape.pl -g 1471881354 -div tic23.divsum > tic23_repeats.html &

# Gene annotation 

python run_BUSCO.py -r -i /home/icruz/data/ensamble_hibrido_teo1/merge/maker/final_genome_teotihuacan.fasta -o busco_finaldraft_genome_teotihuacan -l /home/icruz/data/sec_illumina_tic23/merge/quickmerge2/quickmerge3/quickmerge4/finish/MARKER/busco/solanaceae_odb10 -m geno Augustus v.3.2.2

### MAKER was run following this pipeline 

* http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018

### Running MAKER 

maker -base maker_gff_final -g draffinal_genome_ticuman.fasta -fix_nucleotides

### AED quality filtering was done using the script from MAKER: quality_filter.pl 

* available here: https://groups.google.com/forum/#!searchin/maker-devel/quality_filter.pl%7Csort:relevance/maker-devel/LC4STWWlwgo/XV4nhGiHsfIJ

### Functional annotation was done using the MAKER and names of the genes were edited using the program AHRD, alternative annotation was done using Mercartor with the database MapMan4

* Link for MAKER pipeline: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018 

* Link for AHRD program: https://github.com/asishallab/AHRD https://github.com/asishallab/AHRD

### This script is very useful to obtain proteins from gff3 format files 

* script is part of https://github.com/NBISweden/GAAS

$ gff3_sp_extract_sequences.pl -gff Final_just_genemodels_tic23_AED_0.5 -f draffinal_genome_ticuman.fasta -p -o mkr_snap3_final.all.maker.proteins_extractions.fasta &

### This is for retrieve CDS (also part of GAAS)

$ gff3_sp_extract_sequences.pl -t cds --cfs -g s.pimp.mgm.FINAL.gff -f tomato-scaffolds.abyss.77.fasta -o cds_solpi.fasta &

### InterproScan command to annotate domains of proteins 

$ interproscan.sh -appl TIGRFAM,SUPERFAMILY,Coils,ProSiteProfiles,SMART,PRINTS,ProSitePatterns,Pfam,ProDom -dp -f TSV -goterms -iprlookup -pa -t p -i Capsicum_annuum_cvCM334.fasta Capsicum_annuum_cvCM334.iprscan >& run1.out &

# Comparative analyses

### OrthoFinder program to construct Orthologs and paralogs, you can have proteins families from this program

orthofinder -S diamond -T fasttree -M msa -t 8 -a 8 -f /home/icruz/data/orthofinder_all/all_solanaceas__ortho_correction_CAFE/cafe_headers >& run.out &

### constructing an ultra metric three 

input.args <- "./SpeciesTree_rooted_node_labels.txt"
tree.rooted <- ape::root(read.tree(input.args), outgroup = "Petuniainflata", resolve.root = TRUE)
nto.na.ns.nt.dte.dti.cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Nicotianatomentosiformis","Nicotianaattenuata","Nicotianasylvestris","Nicotianatabacum","DaturastramoniumTeo1","DaturastramoniumTic23","Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
nto.na.ns.nt.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Nicotianatomentosiformis","Nicotianaattenuata","Nicotianasylvestris","Nicotianatabacum")))
dte.dti.cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("DaturastramoniumTeo1","DaturastramoniumTic23","Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
cag.cam.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Capsicumannunmglabriusculum","Capsicumannummorelia")))
dte.dti.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("DaturastramoniumTeo1","DaturastramoniumTic23")))
st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
aet.split <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Petuniainflata","Solanumlycopersicum")))
chr.df <- data.frame(node = c(nto.na.ns.nt.dte.dti.cag.cam.st.sp.sly.spi.mrca, nto.na.ns.nt.mrca, dte.dti.cag.cam.st.sp.sly.spi.mrca, cag.cam.st.sp.sly.spi.mrca, cag.cam.mrca, dte.dti.mrca, st.sp.sly.spi.mrca, sp.sly.spi.mrca, aet.split), age.min = c(31, 10, 30, 19, 0.01, 7.9, 1.5, 35), stringsAsFactors = FALSE)

chr.df$age.max <- chr.df$age.min
chr.df[[which(chr.df$node == aet.split), "age.max"]] <- 35
chr.df$soft.bounds <- FALSE
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)

#### Scale the rooted tree
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)  # Lambda value obtained from references
### Round the edge lengths:
tree.rooted.chrono$edge.length <- round(tree.rooted.chrono$edge.length, digits = 1)
#' Save the results
write.tree(tree.rooted.chrono, "./allsol_forCAFE_concate_ml_gamma_rescaled_ultrametric.newick")

#### Plot the tree:
pdf( "./all_sol_forCAFE_concate_ml_gamma_rescaled_ultrametric.pdf" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

### Protein families expansions and contractions 

$ cafetutorial_clade_and_size_filter.py -i cafe_all_sol.tab -o filtered_cafe_input_prueba.tab -s

$ cafe cafe_script_fil.sh

$ cafetutorial_report_analysis.py -i resultfile_run2.cafe -o summary_run_2_cafe

### Families with signal of expansion and contraction can be retrieve from CAFE outputs, then this script will carry out enrichment test in protein families with significant expansion and contraction 

* https://github.com/asishallab/SlydGeneFamsAnalyses/enrichedAnnosInExpContrFams.R 

### Running MAPP for physicochemical divergence 

##### first we have to retrieve all of the multiple sequence alignments and gene therefore each protein family from Orthofinder results. OrthoFinder create two folders with all of the MSA and genes three. This script will do that and create a folder for each protein family storing both files (MSA and the gene three)

* https://github.com/asishallab/SlydGeneFamsAnalyses/blob/icruz/exec/readAndParseOrthogroupsTxt.R

### Running MAPP en a loop mode for all gene families

for d in *; do
  if [ -d "$d" ]; then         # or:  if test -d "$d"; then
    ( cd "$d" && java -jar /home/icruz/data/MAPPrun/MAPPprogram/MAPP.jar -f ${d}_edited_AA_MSA.fa -t ${d}_phyl_tree.newick -o ${d}_mapp.out )
  fi;
done

### Once MAPP is done use this script to read the results and also obtain a list with enriched proteins and its function 

$ readMappResults.R --workDir /working_dir --outDir resulstMAPP --allInterProTable "here load the Interproannotations of the proteins"/"Also load MapMan4 anotation"

### Positive selection analyses with FUBAR, this pipeline require CDS and the multiple sequences aligments from OrthoFinder. The CDS sequences need to have the same ID name of the alignments, protein alignments are already stored in each folder from previous analyses. Each protein family needs its CDS to construct a CDS alignment in each family. 

### Remove stop codons from proteins alignments in each protein family

$ for f in *.fa; do sed 's/*//g' $f > ${f%\.*}_cleanedAs.fa; done

# Concatenate all the CDS

$ cat cds.sl.fasta cds.ds.fasta. cds.spe.fasta alltheCDS.fasta

### Retrieve all the ID names from each protein aligment in all the folders

$ for d in *; do if [ -d "$d" ]; then ( cd "$d" && grep "^>" *_AA_MSA_orig_gene_ids.fa > ${d}_listIDS.txt ); fi; done &> runningheaderextrac.txt &

### Extract the CDS specific for each protein family in all the folders

$ for d in *; do if [ -d "$d" ]; then  ( cd "$d" && perl seqextr.pl $f Solanumtuberosum.fasta > ${f%\.*}_CDS.fa; done) fi; done

### Running PAL2NAL to generate a CDS alignment for all folders

for d in *; do if [ -d "$d" ]; then ( cd "$d" && perl /home/icruz/apps/pal2nal.v14/pal2nal.pl *_AA_MSA_orig_gene_ids.fa *_CDS.fa -codontable 1 -output fasta > ${d}_pal2nal_CDS.fa ); fi; done

### Running FUBAR for all folders 

$ for d in *; do if [ -d "$d" ]; then ( cd "$d" && hyphy fubar --alignment ${d}_pal2nal_CDS.fa --tree ${d}_phyl_tree_orig_gene_ids.newick ); fi; done &> running_FUBAR.txt &

### Parsing FUBAR outputs

myext = phyphy.Extractor("/home/icruz/data/PRUEBA_MEME_fams/PRUEBA_MEME_withnewpipeline/working_dir/MEME_22001_23000/OG0022001/OG0022001_pal2nal_CDS.fa.MEME.json")
myext.extract_csv("fest.tsv", delim = "\t")

### Parsing .json format for all directories

for d in *; do if [ -d "$d" ]; then ( cd "$d" && python /LUSTRE/Genetica/ivan/prueba_parsing_fubar/parsejason_v2.py -f *_pal2nal_CDS.fa.FUBAR.json ); fi; done &> running_parseJSON.txt &

### Now you will have a txt file in each directory or folder with the sites under positive selection, the next step is to parse that file to just keep the sites with a byes factor > 100, significant posterior probabilities ≥ 0.98. Use this script for the purpose. 

* https://github.com/asishallab/GeneFamilies/blob/master/exec/loadFubarResults.R

### Enrichment tests can be done using InterProscan and Mapman4 annotations with this script 

* https://github.com/asishallab/SlydGeneFamsAnalyses/identifyDomainsAtSelectedSites.R
