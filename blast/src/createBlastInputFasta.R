#
# First written and executed on 03-28-16
# 04-17-16: editing to change to new location for blast in the GeneFusion project
#
# Location of the input file which contains the GeneFusion event and associated peptide sequences
inputpepdir = "/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/03-21-16/FastaFile_Ramila/"
inputpepfile = "TN_fusion_031916.txt"
# Location of generated fasta file that will be used for blast
fastadir = "/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/blast/data/input/pepfasta/" 
#
print(paste("input dir for genefusion events and peptides:", inputpepdir))
print(inputpepfile)
df=read.table(paste(inputpepdir, inputpepfile, sep=""), stringsAsFactors=F, header=T, sep="\t")
str(df)
#set location for output fasta files to be written to
print(paste("dir for writing out fasta files", fastadir))
setwd(fastadir)
#
# generate fasta files, one for each peptide sequence - one fusion event could have multiple peptides
#
for(i in 1:nrow(df)){
 id=df[i,"fusionPeptide_ID"]
 print(id);
 sink(paste(id,".fasta",sep=""));
 cat(paste(">",id,sep=""));
 cat("\n");
 cat(df[i,"Peptide"]);
 cat("\n");
 sink();
}
