df=read.table("/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/03-21-16/FastaFile_Ramila/TN_fusion_031916.txt", stringsAsFactors=F, header=T, sep="\t")
str(df)
for(i in 1:nrow(df)){
 id=df[i,"fusionPeptide_ID"]
 print(id);
 sink(paste(id,".fasta",sep=""));
 cat(paste(">",id,sep=""));
 cat("\n");
 cat(df[1,"Peptide"]);
 cat("\n");
 sink();
}

