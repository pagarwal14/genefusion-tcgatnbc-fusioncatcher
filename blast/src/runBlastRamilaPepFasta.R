library(stringr)
peptides.fasta.dir = "/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/data/input/blast/";
peptides.fasta.files = dir(peptides.fasta.dir, pattern = ".fasta");
output.dir="/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/data/results/blast/nr/";
numFiles = length(peptides.fasta.files);
print(paste("Number of fasta files =", numFiles));
#[1] "Number of fasta files = 81"
for(i in 1:numFiles) {
 in.file = paste(peptides.fasta.dir, peptides.fasta.files[i],sep="");
 #print(in.file);
 in.file.basename = basename(in.file);
 #print(in.file.basename);
 locFileExt = str_locate(in.file.basename, '.fasta');
 fusionEvent = substr(in.file.basename, 0, locFileExt[1,1]-1);
 #print(fusionEvent);
 blastcmd = paste( "blastp -query", in.file, "-db nr -out", paste( output.dir, fusionEvent, "_blastp_nr.out", " -remote", sep="" ) );
 #print(blastcmd);
 results = system( blastcmd );
 print(results);
}

