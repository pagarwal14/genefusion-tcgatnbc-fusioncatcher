library(stringr)
#
# Get the cell line name from cmd line
#
cmdFlag = T # change to F if hard coding cell line name in code.
if(cmdFlag==T){
args = commandArgs(trailingOnly=TRUE)
numargs = length(args)
print(paste("Number of cmd line args",numargs))
if(length(args)==0) {
stop("At least one argument must be supplied (input file).tsv", call.=FALSE)
}else if (length(args)==1) {
cellLineName = args[1]
}else{
stop("Only one argument must be supplied (Cell Line Name)", call.=FALSE)
}
#
}else{#default cell line name
cellLineName = "MDA231"
}
print(paste("cellLineName", cellLineName))
#
#
if(cellLineName=="MDA231"){
	cellLinePepHitFileName = "MDA231-2.txt";
} else if(cellLineName=="MCF7"){
	cellLinePepHitFileName = "MCF7.txt";
} else if(cellLineName=="MCF10"){
        cellLinePepHitFileName = "MCF10.txt";
} else if(cellLineName=="Xenograft"){
        cellLinePepHitFileName = "Xenograft.txt";
}else{
stop("Incorrect cell line name");
}
#
#results dir
resultsdir = "/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/peptideMS/results/"
#
dfAllFusion = read.table("/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/peptideMS/data/TN_fusion_032116.txt", sep="\t", stringsAsFactors=F, header=T)
str(dfAllFusion)
dfCellLine = read.table(paste("/srv/agarw005/projects/GeneFusion/analysis/TCGA_TNBC/rna-seq/147files/FusionCatcher_10-09-15/analyze_fusioncatcher_prod_10-18-15/pipeline/genefusion-tcgatnbc-fusioncatcher/peptideMS/data/pepsearch040516/", cellLinePepHitFileName, sep=""), sep="\t", stringsAsFactors=F, header=T)
str(dfCellLine)
dfCellLine = dfCellLine[,c("Protein.Group.Accessions", "Sequence")]
str(dfCellLine)
npeps = nrow(dfCellLine)
npeps
#
#
getPep = function(fusionid){
 pepstr = dfAllFusion[dfAllFusion$fusionPeptide_ID==fusionid,c("Peptide")]; # this works
 #pepstr = dfAllFusion["fusionPeptide_ID" == fusionid, c("Peptide")]; # this is not correct
 return(pepstr); # return value must be within brackets
}
##
##
dfmaster = data.frame(CellLine = character(), FusionID = character(), HitSequence = character(), FullPepSequence = character(), StartPos = numeric(), EndPos = numeric(), stringsAsFactors=F) #create empty data frame with col header name and data type specify - GOOD
str(dfmaster)
##
##
for(i in 1:npeps){
  #  print(dfCellLine[i,1]);
  #  print(dfCellLine[i,"Sequence"]);
  print("START FOR LOOP NEW ITERATION"); 
  ids = dfCellLine[i,1];
  #print(ids)
  #Convert peptide sequence that was reported as MS hit to upper case in case there are lower case letters ("m" is lower case at least)
  seq = toupper(dfCellLine[i,"Sequence"]);
  qrystrlength = nchar(seq);
  print(paste("qrystrlength", qrystrlength));
  locssep = str_locate_all(pattern=";", ids);
  print(length(locssep[[1]][,1]));
  # if there are more than one rows, means there are more than one id in this string
  lenrowids = length(locssep[[1]][,1])
  print(paste("lenrowids=",lenrowids));
  if( lenrowids > 0 ){
   print(ids);
   print(seq);
   row=c(ids, seq);
   print(row);	
   #print(paste("row=",row));	
   splitids = unlist(strsplit(ids, ";"));#split on token - gives list - unlist to vector
   numids = length(splitids);
   for(j in 1:numids){
     #print(j);
     ids = splitids[j];
     ##
     querystring = seq;
     targetstring = getPep(ids);
     qryloctarget = str_locate_all(pattern=seq, targetstring);
     qryloctargetStartPos = unname(qryloctarget[[1]][1,1]);
     qryloctargetEndPos = unname(qryloctarget[[1]][1,2]);
     lenQueryString = nchar(querystring);
     lenQueryStringInTarget = (qryloctargetEndPos - qryloctargetStartPos)+1;
     ##
     ## 
     #dfmaster = rbind(dfmaster, data.frame(splitids[j], seq, stringsAsFactors=F));# can also specify the header col ("FusionID" = ids, "HitSequence" = seq)
     ###dfmaster = rbind(dfmaster, data.frame(ids, seq, stringsAsFactors=F));# can also specify the header col ("FusionID" = ids, "HitSequence" = seq)
     dfmaster = rbind(dfmaster, data.frame(CellLine=cellLineName, FusionID=ids, HitSequence=seq, FullPepSequence=targetstring, StartPos=qryloctargetStartPos, EndPos=qryloctargetEndPos, stringsAsFactors=F));# can also specify the header col ("FusionID" = ids, "HitSequence" = seq)
   }
 } else {
   #row=c(ids, seq);
   print("\nFUSION EVENT\n");
   print(ids);
   querystring = seq;
   print("\nQUERY\n");
   print(querystring);
   print("\nTARGET\n");
   targetstring = getPep(ids);
   print(targetstring);
   qryloctarget = str_locate_all(pattern=seq, targetstring);
   qryloctargetStartPos = unname(qryloctarget[[1]][1,1]);
   #print(qryloctargetStartPos);
   qryloctargetEndPos = unname(qryloctarget[[1]][1,2]);
   #print(qryloctargetEndPos);
   lenQueryString = nchar(querystring);
   #print(paste("*",lenQueryString));
   lenQueryStringTarget = qryloctargetEndPos - qryloctargetStartPos;
   #print(paste("**",lenQueryStringTarget));
     ##
     ##
	print("####");
	print(ids);
	print(seq);
	print(targetstring);
	print(qryloctargetStartPos);
	print(qryloctargetEndPos);
	print("####");
     ## 
   dfmaster = rbind(dfmaster,data.frame(CellLine=cellLineName, FusionID=ids, HitSequence=seq, FullPepSequence=targetstring, StartPos=qryloctargetStartPos, EndPos=qryloctargetEndPos, stringsAsFactors=F));# can also specify the header col ("FusionID" = ids, "HitSequence" = seq)
 }
  print("END FOR LOOP ITERATION");
}
str(dfmaster)
write.table(dfmaster, paste(resultsdir, cellLineName, "_Hits_FusionEvents.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
