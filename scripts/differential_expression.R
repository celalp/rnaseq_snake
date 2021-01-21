#Differential expression script for DE pipeline

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))

regulation<-function(df, fccol, sigcol){
  #check the regulation direction based on the log2FC sign
  fc <- df[[fccol]]
  issig <- df[[sigcol]]
  regs <- rep(NA, nrow(df))
  regs[fc > 0 & issig] <- 'Up'
  regs[fc < 0 & issig] <- 'Down'
  return(regs)
}

get_counts<-function(file_path, extension, tool_type){
  #search files in the file_path
  #make into file path for tximport input
  files<-dir(path=file_path, pattern = extension)
  file_list<-file.path(file_path,files)
  
  transcript_abundances  <- tximport(file_list, type = tool_type, txIn = FALSE, txOut = FALSE, abundanceCol = 'TPM', lengthCol = "length")
  
  #get the sample name, usually the prefix of the file
  samplenames<-gsub(extension, "", files)

  #subset the count result and rename the columns with sample name
  counts<-transcript_abundances$counts
  colnames(counts)<-samplenames
  return(counts)
}

do_de<-function(counts, sample_meta, formula, comparison, p_cutoff, log2_cutoff){
  #build DESeq object
  dds<-DESeqDataSetFromMatrix(round(counts), colData = sample_meta, design = formula)
  dds<-DESeq(dds)
  
  res<-as.data.frame(results(dds, contrast = comparison))
  res$ensemblID<-rownames(res)#add ensembleID column
  res<- res %>% select(ensemblID, everything())#re-order the dataframe
  
  #add columns of sig and regulation direction
  res$sig<-ifelse(res$padj<p_cutoff & abs(res$log2FoldChange)>log2_cutoff, T, F)
  res$regulation<-regulation(res, fccol="log2FoldChange", sigcol="sig")
  
  return(res)
}

#########################

option_list = list(
  make_option(c('-s','--samples_table'), type="character", help="input sample info table"),
  make_option(c('-f','--formula'), type="character",  help="formula for DESeq2 condition"),
  make_option(c('-t','--type'), type="character",  help="quantification tool type (e.g rsem, sailfish, salmon)"),
  make_option(c('-l','--folder_path'), type="character", help="path to the quantification output folder"),
  make_option(c('-e','--extension'), type="character", help="quantification output files extension,e,g .gene.results"),
  make_option(c('-p','--p_value'), type="numeric", default = 0.01, help="p value cutoff"),
  make_option(c('-g','--logfc'), type="numeric", default = 1.5, help="log2 FC cutoff"),
  make_option(c("-d", '--directory'), type="character",help="directory to find results in"),
  make_option(c("-c", '--contrast'), type="character", help="comparison in DESeq2")
) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#read the sample meta table
samples_table<-read.table(opt$samples_table, header = T, sep = "\t", stringsAsFactors = F)

results_dir <- opt$directory

#the input comes from python. It is a string connected with '-'
#example: 'condition,KO_35,KO_14-condition,KO_35,WT_14'
comparison_list <- unlist(strsplit(opt$contrast,'-'))#[1] "condition,KO_35,KO_14" "condition,KO_35,WT_14" "condition,KO_35,WT_35"

#build the formula
condition_formula <- as.formula(paste0("~",opt$formula))#"~condition"

#get counts information
gene_counts <- get_counts(file_path=opt$folder_path, extension=opt$extension,tool_type=opt$type)

#in case there are multiple comparison
for (c in comparison_list){
  comparison_vector<-unlist(strsplit(c,','))
  result <- do_de(counts=gene_counts, sample_meta=samples_table, formula=condition_formula, 
        comparison=comparison_vector, p_cutoff=opt$p_value, log2_cutoff=opt$logfc)
  outname_prefix<-paste0(strsplit(c,',')[[1]][2],'vs',strsplit(c,',')[[1]][3])#"KO_35vsKO_14"
  outname <- file.path(results_dir, paste0(outname_prefix, '.txt'))
  write.table(result, outname, col.names = T, row.names = F, quote = F, sep = '\t')
}

print('DE Done')

