#Differential expression script for DE pipeline

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(yaml))

regulation<-function(df, fccol, sigcol){
  #check the regulation direction based on the log2FC sign
  fc <- df[[fccol]]
  issig <- df[[sigcol]]
  regs <- rep(NA, nrow(df))
  regs[fc > 0 & issig] <- 'Up'
  regs[fc < 0 & issig] <- 'Down'
  return(regs)
}

make_tximport_parameters_list <- function(default_list, config_tool_setting){
  #create the parameter list for tximport
  
  #get the tool type from the config file
  tool_type <- config_tool_setting[['type']]
  
  p_list <- list()
  if (tool_type %in% names(default_list)){
    #TRUE -> the tool has set default setting
    #check if we are using config parameter or default
    for (i in names(default_list[[tool_type]])){
      if (is.null(config_tool_setting[[i]])){
        #config doesnt have it, use the default
        p_list[[i]] <- default_list[[tool_type]][[i]]
      }
      else{
        #config has it
        p_list[[i]] <- config_tool_setting[[i]]
      }
    }
  }else{
    #use the basic setting combining with the config setting
    p_list <- c(default_list[['basic']],config_tool_setting)
  }
  
  return(p_list)
}

get_counts<-function(file_path, extension, parameter_list, sample_meta, formula){
  #return a DeSeq object for next step comparison
  
  #search files in the file_path
  #make into file path for tximport input
  files<-dir(path=file_path, pattern = extension)
  file_list<-file.path(file_path,files)
  parameter_list[['files']] <- file_list
  
  transcript_abundances  <- do.call(tximport, parameter_list)
  
  #get the sample name, usually the prefix of the file
  samplenames<-gsub(extension, "", files)
  
  #subset the count result and rename the columns with sample name
  counts<-transcript_abundances$counts
  colnames(counts)<-samplenames
  
  #build DESeq object
  dds<-DESeqDataSetFromMatrix(round(counts), colData = sample_meta, design = formula)
  dds<-DESeq(dds)
  
  return(dds)
}

do_de<-function(dds, comparison, p_cutoff, log2_cutoff){
  
  res<-as.data.frame(results(dds, contrast = comparison))
  res$ensemblID<-rownames(res)#add ensembleID column
  res<- res %>% select(ensemblID, everything())#re-order the dataframe
  
  #add columns of sig and regulation direction
  res$sig<-ifelse(res$padj<p_cutoff & abs(res$log2FoldChange)>log2_cutoff, T, F)
  res$regulation<-regulation(res, fccol="log2FoldChange", sigcol="sig")
  
  return(res)
}

#########################

#some default parameters
#if the tool input is not listed in the parameter list, will use the basic setting and the config setting
tximport_default_parameters<-list(rsem=list(type="rsem", abundanceCol="TPM",lengthCol = "length",
                                            txIn = FALSE, txOut = FALSE), basic=list(txIn = FALSE, txOut = FALSE))
#default pvalue and log2fc cutoff
pvalue<-0.01
log2fc<-1.5

option_list = list(
  make_option(c('-s','--samples_table'), type="character", help="input sample info table"),
  make_option(c('-y','--yaml_config'), type="character", help="path to the DE yaml config file"),
  make_option(c('-l','--folder_path'), type="character", help="path to the quantification output folder"),
  make_option(c("-d", '--directory'), type="character",help="directory to find results in")
) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#load yaml file
de_config <- yaml.load_file(opt$yaml_config)

#read the sample meta table
samples_table<-read.table(opt$samples_table, header = T, sep = "\t", stringsAsFactors = F)

#set and create the output directory /result/diff_ex
results_dir <- opt$directory
dir.create(results_dir,mode="0770")

#get the list of comparisons
compare_condition <- names(de_config[['contrast']])
comparison_list<- lapply(de_config[['contrast']][[compare_condition]], function(x) {c(compare_condition,unlist(strsplit(x, ',')))})
#[[1]]
#[1] "condition" "KO_35"     "WT_35"
#[[2]]
#[1] "condition" "KO_21"     "WT_21"    

#build the formula
condition_formula <- as.formula(paste0("~",de_config[['formula']]))#"~condition"

#set up tximport parameters
para_list<-make_tximport_parameters_list(tximport_default_parameters, de_config[["tool_setting"]])

#set up pvalue and log2fc cutoff
if (!is.null(de_config[['pvalue_cutoff']])){
  pvalue<-de_config[['pvalue_cutoff']]
}

if (!is.null(de_config[['log2fc_cutoff']])){
  log2fc<-de_config[['log2fc_cutoff']]
}


#get counts information and create DEseq object
result_dds<-get_counts(file_path=opt$folder_path,extension=de_config[['quantification_output_extension']],
                     parameter_list=para_list,sample_meta=samples_table, formula=condition_formula)

#in case there are multiple comparison
for (c in comparison_list){
  result <- do_de(test_dds,comparison=c, p_cutoff=pvalue, log2_cutoff=log2fc)
  outname_prefix<-paste0(c[2],'vs',c[3])#"KO_35vsKO_14"
  outname <- file.path(results_dir, paste0(outname_prefix, '.txt'))
  write.table(result, outname, col.names = T, row.names = F, quote = F, sep = '\t')
}

print('DE Done')


