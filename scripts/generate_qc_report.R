suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-s", "--samplename"), type="character", help="samplename"),
    make_option(c("-r", '--rmarkdown'), type="character",  help="path to R markdown file"),
      make_option(c('--analysis_path'), type="character", help="path to analysis output directory")
      ) 

      opt_parser <- OptionParser(option_list=option_list)
      opt <- parse_args(opt_parser)

      samplename <- opt$samplename
      out_dir <- opt$analysis_path
      out_file <- paste0(samplename, ".intermediate_qc_report.html")
      db_file <- file.path(out_dir, paste0(samplename, ".intermediate.db"))

      rmarkdown::render(opt$rmarkdown, 
                        params = list(samplename = samplename, 
			                                results_db_file = db_file),
							                  output_dir = out_dir,
									                    output_file = out_file)
