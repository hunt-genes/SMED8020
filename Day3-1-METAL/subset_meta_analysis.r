options(stringsAsFactors=F)
require(data.table)
require(stringr)
require(optparse)
option_list <- list(
            make_option("--input", type="character", default="", help="Input file"),
            make_option("--output", type="character", default="", help="Output file"))
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

##read in file
df<-fread(opt$input)

##how many studies was the marker in?
df$num_study_missing<-str_count(df$Direction,"\\?")

##subset to markers in more than 1 study (max for this example is 2)
subset<-df[df$num_study_missing<=1,]

#write output
write.table(subset,opt$output,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")