options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

option_list <- list(
  make_option("--prefix_inp", type="character", default="",
    help="Prefix of input files"),   
  make_option("--prefix_out", type="character", default="",
    help="Prefix of output files"),   
  make_option("--prefix_pheno", type="character", default="",
    help="Prefix of phenotype files"),   
  make_option("--analysis", type="character", default="BOLTLMM",
    help="Analysis software BOLTLMM or SAIGE [default='BOLTLMM]"),        
  make_option("--minMAC", type="numeric", default=25,
    help="minimal minor allele count threshold [default=25]"),
  make_option("--minMAF", type="numeric", default=0,
    help="minimal minor allele frequency [default = 0]"),
  make_option("--hitregion", type="character", default="",
    help="File with candidate regions, CHROM;START;END;COL;LEGENDTEXT [default='']"),
  make_option("--outputClean", type="logical", default=F,
    help="Output MAF filtered and GC corrected results T/F [default=F]"),
  make_option("--format", type="character", default="dose",
    help="Input format of analysis [default='dose']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# inpput files
file_prefixI <- paste0(opt$prefix_inp,"_",opt$analysis,".imputed")
# output files
file_prefixO <- paste0(opt$prefix_out,"_",opt$analysis,".imputed")

# Get some basic information about phenotypes to determine title and MAF threshold 
# read phenotype summary
phenoSummary <- read.table(paste0(opt$prefix_pheno,"_summary.txt"),sep="\t",comment.char="",
	header=T,check.names=F)
traitname <- phenoSummary[1,1]

binary <- ifelse(phenoSummary[1,"Summary_Controls"] == "" || is.na(phenoSummary[1,"Summary_Controls"]),F,T)

# Read input
file_inp <- paste0(file_prefixI,".results.txt")

if(opt$analysis == "BOLTLMM"){
	pcol <- "P_BOLT_LMM_INF"
	ycol <- "log10P_GC"
	chicol <- "CHISQ_BOLT_LMM_INF"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "BP"
	idcol <- "SNP"
	refcol <- "ALLELE0"
	altcol <- "ALLELE1"
	accol <- ""
	maccol <- ""
	afcol <- "A1FREQ"
	af0col <- ""
	af1col <- ""
} else if (opt$analysis == "SAIGE"){
  if (opt$format == "dose") {
	pcol <- "p.value"
	ycol <- "log10P_GC"
	chicol <- "CHISQ"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "POS"
	idcol <- "SNPID"
	refcol <- "Allele1"	
	altcol <- "Allele0"
	accol <- ""
	maccol <- ""
	afcol <- "AF"
	af0col <- ""
	af1col <- ""
  } else if (opt$format == "vcf") {
	pcol <- "p.value"
	ycol <- "log10P_GC"
	chicol <- "CHISQ"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "POS"
	idcol <- "SNPID"
	refcol <- "Allele1"
	altcol <- "Allele2"
	accol <- ""
	maccol <- ""
	afcol <- "AF_Allele2"
	af0col <- ""
	af1col <- ""
  } else if (opt$format == "bgen") {
	pcol <- "p.value"
	ycol <- "log10P_GC"
	chicol <- "CHISQ"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "POS"
	idcol <- "SNPID"
	refcol <- "Allele1"
	altcol <- "Allele2"
	accol <- ""
	maccol <- ""
	afcol <- "AF_Allele2"
	af0col <- ""
	af1col <- ""
  } else {
    stop("Unknown format")
  }
} else {
	stop("Unknown Analysis")
}

# quantitative trait
if(!binary){
    	cases <- as.integer(gsub("N=([0-9]+); Mean=.+","\\1",phenoSummary[1,"Summary_Cases"]))
	maintitle <- paste0("Phenotype '",traitname,"': ",
			format(cases,big.mark=",",scientific=FALSE),
			" Individuals")
	maf_threshold <- max(opt$minMAF, signif(opt$minMAC / (2*cases),2))
# binary trait
} else {
  	cases <- as.integer(gsub("N=","",phenoSummary[1,"Summary_Cases"]))
	controls <- as.integer(gsub("N=","",phenoSummary[1,"Summary_Controls"]))
	maintitle <- paste0("Phenotype '",traitname,"': ",
			format(cases,big.mark=",",scientific=FALSE),
			" cases versus ",
			format(controls,big.mark=",",scientific=FALSE)," controls")
	maf_threshold <- max(opt$minMAF, signif(opt$minMAC / (2*min(cases,controls)),2))
}

# convert -log10(P) values to as.character(P)
log10toP <- function(log10P){
		log10P <- abs(as.numeric(log10P))
		if(is.na(log10P)) return(NA)
		if(log10P > 300){
			part1 <- log10P%/%100*100
			part2 <- log10P-part1
			if(part2 != 0){
				P <- format(signif(10^-part2,3), scientific = T)
				P <- paste(as.numeric(gsub("e-.+","",P)),"e-",as.numeric(gsub(".+-","",P),sep="")+part1,sep="")
			} else {
				P <- paste("1e-",part1,sep="")
			}
		} else {
			P <- signif(10^-log10P,3)
		}
	return(as.character(P))
}

makeFinal <- function() {
    message(paste0("Read file: ", file_inp))

    # Read input file; add -log10 P values
    resIn <- fread(file_inp)

    # Remove invalid results (SE ==0)
    if(length(which(resIn[[secol]] == 0)) > 0){
	resIn <- resIn[which(resIn[[secol]] > 0),]
    }

    resIn$log10P <- -log10(resIn[[pcol]])
    # If p values are too small for R (=0). Get log10 p-values from CHISQ or SCORE statistics
    noP <- which(resIn[[pcol]] == 0)
    # log10(P) from CHISQ
    message("Pvalue")
    if(length(noP)>0 & chicol %in% colnames(resIn)){
	resIn$log10P[noP] <- signif(-pchisq(resIn[[chicol]][noP],df=1,lower.tail=F,log.p = T)/log(10),5)
	resIn[[pcol]][noP] <- sapply(resIn$log10P[noP],log10toP)
    # log10(P) from SCORE statistics
    } else if (length(noP)>0 & !chicol %in% colnames(resIn)){
	    resIn$log10P[noP] <- rep(-log10(as.numeric(c(unlist(format(.Machine)))[3])),length(noP))
	    #sapply(resIn$Tstat[noP],Score2log10P)
	resIn[[pcol]][noP] <- sapply(resIn$log10P[noP],log10toP)    
    }

    # Add overall MAF
    message("MAF")

    resIn$MAF <- ifelse(resIn[[afcol]] <= 0.5,resIn[[afcol]],1-resIn[[afcol]])
    resIn <- resIn[which(resIn$MAF >= maf_threshold),]

    # ADD CHISQ, calculate from p value;  could be faster with a parallel solution
    message("chisq")
    if(!chicol %in% colnames(resIn)){
	    resIn$CHISQ <- qchisq(-resIn$log10P*log(10), df=1, lower.tail = F,log.p=T)
    }

    # Minimal MAF
    #minMAF <- min(resIn$MAF,na.rm=T)

    # Lambda estimation using SNPs with MAF > 0.01 and GC correction if necessary
    message("lambda (maf>=0.01)")
    denom<-qchisq(0.5, df=1)
    Lambda <- qchisq(10^-median(resIn$log10P[which(resIn$MAF >= 0.01)],na.rm=T), df=1, lower.tail = F) / denom
    print(Lambda)
    message("lambda (maf>=0.05)")
    Lambda5pc <- qchisq(10^-median(resIn$log10P[which(resIn$MAF >= 0.05)],na.rm=T), df=1, lower.tail = F) / denom
    print(Lambda5pc)
    if(Lambda > 1){
	    resIn$SEBETA_GC <- signif(resIn[[secol]]*sqrt(Lambda),5)
    	    resIn$CHISQ_GC <- (resIn[[betacol]] / resIn$SEBETA_GC)^2
	    resIn$log10P_GC <- -pchisq(resIn$CHISQ_GC, df=1, lower.tail = F, log.p = T)/log(10)
    } else {
	    resIn$log10P_GC <- resIn$log10P
	    resIn$SEBETA_GC <- resIn[[secol]]
    	    resIn$CHISQ_GC <- resIn[[chicol]]
    }


    ## Clean GWAS output 

    if(opt$analysis == "BOLTLMM"){
	    resOut <- data.frame(
		    'CHROM'=resIn[[chrcol]],
		    'POS'=resIn[[poscol]],
		    'ID'=resIn[[idcol]],
		    'REF'=resIn[[refcol]],
		    'ALT'=resIn[[altcol]],
    #		'AC'=NA,
    #		'MAC'=NA,
		    'MAF'=resIn$MAF,
		    'AF'=resIn$A1FREQ,
    #		'AF.CTRL'=NA,
    #		'AF.CASE'=NA,
		    'CHISQ'=resIn[[chicol]],
		    'BETA'=resIn[[betacol]], #Approxmation
		    'SEBETA_GC'=resIn$SEBETA_GC, #Approxmation
		    'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		    'LOG10P_GC'=resIn[[ycol]])
    } else if (opt$analysis == "SAIGE"){
      if(binary){
        if ("AF.Controls" %in% colnames(resIn)) {
	    resOut <- data.frame(
		    'CHROM'=resIn[[chrcol]],
		    'POS'=resIn[[poscol]],
		    'ID'=resIn[[idcol]],
		    'REF'=resIn[[refcol]],
		    'ALT'=resIn[[altcol]],
		    'AC'=resIn$AC,
		    'MAC'=ifelse(resIn$AC > resIn$N,resIn$N*2-resIn$AC,resIn$AC),
		    'MAF'=resIn$MAF,
		    'AF'=resIn[[afcol]],
		    'AF.CTRL'=resIn$AF.Controls,
		    'AF.CASE'=resIn$AF.Cases,
		    'CHISQ'=resIn[[chicol]],
		    'BETA'=resIn[[betacol]],
		    'SEBETA_GC'=resIn$SEBETA_GC,
		    'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		    'LOG10P_GC'=resIn[[ycol]])
	} else {
	    resOut <- data.frame(
		    'CHROM'=resIn[[chrcol]],
		    'POS'=resIn[[poscol]],
		    'ID'=resIn[[idcol]],
		    'REF'=resIn[[refcol]],
		    'ALT'=resIn[[altcol]],
		    'AC'=resIn$AC,
		    'MAC'=ifelse(resIn$AC > resIn$N,resIn$N*2-resIn$AC,resIn$AC),
		    'MAF'=resIn$MAF,
		    'AF'=resIn[[afcol]],
		    'CHISQ'=resIn[[chicol]],
		    'BETA'=resIn[[betacol]],
		    'SEBETA_GC'=resIn$SEBETA_GC,
		    'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		    'LOG10P_GC'=resIn[[ycol]])
       }
      } else {
	    resOut <- data.frame(
		    'CHROM'=resIn[[chrcol]],
		    'POS'=resIn[[poscol]],
		    'ID'=resIn[[idcol]],
		    'REF'=resIn[[refcol]],
		    'ALT'=resIn[[altcol]],
		    'AC'=resIn$AC,
		    'MAC'=ifelse(resIn$AC > resIn$N,resIn$N*2-resIn$AC,resIn$AC),
		    'MAF'=resIn$MAF,
		    'AF'=resIn[[afcol]],
		    #'AF.CTRL'=resIn$AF.Controls,
		    #'AF.CASE'=resIn$AF.Cases,
		    'CHISQ'=resIn[[chicol]],
		    'BETA'=resIn[[betacol]],
		    'SEBETA_GC'=resIn$SEBETA_GC,
		    'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		    'LOG10P_GC'=resIn[[ycol]])
      }
    }
    write.table(resOut,paste0(file_prefixO,".finalResults.txt"),
    	    sep="\t",col.names=T,row.names=F,quote=F)
    resOutSuggT = resOut[which(resOut$LOG10P_GC >= 5),]
    resOutSuggT = resOutSuggT[order(resOutSuggT$LOG10P_GC, decreasing=T),]
    write.table(resOutSuggT,paste0(file_prefixO,".finalResults_suggesT.txt"),
	    sep="\t",col.names=T,row.names=F,quote=F)
    res = list(resIn=resIn, Lambda = Lambda)
    return(res)
}

print(paste0("maf_threshold: ", maf_threshold))

res = makeFinal()
resIn = res$resIn
Lambda = round(res$Lambda,5)
maintitle = paste0(maintitle,", Lambda=", Lambda, ", ", opt$analysis)
if (opt$minMAF != 0) maintitle  = paste0(maintitle, ", MAF>",opt$minMAF)

# QQ plot Rscript downloaded from here: https://github.com/ilarsf/gwasTools
cmd = paste(
	"Rscript --vanilla QQplot.r",
	"--input",paste0(file_prefixO,".finalResults.txt"),
	"--prefix",file_prefixO,
	"--maf","MAF",
	"--pvalue","LOG10P_GC",
	"--log10p","T",
	"--maintitle", paste0("\'",gsub("[\"\']","`",maintitle),"\'")
	)
system(cmd)

# Settings for plot
gwthreshold <- 5E-8

# Determine candidate regions if not provided
candidateHits <- which(resIn$log10P_GC > -log10(gwthreshold))
# BOLTLMM
if (length(candidateHits) >0) {
  if(opt$analysis =="BOLTLMM"){
	tophits <- resIn[candidateHits,]

	# merge candidate regions if they are within 0.025 cM
	tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
	x <- as.numeric(tophits$GENPOS)
	y <- tophits$numCHR
	start = c(1, which(diff(y) != 0 | diff(x) <= -0.025 | diff(x) >= 0.025) + 1)
	end = c(start - 1, length(x))
	candidateRegions <- data.frame(
		'CHROM'=tophits$CHR[start],
		'START'=tophits$BP[start] - 50000,
		'END'=tophits$BP[end] + 50000,
		'COL'="blue",
		'LEGENDTEXT'="Candidate Regions [Top Hits ±(0.025cM + 50kb)]"
		)
  # SAIGE	
  } else if (opt$analysis =="SAIGE"){
	tophits <- resIn[candidateHits,]
	tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
	x <- as.numeric(tophits$POS)
	y <- tophits$numCHR

	start = c(1, which(diff(y) != 0 | diff(x) <= -500000 | diff(x) >= 500000) + 1)
	end = c(start - 1, length(x))
	candidateRegions <- data.frame(
		'CHROM'=tophits$CHR[start],
		'START'=tophits$POS[start] - 500000,
		'END'=tophits$POS[end] + 500000,
		'COL'="blue",
		'LEGENDTEXT'="Candidate Regions [±500kb]"
		)
  }
} else {
	candidateRegions <- data.frame(
		'CHROM'=character(0),
		'START'=numeric(0),
		'END'=numeric(0),
		'COL'=character(0),
		'LEGENDTEXT'=character(0)
		)
}

manPrefix = file_prefixO
if (opt$hitregion != ""){
  knownRegions <- read.table(opt$hitregion, sep='\t', header=T)
  candidateRegions = rbind(knownRegions, candidateRegions)
  manPrefix = paste0(file_prefixO,'.known')
}

write.table(candidateRegions[c("CHROM","START","END","COL","LEGENDTEXT")],
	paste0(manPrefix,".hitregions.txt"),sep="\t",col.names=T,row.names=F,quote=F)

thinratio=1
# Manhattan plot Rscript downloaded from here: https://github.com/ilarsf/gwasTools
cmd = paste(
	"Rscript --vanilla ManhattanPlot.r",
	"--input",paste0(file_prefixO,".finalResults.txt"),
	"--prefix",manPrefix,
	#"--hitregion",ifelse(dim(candidateRegions)[1] > 0,paste0(file_prefixO,".hitregions.txt"),""),
	"--hitregion",paste0(manPrefix,".hitregions.txt"),
	"--chr","CHROM",
	"--pos","POS",
	"--pvalue","LOG10P_GC",
	"--log10p","T",
        "--thinratio",thinratio,
	"--maintitle", paste0("\'",gsub("[\"\']","`",maintitle),"\'")
	)

system(cmd)
