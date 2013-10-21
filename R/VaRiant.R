#############################################################
# Define a VCF Class
# 

require("stringr")

VCF <- function(vcf_file, ...) UseMethod("VCF")

read.vcf <- function(vcf_file){
    open(con <- file(vcf_file))
    num_headers = 0
    while(str_detect(readLines(con,n=1L),"##"))
        num_headers = num_headers + 1
    close(con)
    read.table(vcf_file,comment.char="",skip=num_headers,header=T)
}

VCF.default <- function(vcf_file, ... ){
    # Set the default method 
    self <- list()
    class(self) <- 'VCF' 
    # read in the VCF file
    tryCatch(
        { 
            self$tbl <- read.vcf(vcf_file)
            # VCF files have 8 mandatory fields:
            self$CHROM <- self$tbl[,1]
            self$POS   <- self$tbl[,2]
            self$ID    <- self$tbl[,3]
            self$REF   <- self$tbl[,4]
            self$ALT   <- self$tbl[,5]
            self$QUAL  <- self$tbl[,6]
            self$FILTER <- self$tbl[,7]
            self$INFO   <- self$tbl[,8]
            # If Genotype data is present, there will be a FORMAT COLUMN
            if("FORMAT" %in% names(self$tbl)){
                self$FORMAT <- self$tbl[,9]
                # Assign the rest as an allele table
                self$INDS <- names(
                    self$tbl[!names(self$tbl) %in% c("X.CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
                )
                self$ALLELES <- Alleles(self$tbl[,10:ncol(self$tbl)])
            }
            # Finally, set some 
            self$VARIDs <- do.call(paste,c(self$tbl[,c('X.CHROM',"POS")],sep='.'))

        },
        error = function(){stop(paste("Could not read input file: ", vcf_file, sep=""))}
    )
    self
}
# IDS - retrieve ids for variants in the VCF 
ids <- function(VCF,rows=c(), ... ) UseMethod("ids")
ids.VCF <- function(VCF, rows=c(), ... ){
    if(length(rows) > 0)
        VCF$ids[rows]
    else
        VCF$ids
}
# Index - retrieve an index by a variant ID 
index <- function(VCF, ids=c(), ... ) UseMethod('index')
index.VCF <- function(VCF, ids=c(), ...){
    which(VCF$ids %in% ids)
}
# Quality Scores - retrieve scores as a vector 
#   These are phred based quality scores which indicate assertions made in the ALT genotype
#   i.e. the -10log_10(prob ALT call is wrong)
scores <- function(x, ...) UseMethod("scores")
scores.VCF <- function(VCF, ... ){
    VCF$tbl[,c('QUAL')]
}
# Alleles - retrieve alleles of specified individuals, optional specificaltion of which rows
Alleles.VCF <- function(VCF, inds, rows, ... ){
    if(is.na(inds)) stop("you must specify individuals")
    Alleles(VCF$tbl[rows,inds])
}

length.VCF <- function(VCF){
    length(VCF$POS)
}

##################################################################
# Allele Table
#
Alleles <- function(x,... ) UseMethod("Alleles")
Alleles.default <- function(tbl, ...){
    self <- list()
    class(self) <- 'Alleles'
    tryCatch(
        {
            self$tbl <- data.frame(lapply(tbl,as.character),stringsAsFactors=FALSE)
        },
        error = function(){stop("Unable to create allele table")}    
    )       

    self
}
length.Alleles <- function(Alleles,...){
    nrow(Alleles$tbl)
}
# homozygous - determine if an allele table is homozygous for a certain genotype (default=0)
homozygous <- function(x,...) UseMethod("homozygous")
homozygous.Alleles <- function(Alleles, genotype = 0 ){
    apply(Alleles,1, function(row){
        all(grepl(paste(genotype,"/",genotype,".*",sep=''),row))
    })
}

