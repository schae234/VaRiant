\name{VCF.default}
\alias{VCF.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VCF.default(vcf_file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf_file}{
%%     ~~Describe \code{vcf_file} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vcf_file, ...) 
{
    self <- list()
    class(self) <- "VCF"
    tryCatch({
        self$tbl <- read.table(vcf_file, header = T)
        self$CHROM <- self$tbl[, 1]
        self$POS <- self$tbl[, 2]
        self$ID <- self$tbl[, 3]
        self$REF <- self$tbl[, 4]
        self$ALT <- self$tbl[, 5]
        self$QUAL <- self$tbl[, 6]
        self$FILTER <- self$tbl[, 7]
        self$INFO <- self$tbl[, 8]
        if ("FORMAT" \%in\% names(self$tbl)) {
            self$FORMAT <- self$tbl[, 9]
            self$INDS <- names(self$tbl[!names(self$tbl) \%in\% 
                c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                  "FILTER", "INFO", "FORMAT")])
            self$ALLELES <- Alleles(self$tbl[, 10:ncol(self$tbl)])
        }
        self$VARIDs <- do.call(paste, c(self$tbl[, c("CHROM", 
            "POS")], sep = "."))
    }, error = function() {
        stop(paste("Could not read input file: ", vcf_file, sep = ""))
    })
    self
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
