\name{GRangePairs-class}
\docType{class}
% Class:
\alias{class:GRangePairs}
\alias{GRangePairs-class}
% Constructors:
\alias{GRangePairs}

% Accessors:
\alias{names,GRangePairs-method}
\alias{names<-,GRangePairs-method}
\alias{length,GRangePairs-method}
\alias{first}
\alias{first,GRangePairs-method}
\alias{last}
\alias{last,GRangePairs-method}
\alias{seqnames,GRangePairs-method}
\alias{strand,GRangePairs-method}

\title{
GRangePairs objects
}

\description{
  The \code{GRangePairs} class is a container for a pair of 
  \code{GRanges} object that have same lengths.
}

\details{
  A \code{GRangePairs} object is a list-like object where each element
  describes a pair of genomic range.
  They do not necessarily have the same \code{seqinfo}, \emph{i.e.},
  the coordinates from the same assembly.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GRangePairs(first=GRanges(), last=GRanges(), names=NULL)}:
        GRangePairs constructor.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GRangePairs object.
  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of granges pairs in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names on \code{x}.
    }
    \item{}{
      \code{first(x)},
      \code{last(x)}:
      Get the "first" or "last" GRange for each grange pair in
      \code{x}. The result is a \link{GRanges} object of the same length
      as \code{x}.
    }
    \item{}{
      \code{seqnames(x)}:
      Get the seqname of first GRanges and last GRanges and return in a 
      \code{DataFrame} object.
    }
    \item{}{
      \code{strand(x)}:
      Get the strand for each grange pair in \code{x}.
    }
    }
  }



\author{
  Ge Tan
}

\seealso{
  \code{\link{Axt}}
}
\examples{
  library(GenomicRanges)
  first <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr3"),
                   ranges=IRanges(start=c(1, 20, 2, 3),
                                  end=c(10, 25, 10, 10)),
                   strand="+")
  last <- GRanges(seqnames=c("chr1", "chr10", "chr10", "chr20"),
                  ranges=IRanges(start=c(1, 25, 50, 5),
                                 end=c(8, 40, 55, 16)),
                  strand="+")
  namesGRangePairs <- c("a","b","c","d")
  grangesPairs1 <- GRangePairs(first, last, names=namesGRangePairs)
  grangesPairs2 <- GRangePairs(first, last)
  
  ## getters
  names(grangesPairs1)
  length(grangesPairs1)
  first(grangesPairs1)
  last(grangesPairs1)
  seqnames(grangesPairs1)
  strand(grangesPairs1)
  
  ## setters
  names(grangesPairs2) <- namesGRangePairs
  
}
