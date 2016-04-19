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
\alias{seqinfo,GRangePairs-method}

% List methods:
\alias{[[,GRangePairs,ANY,ANY-method}
\alias{unlist,GRangePairs-method}

% Coercion:
\alias{grglist,GRangePairs-method}
\alias{coerce,GRangePairs,GRangesList-method}
\alias{coerce,GRangePairs,GRanges-method}

% Combining:
\alias{c,GRangePairs-method}

% Swaping:
\alias{swap}
\alias{swap,GRangePairs-method}

% "show" method:
\alias{show,GRangePairs-method}

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
    \item{}{
      \code{seqinfo(x)}:
       Get the information about the underlying sequences.
    }
    }
  }


\section{Vector methods}{
  In the code snippets below, \code{x} is a GRangePairs object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new GRangePairs object made of the selected
      genomic ranges pairs.
    }
  }
}

\section{List methods}{
  In the code snippets below, \code{x} is a GRangePairs object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the i-th alignment pair as a \link{GRangePairs} object
      of length 2. As expected \code{x[[i]][1]} and \code{x[[i]][2]} are
      respectively the "first" and "last" granges in the pair.
    }
    \item{}{
      \code{unlist(x, use.names=TRUE)}:
      Return the \link{GRangePairs} object conceptually defined
      by \code{c(x[[1]], x[[2]], ..., x[[length(x)]])}.
      \code{use.names} determines whether \code{x} names should be
      propagated to the result or not.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRangePairs object.

  \describe{
    \item{}{
      \code{grglist(x, use.mcols=FALSE)}:

      Return a \link{GRangesList} object of length \code{length(x)}
      where the i-th element represents the ranges (with respect to the
      reference) of the i-th grange pair in \code{x}.

      Note that this results in the ranges being \emph{always} ordered
      consistently with the original "query template", that is, being in the
      order defined by walking the "query template" from the beginning to
      the end.

      If \code{use.mcols} is TRUE and \code{x} has metadata columns on it
      (accessible with \code{mcols(x)}), they're propagated to the returned
      object.
    }
    \item{}{
      \code{as(x, "GRangesList")}:
      Alternate ways of doing \code{grglist(x, use.mcols=TRUE)}.
    }
    \item{}{
      \code{as(x, "GRanges")}:
      Equivalent of \code{unlist(x, use.names=TRUE)}.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{swap(x)}:
        Swap the first, last GRanges.
    }
    \item{}{
      \code{show(x)}:
      By default the \code{show} method displays 5 head and 5 tail
      elements. This can be changed by setting the global options
      \code{showHeadLines} and \code{showTailLines}. If the object
      length is less than (or equal to) the sum of these 2 options
      plus 1, then the full object is displayed.
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
  ## Constructor
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
  seqinfo(grangesPairs1)
  
  ## setters
  names(grangesPairs2) <- namesGRangePairs
  
  ## Vector methods
  grangesPairs1[1]
  
  ## List methods
  grangesPairs1[[1]]
  unlist(grangesPairs1)
  
  ## Coersion
  grglist(grangesPairs1)
  as(grangesPairs1, "GRangesList")
  as(grangesPairs1, "GRanges")
  as(grangesPairs1, "DataFrame")
  as.data.frame(grangesPairs1)
  
  ## Combining
  c(grangesPairs1, grangesPairs2)
}
