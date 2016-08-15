\name{N50}
\alias{N50}
\alias{N90}

\title{
  Assembly statistics.
}
\description{
  Calculate the N50, N90 values for a fasta or 2bit file.
}
\usage{
N50(fn)
N90(fn)
}

\arguments{
  \item{fn}{
    \code{character}(1): The path of a fasta or 2bit file.
  }
}
\details{
  This function calculates the N50, N90 values for an assembly.
  The N50 value is calculated by first ordering every contig/scaffold by length
  from longest to shortest.
  Next, starting from the longest contig/scaffold, the lengths of each contig
  are summed, until this running sum equals one-half of the total length of 
  all contigs/scaffolds in the assembly.
  Then the length of shortest contig/scaffold in this list is the N50 value.
  Similar procedure is used for N90 but including 90\% of the assembly.
}
\value{
  An integer value of N50 or N90 value.
}

\author{
  Ge Tan
}
\examples{
  twoBitFn <- file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit")
  N50(twoBitFn)
}