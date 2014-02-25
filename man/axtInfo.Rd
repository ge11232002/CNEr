\name{axtInfo}
\alias{axtInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  axtInfo function
}
\description{
  Given the path of axt file, retrieve the alignments' withs information.
}
\usage{
axtInfo(axtFiles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{axtFiles}{
  Object of \code{character}. The length can be one or more.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A vector of \code{integer} is returned. 
  It stores the withds of all the alignments.
}

\author{
  Ge Tan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{link{readAxt}}
}
\examples{
  axtFilesHg19DanRer7 = 
    file.path(system.file("extdata", package="CNEr"), 
              "hg19.danRer7.net.axt")
    ans = axtInfo(axtFilesHg19DanRer7)

}
