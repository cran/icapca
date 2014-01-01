\name{initializations}
\alias{initializations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes number of initializations to be performed by ica_pca
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For a given number of sub-Gaussian and super-Gaussian sources
and Gaussian components, the function ica_pca will initialize
the model multiple times.  This function will compute the number
of initializations that will be performed.
}
\usage{
initializations(subgaussians, supergaussians, gaussians)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subgaussians}{
     the number of sub-Gaussian sources in the model
}
  \item{supergaussians}{
     the number of super-Gaussian sources in the model
}
  \item{gaussians}{
     the number of Gaussian sources in the model
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If the number of initializations is small (less than 50 to 100),
the ica_pca function may fail to identify the optimal model;
models with small numbers of initializations should be run
several times using different values for seed and/or
offset_random.  As the number of sources and components gets 
large (e.g., with totals more than 10) the number of
initializations grows quickly.  To a first approximation, 
computation time is proportional to the number of 
initializations.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
returns the number of initializations
}
\author{
Roger P. Woods, M.D.
}
\examples{
initializations(4,1,2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
