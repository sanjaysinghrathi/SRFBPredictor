\name{SRFBPredictor}
\alias{SRFBPredictor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SRFBPredictor
}
\description{
SRFBPredictor is a package to predict the reponse of a colorectal cancer patient given Chemoradiotherapy.
}
\usage{
SRFBPredictor(training_data_exp_matrix, training_data_phenodata,
              testing_data_exp_matrix, testing_data_phenodata,
              biomarker)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training_data_exp_matrix}{Expression matrix for training dataset having probes as rows and samples as columns}
  \item{training_data_phenodata}{Phenodata (clinical) data for training dataset having samples as rows and various clinical attributes as columns. Response should be first column followed by other important clinical attributes}
  \item{testing_data_exp_matrix}{Expression matrix for test dataset having probes as rows and samples as columns}
  \item{testing_data_phenodata}{Phenodata (clinical) data for test dataset having samples as rows and various clinical attributes as columns. Response should be first column followed by other important clinical attributes. Columns for test phenodata should be same as training dataset}
  \item{biomarker}{Signature to be test.}
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
load("SRFBPredictor.RData")

training_em <- expression_data[,1:100]
training_pData <- pData[1:100,1:3]

testing_em <- expression_data[,101:123]
testing_pData <- pData[101:123, 1:3]

biomarker <- as.character(signature[1:50])

SRFBPredictor(training_data_exp_matrix = training_em, training_data_phenodata = training_pData,
              testing_data_exp_matrix = testing_em, testing_data_phenodata = testing_pData,
              biomarker = biomarker)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
