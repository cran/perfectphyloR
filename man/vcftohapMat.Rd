% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_to_hapMat.R
\name{vcftohapMat}
\alias{vcftohapMat}
\title{Create a \code{hapMat} object from variant call format (\code{vcf}) file.}
\usage{
vcftohapMat(vcf_file_path)
}
\arguments{
\item{vcf_file_path}{File path to the \code{vcf} file.}
}
\value{
An object of class \code{hapMat}.
}
\description{
This function creates a \code{hapMat} object from variant call format (\code{vcf}) file.
}
\examples{
\dontrun{
# Specify the file path.
vcf_file_path <- "C:/vcfData/vcfData.vcf.gz"
# Create a hapMat object from the vcf file.
ex_vcf_hapMat <- vcftohapMat(vcf_file_path) 
}  

}
