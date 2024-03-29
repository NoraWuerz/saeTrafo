% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlme_methods.R
\name{getResponse}
\alias{getResponse}
\alias{getResponse.NER}
\title{Extract response variable from an saeTrafo object}
\usage{
\method{getResponse}{NER}(object, ...)
}
\arguments{
\item{object}{an object of type "NER".}

\item{...}{additional arguments that are not used in this method.}
}
\value{
Vector containing the response variable.
}
\description{
Method \code{getResponse.NER} extracts the response variable from a saeTrafo
object.
}
\details{
The generic function \code{getResponse} is imported from package
\pkg{nlme} and re-exported to make the S3-methods available, even though
the \pkg{nlme} package itself is not loaded or attached. For default
documentation, see \code{\link[nlme]{getResponse}}.
}
\examples{
# Example to extract the response variable

# Load Data
data("eusilcA_smp")
data("pop_area_size")
data("pop_mean")
data("pop_cov")

# Nested error regression model
NER_model <- NER_Trafo(fixed = eqIncome ~ gender + eqsize + cash +
                       self_empl + unempl_ben + age_ben + surv_ben +
                       sick_ben + dis_ben + rent + fam_allow + house_allow +
                       cap_inv + tax_adj,
                       smp_domains = "district",
                       pop_area_size = pop_area_size,
                       pop_mean = pop_mean, pop_cov = pop_cov,
                       smp_data = eusilcA_smp)

getResponse(NER_model)
}
\seealso{
\code{\link{NER_Trafo}}, \code{\link[nlme]{getResponse}}
}
