% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stride_template.R
\docType{data}
\name{stride_template}
\alias{stride_template}
\title{Walking Stride Pattern Templates}
\format{A \code{list} with four named elements:
\itemize{
  \item \code{left_wrist},
  \item \code{left_hip},
  \item \code{left_ankle},
  \item \code{right_ankle}.
}
Each of the above is a five-element \code{list} of \code{matrix} objects.
The matrices are collection of (sub)population-specific stride pattern templates.
For example,
\itemize{
  \item \code{stride_template$left_wrist[[1]]} is a \code{1 x 200} \code{matrix}
  of one population-specific stride template derived from accelerometry data
  collected at left wrist.
  \item \code{stride_template$left_wrist[[2]]} is a \code{2 x 200} \code{matrix}
  of two distinct subpopulation-specific stride templates derived from accelerometry data
  collected at left wrist. Each row is a one subpopulation-specific stride template.
  \item \code{stride_template$right_ankle[[5]]} is a \code{5 x 200} \code{matrix}
  of five distinct subpopulation-specific stride templates derived from accelerometry data
  collected at right ankle.
}}
\usage{
stride_template
}
\description{
Walking stride pattern templates derived from accelerometry
data collected at four body locations: left wrist, left hip, left ankle,
and right ankle.
}
\details{
Raw accelerometry data used to derive walking stride pattern templates
were collected as a part of the study on Identification of Walking,
Stair Climbing, and Driving Using Wearable Accelerometers, sponsored by the Indiana
University CTSI grant and conducted at the Department of Biostatistics,
Fairbanks School of Public Health at Indiana University. The study was led by Dr.
Jaroslaw Harezlak, assisted by Drs. William Fadel and Jacek Urbanek.
It enrolled 32 healthy participants between 23 and 52 years of age.
Participants were asked, among others, to perform self-paced, undisturbed,
outdoor walking on the sidewalk. Accelerometry data were collected
at four body locations: left wrist, left hip, left ankle,
and right ankle.

To derive empirical stride pattern, firstly, from each body location,
642 data segments corresponding to individual walking strides
were manually segmented. Secondly, Vector Magnitude (VM), which is a univariate summary
of three-dimensional time-series of raw accelerometry data, was computed.
Lastly, 642 univariate vectors of VM were interpolated to have the same vector length,
scaled, and clustered into 1-5 clusters via correlation clustering.
The respective
cluster means
were defined to be 1-5 element subpopulation-specific stride templates, respectively.
}
\keyword{datasets}
