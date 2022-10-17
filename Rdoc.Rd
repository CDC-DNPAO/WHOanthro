\name{whoanthro}
\alias{whoanthro}
\alias{age_in_months}
\alias{weight_kg}
\alias{height_cm}
\alias{ref_data}
\alias{nhanes}
\alias{data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GENERATE SEX- AND AGE (or LENGTH)-STANDARDIZED WEIGHT, HEIGHT/LENGHT, AND BMI METRICS FROM THE WHO GROWTH CHARTS

}
\description{
Generate z-scores, percentiles, and other metrics for weight, height, BMI, and weight for length/height based on the WHO growth charts.
Has a single function, 'whoanthro'.  Requires the package data.table to be
installed; library(whoanthro) will also attach data.table.

}
\usage{
whoanthro(data, agedays = age_in_days, wt = weight_kg, lenhei = length or height, headc=head_circumference, bmi = bmi)
# If a variable (e.g, head circumference) is missing for all records or if you're not interested in the z-score, specify as '= NA' in funtion.
# So, if not interested in head circumference, use: whoanthro(data, agedays, wt, lenhei, headc=NA, bmi)

# Do NOT put arguments in quotation marks

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data.frame or data.table}
  \item{agedays}{agedays in days.}
  \item{wt}{weight (kg).}
  \item{lenhei}{recumbent length or standing height (cm).}
  \item{headc}{head circumference (cm).}
  \item{bmi}{BMI, kg/m^2.}
}

\details{
Expects 'sex' to be a variable in the dataset. Can be coded as either 'boys/girls' or 
'male/female' or '1/2'.  Character values can be upper or lower case; only the first character is considered.

Weight is in kg, and lenhei is length or height in cm.  BMI is kg/m^2.

If age is in completed number of months, multiply by 30.4375 and add 15.219 (1/2 of 30.4375).  Then take the floor of this value. Age in the WHO reference data is givent as 0 to 1826 days.
If there is a single 'height' field, it iss likely that recumbent length was measured for children under 2.0 years.  If there are separate height and length fields, let 'lenhei' be equal to length 
for children under 2.0 years and be equal to height for children over 2.0 y.
}

\value{
Returns a data.table containing the original data and various z-scores and percentiles.  If convert a data.table to this to a dataframe with 'setDF(output_data)'.

Variables in output:

waz: weight-for-sex/age z-score
lhaz: length (or height) for sex/age z-score
wflz: weight-for-length (and sex) z-score
bmiz: BMI-for-sex/age z-score
headcz: head circumference for sex/age z-score

  
}
\references{
World Health Organization. Department of Nutrition for Health and Development. WHO Child Growth Standards. Length/height-for-age, weight-for-age, weight-for-length, weight-for-height, and body mass index-for-age. Methods and Development [Internet]. Geneva; 2006.  
Available from: http://www.who.int/childgrowth/standards/technical_report/en/

}
\author{
David Freedman
 
}
\note{
Do NOT put arguments in quotation marks, such as cdcanthro(data,'agedays', 'wt','ht','headc', 'bmi').  Use: cdcanthro(data, age, wt, ht, headc, bmi)

Reference data are the expanded merged LMS data tables (at the bottom of various pages) on the WHO web site.  
For example, length/height for age is at the bottom of 
\url{https://www.who.int/tools/child-growth-standards/standards/length-height-for-age}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}

\examples{
nhanes   # NHANES data (1999/2000 through 2017/18)
data = whoanthro(nhanes, agedays, wt, lenhei, headc, bmi) 
round(data,2)

# if head circumference, for example, is not in data or all its values are missing: 
data = whoanthro(nhanes, agedays, wt, lenhei, headc=NA, bmi) # if head circ. if not in input dataset
round(data, 2)

# setDF(data) to convert to a dataframe


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
