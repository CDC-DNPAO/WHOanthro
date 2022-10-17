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
To install package:
install.packages('https://raw.github.com/CDC-DNPAO/WHOAnthro/master/whoanthro_0.1.1.tar.gz', type='source', repos=NULL)

whoanthro(data, agedays = age_in_days, wt = weight_kg, lenhei = length or height, headc=head_circumference, bmi = bmi)
# If a variable (e.g, head circumference) is missing for all records, specify as '= NA' in funtion.
# So, whoanthro(data, agedays, wt, lenhei, headc=NA, bmi)

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

}

\value{
Returns a data.table containing the original data and various weight, height, and BMI metrics.  Can convert this to a dataframe with 'setDF(output_data)'.

Variables in output:

waz: weight-for-sex/age z-score
lhaz: length (or height) for sex/age z-score
wflz: weight-for-length for sex z-score
bmiz: BMI-for-sex/age z-score
headcz: head circumference for sex/age z-score

}
\references{
World Health Organization. Department of Nutrition for Health and Development. WHO Child Growth Standards. Length/height-for-age, weight-for-age, weight-for-length, weight-for-height, and body mass index-for-age. Methods and Development [Internet]. Geneva; 2006.  

http://www.who.int/childgrowth/standards/technical_report/en/
}

\author{
David Freedman
}
\note{
Do NOT put arguments in quotation marks, such as whoanthro(data,'agedays','wt','lenhei', 'headc','bmi').  Use: whocanthro(data, agedays, wt, lenhei, headc)

Reference data are the EXPANDED LMS data tables (at bottom of page) files at
\url{https://www.who.int/tools/child-growth-standards/standards/weight-for-length-height}
\url{https://www.who.int/tools/child-growth-standards/standards/length-height-for-age}
\url{https://www.who.int/tools/child-growth-standards/standards/weight-for-age}
\url{https://www.who.int/toolkits/child-growth-standards/standards/body-mass-index-for-age-bmi-for-age}
\url{https://www.who.int/tools/child-growth-standards/standards/head-circumference-for-age}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{http://www.who.int/childgrowth/standards/technical_report/en/}
}
\examples{
nhanes   # NHANES data (2015/16 and 2017/18)
data = whoanthro(nhanes, agedays, wt, lenhei, headc, bmi)
round(data,2)

# if head circumference, for example, is not in data or all its values are missing: 
data = whoanthro(nhanes, agedays, wt, lenhei, headc=NA, bmi)
round(data, 2)


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
