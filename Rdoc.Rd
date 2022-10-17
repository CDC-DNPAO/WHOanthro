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

waz, lhaz, bmiz, and headcz: WHOC --for-sex/age z-scores for Weight, Length/Height, BMI, and head circumference. 
  
}
\references{



}
\author{
David Freedman
}
\note{
Do NOT put arguments in quotation marks, such as cdcanthro(data,'age','wt','ht','bmi').  Use: cdcanthro(data, age, wt, ht, bmi)

Reference data are the merged LMS data files at
\url{https://www.cdc.gov/growthcharts/percentile_data_files.htm}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm}
}
\examples{
data = expand.grid(sex=1:2, agem=120.5, wtk=c(30,60), htc=c(135,144));
data$bmi = data$wtk / (data$htc/100)^2;
data = cdcanthro(data, age=agem, wt=wtk, ht=htc, bmi);
# OR data = cdcanthro(data, agem, wtk, htc, bmi);
round(data,2)
# setDF(data) to convert to a dataframe

nhanes   # NHANES data (2015/16 and 2017/18)
nhanes  = nhanes[!is.na(bmi),] # exclude subjects with missing wt/ht
nhanes$agemos = nhanes$agemos + 0.5   # because agemos is completed number of months
data = cdcanthro(nhanes, age=agemos, wt, ht, bmi, all=TRUE)
# OR data = cdcanthro(nhanes, agemos, wt, ht, bmi, all=TRUE)
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
