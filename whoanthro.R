
# need data in wide format for BMIz and WFLz
# need to have all variables in data; agedays,wt,lenhei,headc, and bmi
# but can run as wanthro(data, agedays, wt, lenhei, bmi, headc=NA);

# lenhei is length up to age 2.0 y and height for >2 y
# http://www.who.int/childgrowth/standards/technical_report/en/

# WHO: The WHO length- and height-based BMI-for-age standards do not 
# overlap, i.e. the length-based interval ends at 730 days and the 
# height-based interval starts at 731 days.
# there is a 'big' jump in M between 730 and 731 days

# who_ref_data <- fread('~/Sync/R/Anal/Growth_Charts/Data/WHOref_d.csv'); 

#############################
# this is for children between 0 and <24 months of age
# computes z-score, modified z-score, extended WHZ (and percentiles) 
# and other BMI metrics based on the WHO growth charts

# Note - sex can be coded as 1 (boys) or 2 (girs) OR BOYS/GIRLS or b/g
# Alen_so need wt (kg), length (cm), and age in *** days ***
# The names of wt, length and age can be anything in your data
# Alen_so need the reference data for the len_mS files

# About the WHO reference values for the L M and S parameters:
# refdata_dir is folder that contains 'CDCref_d.csv' - this file can be download from  
# https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm
# (4th paragraph, before 'Instructions for SAS Users)

# lenhei

.c <- function (...) as.character(substitute(c(...))[-1L])

fnames <- function(vars,d){
      vars=unique(vars);
      vars[vars %nin% names(d)]
      }

fdel <- function(x,d){
      x=x[x %nin% fnames(x,d)];
      if (length(x)>0) d[,(x):=NULL]; 
}   

fvalid <- function(x){
   length(x[!is.na(x)])
}

set_cols_first <- function (DT, colen_s, intersection = TRUE) # thanks to hutils
   {
      if (intersection) {
         return(setcolorder(DT, c(intersect(colen_s, names(DT)), 
                                  setdiff(names(DT), colen_s))))
      }
      else {
         return(setcolorder(DT, c(colen_s, setdiff(names(DT), colen_s))))
      }
}

w_zscore=function(var, l, m, s){ 
   ls=l*s; invl=1/l
   z = (((var/m) ^ l) -1) / (ls) # z-score formula
   sdp2 = (m * (1 + 2*ls) ^ (invl)) 
   sdp3 = (m * (1 + 3*ls) ^ (invl))
   sdm2 = (m * (1 - 2*ls) ^ (invl)) 
   sdm3 = (m * (1 - 3*ls) ^ (invl))
   z=fcase(
      z >= -3 & z < 3, z,
      z > 3, 3 + (var - sdp3)/(sdp3 - sdp2),
      z < (-3), -3 - abs((var - sdm3)/(sdm2 - sdm3))
   )
   list(z, sdp2, sdp3, sdm2,sdm3)
}

whoanthro <- function(data,
               agedays = agedays, 
               wt = wt, 
               lenhei = lenhei,
               headc = headc,
               bmi = bmi
               )
{
      wfl_m <- seq_ <- denom <- 
      waz <- sdp2 <- sdp3 <- sdm2 <- sdm3 <- len <- 
      lhaz <- bmiz <- headz <- lhaz <- wfl_l <- wflm <- wfl_s <- wflz <- 
      weight <- height <- bmi_l <- bmi_s <-  bmi_m <- NULL
   
   if (is.data.table(data) == FALSE) data <- as.data.table(data) 
   data$seq_ <- 1L:nrow(data) # for merging back with original data
   set_cols_first(data,'seq_')
   dorig <- copy(data)    
   
   nms <- grep('^sex$',names(data),ignore.case = TRUE, value = TRUE)
   if (length(nms) != 1) {
      stop ("A child's sex MUST be named 'sex' or 'SEX'; this is case insensitive.
             Also, you cannot have both 'sex' and 'SEX' as variables in your data.")
   }
   if (nms!='sex') {names(data)[which(names(data)==nms)] <- 'sex'}
   
   data$agedays <- data[[deparse(substitute(agedays))]]
   data$wt <- data[[deparse(substitute(wt))]]
   data$bmi <- data[[deparse(substitute(bmi))]]
   data$lenhei <- data[[deparse(substitute(lenhei))]]
   data$headc <- data[[deparse(substitute(headc))]]
   
   if (('agedays' %chin% names(data)) == FALSE){
      stop('There must be an variable for age in days in the data')
   }
   
   # sex can be coded almost any way
   data[,sexn:=toupper(substr(sex,1,1))]
   data[,sexn:=fcase(
      sexn %in% c(1,'B','M'), 1L,
      sexn %in% c(2,'G','F'), 2L
   )]
   
   data <- data[agedays<=1856, .(seq_, sex,sexn,agedays,wt,lenhei,headc,bmi)]; 
   # WFL calculations do not go over 110 cm because they're 'for length'
   # (ame as WFH calculation in CRAN anthro)

   dref1 <- who_ref_data[(denom=='forage' & agedays<1857) | denom=='forlen']
   # setkey(data,sexn,agedays); setkey(dref1,sexn,agedays)
   dt1 <- dref1[data, on=c('sexn','agedays'),nomatch=0]; 
   
   # waz
   dt1[,.c(waz, sdp2,sdp3,sdm2,sdm3):= 
          w_zscore(dt1$wt, dt1$wei_l, dt1$wei_m, dt1$wei_s)]
 
   # lhaz length/height for age z-score
   dt1[,.c(lhaz, sdp2,sdp3,sdm2,sdm3):= 
          w_zscore(dt1$lenhei, dt1$len_l, dt1$len_m, dt1$len_s)]
   
   # bmiz
   if (fvalid(dt1$bmi) > 0){
   dt1[,.c(bmiz, sdp2,sdp3,sdm2,sdm3):= 
          w_zscore(dt1$bmi, dt1$bmi_l, dt1$bmi_m, dt1$bmi_s)]
   }
   
   # head circumference z-score
   if (fvalid(dt1$headc) > 0){
   dt1[,.c(headcz, sdp2,sdp3,sdm2,sdm3):= 
          w_zscore(dt1$headc, dt1$headc_l, dt1$headc_m, dt1$headc_s)]
   }
   
   # wflz (weight-for-length z-score)
   dref2 <- who_ref_data[denom=='forlen',.(sexn,len,wfl_l,wfl_m,wfl_s)]
   
   setkey(data,sexn,lenhei); setkey(dref2,sexn,len)
   dt2 <- dref2[data,nomatch=0]; 
   dt2[,.c(wflz, sdp2,sdp3,sdm2,sdm3):= 
          w_zscore(dt2$wt, dt2$wfl_l, dt2$wfl_m, dt2$wfl_s)]
   
   x <- unique(grep('sd[mp]|[incl]_|sexn', names(dt1), value=T)); 
   dt1[,(x):=NULL] 
   x2 <- unique(grep('sd[mp]|wfl_|sexn', names(dt2), value=T)); 
   dt2[,(x2):=NULL] 
 
   dt <- merge(dt1,dt2[,.(seq_,wflz)], by='seq_', all=T)
   dt[,wflz:=fifelse(agedays>730,NA_real_,wflz)] 
   vars<- grep('seq_|z$', names(dt), value=TRUE); 
   dt <- dt[,..vars]
   
   dtot <- dt[dorig, on='seq_']; 
   set_cols_first(dtot,names(dorig))
   dtot[,seq_:=NULL]
   
   dtot[]
}


