# whoanthro: WHO growth-standard z-scores for young children (0 to ~5 y)
#
# Computes WHO Child Growth Standard z-scores, plus biologically-implausible-
# value (BIV) flags, for weight-for-age, length/height-for-age, BMI-for-age,
# head-circumference-for-age, and weight-for-length.
#
# lenhei is recumbent length up to age 2.0 y (<= 730 d) and standing height for
# > 2 y (>= 731 d). The length- and height-based intervals do not overlap, so
# there is a deliberate jump in M between 730 and 731 days.
# Weight-for-length is 'for length' and does not go over 110 cm.
#
# sex can be coded 1/2, BOYS/GIRLS, b/g, m/f (case-insensitive).
# wt in kg, length/height (lenhei) in cm, age in DAYS.
#
# Load the reference table into `who_ref_data` before calling whoanthro(), e.g.:
# who_ref_data <- fread('/Volumes/Samsung_T5/Sync/R/Anal/Growth_Charts/WHO gcs/Data/WHO_ref.csv')
#
# 'BIV' flags (match the CRAN anthro package): f* = 1 if the z-score is outside the
# WHO plausibility range, 0 if inside, NA if the z-score is NA.
#   fwaz    : -6 <= waz    <= 5
#   flhaz   : -6 <= lhaz   <= 6
#   fwflz   : -5 <= wflz   <= 5
#   fbmiz   : -5 <= bmiz   <= 5
#   fheadcz : -5 <= headcz <= 5


 
# ---- helpers ----------------------------------------------------------------
 
# capture bare argument names as a character vector: cc(a, b) -> c("a","b")
cc <- function (...) as.character(substitute(c(...))[-1L])
 
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
 
# who_z should return only z, and compute the SD limits only on the rows 
# that need them. The old version computee four ^ operations for every row 
# even when correct=FALSE.  
who_z <- function(var, l, m, s, correct = FALSE) {
   z <- ((var/m)^l - 1) / (l*s)
   if (correct) {
      i <- which(z > 3)
      if (length(i) > 0) {
         li <- l[i]; mi <- m[i]; si <- s[i]; iv <- 1/li
         sdp2 <- mi*(1 + 2*li*si)^iv;  sdp3 <- mi*(1 + 3*li*si)^iv
         z[i] <- 3 + (var[i] - sdp3)/(sdp3 - sdp2)
      }
      i <- which(z < -3)
      if (length(i) > 0) {
         li <- l[i]; mi <- m[i]; si <- s[i]; iv <- 1/li
         sdm2 <- mi*(1 - 2*li*si)^iv;  sdm3 <- mi*(1 - 3*li*si)^iv
         z[i] <- -3 - abs((var[i] - sdm3)/(sdm2 - sdm3))
      }
   }
   z
}
 
 
whoanthro <- function(data,
               agedays = agedays,  
               wt      = wt,
               lenhei  = lenhei,
               bmi     = bmi,
               headc   = headc
               )
{
 
   # Work on a copy so we never modify the caller's object by reference.
   # as.data.table() copies a data.frame/tibble; a data.table it would return
   # unchanged, so we copy() that case explicitly.
   if (!is.data.table(data)) {
      data <- as.data.table(data)
   } else {
      data <- data.table::copy(data)
   }
   
   set(data, j = "seq_", value = seq_len(nrow(data))) # for merging back with original data
   set_cols_first(data,'seq_')
   dorig <- copy(data)
   
   nms <- grep('^sex$',names(data),ignore.case = TRUE, value = TRUE)
   if (length(nms) != 1) {
      stop ("A child's sex MUST be named 'sex' or 'SEX'; this is case insensitive.
             Also, you cannot have both 'sex' and 'SEX' as variables in your data.")
   }
   if (nms!='sex') {names(data)[which(names(data)==nms)] <- 'sex'}
   
   # ---- resolve the bare column-name arguments to strings, once, up front ----
   # Each *_var holds the NAME of a column (substitute() captures the symbol the
   # caller typed, deparse() makes it a string). headc and bmi are optional --
   # when omitted, missing() is TRUE and *_var is NULL, handled below.
   agedays_var <- deparse(substitute(agedays))
   wt_var      <- deparse(substitute(wt))
   lenhei_var  <- deparse(substitute(lenhei))
   headc_var   <- deparse(substitute(headc))
   bmi_var     <- deparse(substitute(bmi))
   
   # ---- check which measurement variables are available ---------------------
   # A variable is usable only if it was passed in the call AND that column
   # actually exists in the data.
   has_wt     <- !missing(wt)        && wt_var     %in% names(data)
   has_lenhei <- !missing(lenhei)    && lenhei_var %in% names(data)
   has_headc  <- !missing(headc) && headc_var  %in% names(data)
   has_bmi    <- !missing(bmi)   && bmi_var    %in% names(data)
   have_bmi <- has_bmi || (has_wt && has_lenhei)
   
   
   if (!has_wt)
      warning("Weight (wt) not called: waz, wflz, and bmiz will not be calculated.")
   if (!has_lenhei)
      warning("Length/height (lenhei) not called: lhaz and wflz will not be calculated.")
   if (!has_headc)
      warning("Head circumference (headc) not called: headcz will not be calculated.")
   if (!has_bmi && has_wt && has_lenhei)
      warning("BMI not supplied: it will be calculated as wt / (lenhei/100)^2.")
   if (!has_bmi && !(has_wt && has_lenhei))
      warning("BMI not found and cannot be calculated: bmiz will not be calculated.")
   
   # In the following code, data is the table to modify.
   # i is omitted, which means all rows. 
   # j = "wt" — the target column, given as a name. 
   #    If wt already exists it's overwritten; if it doesn't, it's created.
   # value = data[[wt_var]] — what to put there. 
   if (!agedays_var %in% names(data))
      stop("Age variable '", agedays_var, "' not found in data.")
   set(data, j = "agedays", value = data[[agedays_var]])
   
   av <- data$agedays
   if (!is.numeric(av))
      stop("Age variable '", agedays_var, "' must be numeric (age in days).")
   
   n_frac <- sum(abs(av - round(av)) > 1e-8, na.rm = TRUE)
   if (n_frac > 0) {
      warning(format(n_frac, big.mark = ","),
              " row(s) had non-integer agedays; rounded to whole days.")
      set(data, j = "agedays", value = round(av))
   }

   if (has_wt) {
      set(data, j = "wt", value = data[[wt_var]])
   } 
   if (has_lenhei) {
      set(data, j = "lenhei", value = data[[lenhei_var]])
   } 
   if (has_headc) {
      set(data, j = "headc", value = data[[headc_var]])
   } 
   
   if (has_bmi) {
      set(data, j = "bmi", value = data[[bmi_var]])
   } else if (has_wt && has_lenhei) {
      set(data, j = "bmi", value = data$wt / (data$lenhei / 100)^2)
   } 
   
   # sex can be coded almost any way
   # sex -> sexn (1 = boys, 2 = girls), tolerant of many codings.
   # Fast path when sex is already numeric 1/2. Otherwise map only the handful of
   # UNIQUE values and expand back with chmatch(), so substr()/toupper() run on
   # a few elements instead of all N rows.
   sx <- data$sex
   if (is.numeric(sx)) {
      sexn <- fcase(sx == 1, 1L,
                    sx == 2, 2L,
                    default = NA_integer_)
   } else {
      sx    <- as.character(sx)
      usex  <- unique(sx)
      fc    <- toupper(substr(usex, 1L, 1L))
      usexn <- fifelse(fc %chin% c("1", "B", "M"), 1L,
                       fifelse(fc %chin% c("2", "G", "F"), 2L, NA_integer_))
      sexn  <- usexn[chmatch(sx, usex)]
   }
   set(data, j = "sexn", value = sexn)
   
   # print(qsu(data))
   keep <- intersect(cc(seq_, sex, sexn, agedays, wt, lenhei, headc, bmi), 
                     names(data))
   data <- data[agedays <= 1856, ..keep]
   
   if (has_lenhei) {
      data[, len := NA_real_]
      data[agedays < 731, len := lenhei]
   }
   
   # ---- age-based indicators (waz, lhaz, bmiz, headcz) -----------------------
   # Only 'forage' rows can match on agedays, and we carry just the LMS columns
   # we actually use -- anything else would be duplicated across every data row.
   lmscols <- cc(wei_l, wei_m, wei_s, len_l, len_m, len_s,
                 bmi_l, bmi_m, bmi_s, headc_l, headc_m, headc_s)
   lmscols <- intersect(lmscols, names(who_ref_data))
   dref1 <- who_ref_data[denom=='forage' & agedays <= 1856,
                         c('sexn','agedays', lmscols), with = FALSE]
   dt1 <- dref1[data, on=c('sexn','agedays'),nomatch=0];
   
   # Calculate z-scores
   if (has_wt){
      dt1[, waz := who_z(wt, wei_l, wei_m, wei_s, correct = TRUE)]
   }
   
   # lhaz length/height for age z-score  (correct = FALSE)
   if (has_lenhei){
      dt1[, lhaz := who_z(lenhei, len_l, len_m, len_s, correct = FALSE)]
   }
   
   # bmiz  (weight-based -> correct = TRUE)
   if (have_bmi){
      dt1[, bmiz := who_z(bmi, bmi_l, bmi_m, bmi_s, correct = TRUE)]
   }
   
   # head circumference z-score  (correct = FALSE)
   if (has_headc){
      dt1[, headcz := who_z(headc, headc_l, headc_m, headc_s, correct = FALSE)]
   }
   
   # ---- weight-for-length (with reference interpolation) ---------------------
   # WFL needs BOTH weight and length; skipping the block also avoids building
   # dref2 and running a join that could not produce anything. wflz is created
   # INSIDE the block, so no all-NA wflz/fwflz columns appear when WFL can't
   # be computed.
   if (has_wt && has_lenhei) {
      dref2 <- who_ref_data[denom=='forlen' & len>=45 & len<=110,
                            .(sexn,len,wfl_l,wfl_m,wfl_s)]
      
      # get ulens in 'data'
      ulens <- unique(data$len[!is.na(data$len)])
      # interpolate reference data to match each length in the input data
      dlen  <- length(setdiff(ulens, dref2$len))
      
      # --- REFERENCE CHART INTERPOLATION ---
      # If the input contains lengths not mapped directly inside the standard WHO
      # tables, we execute linear interpolation via stats::approx across
      # LMS parameters to generate custom growth metrics matching the user's data
      interp_vars <- cc(wfl_l, wfl_m, wfl_s)
      
      if (dlen > 0) {
         dref2 <- dref2[, c(
            list(len = ulens),
            lapply(.SD, function(col) approx(len, col, xout = ulens)$y)
         ), by = sexn, .SDcols = interp_vars]
      }
      
      # on= avoids setkey(), which would physically sort all N rows of `data`
      dt2 <- dref2[data, on = c('sexn','len'), nomatch = 0]
      dt2[, wflz := who_z(wt, wfl_l, wfl_m, wfl_s, correct = TRUE)]
      
      dt1[, wflz := NA_real_]
      dt1[dt2, wflz := i.wflz, on = 'seq_']
      dt1[agedays > 730, wflz := NA_real_]   # WFL applies only to length (<= 2 y)
   }
   
   # blank z-scores outside the valid age ranges while agedays is still here.
   # The BIV flags below then come out NA for these rows automatically.
   for (cl in intersect(cc(waz, lhaz, bmiz, headcz), names(dt1)))
      dt1[agedays > 1856, (cl) := NA_real_]
   zs <- intersect(cc(waz, lhaz, wflz, bmiz, headcz), names(dt1))
   dt <- dt1[, c('seq_', zs), with = FALSE]
   
   # 'biologically implausible value (BIV)' flags (see header for ranges).
   if ('waz'    %in% names(dt)) dt[, fwaz    := as.integer(waz    < -6 | waz    > 5)]
   if ('lhaz'   %in% names(dt)) dt[, flhaz   := as.integer(lhaz   < -6 | lhaz   > 6)]
   if ('wflz'   %in% names(dt)) dt[, fwflz   := as.integer(wflz   < -5 | wflz   > 5)]
   if ('bmiz'   %in% names(dt)) dt[, fbmiz   := as.integer(bmiz   < -5 | bmiz   > 5)]
   if ('headcz' %in% names(dt)) dt[, fheadcz := as.integer(headcz < -5 | headcz > 5)]
   
   # Rename any original columns whose names collide with the computed output
   # columns by appending '0'. Without this the join produces duplicate names
   # (i.waz, etc.) and the user's original values can mask the computed ones.
   overlap <- setdiff(intersect(names(dorig), names(dt)), 'seq_')
   # subtracts 'seq_' from the intersection. it's the 'y' in setdiff(x,y)
   if (length(overlap) > 0) {
      setnames(dorig, overlap, paste0(overlap, '0'))
      warning(paste(overlap, collapse = ', '),
              " in original data renamed to ",
              paste(paste0(overlap, '0'), collapse = ', '),
              " to avoid conflict with computed values.")
   }
   
   message(
      "Z-scores for weight, lenhei (length or height), BMI, and headc are ",
      "calculated for children <= 1856 days.\n",
      "Z-scores for WFL (weight-for-length) are calculated for children ",
      "<= 730 days who have a length between 45 and 110 cm (inclusive)."
   )
   
   dt <- dt[dorig, on='seq_']
   set_cols_first(dt, names(dorig))
   dt[, seq_ := NULL]
   dt[]
   }
