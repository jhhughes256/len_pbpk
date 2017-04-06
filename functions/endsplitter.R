# This function is designed specifically for use with the spreadsheet
# "All_Tissue_PK_Data_Summary.xls" after the data has been processed by
# datacheck_nonclin.R up to:
# + line 152 for IV datasets
# + line ### for PO datasets
# + line ### for PR datasets

# Can be debugged with the test vector below
x <- c(
   "5hr_1", "45min_2", "1_5mg_1h_4", "1_5mg_45m_3", "1_5mg_45m_3_Repeat",
   "1_5mg_20m_5_090925052635", "300min_1", "1_5hr_5", "0_5mg_Gav_1_5h_1",
   "Gavage_16h_1", ""
 )
# Output should be:
#Debug Output
#  UID ID TADNOM
#1  1  1   300
#2  2  2    45
#3  3  4    60
#4  4  3    45
#5  4  3    45
#6  5  4    20
#7  6  1   300
#8  7  5    90
#9  8  1    90
#10  9  1   960
#11 10  2   960

end.splitter <- function(x) {
  xsplit.t <- 0  # placeholder
	tsplit.t <- 0
	uid.t <- 0
	uid.mark <- 0	 # to be subtracted from i to determine unique ID
  for (i in 1:length(x)) {
    xstr <- str_sub(x[i], -1)  # final character in string
		ystr <- str_sub(x[i], -2, -2)  # penultimate character in string
		zstr <- str_sub(x[i], -3, -3)  # maybe final character of time string
    if (length(unique(0:9 %in% xstr)) == 2) {  # separate numbers from text
			if (match("_", ystr, nomatch = FALSE) == 1) {  # separates ID from junk
        # data here is normal ID
			  xsplit.t <- c(xsplit.t, xstr)  # save ID to vector
			  if (match("r", zstr, nomatch = FALSE) == 0) {  # remove hr text
			    if (match("n", zstr, nomatch = FALSE) == 0) {  # remove min text
				    suppressWarnings(avec <- as.numeric(str_sub(x[i], -5, -4)))
					  tstr <- str_sub(x[i], -4, -4)
				    if (match("h", zstr, nomatch = FALSE) == 0) {	# double digit m
						  if (!is.na(avec) == TRUE) {
								tsplit.t <- c(tsplit.t, as.character(avec/60))
							} else {  # single digit m
								tvec <- as.numeric(tstr)
								tsplit.t <- c(tsplit.t, as.character(tvec/60))
							}
						} else {  #non split and split h
							astr <- str_sub(x[i], -6, -6)
							if (length(unique(0:9 %in% astr)) == 2) {	# split time h
								tsplit.t <- c(tsplit.t, paste(astr, tstr, sep = "."))
							} else {  # single or double digit h
								if (!is.na(avec) == TRUE) {  # double digit h
									tsplit.t <- c(tsplit.t, avec)
								} else {  # single digit h
									tsplit.t <- c(tsplit.t, tstr)
								}
							}
						}
					} else {  # handles min from here
						suppressWarnings(avec <- as.numeric(str_sub(x[i], -7, -6)))
						if (!is.na(avec) == TRUE){  # double or triple digit min
							suppressWarnings(tvec <- as.numeric(str_sub(x[i], -8, -6)))
							if (!is.na(tvec) == TRUE){  # triple digit min
								tsplit.t <- c(tsplit.t, as.character(tvec/60))
							} else {  # double digit min
								tsplit.t <- c(tsplit.t, as.character(avec/60))
							}
						} else {  # single digit min
							tvec <- as.numeric(str_sub(x[i], -6, -6))
							tsplit.t <- c(tsplit.t, as.character(tvec/60))
						}
					}
				} else {  # handles hr from here
					tstr <- str_sub(x[i],-5,-5)
					astr <- str_sub(x[i],-7,-7)
					if (length(unique(0:9 %in% astr)) == 2){
						tsplit.t <- c(tsplit.t,paste(astr, tstr, sep="."))  #split time hr
					} else {
						tsplit.t <- c(tsplit.t,tstr)  #single digit hr
					}
				}
			} else {
        # Data here has junk ID attached to end of normal ID
        # The aim is to set new number to replace junk with last value + 1
        # Time must be of format "double digit min"
				temp <- as.numeric(tail(xsplit.t,1)) + 1
				xsplit.t <- c(xsplit.t, temp)
				tvec <- as.numeric(str_sub(x[i], -18, -17))
				tsplit.t <- c(tsplit.t, as.character(tvec/60))
        }
		} else {
		  if (str_detect(x[i], "_")) {  # data here is for the repeated sample
				xsplit.t <- c(xsplit.t,tail(xsplit.t,1))  # use last number
				uid.mark <- uid.mark + 1
				tsplit.t <- c(tsplit.t,tail(tsplit.t,1))
			} else {  # data here is missing sample names
				temp <- as.numeric(tail(xsplit.t, 1)) + 1
				xsplit.t <- c(xsplit.t, temp)  # replace junk with last value + 1
				temp2 <- tail(tsplit.t, 1)
				tsplit.t <- c(tsplit.t, temp2)  # use last value
			}
    }
		uid.t <- c(uid.t, i - uid.mark)
  }
	xsplit <- xsplit.t[-1]	  #strip first value from vector to leave desired output
	tsplit <- tsplit.t[-1]
	uid <- uid.t[-1]
	output <- data.frame(uid, xsplit, as.numeric(tsplit)*60)
	colnames(output) <- c("UID", "ID", "TADNOM")
	output
}
