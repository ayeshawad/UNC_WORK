# Filter samples to use in this run                                                     ## These are just options from Terry's template
#coldata = coldata[coldata$TIN_median > 68,] # Filter low TIN values                    ## I don't have TIN values in coldata

##coldata <- coldata %>% # need this an any of the following options. 
#dplyr::filter(!grepl("IL", subtype)) %>% # Filter CD IL samples                      ## I don't have IL data
#dplyr::filter(!grepl("UC", disease)) %>% # Filter UC samples                         ## More complicated to reformat
#dplyr::filter(!sample %in% c("64", "85")) # Filter particular samples                ## An easier way to eliminate UC Pt, reformat

## I am removing these things so that there are no dups                                 ## The above filters are not working for me. I will try the which fx
x <-  coldata[which(coldata[,"disease"] != "UC"),]
y <-  x[which(x[,"inflammation"] != "inflamed"),]
z <-  y[which(y[,"tissue"] != "ileum"),]
coldata <- z