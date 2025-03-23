# PLSI-SCCS

To fit the PLSI-SCCS model, try the following codes (it may take hours to converge):

```R
library(splines)
library(tidyverse)
library(Matrix)
library(SCCS)
source("./SCCS_PLSI.R")
path_name <- "~/"
M_1 <- 10
M_2 <- c(4,4)
gamma <- exp(-2)

astart <- c(condat$sta)
aend <- c(condat$end)
case <- c(condat$case)
adrug <- cbind(condat$hib,condat$mmr)
adrug[is.na(adrug)] <-  999
aedrug <- cbind(adrug[,1]+49,adrug[,2]+49)
aevent <- cbind(condat$conv)

agebreaks <- seq(366,730,by=28)
agebreaks <- agebreaks[-c(1,length(agebreaks))]

model_fit <- SCCSFit_nopar(astart,
                           aend, 
                           case, 
                           aevent,
                           adrug, 
                           aedrug,
                           agebreaks,
                           M_1 = M_1,
                           M_2 = M_2,
                           gamma = gamma,
                           same_expo = FALSE,
                           ninitials = 1,
                           warmstart = NA)


saveRDS(model_fit,paste0(path_name,"model_fit.rds"))

```

