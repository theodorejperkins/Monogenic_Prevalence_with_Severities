# This is code for estimating the prevalence of monogenic autosomal recessive diseases, 
# based on a list of known pathogenic variants and their observed frequencies from some
# population database. It accompanies the paper "Estimating the Prevalence of GNE Myopathy
# using Population Genetic Databases" by Alexa Derksen et al. This code was based in part
# on the code at https://github.com/leklab/LAMA2-CMD_prevalence. However the code is
# considerably revised to include three different prevalence estimates, to account for
# three different severities of pathogenic variants, and to estimate confidence intervals
# based on resampling. The code was developed by Madeeha Shaikh and Theodore J. Perkins at 
# the University of Ottawa / Ottawa Hospital Research Institute in 2024.

# The code depends on three inputs which are specified below: 
# 1) One (or more) variant severity and frequency files,
# 2) A file stating Bayesian prior parameters for different types of variants,
# 3) Whether or not to calculate confidence intervals
# The variant severity/frequency files are similar to spreadsheets one can obtain from 
# the gnomad (https://gnomad.broadinstitute.org/) website. The necessary columns are the 
# POP_AC and POP_AN columns for different populations POP, as well as an additional column
# "Disease association" which gives a value of "Homozygous non disease causing", or 
# "Homozygous disease causing" or "Homozygous embryonic lethal" for each variant.
# The code detects which populations POP are present in the variants file, and will
# output six estimates for prevalence in each population: maximum likelihood and two
# different Bayesian estimates (see manuscript), each with or without the assumption of
# different variant severities. If requested, it will use 1000 resamples of the data to
# generate a 95% confidence interval for every estimate. But be warned that this can take
# considerable computation time.

VariantFiles <- c('pathogenic_hg37_0129.csv', 'Pathogenic_+CADD23_hg37_0129.csv', 'Pathogenic_+CADDgreaterthan13_insilicopath_hg37_0129.csv', 'Pathogenic+all CADDgreaterthan13_hg37_0129.csv', 'hg38_pathogenic_0129.csv')
PriorsFile <- 'beta_parameter_prior_ExAC_20240218.txt'
CalculateConfidenceIntervals = TRUE

library(dplyr)

#posterior distribution
get_posterior_w_v <- function(x,n,w,v){
  a <- x+v
  b <- n-x+w
  return(list(a,b))
}

#mean 
get_mean_normal <- function(vs,ws){
  tmp <- vs/(vs+ws)
  return(sum(tmp))
}

#variance
get_var_normal <- function(vs,ws){
  vars <- vs*ws/((vs+ws+1)*(vs+ws)^2)
  return(sum(vars))
}

#input dataframe and population to examine and get variables needed for calculation: mean normal, variance, and mean allele frequency
analyze_data <- function(data, params, pop) {
  pop_ac_col = paste(pop,"AC",sep="_")
  pop_an_col = paste(pop,"AN",sep="_")
  
  posterior_param <- data.frame(v=rep(0,nrow(data)),w=rep(0,nrow(data)))
  rows = nrow(data)
  max_freq = 0
  
  for (i in 1:rows) {
    v <- params$v[which(params$type == data$VEP.Annotation[i])]
    w <- params$w[which(params$type == data$VEP.Annotation[i])]
    
    AC_interested <- pull(data, pop_ac_col)[i]
    AN_interested <- pull(data, pop_an_col)[i]
    
    post_params <- get_posterior_w_v(AC_interested, AN_interested, w, v)
    
    posterior_param[i, 1] <- post_params[[1]]
    posterior_param[i, 2] <- post_params[[2]]
    
    max_freq = max_freq + (AC_interested/AN_interested)  
  }
  
  mu <- get_mean_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])
  sigma.2 <- get_var_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])
  
  return(list(mu,sigma.2, max_freq))
}

# calculate_all_estiates: 
# Input gnomad data sheet, bayesian prior parameters, and population name/type to be analyzed. 
# Will call "analyze_data" function to get parameters needed to run all 6 calculations.
# Returns a list of the six esimates.
calculate_all_estimates <- function(data,params,pop) {
  
  NC_df <- data %>% filter(`Disease.association` == "Homozygous non disease causing")
  DC_df <- data %>% filter(`Disease.association` == "Homozygous disease causing")
  EL_df <- data %>% filter(`Disease.association` == "Homozygous embryonic lethal")
  
  NC_val <- analyze_data(NC_df, params, pop)
  DC_val <- analyze_data(DC_df, params, pop)
  EL_val <- analyze_data(EL_df, params, pop)
  all_val <- analyze_data(data, params, pop)
  
  f = NC_val[[1]]
  g = DC_val[[1]]
  h = EL_val[[1]]
  var = DC_val[[2]]
  
  f_max = NC_val[[3]]
  g_max = DC_val[[3]]
  h_max = EL_val[[3]]
  
  mu <- all_val[[1]]
  sigma.2 <- all_val[[2]]
  q <- all_val[[3]]
  
  Bay_prob_disease_noVar = 2*f*g + 2*f*h + 2*g*h + g^2
  Bay_prob_disease_wVar = 2*f*g + 2*f*h + 2*g*h + g^2 + var
  prob_disease_max = 2*f_max*g_max + 2*f_max*h_max + 2*g_max*h_max + g_max^2
  
  NoASSOC_Bay_prob_disease_noVar = mu^2
  NoASSOC_Bay_prob_disease_wVar = mu^2 + sigma.2
  NoASSOC_HW = q^2
  
  all_estimates = c(prob_disease_max,
                    Bay_prob_disease_noVar,
                    Bay_prob_disease_wVar,
                    NoASSOC_HW,
                    NoASSOC_Bay_prob_disease_noVar,
                    NoASSOC_Bay_prob_disease_wVar
                    )
  
  return(all_estimates)
}


# A function that resamplings the AC values for a specified population,
# by assuming they are binomially distributed with success probability
# AC/AN and number of tries AN.
resample_pop <- function(data,pop) {

  # Resampled dataframe
  resdata <- data
  
  # Find column of data that is to be resampled
  ACName = paste(pop,"_AC",sep="")
  ACi = which(names(data)==ACName)

  # Loop through rows of matrix resampling the values
  for (rowi in 1:nrow(data)) {
    AC <- data[rowi,ACi]
    AN <- data[rowi,ACi+1]
    resdata[rowi,ACi] <- rbinom(1,round(AN),AC/AN)
  }

  # Return
  return(resdata)
}

# Calculate confidence intervals in each estimate by resampling
# the data, and re-running the estimators, many times!
all_estimates_confidence <- function(data,params,pop) {
  NTimes <- 1000
  Lo <- 25
  Hi <- 975
  all_est1 <- c()
  all_est2 <- c()
  all_est3 <- c()
  all_est4 <- c()
  all_est5 <- c()
  all_est6 <- c()
  
  # Resample many times, and accumulate estimates in separate lists
  for (i in 1:NTimes) {
    set.seed(i)
    resdata <- resample_pop(data,pop)
    resest <- calculate_all_estimates(resdata,params,pop)
    all_est1 <- c(all_est1,resest[1])
    all_est2 <- c(all_est2,resest[2])
    all_est3 <- c(all_est3,resest[3])
    all_est4 <- c(all_est4,resest[4])
    all_est5 <- c(all_est5,resest[5])
    all_est6 <- c(all_est6,resest[6])
  }

  # Sort the lists
  all_est1 <- sort(all_est1)
  all_est2 <- sort(all_est2)
  all_est3 <- sort(all_est3)
  all_est4 <- sort(all_est4)
  all_est5 <- sort(all_est5)
  all_est6 <- sort(all_est6)
  
  # Output confidence bounds
  ci <- c(all_est1[Lo],all_est1[Hi],
          all_est2[Lo],all_est2[Hi],
          all_est3[Lo],all_est3[Hi],
          all_est4[Lo],all_est4[Hi],
          all_est5[Lo],all_est5[Hi],
          all_est6[Lo],all_est6[Hi])
}

# Find populations -- look through the data (gnomad dataframe) to identify
# pairs of columns ending in _AC and _AC, which are distinct populations
# to analyze.
get_populations <- function(data) {
  
  # What are all the column names?
  ColNames <- names(data)
  
  # Accumulate population names
  PopNames = c()
  for (coli in 1:(length(ColNames)-1)) {
    ThisCol <- ColNames[coli]
    ThisLen <- nchar(ThisCol)
    NextCol <- ColNames[coli+1]
    NextLen <- nchar(NextCol)
    if (ThisLen>=3 & NextLen>=3) {
      ThisEnd <- substr(ThisCol,ThisLen-2,ThisLen)
      NextEnd <- substr(NextCol,ThisLen-2,ThisLen)
      if ((ThisEnd=="_AC") & (NextEnd=="_AN")) {
        # We found a column pair!
        ThisPop <- substr(ThisCol,1,ThisLen-3)
        PopNames = c(PopNames,ThisPop)
      }
    }
  }
  
  return(PopNames)
}


# A function for printing out the six estimates, and confidence
# intervals, if provided.
pretty_print_all_estimates <- function(all_estimates,pop,conf=FALSE) {

  if (length(conf)==1) {
    cat("\n")
    cat(pop,"Maximum likelihood estimated prevalence with variant severity (per million):",all_estimates[1]*1e6,"\n")
    cat(pop,"Bayesian estimated prevalence without variance with variant severity (per million):",all_estimates[2]*1e6,"\n")
    cat(pop,"Bayesian estimated prevalence with variance with variant severity (per million):",all_estimates[3]*1e6,"\n")
    cat(pop,"Maximum likelihood estimated prevalence without variant severity (per million): ",all_estimates[4]*1e6,"\n")
    cat(pop,"Bayesian estimated prevalence without variance without variant severity (per million):",all_estimates[5]*1e6,"\n")
    cat(pop,"Bayesian estimated prevalence with variance without variant severity (per million):",all_estimates[6]*1e6,"\n")
  } else {
    cat("\n")
    cat(pop,"Maximum likelihood estimated prevalence with variant severity (per million):",all_estimates[1]*1e6," (",conf[1]*1e6,",",conf[2]*1e6,")\n")
    cat(pop,"Bayesian estimated prevalence without variance with variant severity (per million):",all_estimates[2]*1e6," (",conf[3]*1e6,",",conf[4]*1e6,")\n")
    cat(pop,"Bayesian estimated prevalence with variance with variant severity (per million):",all_estimates[3]*1e6," (",conf[5]*1e6,",",conf[6]*1e6,")\n")
    cat(pop,"Maximum likelihood estimated prevalence without variant severity (per million):",all_estimates[4]*1e6," (",conf[7]*1e6,",",conf[8]*1e6,")\n")
    cat(pop,"Bayesian estimated prevalence without variance without variant severity (per million):",all_estimates[5]*1e6," (",conf[9]*1e6,",",conf[10]*1e6,")\n")
    cat(pop,"Bayesian estimated prevalence with variance without variant severity (per million):",all_estimates[6]*1e6," (",conf[11]*1e6,",",conf[12]*1e6,")\n")
  }
}



# "Main" code starts here!

params = read.csv(PriorsFile, header = TRUE, sep = "\t")
# Loop through each data file
for (DataFile in VariantFiles) {
  data = read.csv(DataFile, header = TRUE, sep = ",")
  cat("\n\n\nFile:", DataFile, "\n")

  # What populations are present in this data?
  popnames <- get_populations(data)
  cat("Found populations:", popnames,"\n")

  # Loop through detected populations, making estimate for each
  for (pop in popnames)  {
    all_est <- calculate_all_estimates(data,params,pop)
    if (CalculateConfidenceIntervals) {
      conf_int <- all_estimates_confidence(data,params,pop)
    } else {
      conf_int <- FALSE
    }
    pretty_print_all_estimates(all_est,pop,conf_int)
    flush.console()
  }
}
