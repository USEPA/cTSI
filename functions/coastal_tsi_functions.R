
## Coastal TSI Functions

# J. Hagy
# 7/15/2022
# 
# Functions supporting data analysis for
# 
# Quantifying coastal ecosystem condition and a trophic state index with a Bayesian analytical framework
# Hagy, JD, B. Kreakie, M. Pelletier, F. Nojavan, J. Kiddon, and A. Oczkowski

# Replace region names with better variable names
# and change to factors
fixRegions <- function(inData){
  outData <- inData %>% 
  mutate(REGION=as.factor(REGION),
           SUBREGIONS=as.factor(SUBREGIONS)) %>% 
    mutate(levelIII=as.factor(levelIII),
           AggEco=as.factor(AggEco),
           EMAP=as.factor(EMAP_BioGeoProv),
           BioGeo=as.factor(Biogeo_subregion)) %>% 
    select(-Biogeo_subregion,-EMAP_BioGeoProv)
  return(outData)
}

# Add secchi data for 2015 NCA data
add2015secchi <- function(inData) {
# Read 2015 secchi data
secchi2015 <- read.csv("Raw/Secchi2015.csv", stringsAsFactors = FALSE)
names(secchi2015) <- c("SITE_ID","VISNUM","SECCHI_MEAN..m.","Visible_on_Bottom")

# A few of the secchi depths have more than one value for SITE_ID and VISNUM.  
# I averaged them
secchi2015 <- secchi2015 %>% 
  group_by(SITE_ID,VISNUM) %>% 
  summarize(secchi_update=mean(SECCHI_MEAN..m.)) %>% ungroup() %>% 
  mutate(SAMPYEAR=2015)

# Join 2015 secchi data to coastal NCA data
tmp <- left_join(as_tibble(inData),as_tibble(secchi2015),
                 by=c("SITE_ID"="SITE_ID","VISNUM"="VISNUM","SAMPYEAR")) %>%  
  as.data.frame() %>% 
  filter(VISNUM==1)

# Add secchi depth to original data
tmp$SECCHI_MEAN..m.[is.na(tmp$SECCHI_MEAN..m.)] <- tmp$secchi_update[is.na(tmp$SECCHI_MEAN..m.)]
tmp <- select(tmp,-secchi_update)
return(tmp)
}

# Determine how many non-detects there are in the complete cases.
# only consider complete cases because these are the data used in the 
# analysis.
evaluateNonDetects <- function(inData,mdl) {
    
  tmp <- inData[complete.cases(inData[,c("SECCHI_MEAN..m.","DIP..mgP.L.","DIN..mgN.L.",
                                         "TN..mgN.L.","TP..mgP.L.","SUBREGIONS","CHLA..ug.L.")]),]
  
  tmp$TP[tmp$TP..mgP.L.<=mdl["TP"]] <- "nondetect"
  tmp$TP[tmp$TP..mgP.L.>mdl["TP"]] <- "detect"
  tmp$TN[tmp$TN..mgN.L.<=mdl["TN"]] <- "nondetect"
  tmp$TN[tmp$TN..mgN.L.>mdl["TN"]] <- "detect"
  tmp$DIN[tmp$DIN..mgN.L.<=mdl["DIN"]] <- "nondetect"
  tmp$DIN[tmp$DIN..mgN.L.>mdl["DIN"]] <- "detect"
  tmp$DIP[tmp$DIP..mgP.L.<=mdl["DIP"]] <- "nondetect"
  tmp$DIP[tmp$DIP..mgP.L.>mdl["DIP"]] <- "detect"
  
  tmp2 <- tmp %>% dplyr::select(TP,TN,DIN,DIP) %>% 
    pivot_longer(cols = c("TP","TN","DIN","DIP"),names_to="Var",values_to = "Flag") %>% 
    group_by(Var,Flag) %>% 
    summarize(n=n()) %>% 
    pivot_wider(names_from="Flag",values_from = "n",values_fill = list(n=0)) %>% 
    mutate(percent_non_detect = 100*nondetect/detect ) %>% 
    dplyr::select(Var,percent_non_detect)         
  
  return(tmp2)
  
}

# Uses saved mean and standard deviation from 2010 NCA data to centering variables
# in any data set.
centerParameters <- function(inData,centeringParameters) {
    cp <- as.matrix(centeringParameters[,c(2,3)])
    rownames(cp) <- centeringParameters[,1] %>% as.matrix()
    SDD.C <<- as.numeric(scale(log(inData$SECCHI_MEAN..m.), 
                               center = cp["logSD","mean"], scale = cp["logSD","sd"]))
    TN.C <<- as.numeric(scale(log(inData$TN..mgN.L.), 
                             center = cp["logTN","mean"], scale = cp["logTN","sd"]))
    TP.C <<- as.numeric(scale(log(inData$TP..mgP.L.), 
                             center = cp["logTP","mean"], scale = cp["logTP","sd"]))
    DIN.C <<- as.numeric(scale(log(inData$DIN..mgN.L.), 
                              center = cp["logDIN","mean"], scale = cp["logDIN","sd"]))
    DIP.C <<- as.numeric(scale(log(inData$DIP..mgP.L.), 
                              center = cp["logDIP","mean"], scale = cp["logDIP","sd"]))
}

# Create variables with parameters for Bayesian POLR model in JAGS
# pass NULL or "test" for testing, "final" for final fit
# returns parameters to global environment
# for "final" nChains is 3 instead of 1 for testing and
# numSavedSteps is 50000 instead of 10000.
jags_run_parameters <- function(phase=NULL) {
  
  # Number of steps to "tune" the samplers.
  adaptSteps <<- 3000          
  
  # Number of steps to "burn-in" the samplers.
  #changes from 1000 the initail setting
  burnInSteps <<- 5000    
  
  # Number of chains to run. 3 for final, 1 for anything else
  if (phase=="final") {
    nChains <<- 3 # nChains is 3 for final run
    numSavedSteps <<- 50000 # 50000 for the final run
  } else {
    nChains <<- 1
    numSavedSteps <<- 10000
  }
  
  # Number of steps to "thin" (1=keep every step).
  # records every 5th step in the MC 
  thinSteps <<- 5
  
  # Steps per chain.
  nIter <<- ceiling( ( numSavedSteps * thinSteps ) / nChains )
  
}

# Create variable importance plot from Random Forest Models as column graph using ggplot
plot_vImp <- function(inRF,title) {
  plt <- data.frame(Variables=as.character(variableLabels),
                    Importance=as.numeric(importance(inRF,type=1))) %>% 
    ggplot(.,aes(y=reorder(Variables,Importance),x=Importance/100,label=round(Importance,0))) + 
    geom_col(orientation="y",fill="light blue") +
    geom_text(nudge_x = 0.06) + 
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    labs(y="Mean Increase in Prediction Error if Permuted",x="Variable")+
    ggtitle(title)
  
  return(plt)
}

# Generate random subsets of data, returning variables "Evaluation" and "Model"
# to global environment.  By default, the evaluation fraction is 10% and 
# the model data has 90% of the data.
randomSubset <- function(inputData,evalFraction=0.1,seed=100) {
  set.seed(seed)
  # Generate a random list of row indices that specify a random sample, without r   
  # replacement, that is evalFraction (e.g., 10%) of the size of the input data set.  This 
  # is used for model evaluation, while the non-sampled data are used for model fit
  Sample <- sample(nrow(inputData),size=round(evalFraction*dim(inputData)[1]),
                   replace=FALSE)
  Evaluation <<- inputData[Sample,]  # select evaluation data
  Model <<- inputData[-Sample,] # use the data that isn't for evaluation to fit the model
}

# Build the naive or uninformed prior with specified number of subregions
build_naive_prior <- function(numRegions) {
  naive_prior <- list(mn_cut_pts = rep(0,3),
                      tau_cut_pts = rep(0.0001,3),
                      mn_alpha_SD = 0,
                      tau_alpha_SD = 0.0001,
                      mn_alpha_N = 0,
                      tau_alpha_N = 0.0001,
                      mn_alpha_P = 0,
                      tau_alpha_P = 0.0001,
                      mn_alpha_DIN = 0,
                      tau_alpha_DIN = 0.0001,
                      mn_alpha_DIP = 0,
                      tau_alpha_DIP = 0.0001,
                      mn_alpha_SubR = rep(0,numRegions),
                      tau_alpha_SubR = rep(0.0001,numRegions),
                      mn_s = 0,
                      tau_log_s=0.0001
  )
  return(naive_prior)
}

build_informed_prior <- function(posterior,numRegions) {
  updatedPrior <- list(mn_cut_pts = posterior[1:3,1],
                       tau_cut_pts = posterior[1:3,2]^-2,
                       mn_alpha_SD = posterior[8,1],
                       tau_alpha_SD = posterior[8,2]^-2,
                       mn_alpha_N = posterior[6,1],
                       tau_alpha_N = posterior[6,2]^-2,
                       mn_alpha_P = posterior[7,1],
                       tau_alpha_P = posterior[7,2]^-2,
                       mn_alpha_DIN = posterior[4,1],
                       tau_alpha_DIN = posterior[4,2]^-2,
                       mn_alpha_DIP = posterior[5,1],
                       tau_alpha_DIP = posterior[5,2]^-2,
                       mn_alpha_SubR = posterior[9:(8+numSubregions),1],
                       tau_alpha_SubR = posterior[9:(8+numSubregions),2]^-2,
                       mn_s = log(posterior["s","mean"]),
                       tau_log_s = posterior["s",2]^-2
  )
  return(updatedPrior)
}

# This function is not used in the final paper, but was 
# part of our development process.  It allowed us to specify
# some of the variables with informed priors, and others 
# that remained uninformed.
build_mixed_prior <- function(posterior,numRegions,useNaive) {
  updatedPrior <- build_informed_prior(posterior,numRegions)
  if (useNaive["cutpts"]) {
    updatedPrior$mn_cut_pts <- rep(0,3)
    updatedPrior$tau_cut_pts  <-  rep(0.0001,3)
  }
  if (useNaive["sd"]) {
    updatedPrior$mn_alpha_SD <- 0
    updatedPrior$tau_alpha_SD  <-  0.0001
  }
  if (useNaive["N"]) {
    updatedPrior$mn_alpha_N <- 0
    updatedPrior$tau_alpha_N  <-  0.0001
  }
  if (useNaive["P"]) {
    updatedPrior$mn_alpha_P <- 0
    updatedPrior$tau_alpha_P  <-  0.0001
  }
  if (useNaive["DIN"]) {
    updatedPrior$mn_alpha_DIN <- 0
    updatedPrior$tau_alpha_DIN  <-  0.0001
  }
  if (useNaive["DIP"]) {
    updatedPrior$mn_alpha_DIP <- 0
    updatedPrior$tau_alpha_DIP  <-  0.0001
  }
  if (useNaive["SubR"]) {
    updatedPrior$mn_alpha_SubR <- rep(0,numRegions)
    updatedPrior$tau_alpha_SubR  <-  rep(0.0001,numRegions)
  }

return(updatedPrior)
}
 
# 
#  Functions from original program created by Farnaz
expected <- function(x, c1.5, c2.5, c.3.5, sigma){
  p1.5 <- invlogit((x-c1.5)/sigma)
  p2.5 <- invlogit((x-c2.5)/sigma)
  p3.5 <- invlogit((x-c3.5)/sigma)
  return((1*(1-p1.5)+2*(p1.5-p2.5)+3*(p2.5-p3.5)+4*p3.5))
  # return((1*(1-p1.5)+2*(p1.5-p2.5)+3*p2.5))
}

# for plotting logistic regression model
jitter.binary <- function(a, jitt=.05, up=1){
  up*(a + (1-2*a)*runif(length(a),0,jitt))
}
logit <- function(x) return(log(x/(1-x)))
invlogit <- function(x) return(1/(1+exp(-x)))

# Sample Posterior distribution using JAGS

# coda.samples is a wrapper function for jags.samples which sets a trace monitor for all requested nodes, updates the model, and coerces the output to a single mcmc.list object.
#  objects is called "coastal_jags" (i.e., no ".R")
#
# jags.samples extracts random samples from the posterior distribution of the parameters of # a jags model
#
# coastal_jags = jags model object.
# parameters = the names of the variables that are monitored
# n.iter = the number of iterations to monitor.  Use 10000 during development, 
# change to 50,000 for the final run
sample_posterior_distribution <- function(inputModel,parameters,iterations) {
  coda_output <- coda.samples(inputModel, parameters, n.iter=iterations)
  return(coda_output)
}

# Combine Markov chains for model evaluation
# pass the coda object with samples of posterior distribution and the 
# number of chains.  coda_input is a MCMC list object.  Returned object
# is an array
combine_chains <- function(coda_input,nChains) {
  # 3 Chains Combined
  combined_coda <- NULL
  for (i in 1:nChains) {
    combined_coda <- rbind(combined_coda,          
                           coda_input[[i]])
  }
  return(combined_coda)
}

# Build a summary of the coefficients using the combined Markov Chains
build_coefficient_summary <- function(combined_coda) {
  
  # Build a summary of the coefficients
  numVars <- ncol(combined_coda)
  Coeff.Summary <- matrix(NA, numVars, 4)
  for (i in 1:numVars){ 
    Coeff.Summary[i,] <- cbind(mean(combined_coda[,i])
                               , sd(combined_coda[,i])
                               , quantile(combined_coda[,i], c(0.025), type = 1)
                               , quantile(combined_coda[,i], c(0.975), type = 1))
  }
  colnames(Coeff.Summary) <- cbind("mean", "sd", "2.5%", "97.5%")
  rownames(Coeff.Summary ) <-colnames(as.mcmc(combined_coda))
  return(Coeff.Summary)
}

# Build Coefficient (Alpha) Matrix
build_coefficient_matrix <- function(Coeff.Summary) {
  # check number of rows in Coeff.Coastal.Summary ... depends on subregions
  numRegions <- nrow(Coeff.Summary)-9
  Alpha <- rbind(Coeff.Summary["alpha_SD",]
                 , Coeff.Summary["alpha_N",]
                 , Coeff.Summary["alpha_P",]
                 , Coeff.Summary["alpha_DIN",]
                 , Coeff.Summary["alpha_DIP",]
                 , Coeff.Summary[seq(9,8+numRegions),])
  return(Alpha)
}

# Build region matrix (part of design matrix) using data and
# provided region variable.
build_region_matrix <- function(inputData,regionVariable) {
  regionMatrix <- matrix(0, dim(inputData)[1],   
                         length(levels(inputData[,regionVariable])))
  for(j in 1:dim(inputData)[1]){
    for(i in 1:length(levels(inputData[,regionVariable]))){
      if (factor(inputData[j, regionVariable])==levels(inputData[,regionVariable])[i]) 
      { regionMatrix[j,i] <- 1 }
    }
  }
  return(regionMatrix)
}

# Calculate Predicted Class if TSI is already calculated
findPredictedClass.TSI <- function(inData,Coeff.Summary) {
  Pred.CatAll <- cut(inData$TSI,
             breaks=c(-Inf,Coeff.Summary[1:3,"mean"]/Coeff.Summary["s","mean"],Inf),
             labels=c("Oligo", "Meso", "Eu", "Hyper"))
  return(Pred.CatAll)  
}

# Calculate Predicted Class starting with data
# this function is not used in the final analysis
findPredictedClass <- function(inData,regionMatrix,Coeff.Summary,returnTSI=FALSE) {
  # Center, scale, and log transform the evaluation data
  Eval.SDD.C <- as.numeric(scale(log(inData$SECCHI_MEAN..m.), 
                                 center = TRUE, scale =   TRUE))
  Eval.TN.C <- as.numeric(scale(log(inData$TN..mgN.L.),
                                center = TRUE, scale = TRUE))
  Eval.TP.C <- as.numeric(scale(log(inData$TP..mgP.L.), 
                                center = TRUE, scale = TRUE))
  Eval.DIN.C <- as.numeric(scale(log(inData$DIN..mgN.L.), 
                                 center = TRUE, scale = TRUE))
  Eval.DIP.C <- as.numeric(scale(log(inData$DIP..mgP.L.), 
                                 center = TRUE, scale = TRUE))
  # Evaluation predictors
  Eval.Predictors <- cbind(Eval.SDD.C, Eval.TN.C, Eval.TP.C, 
                           Eval.DIN.C, Eval.DIP.C, regionMatrix)
  
  # Build Coefficient Matrix
  Alpha <- build_coefficient_matrix(Coeff.Summary)
  # Calculate matrix product of Eval.Predictors and mean model Coefficients.
  predict.EvaluationAll <- Eval.Predictors %*% Alpha[,"mean"] 
  # Remove NA's
  predict.EvaluationAll <- predict.EvaluationAll[!is.na(predict.EvaluationAll)]
  
  Predict.CatAll <-  vector(length = length(predict.EvaluationAll))
  
  C <- Coeff.Summary[c("C[1]","C[2]","C[3]"),"mean"]
  
  for (i in 1:length(predict.EvaluationAll)) {
    if (predict.EvaluationAll[i]< C[1]) { Predict.CatAll[i] <- "Oligo" }
    else if (predict.EvaluationAll[i] >= C[1] && predict.EvaluationAll[i] < C[2])
    { Predict.CatAll[i] <- "Meso" }
    else if (predict.EvaluationAll[i] >= C[2] && predict.EvaluationAll[i] < C[3]) 
    {Predict.CatAll[i] <- "Eu"}
    else if (predict.EvaluationAll[i] >= C[3]) { Predict.CatAll[i] <- "Hyper"}
    else {Predict.CatAll[i] <- NA}
  }
  
  Pred.CatAll <- factor(Predict.CatAll, 
                        levels=c("Oligo", "Meso", "Eu", "Hyper"), ordered=TRUE)
  if (returnTSI) {
    return(predict.EvaluationAll)
  } else {
  return(Pred.CatAll)
  }
}

## Several functions are used to create plots for visualizing results of POLR model
# Generate POLR scatterplot
plot_POLR_results <- function(Model,Alpha,Coeff.Summary,regionVariable) {
  
  # build region matrix for Model Data
  Model_SubRegion <- build_region_matrix(Model,regionVariable)
  
  beta <- Alpha[,"mean"]
  C <- Coeff.Summary[c("C[1]","C[2]","C[3]"),"mean"]
  c1.5 <- C[1]
  c2.5 <- C[2]
  c3.5 <- C[3]
  sigma <- Coeff.Summary["s","mean"]
  
  par(mar=c(3,3,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
  
  plot(0, 0, xlim=c(-600,600), ylim=c(1,4), xlab="TSI", ylab="TS",
       type="n", axes=F)
  
  axis(1)
  axis(2, at=1:4, labels=c("Oligo","Meso","Eutro", "Hyper"), las=1)
  lines(rep(c1.5, 2), c(1,2))
  lines(rep(c2.5, 2), c(2,3))
  lines(rep(c3.5, 2), c(3,4))
  curve(expected(x, c1.5, c2.5, c3.5, sigma), add=TRUE)
  
  with(Model, points(cbind(SDD.C, TN.C, TP.C, DIN.C, DIP.C, Model_SubRegion)%*%beta,
                     jitter.binary(as.numeric(ordered(Model[,"TS_Chla_Q"]))),  col="azure4"))
  
}

generate_graphic_model <- function(inData,Coeff.Summary,regionVariable) {
  
  Alpha <- build_coefficient_matrix(Coeff.Summary)
  beta_AllVar <-  Alpha[,"mean"]
  kappa_AllVar <- Coeff.Summary[c("C[1]","C[2]","C[3]"),"mean"]
  c <- kappa_AllVar
  
  c1.5_AllVar <- kappa_AllVar[1]
  c2.5_AllVar <- kappa_AllVar[2]
  c3.5_AllVar <- kappa_AllVar[3]
  sigma_AllVar <- Coeff.Summary["s","mean"]
  
  SDD.C <- as.numeric(scale(log(inData$SECCHI_MEAN..m.), center = TRUE, scale = TRUE))
  TN.C <- as.numeric(scale(log(inData$TN..mgN.L.), center = TRUE, scale = TRUE))
  TP.C <- as.numeric(scale(log(inData$TP..mgP.L.), center = TRUE, scale = TRUE))
  DIN.C <- as.numeric(scale(log(inData$DIN..mgN.L.), center = TRUE, scale = TRUE))
  DIP.C <- as.numeric(scale(log(inData$DIP..mgP.L.), center = TRUE, scale = TRUE))
  
  regionMatrix <- build_region_matrix(inData,regionVariable)
  
  X <- cbind(SDD.C, TN.C, TP.C, DIN.C, DIP.C, regionMatrix)
  
  TSI <- X %*% beta_AllVar
  TSI <- TSI[!is.na(TSI)]
  
  c <- Coeff.Summary[c("C[1]","C[2]","C[3]"),"mean"]
  # se of kappas
  se.c <-  Coeff.Summary[c("C[1]","C[2]","C[3]"),"sd"]
  Ibcg <- seq(range(TSI)[1],range(TSI)[2], length.out = 100)
  pA <- invlogit((kappa_AllVar[1] - Ibcg)/sigma_AllVar)
  pB <- invlogit((kappa_AllVar[2] - Ibcg)/sigma_AllVar) -  invlogit((kappa_AllVar[1] -
                                                                       Ibcg)/sigma_AllVar)
  pC <- invlogit((kappa_AllVar[3] - Ibcg)/sigma_AllVar) -  invlogit((kappa_AllVar[2] -
                                                                       Ibcg)/sigma_AllVar)
  pNA <- 1.0 - invlogit((kappa_AllVar[3] - Ibcg)/sigma_AllVar)
  
  # Plot the same graph using ggplot
  alllines <- rbind(data.frame(TSI=Ibcg,Prob=pA,class="Oligotrophic"),
                    data.frame(TSI=Ibcg,Prob=pB,class="Mesotrophic"),
                    data.frame(TSI=Ibcg,Prob=pC,class="Eutrophic"),
                    data.frame(TSI=Ibcg,Prob=pNA,class="Hypertrophic"))
  
  plt <- ggplot()+
    geom_rect(aes(xmin=c(c1.5_AllVar-se.c[1],c2.5_AllVar-se.c[2],c3.5_AllVar-se.c[3]),
                  xmax=c(c1.5_AllVar+se.c[1],c2.5_AllVar+se.c[2],c3.5_AllVar+se.c[3]),
                  ymin=rep(0,3),ymax=rep(1,3)),fill="lightblue1")+
    geom_vline(xintercept=c)+
    geom_line(data=alllines,aes(x=TSI,y=Prob,color=class,linetype=class),size=1)+
    scale_color_brewer(name="Trophic Class", palette="Set2") +
    guides(linetype=FALSE)+
    theme_classic()+
    theme(legend.position="top")+
    scale_y_continuous(expand=c(0,0))+
    labs(x="Trophic State Index",y="Trophic Class Probability")
  return(plt)
}

# Calculates the probability of a new observation being in one of the trophic 
# state categories given its TSI value and coefficients of the POLR model.
# This was used to apply the analysis to Boston Harbor
# can could be used for other data.
calculate_ts_probabilities <- function(Coeff.Summary,TSIinput) {
  
  Alpha <- build_coefficient_matrix(Coeff.Summary)
  sigma <- Coeff.Summary["s","mean"] 
  c <- Coeff.Summary[c("C[1]","C[2]","C[3]"),"mean"]/sigma
  se.c <-  Coeff.Summary[c("C[1]","C[2]","C[3]"),"sd"]
  Ibcg <- TSIinput
  pOligotrophic <- invlogit((c[1] - Ibcg))
  pMesotrophic <- invlogit((c[2] - Ibcg)) - invlogit((c[1] - Ibcg))
  pEutrophic <- invlogit((c[3] - Ibcg)) - invlogit((c[2] - Ibcg))
  pHypertrophic <- 1.0 - invlogit((c[3] - Ibcg))
  
  probs <- data.frame(pOligotrophic,pMesotrophic,pEutrophic,pHypertrophic)
  return(probs)
}

# Used to apply ROS method to Boston Harbor data
impute_non_detects <- function(inData,concVar,flagVar) {
  tmp <- inData[,c(concVar,flagVar)]
  names(tmp) <- c("conc","flag")
  tmp$censored <- tmp$flag=="<"
  tmp$order <- seq.int(nrow(tmp))
  tmp <- tmp[order(tmp$conc),]
  
  myROS <- ros(tmp$conc,tmp$censored)
  myROS.df <- as.data.frame(myROS)
  
  tmp <- cbind(tmp %>% select(-censored),myROS.df)
  
  tmp <- tmp[order(tmp$order),]
  
  return(tmp)
}

impute_time_window <- function(inData,startYear,endYear){
  #subset the data
  windowData <- inData %>% filter(year(TimeStamp)>=startYear & year(TimeStamp)<=endYear)
  # fit modeled TDN values using ROS
  tmp <- impute_non_detects(windowData,"TDN","TDN.Flag")
  windowData$TDN.modeled <- tmp$modeled
  # fit modeled PN values using ROS
  tmp <- impute_non_detects(windowData,"PN","PN.Flag")
  windowData$PN.modeled <- tmp$modeled
  # fit modeled NH4 values using ROS
  tmp <- impute_non_detects(windowData,"NH4","NH4.Flag")
  windowData$NH4.modeled <- tmp$modeled
  # fit modeled NOx values using ROS
  tmp <- impute_non_detects(windowData,"NOx","NOx.Flag")
  windowData$NOx.modeled <- tmp$modeled
  # fit modeled TDP values using ROS
  tmp <- impute_non_detects(windowData,"TDP","TDP.Flag")
  windowData$TDP.modeled <- tmp$modeled
  # fit modeled PP values using ROS
  tmp <- impute_non_detects(windowData,"PP","PP.Flag")
  windowData$PP.modeled <- tmp$modeled
  # fit modeled PO4 values using ROS
  tmp <- impute_non_detects(windowData,"PO4","PO4.Flag")
  windowData$PO4.modeled <- tmp$modeled
  return(windowData)
}
