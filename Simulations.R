## ------------------------ WORKSPACE PREPARATION AND LIBRARIES--------------------------------

# Clean the environment
rm(list=ls())

# Change the working directory path and loading the required packages
setwd("C:/Users/nicolas/Dropbox/Utrecht University/Thesis/Simulated Data")
library(MASS)
library(robust)
library(Matrix)
library(svMisc)


## ------------------ FUNCTION TO GENERATE LINEAR REGRESSION DATASETS -------------------------

# Function that generates the datasets
genNormData<-function(betas,xCorrels,numPredictors,numObservations,numDatasets,rSquared)
{
  datasets<-array(0,dim=c(numObservations,numPredictors+1,numDatasets))
  x_means<-rep(0,numPredictors)
  for(i in 1:numDatasets)
  {
    x_matrix<-matrix(0,nrow=numObservations,ncol=numPredictors)
    x_matrix<-mvrnorm(numObservations,mu=x_means,Sigma=xCorrels,empirical=FALSE)
    pred_y<-apply(t(t(cbind(rep(1,numObservations),x_matrix))*betas),MARGIN=1,FUN=sum)
    res_var<-sum((betas[-1]%*%t(betas[-1])*xCorrels))*(1-r_squared)/r_squared
    samp_Y<-pred_y+rnorm(numObservations,0,sqrt(res_var))
    datasets[,,i]<-cbind(x_matrix,samp_Y)
  }
  colnames(datasets)<-c(paste("x",seq(1,numPredictors),sep=""),"y")
  rownames(datasets)<-c(1:numObservations)
  
  return(datasets)
  
}

## ----------- FUNCTION TO GENERATE LR DATASETS WITH HETEROSCEDASTIC RESIDUALS ----------------

# Function that generates the datasets with heteroscedastic residuals depending on X1, type is 1 if it is increasing and -1 if it is decreasing
genHetData<-function(alphas,betas,xCorrels,numPredictors,numObservations,numDatasets,type)
{
  datasets<-array(0,dim=c(numObservations,numPredictors+1,numDatasets))
  x_means<-rep(0,numPredictors)
  for(i in 1:numDatasets)
  {
    x_matrix<-matrix(0,nrow=numObservations,ncol=numPredictors)
    x_matrix<-mvrnorm(numObservations,mu=x_means,Sigma=xCorrels,empirical=FALSE)
    pred_y<-apply(t(t(cbind(rep(1,numObservations),x_matrix))*betas),MARGIN=1,FUN=sum)
    hetError<-rnorm(numObservations,mean=0,sd=exp(alphas[1]+type*alphas[2]*x_matrix[,1]))
    samp_Y<-pred_y+hetError
    datasets[,,i]<-cbind(x_matrix,samp_Y)
  }
  colnames(datasets)<-c(paste("x",seq(1,numPredictors),sep=""),"y")
  rownames(datasets)<-c(1:numObservations)
  
  return(datasets)
  
}

## -------------- FUNCTION TO GENERATE OUTLIERS IN THE X SPACE USING REPLACEMENT ------------------

# cofefsVec is a vetcor that indicates which X's should be considered
outliers_X<-function(data,nOutliers,coeffsVec)
{
  n_datasets<-dim(data)[3]
  n_var<-dim(data)[2]
  n_obs<-dim(data)[1]
  returnDatasets<-array(0,dim=c(n_obs,n_var,n_datasets))
  
  for(i in 1:n_datasets)
  {
    activeData<-data[,,i]
    quant_25<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.25))
    quant_75<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.75))
    IQR<-quant_75-quant_25
    
    assignedCols<-sample(coeffsVec,nOutliers,replace=TRUE)
    
    for(k in 1:n_var-1)
    {
      numOutsPred<-sum(assignedCols==k)
      if(numOutsPred>0)
      {
        sampleRows<-sample(c(1:n_obs),numOutsPred,replace = FALSE)
        imputedVals<-runif(numOutsPred,quant_25[k]-3*IQR[k],quant_25[k]-1.5*IQR[k])
        activeData[sampleRows,k]<-imputedVals
      }
    }
    
    returnDatasets[,,i]<-activeData
  }
  
  colnames(returnDatasets)<-c(paste("x",seq(1,n_var-1),sep=""),"y")
  rownames(returnDatasets)<-c(1:n_obs)
  
  return(returnDatasets)
}

## -------------- FUNCTION TO GENERATE OUTLIERS IN THE X SPACE USING REPLACEMENT ------------------

outliers_Y<-function(data,nOutliers)
{
  n_datasets<-dim(data)[3]
  n_var<-dim(data)[2]
  n_obs<-dim(data)[1]
  returnDatasets<-array(0,dim=c(n_obs,n_var,n_datasets))
  
  for(i in 1:n_datasets)
  {
    activeData<-data[,,i]
    quant_25<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.25))
    quant_75<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.75))
    IQR<-quant_75-quant_25
    
    sampleRows<-sample(c(1:n_obs),nOutliers,replace = FALSE)
    imputedVals<-runif(nOutliers,quant_75[n_var]+1.5*IQR[n_var],quant_75[n_var]+3*IQR[n_var])
    activeData[sampleRows,n_var]<-imputedVals
    
    returnDatasets[,,i]<-activeData
  }
  
  colnames(returnDatasets)<-c(paste("x",seq(1,n_var-1),sep=""),"y")
  rownames(returnDatasets)<-c(1:n_obs)
  
  return(returnDatasets)
}

## --------------- FUNCTION TO GENERATE OUTLIERS IN THE XY SPACE USING REPLACEMENT ------------------

# coeffsVec is a vetcor that indicates which X's should be considered, if "lowHigh" is true, the outliers are created using lo Xs and high Ys, otherwise thy are created with high Xs and high Ys
outliers_XY<-function(data,nOutliers,coeffsVec,lowHigh)
{
  n_datasets<-dim(data)[3]
  n_var<-dim(data)[2]
  n_obs<-dim(data)[1]
  returnDatasets<-array(0,dim=c(n_obs,n_var,n_datasets))
  
  for(i in 1:n_datasets)
  {
    activeData<-data[,,i]
    quant_25<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.25))
    quant_75<-apply(X=activeData,MARGIN=2,FUN=quantile, probs=c(0.75))
    IQR<-quant_75-quant_25
    
    assignedCols<-sample(coeffsVec,nOutliers,replace=TRUE)
    
    for(k in 1:n_var-1)
    {
      numOutsPred<-sum(assignedCols==k)
      if(numOutsPred>0)
      { 
        sampledRows<-sample(c(1:n_obs),numOutsPred,replace = FALSE)
        
        if(lowHigh==TRUE)
        {
          imputed_X<-runif(numOutsPred,quant_25[k]-3*IQR[k],quant_25[k]-1.5*IQR[k])  
        }
        else
        {
          imputed_X<-runif(numOutsPred,quant_75[k]+1.5*IQR[k],quant_75[k]+3*IQR[k])  
        }
        
        imputed_Y<-runif(numOutsPred,quant_75[n_var]+1.5*IQR[n_var],quant_75[n_var]+3*IQR[n_var])
        activeData[sampledRows,k]<-imputed_X
        activeData[sampledRows,n_var]<-imputed_Y
      }
    }
    
    returnDatasets[,,i]<-activeData
  }
  
  colnames(returnDatasets)<-c(paste("x",seq(1,n_var-1),sep=""),"y")
  rownames(returnDatasets)<-c(1:n_obs)
  
  return(returnDatasets)
}

## ----------------------- PARAMETERS SETUP AND DATA GENERATION -------------------------------

# Parameters for the data generation
n_obs<-100
n_datasets<-1000
n_predictors<-3
betas<-c(2,3,4,5)
r_squared<-0.5
alphas<-c(1.4,0.3)
x_correlations<-matrix(c(1,0.11,0.08,0.11,1,0.14,0.08,0.14,1),nrow=3,ncol=3)

# Generating the datasets and performing a lm for a specific dataset
normalDatasets<-genNormData(betas,x_correlations,n_predictors,n_obs,n_datasets,r_squared)
heteroscDatasets<-genHetData(alphas,betas,x_correlations,n_predictors,n_obs,n_datasets,1)

# Generating Outlier Dataset Low X and high Y
coeffs<-c(1)
nOutliers<-1

lowXhighY<-outliers_XY(normalDatasets,nOutliers,coeffs,lowHigh = TRUE)

## ------------------------ FUNCTIONS TO PERFORM LIN REG AND ROB REG --------------------------

# Function that performs linear regression on datasets and extracts the coefficients
resultsLR<-function(datasets)
{
  n_datasets<-dim(datasets)[3]
  n_var<-dim(datasets)[2]
  resultsCoeff<-matrix(0,nrow=n_datasets,ncol=n_var)
  resultsSE<-matrix(0,nrow=n_datasets,ncol=n_var)
  for(i in 1:n_datasets)
  {
    activeData<-data.frame(datasets[,,i])
    for(j in 1:n_var)
    {
      model<-lm(activeData$y~activeData$x1+activeData$x2+activeData$x3)
      resultsCoeff[i,j]<-coefficients(model)[j] 
      resultsSE[i,j]<-coefficients(summary(model))[j,2]
    }
    
  }
  colnames(resultsCoeff)<-c("Intercept","Beta1","Beta2","Beta3")
  rownames(resultsCoeff)<-c(1:n_datasets)

  colnames(resultsSE)<-c("SE_Int","SE_B1","SE_B2","SE_B3")
  rownames(resultsSE)<-c(1:n_datasets)

  return(list(resultsCoeff,resultsSE))
}

# Function that performs robust regression (lmRob) on datasets and extracts the coefficients
resultsRobReg<-function(datasets)
{
  n_datasets<-dim(datasets)[3]
  n_var<-dim(datasets)[2]
  resultsCoeff<-matrix(0,nrow=n_datasets,ncol=n_var)
  resultsSE<-matrix(0,nrow=n_datasets,ncol=n_var)
  for(i in 1:n_datasets)
  {
    activeData<-data.frame(datasets[,,i])
    for(j in 1:n_var)
    {
      model<-lmRob(activeData$y~activeData$x1+activeData$x2+activeData$x3)
      resultsCoeff[i,j]<-coefficients(model)[j] 
      resultsSE[i,j]<-coefficients(summary(model))[j,2]
    }
    
  }
  colnames(resultsCoeff)<-c("Intercept","Beta1","Beta2","Beta3")
  rownames(resultsCoeff)<-c(1:n_datasets)
  
  colnames(resultsSE)<-c("SE_Int","SE_B1","SE_B2","SE_B3")
  rownames(resultsSE)<-c(1:n_datasets)
  
  return(list(resultsCoeff,resultsSE))
}



## -------------------- FUNCTION TO COMPUTE COVERAGE OF ESTIMATES------------------------------

estimatesCoverage<-function(coefResults, SDResults, popCoeffs)
{
  resultsCoverage<-data.frame(matrix(0,nrow=2,ncol=4))
  colnames(resultsCoverage)=c("Intercept","Beta1","Beta2","Beta3")
  rownames(resultsCoverage)=c("Cases","Cov. Prob.")
  for(i in 1:4)
  {
    meanVal<-popCoeffs[i]
    lowB<-coefResults[,i]+qt(0.025,df=n_obs-4)*SDResults[,i]
    upB<-coefResults[,i]+qt(0.975,df=n_obs-4)*SDResults[,i]
    
    resultsCoverage[1,i]<-sum((meanVal>=lowB)&(meanVal<=upB))
    resultsCoverage[2,i]<-resultsCoverage[1,i]/nrow(coefResults)
    
  }
  
  return(resultsCoverage)
}


## ------------------- APPLYING MODELS TO THE SIMULATED DATA ----------------------------------

# Executing linear regression on normal datasets and extracting the coefficent results and SE's
resultsLR_NormalData<-resultsLR(lowXhighY)
lrCoeffs_normalData<-data.frame(resultsLR_NormalData[1])
lrCoeffsSE_normalData<-data.frame(resultsLR_NormalData[2])

# Executing robust regression on normal datasets and extracting the coefficent results and SE's
resultsRR_NormalData<-resultsRobReg(lowXhighY)
robRegCoeffs_normalData<-data.frame(resultsRR_NormalData[1])
robRegCoeffsSE_normalData<-data.frame(resultsRR_NormalData[2])

# Computing coverage of the CI's for the 2 types of regression
estimatesCoverage(lrCoeffs_normalData,lrCoeffsSE_normalData,betas)
estimatesCoverage(robRegCoeffs_normalData,robRegCoeffsSE_normalData,betas)

## ------------------------------------- PLOTS ------------------------------------------------

require(ggplot2)

allcoeffs<-data.frame(cbind(lrCoeffs_normalData,robRegCoeffs_normalData))
names(allcoeffs)<-c("b0ols","b1ols","b2ols","b3ols","b0rr","b1rr","b2rr","b3rr")


a<-ggplot(data=allcoeffs, aes(x=b0ols))+
  geom_histogram(aes(x=b3ols),color="tomato",fill="tomato3",alpha=0.5)+
  geom_histogram(aes(x=b3rr),color="mediumseagreen",fill="lightseagreen",alpha=0.5)+
  labs(x=expression(Beta[3]))
  
setwd("C:/Users/nicolas/Dropbox/Utrecht University/Semester III - 201702/Research Seminar/Research Project/Graphs")  

pdf(file="beta3outliers.pdf", width=6, height=6)
a
dev.off()  


quantiles<-quantile(allcoeffs$b3ols,c(0.025,0.55,0.975))
b<-ggplot(data=allcoeffs, aes(x=b0ols))+
  geom_histogram(aes(x=b3ols),color="darkblue",fill="blue",alpha=0.5)+
  geom_vline(xintercept=quantiles[1],color="red",linetype="dashed")+
  geom_vline(xintercept=quantiles[3],color="red",linetype="dashed")+
  geom_vline(xintercept=quantiles[2],color="black",linetype="dashed",size=1.2)+
  labs(x=expression(Beta[3]))

setwd("C:/Users/nicolas/Dropbox/Utrecht University/Semester III - 201702/Research Seminar/Research Project/Graphs")  

bmp(file="coverage.bmp", width=250, height=350)
b
dev.off()  


vals_x<-c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
weighted_res<-c(0,0,0.3,1,3,6,3,1,0.3,0,0)

residuals_fun<-data.frame(cbind(vals_x,weighted_res))

ggplot(data=residuals_fun, aes(x=vals_x))+
  geom_line(aes(x=vals_x,y=weighted_res),color="tomato3",size=1.2)+
  labs(x="residual",y="weight")


