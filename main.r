
## Step.0
## These packages should be initialized at first.
## In addition , they only need to be initialized during the first execution of this program.
## If these packages do not exist,please copy these files to R's library manually.

# install.packages('Rgraphviz', dependencies = TRUE) 
# install.packages('graph',depend=TRUE) 
# install.packages('arules', depend=TRUE)
# install.packages('bnlearn', depend=TRUE)
# install.packages('igraph', depend=TRUE)
# install.packages('BiocGenerics', depend=TRUE)
# install.packages('grid', depend=TRUE)
 

## Set the current working directory.
setwd(getwd())



## Step.1 
## In this step, we completed the following works:
## > sampling from source data,
## > pretreatment for the aprori, 
## > get association rules,
## > re-process the results.
source("Sampling&Aprori.R")




## Step.2
## To confirm whether the source data need to sort by the association rules, 
## if not, skip to the third step.
Sort1_Or_NotSort2 <- TRUE
if(Sort1_Or_NotSort2){
    source("SortbyRules.R")
}




## Step.3
## Implementation of TPDA algorithm, the algorithm input in two cases:
## > One is the result of direct sampling of source data.
## > The other is the result sorted by the association rule after sampling.
## Then make calculations.
source("TPDA_Function_5_4.R")




## Step.4
## Use nine algorithms to get the results.
source("Nine_Algorithms.R")




## Step.5[1]
## Call the function :TPDA_Algorithm
## For the specified experimental data, 
## here we use different threshold, dispersion, discrete method to get a set of experimental results of TPDA.
## The result is stored in : OUT/Multi_Experiment_Result.csv by default.
TPDA_Output <- TPDA_Algorithm(4 ,0.025,0.025,0.025,'interval',"IN/Gene100_6.csv" )
write.csv(TPDA_Output,"OUT/Multi_Experiment_Result.csv",row.names = F)



## Step.5[2]
## Call the function :NineAlgorithmTest.
## The nine algorithms are used to experiment with the data, 
## including gs,hc,iamb,mmpc,rsmax2,tabu,fastiamb,interiamb,mmhc.
## The result is stored in : OUT/Nine_Result.csv by default.
NineAlgorithmTest('interval',4,"IN/Gene100_6.csv" ,"OUT/Nine_Result.csv")



