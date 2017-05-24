## Step.1 
## In this step, we completed the following works:
## > sampling from source data,
## > pretreatment for the aprori, 
## > get association rules,
## > re-process the results.



library("arules")

## This part is the sampling from source data.
print('Waiting ...')
AllG <- read.csv(
  "IN/weishengsu_expression_change.csv",header = F,stringsAsFactors = FALSE)
AllG <- data.frame(t(AllG))
AllG <- AllG[,-1]
G100 <- sample(AllG,100)
G100 <- data.frame(t(G100))
write.csv(G100,TEST_DATA,row.names = F)


## This part is the pretreatment of association rules.

jRow <- length(rownames(AllG))
jDim <- 2
jFin <- paste(getwd(),TEST_DATA,sep = '/')
jFout <- paste(getwd(),"EX/rules/Mult_data_100.txt",sep = '/')
system("javac preTreatment.java")
system(paste("java preTreatment",jRow,jDim,jFin,jFout,sep = ' '))




##  This part use the encapsulated function 'apriori' to complete the analysis of association rules.

tr <-  read.transactions("EX/rules/Mult_data_100.txt", format="basket", sep = ",")
rules = apriori(tr,parameter = list(support = 0.15,confidence = 0.45,minlen = 2,maxlen = 2)) #pre 0.3  0.6
summary(rules)     
outdata <- as(rules,"data.frame")   
write.table(outdata,"EX/rules/Result.csv",sep = ",")   


# This part is the further processing of the results after association rules.
jRules <- length(rules)
jFin <- paste(getwd(),"EX/rules/Result.csv",sep = '/')
jFout <- paste(getwd(),"EX/rules/sum_data.csv",sep = '/')
system("javac furTreatment.java")
system(paste("java furTreatment",jRules,jFin,jFout,sep = ' '))


 