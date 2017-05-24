
## Step.4
## Use nine algorithms to get the results.
## The nine algorithms are used to experiment with the data, 
## including gs,hc,iamb,mmpc,rsmax2,tabu,fastiamb,interiamb,mmhc.

library(igraph)
library(bnlearn)
library(Rgraphviz) 

## Function : NineAlgorithmTest.
## methodd : discrete method.
## discret : discrete value.
## sourcee : input data / Source data.
## outfile : the path of output file.
NineAlgorithmTest <- function(methodd,discret,sourcee,outfile){
  
    methodd <- methodd
    sourcee <- sourcee
    #Sort1_Or_NotSort2 <- FALSE    
    LSdim <-discret  
    
    NotSort_File <-sourcee                              # without sorting
    Sort_File <- 'EX/rules/express100_after.csv'        # with sorting
    Nin_1 <- "IN/weishengsu_add_5_and_4.csv"
    Nout_1 <- outfile # "OUT/Nine_Result.csv"
    Nout_3 <- "EX/SampleNine.csv" 
    #Nout_2 <- "OUT/Nine_Result.txt"
 
    ## Determine the input source
    ## If sorted by association rules,here we use the data after sorted.
    ## If not,we directly use the sampling data.
    
    Vitamin9 <- read.csv(Nin_1,header = F,stringsAsFactors = FALSE)[-1,]
    if(Sort1_Or_NotSort2 == TRUE){
        Data100 <- read.csv( Sort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    }else{
        Data100 <- read.csv( NotSort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    }
    SampleA <- rbind(Vitamin9,Data100) #100+9
    write.csv(SampleA,Nout_3,row.names = F)
   
    
    ## Start the TPDA.
    
    zp <- read.csv(Nout_3,header = F,stringsAsFactors = FALSE)[-1,]
    row.names(zp) <- zp[,1]
    zp <- zp[-1]
    mydata <- data.frame(zp)  
    tq <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
    rownames(tq) <- rownames(mydata)
    colnames(tq) <- colnames(mydata)
    mydata <- tq
    mydata <- data.frame(t(mydata))
    
    tpN <- discretize(mydata,method = methodd,breaks = discret ) 
    rownames(tpN) <- rownames(mydata)
    colnames(tpN) <- colnames(mydata)
    mydata <- data.frame((tpN))
    #write.csv(mydata,Nout_1,row.names = F,append=TRUE,col.names = F,sep=' ')
     
    
    ## Different from the TPDA's comparison method, here directly use the results of 
    ## 9 algorithms to make comparison to find the amount of relationships to meet the conditions.
    
    ## Specific phenotype and gene.
    TDBX5 <- c('LUT','ZEA','Bcry','AC','BC')
    TDJY4 <- c('GRMZM2G012966','GRMZM2G152135','GRMZM2G300348','GRMZM2G108457')
    CheckInArr <- function(nm,arr){
        cnt <- length(arr)
        for(ci in 1:cnt){
            if(nm[,1] == arr[ci]||nm[,2] == arr[ci]){
                return (TRUE)
            }
        }
        return (FALSE)
    }
    
    ToArrayStr <- function(df){
        sstr <- c()
        longstr <- c()
        for(ssi in 1:length(df[,1])){
            sstr <- c(sstr,paste('[',df[ssi,1],',',df[ssi,2],']',sep=''))
            longstr <- paste(longstr,sstr[ssi])
        }
        return (longstr)
    }
    
    
   
    ## Function : OPinfo.
    ## Implements the execution of nine algorithms, 
    ## as well as the analysis and comparison of the results of each algorithm.
    OPinfo <- function(mtd,mdata){
        dataa <- data.frame()
        timestartNine <- Sys.time()
        if(mtd == 'gs'){
            dataa <- gs(mdata)$arcs
        }else if(mtd == 'hc'){
            dataa <- hc(mdata)$arcs 
        }else if(mtd == 'iamb'){
            dataa <- iamb(mdata)$arcs 
        }else if(mtd == 'mmpc'){
            dataa <- mmpc(mdata)$arcs 
        }else if(mtd == 'rsmax2'){
            dataa <- rsmax2(mdata)$arcs 
        }else if(mtd == 'tabu'){
            dataa <- tabu(mdata)$arcs
        }else if(mtd == 'fastiamb'){
            dataa <- fast.iamb(mdata)$arcs
        }else if(mtd == 'interiamb'){
            dataa <- inter.iamb(mdata)$arcs 
        }else if(mtd == 'mmhc'){
            dataa <- mmhc(mdata)$arcs 
        }
        print(paste(mtd,'Functions calling is OK...'))
        timeendNine <- Sys.time()
        timeNine <- c(timeendNine-timestartNine)#timeNine
        
        resultw <- data.frame(dataa)
        #print(resultw)
        tempBX <- data.frame()
        count <- 0
        for(k in 1:length(resultw[,1])){
            if(CheckInArr(resultw[k,],TDBX5) == TRUE){
                count <- count+1
                tempBX[count,1] <- resultw[k,1]
                tempBX[count,2] <- resultw[k,2]
            }
        }
        NineInfo <- data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
        if(count == 0){
    
           print(paste(mtd,'Comparison is OK...'))
           
           NineInfo[,1] <- mtd
           NineInfo[,6] <- length(resultw[,1])
           NineInfo[,7] <- timeNine
           NineInfo[,8] <- 0
           NineInfo[,9] <- 0
           NineInfo[,10] <- paste('p`=',0)
           NineInfo[,11] <- paste('P =',0)
           NineInfo[,12] <- as.character(Sys.time())
           NineInfo[,13] <- 'NULL'
         
        }else{
           tempG <- data.frame()
           countYes <- 0
           for(k in 1:length(tempBX[,1])){
              if(CheckInArr(tempBX[k,],TDJY4) == TRUE){
                  countYes <- countYes+1
                  
                  tempG[countYes,1] <- tempBX[k,1]
                  tempG[countYes,2] <- tempBX[k,2]
              }
          }
          print(paste(mtd,'Comparison is OK...'))
          strr <- c()
           if(countYes == 0)strr <- 'NULL'
          else strr <- ToArrayStr(tempG)
          NineInfo[,1] <- mtd
          NineInfo[,6] <- length(resultw[,1])
          NineInfo[,7] <- timeNine
          NineInfo[,8] <- length(tempBX[,1])
          NineInfo[,9] <- countYes
          NineInfo[,10] <- paste('p`=',(length(tempBX[,1])-countYes)/100)
          NineInfo[,11] <- paste('P =',countYes/4)
          NineInfo[,12] <- as.character(Sys.time())
          NineInfo[,13] <- strr
       }
       return(NineInfo)
    }
    ## end of OPinfo.
    
    
    ## Formalize the results.
    
    outputN <- data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
    colnames(outputN) = c("which_method",
                          "discrete_method",
                          "dispersion ",
                          "source_data",
                          "is_sorted",
                          "whole_relationships",
                          "whole_time",
                          "unique_phenotypes_amount",
                          "unique_genes_amount",
                          "irrelevant_genes_proporition",
                          "unique_genes_proporition",
                          "nine_runtime",
                          "satisfied_relationships")
    
    outputN[1,] <- OPinfo('gs',mydata)
    outputN[2,] <- OPinfo('hc',mydata)
    outputN[3,] <- OPinfo('iamb',mydata)
    outputN[4,] <- OPinfo('mmpc',mydata)
    outputN[5,] <- OPinfo('rsmax2',mydata)
    outputN[6,] <- OPinfo('tabu',mydata)
    outputN[7,] <- OPinfo('fastiamb',mydata)
    outputN[8,] <- OPinfo('interiamb',mydata)
    outputN[9,] <- OPinfo('mmhc',mydata)
    
    outputN[,2] <- methodd
    outputN[,3] <- LSdim
    outputN[,4] <- sourcee
    outputN[,5] <- as.character(Sort1_Or_NotSort2)   
    
    write.csv(outputN,Nout_1,row.names = F,append = TRUE,col.names = F,sep = ' ')
    #write.table(outputN,Nout_2,row.names = F,append = TRUE,col.names = F,sep = ' ')
  
}


