## Step.3
## Implementation of TPDA algorithm, the algorithm input in two cases:
## > One is the result of direct sampling of source data.
## > The other is the result sorted by the association rule after sampling.


library(igraph)
library(bnlearn)
library(Rgraphviz) 

CntTdpa <- 0
TPDA_Output <- data.frame()
    

## Function : TPDA_Algorithm.
## dim : discrete value.
## wt1 : threshold of the first stage.
## wt2 : threshold of the second stage.
## wt3 : threshold of the third stage.
## methodd : discrete method.
## sourcee : input data / source data.
TPDA_Algorithm <- function(dim,wt1,wt2,wt3,methodd,sourcee){
    
    LSdim <- dim         
    weight_1 <- wt1     
    weight_2 <- wt2     
    weight_3 <- wt3    
    MethodName <- methodd
    
    NotSort_File <- sourcee                         #The input data without sort 
    Sort_File <- 'EX/rules/express100_after.csv'    #The input data with sort
    ftin_1 = "IN/weishengsu_add_5_and_4.csv"       
  
    ftout_1 = "EX/Sample.csv"           
    ftout_2 = "EX/SampleAfterChangeNB.csv"            
    ftout_3 = "EX/c2links.csv"           
    ftout_4 = "EX/TPDA_Result.csv"   
    ftout_5 = "OUT/Experiment_Result.csv"  
    
    
    ## Determine the input source
    ## If sorted by association rules,here we use the data after sorted.
    ## If not,we directly use the sampling data.
    Vitamin9 <- read.csv(ftin_1,header = F,stringsAsFactors = FALSE)[-1,]
    
    if(Sort1_Or_NotSort2 ==  TRUE){
        Data100 <- read.csv( Sort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    }else{
        Data100 <- read.csv( NotSort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    }
    SampleA <- rbind(Vitamin9,Data100) #100+9
    write.csv(SampleA,ftout_1,row.names = F)
    
    
    ## Start the TPDA.
    
    z <- read.csv(ftout_1,header = F,stringsAsFactors = FALSE)[-1,]
    row.names(z) <- z[,1]
    z <- z[-1]
    
    
    ## Function : ExchangeSplit.
    ## The implement of replacement,called by the function 'ChangeNumeric'.
    ExchangeSplit <- function(complexx){
        countPt <- length(complexx) 
        sttr <- complexx[,1]
        r <- data.frame()
        ## Map different ranges to different values.
        for(index in 1:countPt){
            ta <- strsplit(sttr[index],split = ",")
            ta <- unlist(ta)
            r[index,1] <- as.numeric(chartr(old = '(',new = ' ',chartr(old = '[',new = ' ',ta[1])))
            r[index,2] <- as.numeric(chartr(old = ')',new = ' ',chartr(old = ']',new = ' ',ta[2])))
            ## Compare the median of these interval.
            r[index,3] <- 0.5*(r[index,1]+r[index,2]) 
            r[index,4] <- 0
            
        }
        ## Suppose the minimum value of the comparison is 255.
        for(ValueChge in 1:length(r[,1])){
            minn <- 255  
            for(pp in 1:length(r[,1])){
                if(r[pp,4] ==  0&&r[pp,3]<minn){
                    minn <- r[pp,3]
                }
            }
            if(minn != 255){
                for(kk in 1:length(r[,1])){
                    if(r[kk,3] ==  minn)r[kk,4] = ValueChge
                }
            }else{}
        }
        ## Return to the fourth column (value column).
        return (r[,4]) 
    }
    
    
    ## Function : ChangeNumeric.
    ## Discrete the source data and replace the results by interval,
    ## After discreted by n values, the source data of each column divided into n ranges, then we regard these ranges as 1 to n.
    ## For example,[2:4], [0:2], [4:6],3 intervals (disordered),here we changed [0:2] to 1, [2:4] to 2, [4:6] to 3.
    ChangeNumeric <- function(tempz,dls){
        tempz <- data.frame(tempz)
        t  <-  data.frame(matrix(as.numeric(unlist(tempz)),ncol = length(tempz[1,])))
        rownames(t) <- rownames(tempz)
        colnames(t) <- colnames(tempz)
        tempz <- t
        
        write.csv(tempz,'D:\\RST\\SampleAfterNumeric.csv',row.names = F) 
        tempz <- data.frame(tempz)
        
        ## Discreted by specific method and given dispersion
        tempName <- discretize(tempz,method = MethodName,breaks = dls)  
        rownames(tempName) <- rownames(tempz)
        colnames(tempName) <- colnames(tempz)
        tempz <- tempName
        write.csv(tempz,'D:\\RST\\SampleAfterLisan.csv',row.names = F) 
        tempz <- as.matrix(tempz)
        for(p in 1:length(tempz[1,])){  
            #print(paste('row:',p,' complete'))    
            splitt <- data.frame(unique(tempz[,p]))
            tp <- as.matrix(splitt)
            
            ## Call the func 'ExchangeSplit' to compete the exchange.
            ch <- ExchangeSplit(tp) 
            
            for(k in 1:length(tempz[,1])){  
                for(kp in 1:length(tp[,1])){
                    if(tempz[,p][k] ==  tp[,1][kp]){ 
                        tempz[,p][k] = ch[kp]   
                    }else{}  
                }    
            }  
        }  
        return(tempz)
    }
    ## end of ChangeNumeric
     
    
    ## Determine Whether the data needs to be discrete.
    if(LSdim != 0){
      print('Discreting...')
      z <- data.frame(ChangeNumeric(z,LSdim))
      write.csv(z,ftout_2,row.names = F) 
      print('Discrete OK.')
    }
    print('Mutual information calculating ...') 
 
    
    ## Function : getDoubleHang_value.
    ## Calculate the value of the determinant of the covariance matrix of the column.
    getDoubleHang_value <- function(mylist1,mylist2){
        Dx1 = var(mylist1)
        Dx2 = var(mylist2)
        double_h_value = Dx1*Dx2-cov(mylist1,mylist2)*cov(mylist1,mylist2)
    }
    
    ## Calculate the initial mutual information, 
    ## and the final result 'relation' is a tri-form of [from,to,value].
    mylist_a <- c()
    mylist_b <- c()
    mylist_info <- c()
    for(xi in 1:nrow(z)){
        mylist_a <- as.double(c(z[xi,]))
        Dxa = var(mylist_a)
        for(yi in 1:nrow(z)){
            #print(yi)
            mylist_b <- as.double(c(z[yi,]))
            Dxb = var(mylist_b)
            cov_ab = getDoubleHang_value(mylist_a,mylist_b)
            mi = 0.5*log(Dxa*Dxb/cov_ab)
            mylist_info <- c(mylist_info,mi)
        }
    }
    info <- matrix(mylist_info,nrow(z),nrow(z),byrow = TRUE)
    rname <- c(row.names(z))
    row.names(info) <- row.names(z)
    colnames(info) <- row.names(z)
    
    ## Find relationships.
    relation <- data.frame()
    for(xi in 1:(nrow(info)-1)){
        yi = xi
        while(TRUE){
            yi = yi+1
            if(nrow(relation) ==  0){
                from <- c(rname[xi])
                to <- c(rname[yi])
                value <- c(info[xi,yi])
                relation <- data.frame(from,to,value,stringsAsFactors = FALSE)
                relation <- data.frame(xi = c(t(relation)),stringsAsFactors = FALSE) 
            }else
                relation <- data.frame(relation,xi = c(rname[xi],rname[yi],info[xi,yi]),stringsAsFactors = FALSE)
            if(yi ==  nrow(info)) break
        }
    }
    row.names(relation) <- c("from","to","value")
    relation <- data.frame(t(relation),stringsAsFactors = FALSE)

        
    ## Two choices for calculating relationships.
    ## The first choice: all relationship.
    AllResult <- function(){
        from <- relation[,1][order(as.numeric(relation$value),decreasing = T)] 
        to <- relation[,2][order(as.numeric(relation$value),decreasing = T)]
        value <- relation[,3][order(as.numeric(relation$value),decreasing = T)]
        relation2 <- data.frame(from,to,value)
        mark1 <- 1
        for(wx in 1:(nrow(relation2))){
            if( as.numeric(as.character(relation2[wx,3])) > weight_1 ){
                #print(paste(as.character(relation2[wx,3]),weight_1))
                mark1 = mark1 +1 
            }
        }
        result <- data.frame()
        result <- relation2[1:mark1,]
        return(result)
    } 
    ## The second choice: discard the "gene-gene" relationships during the calculation.
    PartlyResult <- function(){
        ## step1. discard the "gene-gene" relationships.
        Gx <- c('LUT','ZEA','Bcry','AC','BC')
        Gx <- data.frame((Gx))
        rlt <- data.frame(1,1,1)
        cp <- 0
        for(i in 1:length(relation[,1])){
            cp = cp+1
            for(j in 1:length(Gx[,1])){
              
                if( (as.character(Gx[j,]) ==  as.character(relation$to[i])) || (as.character(Gx[j,]) ==  as.character(relation$from[i])) ){
                    rlt <- data.frame(rlt,i = c((relation$from[i]),(relation$to[i]),(relation$value[i])),stringsAsFactors = FALSE) 
                }else{
                    #print(paste(as.character(Gx[j,]),"  ",as.character(relation$from[i])))
                    cp = cp-1
                }
            }
        }
        ## step2. get the structured data (horizontal structure).
        row.names(rlt) <- c("from","to","value")
        rlt <- data.frame(t(rlt),stringsAsFactors = FALSE)
        rlt <- data.frame((rlt))
        rlt <- rlt[-(1:3),]
        ## step3. get the sorted data (vertical structure).
        from <- rlt[,1][order(as.numeric(rlt$value),decreasing = T)] ##pai xu (jiang xu)
        to <- rlt[,2][order(as.numeric(rlt$value),decreasing = T)]
        value <- rlt[,3][order(as.numeric(rlt$value),decreasing = T)]
        rlt2 <- data.frame(from,to,value)
        
        mark1 <- 1
        for(wx in 1:(nrow(rlt2))){
            if( as.numeric(as.character(rlt2[wx,3])) > weight_1 ){
                #print(paste(as.character(rlt2[wx,3]),weight_1))
                mark1 = mark1 +1 
            }
        }
        result <- data.frame()
        result <- rlt2[1:mark1,]
        return(result)
    }
 
    
    ## The above mentioned 'relation' and 'rlt2' are the same structure, 
    ## the former is the whole relationship, the latter is the relationship after the filter.
    ## These two choices are implemented by calling different functions:
    ## > PartlyResult() returns part of the data.
    ## > AllResult() returns whole data.
    Time_TPDA_Start <- Sys.time()
    tempRst <- AllResult()  
    write.csv(tempRst,ftout_3,row.names = F) 
    print('Mutual information is OK.')
   
    
    
    #weight_2  <-  0.025    
    #weight_3  <-  0.025   
    
    ## The variable 'gene' is the source data, 
    ## and the variable 'build_data' is the ranking based on the mutual information.
    gene <- read.csv(ftout_1,header = FALSE,stringsAsFactors = FALSE)
    build_data <- read.csv(ftout_3,header = T,stringsAsFactors = F)
    ck <- union(build_data$from,build_data$to)
    
    lname <- gene[,1]
    gene <- gene[-1,]
    
    labe <- intersect(lname,ck)
    
    ## Function : Translat.
    ## Implement the transformation of form 'n.x' into form 'GRMZMxxx'.
    Translat <- function(data_nx){
        retdata <- data_nx
        for(i in 1:length(data_nx)){
            retdata[i] <- labe[which( trans_labe == data_nx[i])]
        }
        return(retdata)
    }
    
    ## Implement the transformation of form 'GRMZMxxx' into form 'n.x'.
    trans_labe <- c()
    for(i in 1:length(labe)){ 
        trans_labe <- c(trans_labe,paste("n.",i,sep = "")) 
    }
    
    ct <- c()
    snum <- length(build_data$from)
    for(i in 1:snum){ 
        ct <- c(ct,build_data[i,1],build_data[i,2]) 
    }
    
    xanum <- ct
    ctnum <- length(data.frame(xanum)[,1])
    
    for(i in 1:length(ct)){ 
        #print(paste(i,ct[i]))
        ct[i] <- paste("n.",which(labe == ct[i]),sep = "") 
    }
    print('Executing the Step I of TPDA...')
 
    
    ## Start the first stage [I] of TPDA algorithm.
    ## Through the previous association rules mining, we can get the current node of the extension of the order, 
    ## in which we can create the initial Bayesian network structure in here.
    graE <- c()
    graR <- c()
    for(i in seq(1,ctnum,2)){
        if( length(union(graE,graE)) != length(union(graE,c(ct[i],ct[i+1])))  ){
            graE <- c(graE,ct[i],ct[i+1])
            g <-  graph(graE, directed = T)
        }else{
            #print(paste(ct[i],"+",ct[i+1]))
            if(edge_connectivity(g, source = ct[i], target = ct[i+1], checks = TRUE) == 0){
                graE <- c(graE,ct[i],ct[i+1])
                g <- graph(graE, directed = T)
            }
            else{
                graR <- c(graR,ct[i],ct[i+1])
            }
        }
    }
    print('Step I is OK.')
 
    
    #plot(g,layout=layout.fruchterman.reingold, vertex.size=10,vertex.color="green")
    
    ## Function : Info
    ## Calculate mutual information,
    ## Note:'gxi' is the origin node,'gyi' is the target node,
    ## 'cutset' is the cut set, 'gen' is the whole node.
     
    Info <- function(gxi,gyi,cutset,gen){
        cutset <- data.frame(cutset,stringsAsFactors = F)
        x <- c()
        y <- c()
        cutset_table <- data.frame()
        for(xi in 1:nrow(gen)){
            if(gxi ==  gen[xi,1]){
                x <- as.numeric(c(gen[xi,2:length(gen)]))
            }
            if(gyi ==  gen[xi,]){
                y <- as.numeric(gen[xi,2:length(gen)])
            }
            for(yi in 1:nrow(cutset)){
                if(gen[xi,1] ==  cutset[yi,1]){
                    if(nrow(cutset_table) ==  0)
                        cutset_table <- data.frame(t(gen[xi,]),stringsAsFactors = F)
                    else
                        cutset_table <- data.frame(cutset_table,t(gen[xi,]),stringsAsFactors = F)
                }
            }
        }
        cutset_table <- data.frame(t(cutset_table),stringsAsFactors = F)
        cutset_table <- cutset_table[,-1]
        to_list <- c()
        cut <- data.frame()
        for(xi in 1:nrow(cutset_table)){
            to_list <- as.numeric(c(cutset_table[xi,]))
            if(nrow(cut) ==  0)
                cut <- data.frame(to_list)
            else
                cut <- data.frame(cut,to_list)
        }
        ## Exclude zero items.
        if(length(x) ==  0){
            cut_x <- data.frame(cut)
        }else{
            cut_x <- data.frame(cut,x)
        }
        if(length(y) ==  0){
            cut_y <- data.frame(cut)
        }else{
            cut_y <- data.frame(cut,y)
        }
        if(length(x) ==  0 && length(y) ==  0 ){
            cut_x_y <- data.frame(cut)
        }else if(length(x) ==  0){
            cut_x_y <- data.frame(cut,y)
        }else if(length(y) ==  0){
            cut_x_y <- data.frame(cut,x)
        }else{
            cut_x_y <- data.frame(cut,x,y)
        }
        t = det(cov(cut_x))*det(cov(cut_y))/(det(cov(cut))*det(cov(cut_x_y)))
        cmi = 0.5*log(t)
    }
    ## end of Info.
    
    
    ## Start the second stage [II] of TPDA algorithm.
    ## Use a more complex conditional independence test to determine 
    ## which pairs of nodes should also be added to the network.
    
    print('Executing the Step II of TPDA...')
    g <-  graph(graE, directed = F)
  
    graRnum <- length(graR)

    for(i in  seq(1,graRnum,2)){
      
        ## Find the path and storage it to 'one_path' with form 'n.x'.
        shortpa <- shortest_paths(g, from = graR[i], to = graR[i+1], mode = c("all"))$vpath
        one_path <- names(V(g))[as.integer(shortpa[[1]])]
         
        ## Calculate cut points and mutual information.
        brek <- one_path[-c(1,length(one_path))]
        if(length(brek) == 0){
            info = 1                     
        }else{                         
            # error mark.         
            info <- Info(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),gene)
             
        }
        ## If the mutual information is greater than 'weight_2' then add it to 'graE'.
        if(info > weight_2){ 
            graE <- c(graE,graR[i],graR[i+1]) 
        }
    }
    rm("one_path","i","shortpa","info","brek","graR")
    print('Step II is OK.')
    
    
    ## Start the second stage [III] of TPDA algorithm.
    ## Each edge is checked and the second stage of the formula is used for the independence test to 
    ## determine whether the node is independent of the condition.If the two nodes are conditional, their edges will be deleted.
    print('Executing the Step III of TPDA...')
    g <- graph(graE, directed = F)
    yanum <- graE
    graEnum <- length(data.frame(yanum)[,1])
    
    ## Delete the current edge and store it in 'g_d',and find whether there is still a path in 'g_d',
    ## if there is a path then calculate the mutual information,if not,then skip.
    for(i in  seq(1,graEnum,2)){
        g_d <- g - edge(paste(graE[i],"|",graE[i+1],sep = ""))  
        
        if( edge_connectivity(g_d, source = graE[i], target = graE[i+1], checks = TRUE) > 0){
            shortpa <- shortest_paths(g_d, from = graE[i], to = graE[i+1], mode = c("all"))$vpath
            one_path <- names(V(g))[as.integer(shortpa[[1]])]
            brek <- one_path[-c(1,length(one_path))]
            
            if(length(brek) > 0){
                info <- Info(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),gene)
                
                ## If the mutual information is less than 'weight_3', then delete it.
                if(info < weight_3){ g <- g_d } 
            }
        }
    }
    print('Step III is OK.')
    
    ## Store the final result of the algorithm.
    result <- as_edgelist(g, names = TRUE)
    for(i in 1:length(result[,1])){
      
        result[i,] <- Translat(result[i,])
    }
    write.csv(result,ftout_4,row.names = F)
    Time_TPDA_End <- Sys.time()
    
  
    

          
    ## For the unique experiment, the results were further analyzed here.
    
    resultw <- read.csv(ftout_4,header = FALSE,stringsAsFactors = FALSE)
    resultw <- data.frame(resultw)
    
    TDBX5 <- c('LUT','ZEA','Bcry','AC','BC') # Specific phenotype.
    TDJY4 <- c('GRMZM2G012966','GRMZM2G152135','GRMZM2G300348','GRMZM2G108457') # Specific gene. 
    
    ## Form check for storage.
    CheckInArr <- function(nm,arr){
       
        cnt <- length(arr)
        for(ci in 1:cnt){
            if(nm[,1] ==  arr[ci]||nm[,2] ==  arr[ci]){
                return (TRUE)
            }
        }
        return (FALSE)
    }
    tempBX <- data.frame()
    countYes <- 0
    for(k in 1:length(resultw[,1])){
        if(CheckInArr(resultw[k,],TDBX5) ==  TRUE){
            countYes <- countYes+1
            #print(resultw[k,])
            tempBX[countYes,1] <- resultw[k,1]
            tempBX[countYes,2] <- resultw[k,2]
        }
    }
    tempG <- data.frame()
    countYes <- 0
    for(k in 1:length(tempBX[,1])){
        if(CheckInArr(tempBX[k,],TDJY4) ==  TRUE){
            countYes <- countYes+1
            tempG[countYes,1] <- tempBX[k,1]
            tempG[countYes,2] <- tempBX[k,2]
        }
      
    }
    tempG <- tempG[order(tempG[,1],decreasing = T),] #sort
     
    cat(paste(' -   is_sorted  ',Sort1_Or_NotSort2,
              ' -   phenotypes_amount is   Z =',length(TDBX5),
              ' -   relevant_genes_amount is   Y =',4,
              ' -   irrelevant_genes_amount is   X =',100,
              ' -   relationships_amount is  ',length(resultw[,1]),
              ' -   tp_threshold is  ',weight_1,' ',weight_2,' ',weight_3,' ',
              ' -   dispersion is  ',LSdim,
              ' -   include_BC is   M =',length(tempBX[,1]),
              ' -   unique_genes_amount is   N =',countYes,
              ' -   irrelevant_genes_proporition is   p`=',(length(tempBX[,1])-countYes)/100,
              ' -   unique_genes_proporition is   P =',countYes/4,
              ' -   whole_time is ',c(Time_TPDA_End-Time_TPDA_Start)
    ))
    
    ## Form exchange for storage.
    ToArrayStr <- function(df){
        sstr <- c()
        longstr <- c()
        for(ssi in 1:length(df[,1])){
            sstr <- c(sstr,paste('[',df[ssi,1],',',df[ssi,2],']',sep = ''))
            longstr <- paste(longstr,sstr[ssi])
        }
        return (longstr)
    }
     
    print('Calculation and Storage...')
    
    OpInfo <- data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
    OpInfo[,1] <- paste(weight_1,weight_2,weight_3)
    OpInfo[,2] <- methodd
    OpInfo[,3] <- LSdim
    OpInfo[,4] <- sourcee
    OpInfo[,5] <- as.character(Sort1_Or_NotSort2)
    OpInfo[,6] <- length(resultw[,1])
    OpInfo[,7] <- c(Time_TPDA_End-Time_TPDA_Start)
    OpInfo[,8] <- length(tempBX[,1])
    OpInfo[,9] <- countYes
    OpInfo[,10] <- paste('p`=',(length(tempBX[,1])-countYes)/100)
    OpInfo[,11] <- paste('P =',countYes/4)
    OpInfo[,12] <- Sys.time()
    OpInfo[,13] <- ToArrayStr(tempG)
    
    write.table(OpInfo,ftout_5,row.names = F,append = TRUE,sep = ' ')
    print('TPDA Competed!')
    cat(paste('\n\n'))
    return (OpInfo)
}




## Call the function 'TPDA_Algorithm'.

CntTdpa<-0
TPDA_Output <- data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
colnames(TPDA_Output) <- c('tp_threshold',
                           'discrete_method',
                           'dispersion',
                           'source_data',
                           'is_sorted',
                           'relationships_amount',
                           'whole_time',
                           'unique_phenotypes_amount',
                           'unique_genes_amount',
                           'irrelevant_genes_proporition',
                           'tpda_runtime',                                
                           'unique_genes_proporition',
                           'satisfied_relationships')






