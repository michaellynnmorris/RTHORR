#' randmf function
#'
#' comparison of correlation matrices using randomization test of hypothesized order relations
#'     outputs two objects in a list: first the RTHOR results for each model, second the comparisons
#' @param n number of variables (rows in the matrix)
#' @param nmat number of matrices to be analyzed
#' @param ord prediction ordering (default "circular6").
#'     For circular models there are two preset inputs for 6 and 8 variables, "circular6" and "circular8".
#'     Also accepts a vector of prediction ordering, for example circular6 = c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)
#' @param input name of the input file
#'
#' @return Data frame of RTHOR model results, one row per matrix
#'
#' @examples
#' randmf_results <- RTHORR::randmf(n=6,
#'                                  nmat=3,
#'                                  ord="circular6",
#'                                  input=system.file("extdata", "input.txt", package = "RTHORR"))
#'
#' @import permute
#'
#' @export






randmf<-function(n,nmat,ord="circular6",input){

  #process ord input
  if (ord == "circular6"){
    ord = c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)
  } else if (ord == "circular8"){
    ord=c(1,2,3,4,3,2,1,1,2,3,4,3,2,1,2,3,4,3,1,2,3,4,1,2,3,1,2,1)
  } else {
    ord = ord
  }


  # library("permute")
  setMaxperm<-(50000)
  requireNamespace("utils")
  np<- ((n*n)-n)/2

  #read input file


  #read in diagonal matrix as vector and then fill in
  za<-scan(input)  #goodf read in vector
  dmatm<-array(dim=c(n,n,nmat))
  for (m in 1:nmat){
    ii<-(m-1)*(np+n)

    for(j in 1:n){
      for(i in 1:n){
        if (i > j) next
        ii<-ii+1
        dmatm[i,j,m]<-za[ii]}}  #set up matrix


    ii<-(m-1)*(np+n)
    for(i in 1:n){
      for(j in 1:n){
        if (i  < j) next
        ii<-ii+1
        dmatm[i,j,m]<-za[ii]}}  #fill in
  } #end nmat loop


  #mathhyp generation
  mathyp<-matrix(nrow=np,ncol=np)
  for(i in 1:np){
    for(j in 1:np){
      if(ord[j]<ord[i]) mathyp[i,j]<-1 else mathyp[i,j]<-0}}


  #count hyp
  nhyp<-0
  for(i in mathyp){
    if(i==1)nhyp<-nhyp+1}

  pair<-nmat*(nmat-1)/2

  zz<-matrix(nrow=1,ncol=n)
  out<- matrix(nrow=nmat,ncol=6)
  dmatp1<-matrix(nrow=n,ncol=n)
  dmat1<-matrix(nrow=n,ncol=n)
  dmatp2<-matrix(nrow=n,ncol=n)
  dmat2<-matrix(nrow=n,ncol=n)
  out2<-matrix(nrow=pair,ncol=9)
  pp<-vector(length=n)

  #set up permutation file
  nper<-1
  for(i in 1:n){
    nper<- i*nper}  ##good n permutation

  if(nper> 50000)nper<-50000
  permat<-matrix(nrow=nper,ncol=n)

  zz<-matrix(nrow=1,ncol=n)
  f<-function (zz){
    zz<-sample.int(n,n,replace=FALSE,prob=NULL)
    return(zz)}

  if(nper<50000)permat<-permute::allPerms(n, control = how(maxperm = 500000), check = TRUE)

  if(nper>=50000) for(ll in 1:nper){permat[ll,]<-t(apply(zz,1,f))}

  cyc<-0

  #do big loop over nmat pairs***********************************
  for(kk1 in 1:nmat){
    for(kk2 in 2:nmat){
      if(kk1>=kk2) next

      dmat1<-dmatm[,,kk1]
      dmat2<-dmatm[,,kk2]

      #run on original data prior to permutations
      #do data 1,0 matc gt=1 eq =2

      ii=0
      scal1<-vector(length=np)
      for(i in 1:n){
        for(j in 1:n){
          if (i >= j) next
          ii<-ii+1
          scal1[ii]<-dmat1[i,j]}}
      ii=0
      scal2<-vector(length=np)
      for(i in 1:n){
        for(j in 1:n){
          if (i >= j) next
          ii<-ii+1
          scal2[ii]<-dmat2[i,j]}}

      #match data matc with mathyp  1=conf 2=tie, 3=less

      matc1 <-matrix(nrow=np,ncol=np)
      for(i in 1:np){
        for(j in 1:np){
          if(scal1[j] > scal1[i]) matc1[i,j]<-1
          if(scal1[j] == scal1[i]) matc1[i,j]<-2
          if(scal1[j] < scal1[i]) matc1[i,j]<-0}}


      # matc2 <-matrix(nrow=np,ncol=np)
      # for(i in 1:np){
      #   for(j in 1:np){
      #     if(scal2[j] > scal2[i]) matc2[i,j]<-1
      #     if(scal2[j] == scal2[i]) matc2[i,j]<-2
      #     if(scal2[j] < scal2[i]) matc2[i,j]<-0}}

      # Replace first pair of nested for loops with vectorized operations
      matc2 <- matrix(nrow = np, ncol = np)
      matc2[] <- as.numeric(scal2[rep(1:np, each = np)] > scal2[rep(1:np, np)])
      matc2[scal2[rep(1:np, each = np)] == scal2[rep(1:np, np)]] <- 2

      nagr1<-0
      nntie1<-0
      nagr2<-0
      nntie2<-0
      bboth<-0
      nn1y2<-0
      yy1n2<-0
      nn1n2<-0
      vviol1<-0
      vviol2<-0

      for(i in 1:np){
        for(j in 1:np){
          m1<-0
          m2<-0
          t1<-0
          t2<-0
          v1<-0
          v2<-0
          if(matc1[i,j]==1 & mathyp[i,j]==1) nagr1<-nagr1+1
          if(matc1[i,j]==2 & mathyp[i,j]==1) nntie1<-nntie1+1

          if(matc2[i,j]==1 & mathyp[i,j]==1) nagr2<-nagr2+1
          if(matc2[i,j]==2 & mathyp[i,j]==1) nntie2<-nntie2+1


          if(matc1[i,j]==1 & mathyp[i,j]==1) m1<-1
          if(matc1[i,j]==2 & mathyp[i,j]==1) t1<-1

          if(matc2[i,j]==1 & mathyp[i,j]==1) m2<-1
          if(matc2[i,j]==2 & mathyp[i,j]==1) t2<-1

          if(matc1[i,j]==0 & mathyp[i,j]==1) v1<-1
          if(matc2[i,j]==0 & mathyp[i,j]==1) v2<-1

          if(m1==1 & m2==1) bboth<- bboth+1
          if(m1==1 & v2==1) yy1n2<- yy1n2+1
          if(v1==1 & m2==1) nn1y2<- nn1y2+1
          if(v1==1 & v2==1) nn1n2<- nn1n2+1
        }} #finish counting fit

      ci1<- (nagr1-(nhyp-(nagr1+nntie1)))/nhyp
      ci2<- (nagr2-(nhyp-(nagr2+nntie2)))/nhyp

      ci3<- (nn1y2-yy1n2)/(bboth+yy1n2+nn1y2+nn1n2)


      count1<- 1
      count2<-1
      count3<-1 #counter for number exceed values in original permuation
      nperx<- nper - 1

      #do small loop over nper ***********************************************
      for(k in 1:nperx){   #minus 1
        #set up new matrix
        pp<-permat[k,]

        for(i in 1:n){
          for(j in 1:n){
            dmatp1[i,j]<-dmat1[pp[i],pp[j]]
            dmatp2[i,j]<-dmat2[pp[i],pp[j]]
          }}

        #do data 1,0 matc gt=1 eq =2

        ii=0
        scal2<-vector(length=np)
        for(i in 1:n){
          for(j in 1:n){
            if (i >= j) next
            ii<-ii+1
            scal2[ii]<-dmatp2[i,j]}}
        ii=0
        scal1<-vector(length=np)
        for(i in 1:n){
          for(j in 1:n){
            if (i >= j) next
            ii<-ii+1
            scal1[ii]<-dmatp1[i,j]}}



        #match data matc with mathyp  1=conf 2=tie, 3=less
        matc1 <-matrix(nrow=np,ncol=np)
        matc2 <-matrix(nrow=np,ncol=np)
        for(i in 1:np){
          for(j in 1:np){
            if(scal2[j] > scal2[i]) matc2[i,j] <- 1
            if(scal2[j] == scal2[i]) matc2[i,j] <- 2
            if(scal2[j] < scal2[i]) matc2[i,j] <- 0
            if(scal1[j] > scal1[i]) matc1[i,j] <- 1
            if(scal1[j] == scal1[i]) matc1[i,j] <- 2
            if(scal1[j] < scal1[i]) matc1[i,j] <- 0
          }}
        nsup1<-0
        ntie1<-0
        nsup2<-0
        ntie2<-0
        both<-0
        n1y2<-0
        y1n2<-0
        n1n2<-0
        viol1<-0
        viol2<-0

        for(i in 1:np){
          for(j in 1:np){
            m1<-0
            m2<-0
            t1<-0
            t2<-0
            v1<-0
            v2<-0

            if(matc1[i,j]==1 & mathyp[i,j]==1) nsup1<-nsup1+1
            if(matc1[i,j]==2 & mathyp[i,j]==1) ntie1<-ntie1+1

            if(matc2[i,j]==1 & mathyp[i,j]==1) nsup2<-nsup2+1
            if(matc2[i,j]==2 & mathyp[i,j]==1) ntie2<-ntie2+1


            if(matc1[i,j]==1 & mathyp[i,j]==1) m1<-1
            if(matc1[i,j]==2 & mathyp[i,j]==1) t1<-1

            if(matc2[i,j]==1 & mathyp[i,j]==1) m2<-1
            if(matc2[i,j]==2 & mathyp[i,j]==1) t2<-1

            if(matc1[i,j]==0 & mathyp[i,j]==1) v1<-1
            if(matc2[i,j]==0 & mathyp[i,j]==1) v2<-1

            if(m1==1 & m2==1) both<- both+1
            if(m1==1 & v2==1) y1n2<- y1n2+1
            if(v1==1 & m2==1) n1y2<- n1y2+1
            if(v1==1 & v2==1) n1n2<- n1n2+1
          }}

        cip<- (n1y2-y1n2)/(both+y1n2+n1y2+n1n2)
        if(cip>=ci3)count3<-count3+1


        if(nsup1 >= nagr1) count1 <- count1+1   #count number of cases where fit is equal or greater

        if(nsup2 >= nagr2) count2 <- count2+1

      } #end first loop pairs**************************************

      prob1<-count1/nper
      prob2<-count2/nper
      prob3<-count3/nper
      cyc<-cyc+1 #counter for out2


      if (kk1==1) out[kk1,1]<-kk1
      if (kk1==1)out[kk1,2]<-nhyp
      if (kk1==1)out[kk1,3]<-nagr1
      if (kk1==1)out[kk1,4]<-nntie1
      if (kk1==1)out[kk1,5]<-ci1
      if (kk1==1)out[kk1,6]<-prob1

      kkm<-kk2

      out[kkm,1]<-kk2
      out[kkm,2]<-nhyp
      out[kkm,3]<-nagr2
      out[kkm,4]<-nntie2
      out[kkm,5]<-ci2
      out[kkm,6]<-prob2


      out2[cyc,1]<-kk1
      out2[cyc,2]<-" vs. "
      out2[cyc,3]<-kk2
      out2[cyc,4]<-bboth
      out2[cyc,5]<-yy1n2
      out2[cyc,6]<-nn1y2
      out2[cyc,7]<-nn1n2
      out2[cyc,8]<-ci3
      out2[cyc,9]<-prob3



    }}#end loop kk1 and kk2 different matrices**************************************************

  colnames(out)<-c("mat","pred","met","tie","CI","p")
  rownames(out)<-c(1:nmat)
  # print(out,quote=FALSE)
  colnames(out2)<-c("mat1","  ", "mat2","bothmet","1met2not", "2met1not","neither","CI","p")
  rownames(out2)<-c(1:pair)
  # print(out2,quote=FALSE)

  #make output dataframes
  out <- data.frame(out)
  out2 <- data.frame(out2)

  #remove unneeded columns and change type to numeric
  out2 <- subset(out2, select = -c(X..) )
  out2[,1:8] <- sapply(out2[,1:8],as.numeric)

  #change column names
  colnames(out2)[4] <- "1met2not"
  colnames(out2)[5] <- "2met1not"

  #create list of two dataframes for output
  output_list <- list(out, out2)

  #name the two list elements of the output_list
  names(output_list) <- c("RTHOR", "comparisons")

  return(output_list)

}  #end randmf


