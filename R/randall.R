#' randall function
#'
#' Randomization test of hypothesized order relations
#' @param n number of variables (rows in the matrix)
#' @param nmat number of matrices to be analyzed
#' @param input name of the input file
#' @param description any information you wish to enter describing each sample
#' @param ord prediction ordering (default 6 type RIASEC circumplex)

#'
#' @return Data frame of RTHOR model results, one row per matrix
#'
#' @examples
#' randall_results <- randall(n=6,nmat=3,input="h:/my documents/all/r/sifbend.txt", description=c("sam1", "sam2", "sam3"), ord=c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1))
#'
#' @import permute
#'
#' @export



# randall(n=6,nmat=3,ord=(c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)),input="h:/my documents/all/r/sifbend.txt")
# randall(n=6,nmat=6,ord=(c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)),input="h:/my documents/all/r/input.txt",samp=(c("total interest","female interest","male interest","total comp","female comp","male comp")))
#
# randall(n=6,nmat=2,ord=(c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)),input="h:/my documents/all/r/sifBend.txt",samp=(c("sample1","sample2")))
#
# randall(n=8,nmat=3,ord=(c(1,2,3,4,3,2,1,1,2,3,4,3,2,1,2,3,4,3,1,2,3,4,1,2,3,1,2,1)),input="h:/my documents/all/r/8type.txt",samp=(c("sample1","sample2")))


#RANDALL program to conduct the Randomization Test of Hypothesized Order Relations.
#Copyright (C) 2016 Terence J. G. Tracey

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


randall<-function(n, nmat, input, samp, ord = c(1,2,3,2,1,1,2,3,2,1,2,3,1,2,1)){
  library("permute")
  setMaxperm<-(50000)
  library("utils")
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

  #print header and nhyp   ********************************
  #out1<-matrix(nrow=2,ncol=2)

  zz<-matrix(nrow=1,ncol=n)
  out<- matrix(nrow=nmat,ncol=7)
  dmatp<-matrix(nrow=n,ncol=n)
  dmat<-matrix(nrow=n,ncol=n)

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

  if(nper<50000)permat<-allPerms(n, control = how(maxperm = 50000), check = TRUE)

  if(nper>=50000) for(ll in 1:nper){permat[ll,]<-t(apply(zz,1,f))}



  #do big loop over nmat  kk
  for(kk in 1:nmat){

    #select matrix from array

    dmat<-dmatm[,,kk]

    #run on original data prior to permutations
    #do data 1,0 matc gt=1 eq =2


    nagr<-0
    ntie<-0

    ii=0
    scal<-vector(length=np)
    for(i in 1:n){
      for(j in 1:n){
        if (i >= j) next
        ii<-ii+1
        scal[ii]<-dmat[i,j]}}

    #match data matc with mathyp  1=conf 2=tie, 3=less

    matc <-matrix(nrow=np,ncol=np)
    for(i in 1:np){
      for(j in 1:np){
        if(scal[j] > scal[i]) matc[i,j]<-1
        if(scal[j] == scal[i]) matc[i,j]<-2
        if(scal[j] < scal[i]) matc[i,j]<-0}}
    nsup<-0
    nntie<-0
    for(i in 1:np){
      for(j in 1:np){
        if(matc[i,j]==1 & mathyp[i,j]==1) nagr<-nagr+1
        if(matc[i,j]==2 & mathyp[i,j]==1) ntie<-ntie+1}}

    ci<- (nagr-(nhyp-(nagr+ntie)))/nhyp


    count<- 1 #counter for number exceed values in original permuation

    nperx<-nper-1  ##check on this correction

    #do small loop over nper
    for(k in 1:nperx){   #minus 1
      #  if(k >50000) break  #stop if greater than 50,000

      #set up new matrix
      pp<-permat[k,]

      for(i in 1:n){
        for(j in 1:n){
          dmatp[i,j]<-dmat[pp[i],pp[j]]
        }}

      #do data 1,0 matc gt=1 eq =2

      ii=0
      scal2<-vector(length=np)
      for(i in 1:n){
        for(j in 1:n){
          if (i >= j) next
          ii<-ii+1
          scal2[ii]<-dmatp[i,j]}}



      #match data matc with mathyp  1=conf 2=tie, 3=less

      matc2 <-matrix(nrow=np,ncol=np)
      for(i in 1:np){
        for(j in 1:np){
          if(scal2[j] > scal2[i]) matc2[i,j] <- 1
          if(scal2[j] == scal2[i]) matc2[i,j] <- 2
          if(scal2[j] < scal2[i]) matc2[i,j] <- 0}}

      nsup<-0
      nntie<-0
      for(i in 1:np){
        for(j in 1:np){
          if(matc2[i,j]==1 & mathyp[i,j]==1) nsup<-nsup+1
          if(matc2[i,j]==2 & mathyp[i,j]==1) nntie<-nntie+1}}

      if(nsup >= nagr) count <- count+1   #count number of cases where fit is equal or greater


    } #end first loop kz
    prob<-count/nper

    out[kk,1]<-kk
    out[kk,2]<-nhyp
    out[kk,3]<-nagr
    out[kk,4]<-ntie
    out[kk,5]<-ci
    out[kk,6]<-prob
    out[kk,7]<-samp[kk]


  }#end loop kk different matrices

  colnames(out)<-c("mat","pred","met","tie","CI","p","description")
  rownames(out)<-c(1:nmat)
  # print(out,quote=FALSE)


  out <- data.frame(out)
  out[,1:6] <- sapply(out[,1:6],as.numeric)

  # out %>% mutate_at('matrixnumber', as.numeric)

  return(out)
}  #end randall

