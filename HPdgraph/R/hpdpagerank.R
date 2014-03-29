# Copyright [2013] Hewlett-Packard Development Company, L.P.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#########################################################
#  File hpdpagerank.R
#  Distributed pagerank 
#
#
#  
#  Aarsh Fard (afard@vertica.com)
#
#########################################################

#   dgraph: A darray which represents an adjacency matrix of a graph. It will be altered by this function.
#   niter: The maximum number of iterations
#   eps: The calculation is considered as complete if the difference of PageRank values between iterations change less than this value for every vertex.
#   damping: The damping factor
#   personalized: Optional personalization vector (of type darray). When it is NULL, a constant value of 1/N will be used where N is the number of vertices.
#   weights: Optional edge weights (of type darray). When it is NULL, a constant value of 1 will be used.
#   trace: boolean, to show the progress
#   na_action: the desired behaviour in the case of infinite data in dgraph
hpdpagerank <- function(dgraph, niter = 1000, eps = 0.001, damping=0.85, personalized=NULL, weights=NULL, 
                        trace=FALSE, na_action = c("pass","exclude","fail")) {
  startTotalTime <- proc.time()
  directed = TRUE # the first version treats any graph as a directed graph
  ### Argument checks
  if (!is.darray(dgraph)) { stop("dgraph argument must be of type darray") }
  nVertices <- nrow(dgraph) # number of vertices in the graph
  if (nVertices != ncol(dgraph)) { stop("The input adjacency matrix must be square") }
  if (nVertices != dgraph@blocks[1]) { stop("The input adjacency matrix must be partitioned column-wise") }
  # Checking missed value
  switch(match.arg(na_action),
      "fail" = {
          anyMiss <- .naCheckPageRank(dgraph, trace)
          if(anyMiss > 0)
              stop("missing values in dgraph")
      },
      "exclude" = {
          anyMiss <- .naCheckPageRank(dgraph, trace, cover=TRUE)
          if(anyMiss > 0)
              warning(anyMiss, " edges are excluded because of missed value")
      },
      "pass" = {
          # do nothing            
      },
      {
          stop("only 'pass', 'exclude', and 'fail' are valid values for na_action")
      }
  )
  # Checking personalized
  if (!is.null(personalized)) {
      if(!is.darray(personalized)) stop("personalized argument must be a darray")
      if(any(dim(personalized) != c(1,nVertices)) || personalized@blocks[2] != dgraph@blocks[2])
        stop("dimensions and partitions in personalized should be compatible with dgraph")
      if(.naCheckPageRank(personalized, trace) > 0) stop("missed values in personalized darray!") 
  }
  # Checking weights
  if (!is.null(weights)) {
    if(!is.darray(weights)) stop("weights argument must be a darray")
    if(weights@sparse != dgraph@sparse) stop("weights should be similar to dgraph for sparse feature")
    if(any(dim(weights) != dim(dgraph)) || any(weights@blocks != dgraph@blocks))
      stop("weights should have the same dimension and partitioning as dgraph")
    if(.naCheckPageRank(weights, trace) > 0) stop("missed values in weights darray!") 
  }

  nparts <- npartitions(dgraph)
  blockSize <- dgraph@blocks[2]

  niter <- as.numeric(niter)
  if (! niter > 0 )
    stop("Invalid iteration count, it should be a positive number")
  eps <- as.numeric(eps)
  if (! eps > 0 )
    stop("Invalid value for 'eps', it should be a positive number")
  damping <- as.numeric(damping)
  if (damping <= 0 || damping >= 1)
    stop("Invalid damping factor, it should be between 0 and 1")

  ### Initialization
  if (trace) {
    cat("Initialization\n")
    starttime<-proc.time()
  }
  maxdiff <- eps
  # initial values for PR is 1/nVertices like most literature (it is 1-damping in igraph package)
  PR <- darray(dim=c(1,nVertices), blocks=c(1,blockSize), sparse=FALSE, data=1-damping)

  #Number of outgoing edges calculated in each partition
  NE<- darray(dim=c(nVertices, nparts), blocks=c(nVertices,1), sparse=TRUE)

  #Transition probabilites. Needs to be broadcast
  TP<- darray(dim=c(nVertices, 1), blocks=c(nVertices, 1), sparse=FALSE)

  #Let's get the number of outgoing edges in each partition
  if (is.null(weights)) { # when there is no weight on the edges
    foreach(i, 1:nparts, progress=trace, function(dg = splits(dgraph,i),sNE = splits(NE, i)) {
      val <- rowSums(dg)
      index <- which(val!=0)
      sNE<-sparseMatrix(i=index,j=rep(1,length(index)),x=val[index],dims=c(nrow(dg),1))	
      update(sNE)
    })
  } else {  # when there are weights on the edges
    foreach(i, 1:nparts, progress=trace, function(dg = splits(dgraph,i),sNE = splits(NE, i), wi=splits(weights,i)) {
      val <- rowSums(dg * wi)
      index <- which(val!=0)
      sNE<-sparseMatrix(i=index,j=rep(1,length(index)),x=val[index],dims=c(nrow(dg),1))	
      update(sNE)
    })
  }  

  #Now use a reducer to get the total edge list values and transition probabilities
  foreach(i, 1:1, progress=trace, function(sNE = splits(NE), sTP = splits(TP)) {
    sTP <- rowSums(sNE)
    sTP <- ifelse(sTP > 0, 1/sTP, 0)
    sTP <- as.matrix(sTP,ncol=1)
    update(sTP)
  })

  # Calculating the transition matrix
  if (is.null(weights)) { # when there is no weight on the edges
    foreach(i, 1:nparts, progress=trace, function(dg=splits(dgraph,i), sTP = splits(TP), damping=damping){
      n <- nrow(dg)
      if(class(dg) == "matrix")
        D <- diag(as.vector(sTP), n)
      else
        D <- .sparseDiagonal(n, as.vector(sTP))
      dg <- damping * (D %*% dg)   # This is the transition matrix
      update(dg)
    })
  } else {  # when there are weights on the edges
    foreach(i, 1:nparts, progress=trace, function(dg=splits(dgraph,i), wi=splits(weights,i), sTP = splits(TP), damping=damping){
      n <- nrow(dg)
      if(class(dg) == "matrix")
        D <- diag(as.vector(sTP), n)
      else
        D <- .sparseDiagonal(n, as.vector(sTP))
      dg <- damping * wi * (D %*% dg)   # This is the transition matrix
      update(dg)
    })    
  }

  # the array that keeps the track of the convergence
  diffArray <- darray(dim=c(1,nparts), blocks=c(1,1), sparse=FALSE)
  if (trace) {    # end of timing step
    endtime <- proc.time()
    spentTime <- endtime-starttime
    cat("Spent time:",(spentTime)[3],"sec\n")
  }

  ### iterations
  iteration_counter <- 1
  if (trace)  iterations_start <- proc.time()
  while (niter >= iteration_counter && maxdiff >= eps) {
    if (trace) {
      cat("Iteration: ",iteration_counter,"\n")
      starttime <- proc.time()
    }
    iteration_counter <- iteration_counter + 1
    # the new PR (PageRank vector)
    PR_new <- PR %*% dgraph

    if(is.null(personalized)) {
      foreach(i, 1:nparts, progress=trace, function(PR_newi=splits(PR_new,i), PRi=splits(PR,i), diff=splits(diffArray,i),
               damping=damping, nVertices=nVertices){
        PR_newi <- PR_newi + (1 - damping) / nVertices
        diff <- max(abs(PR_newi - PRi))
        update(PR_newi)
        update(diff)
      })
    } else {
      foreach(i, 1:nparts, progress=trace, function(PR_newi=splits(PR_new,i), PRi=splits(PR,i), diff=splits(diffArray,i), 
                damping=damping, peri=splits(personalized,i)){
        if(class(peri) == "matrix")
          PR_newi <- PR_newi + (1 - damping) * peri
        else
          PR_newi <- as.matrix(PR_newi + (1 - damping) * peri)
        diff <- max(abs(PR_newi - PRi))
        update(PR_newi)
        update(diff)
      })
    }

    maxdiff <- max(diffArray)
    tempPR <- PR
    PR <- PR_new
    PR_new <- tempPR
    if (trace) {    # end of timing step
      endtime <- proc.time()
      spentTime <- endtime-starttime
      cat("Spent time:",(spentTime)[3],"sec\n")
    }
    
  } # while

  if (trace) {
    iterations_finish <- proc.time()
    iterations_totalTime <- iterations_finish[3] - iterations_start[3]
  }

  if (trace) {
    cat("Normalizing the result\n")
    starttime <- proc.time()
  }
  sumPR <- sum(PR)
  foreach(i, 1:nparts, progress=trace, function(PRi=splits(PR,i), sumPR=sumPR){
    PRi <- PRi / sumPR
    update(PRi)
  })
  if (trace) {    # end of timing step
    endtime <- proc.time()
    spentTime <- endtime-starttime
    cat("Spent time:",(spentTime)[3],"sec\n")
  }

  if (trace) {
    endTotalTime <- proc.time()
    totalTime <- endTotalTime - startTotalTime
    cat("*****************************\n")
    cat("Total running time:",(totalTime)[3],"sec\n")
    iterationTime = iterations_totalTime / (iteration_counter -1)
    cat("Running time of each iteration on average:",iterationTime,"sec\n")
  }

  PR
}

## .naCheckPageRank checks any missed value (NA, NaN, Inf) in X
#   X: the input darray
#   trace: boolean, to show the progress 
#   cover: when it is TRUE, the missed values will be replaced with 0
#   it returns the number of missed values
.naCheckPageRank <- function(X, trace, cover = FALSE) {
  if (trace) {
    cat("Checking for missed values\n")
    starttime<-proc.time()
  }
  nparts <- npartitions(X)
  tempArray = darray(dim=c(1,nparts), blocks=c(1,1), data=0)

  foreach(i, 1:nparts, progress=trace, function(tempArrayi=splits(tempArray,i), xi=splits(X,i), cover=cover){
    missed <- !is.finite(xi)
    if(any(missed)) {
      tempArrayi <- matrix(sum(missed))
      update(tempArrayi)
      if(cover) {
        xi[missed] <- 0
        update(xi)
      }
    }
  })
  found <- sum(tempArray)
  if (trace) {    # end of timing step
    endtime <- proc.time()
    spentTime <- endtime-starttime
    cat("Spent time:",(spentTime)[3],"sec\n")
  }

  found
}

## finding the top ranked page
hpdwhich.max <- function(PR, trace=FALSE) {
  if(!is.darray(PR))  stop("PR must be of type darray")
  if(nrow(PR) != 1) stop("PR must have a single row")

  nVertices <- ncol(PR)
  nparts <- npartitions(PR)
  blockSize <- PR@blocks[2]

  dntopValue <- darray(c(nparts,1),c(1,1))  # top value from each partition of PR
  dntopIndex <- darray(c(nparts,1),c(1,1))  # the index (vertex ID) of top value from each partition of PR
  
  if (trace) {
    cat("Finding the top of each partition\n")
    starttime<-proc.time()
  }
  foreach(i, 1:nparts, progress=trace, function(pri=splits(PR,i), dntkV=splits(dntopValue,i), dntkI=splits(dntopIndex,i), idx=i, blockSize=blockSize){
    offset <- blockSize * (idx-1)
    dntkV <- max(pri)
    dntkI <- which.max(pri) + offset
    update(dntkV)
    update(dntkI)
  })
  if (trace) {    # end of timing step
    endtime <- proc.time()
    spentTime <- endtime-starttime
    cat("Spent time:",(spentTime)[3],"sec\n")
  }

  indices <- getpartition(dntopIndex)
  values <-  getpartition(dntopValue)

  indices[which.max(values)]
}

