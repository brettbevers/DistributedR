\name{hpdpagerank}
\alias{hpdpagerank}
\title{Distributed PageRank}
\description{
hpdpagerank function is a distributed implementation of pagerank algorithm.
}
\usage{
hpdpagerank(dgraph, niter = 1000, eps = 0.001, damping=0.85, 
            personalized=NULL, weights=NULL, trace=FALSE, 
            na_action = c("pass","exclude","fail"))
}

\arguments{
  \item{dgraph}{
    a darray (dense or sparse) which contains the adjacency matrix of the graph. A sparse darray is highly recommended for the sake of efficiency. The darray should be column-wise partitioned. It should be noticed that values of this darray will be altered after running hpdpagerank function.
  }
  \item{niter}{
    the maximum number of iterations
  }
  \item{eps}{
    the calculation is considered as complete if the difference of PageRank values between iterations change less than this value for every vertex.
  }
  \item{damping}{
    the damping factor
  }
  \item{personalized}{
    Optional personalization vector (of type darray). When it is NULL, a constant value of 1/N will be used where N is the number of vertices. This darray should have a single row and the number of its columns should be equal to the number of vertices. Its number of partitions should be the same as dgraph.
  }
  \item{weights}{
    Optional edge weights (of type darray). When it is NULL, a constant value of 1 will be used. The dimensions, sparsity, and partitioning of this darray should be the same as dgrapg.
  }
  \item{trace}{
    when this argument is TRUE, intermediate steps of the progress are displayed.
  }
  \item{na_action}{
    it indicates what should happen when the dgraph contains missed values. Values of NA, NaN, and Inf in the adjacency matrix are treated as missed values. There are three options for this argument 'pass', 'exclude', and 'fail'. The default value is 'pass' which means the missed value will not be checked. When 'exclude' is selected, any edge with missed value will be replaced with zero. When 'fail' is selected, the function will stop in the case of any missed value in the input adjacency matrix.
  }
}

\value{
hpdpagerank returns a darray which contains the PageRank vector.
}
\references{
Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual Web Search Engine. Proceedings of the 7th World-Wide Web Conference, Brisbane, Australia, April 1998.

\url{http://www-db.stanford.edu/~backrub/google.html}
}
\author{
Arash Fard <afard@vertica.com>
}

\examples{
 \dontrun{
    
    library(HPdgraph)
    distributedR_start()
    
    graph <- matrix(0, 6,6)
    graph[2,1] <- 1; graph[2,3] <- 1; graph[3,1] <- 1; graph[3,2] <- 1;
    graph[3,4] <- 1; graph[4,5] <- 1; graph[4,6] <- 1; graph[5,4] <- 1; 
    graph[5,6] <- 1; graph[6,4] <- 1

    dgraph <- as.darray(graph, c(6,3))
    pr <- hpdpagerank(dgraph)
 }
}

\keyword{distributed R}
\keyword{Big Data Analytics}
\keyword{distributed pagerank}
