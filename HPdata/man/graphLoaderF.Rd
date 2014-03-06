\name{graphLoaderF}
\alias{graphLoaderF}
\title{A simple graph loader from a set of files}
\description{
graphLoaderF function is an example for loading an adjacency matrix from an edge-list which is split in a set of files. It is assumed that edge-list is properly split among nfiles. Split files should correspond to column-wise partitioned adjacency matrix; i.e., each file will be used to load one partition of such a matrix. The files should be available from the same path from all the nodes; i.e., they should be located either on NFS or on the same path in different nodes of the cluster. Moreover, the input argument should be set properly according to the graph stored in the files.
}
\usage{
graphLoaderF(pathPrefix, nVertices, nfiles)
}

\arguments{
  \item{pathPrefix}{
    the path and prefix to the files. It should be the same on all nodes of the cluster.
  }
  \item{nVertices}{
    the total number of vertices in the graph
  }
  \item{nfiles}{
    the number of split files which contain the graph
  }
}

\value{
  \item{X}{the sparse darray of adjacency matrix}
}

\author{
Arash Fard <afard@vertica.com>
}

\examples{
 \dontrun{
    require(HPdgraph)
    distributedR_start()

    # Assuming that the graph has 14 vertices and its edge-list is split 
    # column-wise among 7 files.It is also assumed that they are all 
    # located in "~/temp/test/" path which is available from all workers,
    # and their names (before identifying number) starts with "small".
    dgraphF <- graphLoaderF("~/temp/test/small", nVertices=14, nfiles=7)
 }
}

\keyword{ Distributed R }
\keyword{ Graph Analytics }
