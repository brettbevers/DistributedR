\name{graphLoaderDB}
\alias{graphLoaderDB}
\title{A simple graph loader from a database}
\description{
graphLoaderDB function is an example for loading an adjacency matrix from a table of edge-list stored in a database. It is assumed that edge-list is stored in a single table. It is also assumed that the ID of vertices in the table starts from zero. As the returned darray will be column-wise partitioned, it would be more efficient when there is projection (index) on the 'to' column of the table.
}
\usage{
graphLoaderDB(tableName, from, to, conf)
}

\arguments{
  \item{tableName}{
    it is the name of the table in the database in string format.
  }
  \item{from}{
    the name of the column in the table which stores the source vertices
  }
  \item{to}{
    the name of the column in the table which stores the target vertices
  }
  \item{conf}{
    the name of configuration in ODBC.INI file for connecting to the database
  }
}

\value{
  \item{X}{the darray of adjacency matrix}
}

\author{
Arash Fard <afard@vertica.com>
}

\examples{
 \dontrun{
    require(HPdgraph)
    distributedR_start()
    # Assuming that edge-list is stored in a table named "graph",
    # and the names of the columns are "x", "y". 
    dgraphDB <- graphLoaderDB("graph", from="x", to="y", conf="RDev")
 }
}

\keyword{ Database }
\keyword{ Distributed R }
\keyword{ Graph Analytics }
