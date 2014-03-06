\name{centerLoader}
\alias{centerLoader}
\title{A simple loader for initial centers of a cluster}
\description{
centerLoader function is an example for loading a set of centers stored in a table to a matrix. It is assumed that centers are stored in a single table. All the rows of the table will be read, and each row will be a center.
}
\usage{
centerLoader(tableName, features = list(...), conf)
}

\arguments{
  \item{tableName}{
    it is the name of the table in the database in string format.
  }
  \item{features}{
    the list of the name of columns corresponding to the features of a center.
  }
  \item{conf}{
    the name of configuration in ODBC.INI file for connecting to the database.
  }
}

\value{
  \item{X}{the matrix of centers}
}

\author{
Arash Fard <afard@vertica.com>
}

\examples{
 \dontrun{
    # Assuming that centers are stored in a table named "centers",
    # and the names of the columns are "col1", "col2", "col3", and "col4". 
    loadedCenters <- centerLoader ("centers", list(col1", "col2", 
    "col3", "col4"), conf="RDev")
 }
}

\keyword{ Database }
\keyword{ K-means }
