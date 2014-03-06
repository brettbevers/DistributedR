\name{sampleLoader}
\alias{sampleLoader}
\title{A simple sample loader}
\description{
sampleLoader function is an example for loading a set of samples stored in a table to a darray. It is assumed that samples are stored in a single table, and the table contains a column called 'rowid'. It is also assumed that 'rowid' starts from 0, and there is no missed 'rowid'.
}
\usage{
sampleLoader(tableName, features = list(...), conf)
}

\arguments{
  \item{tableName}{
    it is the name of the table in the database in string format.
  }
  \item{features}{
    the list of the name of columns corresponding to the features of a sample.
  }
  \item{conf}{
    the name of configuration in ODBC.INI file for connecting to the database.
  }
}

\value{
  \item{X}{the darray of samples}
}

\author{
Arash Fard <afard@vertica.com>
}

\examples{
 \dontrun{
    require(distributedR)
    distributedR_start()
    # Assuming that samples are stored in a table named "samples",
    # and the names of the columns are "col1", "col2", "col3", and "col4". 
    loadedSamples <- sampleLoader ("samples", list(col1", "col2", 
    "col3", "col4"), conf="RDev")
 }
}

\keyword{ Database }
\keyword{ Distributed R }
