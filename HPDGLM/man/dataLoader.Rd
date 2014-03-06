\name{dataLoader}
\alias{dataLoader}
\title{
    A Simple Data Loader from Database
}
\description{
    dataLoader function is a simple function to load a dataset (a set of samples) from a table in Vertica database to a pair of darrays which correspond to responses and predictors of a predictive model. It is assumed that samples (including responses and predictors) are stored in a single table, and the table contains a column called 'rowid'. It is also assumed that 'rowid' starts from 0, and there is no missed 'rowid'.
}
\usage{
dataLoader(tableName, resp = list(...), pred = list(...), conf)
}
\arguments{
  \item{tableName}{
    it is the name of the table in the database in string format.
  }
  \item{resp}{
    the list of the name of columns corresponding to responses.
  }
  \item{pred}{
    the list of the name of columns corresponding to predictors.
  }
  \item{conf}{
    the name of configuration in ODBC.INI file for connecting to the database.
  }
}

\value{
  \item{Y }{the darray of responses}
  \item{X }{the darray of predictors}
}

\author{
    Arash Fard
}

\examples{
 \dontrun{
    require(distributedR)
    distributedR_start()
    # Assuming that samples are stored in a table named mortgage.
    # The name of response column is def,
    # and the names of predictive columns are mltvspline1, mltvspline2, 
    # agespline1, agespline2, hpichgspline, and ficospline.
    loadedData <- dataLoader ("mortgage", list("def"), list("mltvspline1", 
    "mltvspline2", "agespline1", "agespline2", "hpichgspline",
    "ficospline"), conf="RDev")
 }
}

\keyword{ Database }
\keyword{ Distributed R }
