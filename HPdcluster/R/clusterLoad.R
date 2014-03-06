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

sampleLoader <- function(tableName, features = list(...), conf) {

    norowid = FALSE # for further improvement in the future
    if(!is.character(tableName))
        stop("The name of the table should be specified")
    if(is.null(conf))
        stop("The ODBC configuration should be specified")

    nFeatures <- length(features)  # number of features
    ps <- distributedR_status()
    nExecuters <- sum(ps$Inst)   # number of executors asked from distributedR
    noBlock2Exc <- 1             # ratio of block numbers to executor numbers
    nparts <- nExecuters * noBlock2Exc  # number of partitions

    # loading vRODBC or RODBC library for master
    if (! require(vRODBC) )
        library(RODBC)

    # the columns of the tabel
    columns <- ""
    if(nFeatures > 1) {
        for(i in 1:(nFeatures-1)) {
            columns <- paste(columns, features[i], ',')
        }
    }
    columns <- paste(columns, features[nFeatures])

    # connecting to Vertica
    connect <- odbcConnect(conf)
    # checking availabilty of all columns
    qryString <- paste("select", columns, "from", tableName, "limit 1")
    oneLine <- sqlQuery(connect, qryString)
    # check valid response from the database
    if (! is.data.frame(oneLine) ) {
        odbcClose(connect)
        stop(oneLine)
    }
    # reading the number of observations in the table
    qryString <- paste("select count(*) from", tableName)
    nobs <- sqlQuery(connect, qryString)
    # check valid response from the database
    if (! is.data.frame(nobs) ) {
        odbcClose(connect)
        stop(nobs)
    }
    if (!norowid) {
        # check valid number of rows based on rowid assumptions    
        qryString <- paste("select count(distinct rowid) from", tableName, "where rowid >=0 and rowid <", nobs)
        distinct_nobs <- sqlQuery(connect, qryString)
        if( nobs != distinct_nobs ) {
            odbcClose(connect)
            stop("There is something wrong with rowid. Check the assumptions about rowid column in the manual.")
        }
    }
    odbcClose(connect)
    
    nobs <- as.numeric(nobs)
    rowsInBlock <- ceiling(nobs/nparts) # number of rows in each partition
  
    # creating darray for features
    X <- darray(dim=c(nobs, nFeatures), blocks=c(rowsInBlock,nFeatures))

    #Load data from Vertica to darray
    foreach(i, 1:npartitions(X), initArrays <- function(x = splits(X,i), myIdx = i, rowsInBlock = rowsInBlock, 
            nFeatures=nFeatures, tableName=tableName, features=features, conf=conf, norowid=norowid, columns=columns) {

        # loading RODBC for each worker
        if (! require(vRODBC) )
            library(RODBC)
        
        if (norowid) {
            start <- (myIdx-1) * rowsInBlock + 1 # start row in the block from 1
            end <- nrow(x) + start   # end row in the block
        } else {
            start <- (myIdx-1) * rowsInBlock # start row in the block, rowid starts from 0
            end <- nrow(x) + start   # end row in the block
        }

        if (norowid) {
            qryString <- paste("select", columns, " from (select", columns, ", row_number() over(order by", columns, ") as rowid from", tableName, ") T where T.rowid >=", start,"and T.rowid <", end)
        } else {
            qryString <- paste("select",columns, " from", tableName, "where rowid >=", start,"and rowid <", end)
        }
        # each worker connects to Vertica to load its partition of the darray 
        connect <- odbcConnect(conf)
        segment<-sqlQuery(connect, qryString, buffsize= end-start)
        odbcClose(connect)

        x <- NULL
        for (i in 1:nFeatures) {
            x <- cbind(x, segment[[i]])
        }

        update(x)
    })
    
    colnames(X) <- features

    X
}

centerLoader <- function(tableName, features = list(...), conf) {

    if(!is.character(tableName))
        stop("The name of the table should be specified")
    if(is.null(conf))
        stop("The ODBC configuration should be specified")

    nFeatures <- length(features)  # number of features
    # loading vRODBC or RODBC library for master
    if (! require(vRODBC) )
        library(RODBC)
    # connecting to Vertica and reading the number of observations in the table
    qryString <- "select"
    if(nFeatures > 1) {
        for(i in 1:(nFeatures-1)) {
            qryString <- paste(qryString, features[i], ',')
        }
    }
    qryString <- paste(qryString, features[nFeatures])
    qryString <- paste(qryString, " from", tableName)

    connect <- odbcConnect(conf)
    segment<-sqlQuery(connect, qryString)
    odbcClose(connect)
    # check valid response from the database
    if (! is.data.frame(segment) )
        stop(segment)

    as.matrix(segment)
}
# Example:
# loadedSamples <- sampleLoader("mortgage2", list("mltvspline1", "mltvspline2", "agespline1", "agespline2", "hpichgspline", "ficospline"), conf="RDev")
# centers <- centerLoader("mortgage", list("mltvspline1", "mltvspline2", "agespline1", "agespline2", "hpichgspline", "ficospline"), conf="RDev")
