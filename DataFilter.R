# R script for filtering rows and columns of a data frame ,
# or data of multiple dimensions
# Filters for row: Numeric ranges
# Filters for column: Column names

filterDataRow <- function(data, filterVarList=list(list()))
{
    indexes <- rep(TRUE, nrow(data))
    if (length(filterVarList) < 1)
        return(data)
    for (filterVar in filterVarList)
    {
        if (is.null(filterVar[['name']]))
            next
        if (is.null(filterVar[['format']]))
            dataColumn <- data[, filterVar$name]
        else
            dataColumn <- format(data[, filterVar$name], filterVar$format)
        
        ranges <- filterVar$range
        dim(ranges) <- c(2, length(ranges)/2)
        
        tempIndexes <- rep(FALSE, nrow(data))
        for (i in seq(1, ncol(ranges)))
        {
            tempIndexes <- tempIndexes | 
                           dataColumn >= ranges[1, i] &
                           dataColumn <= ranges[2, i]
        }
        indexes <- indexes & tempIndexes
    }
    return(data[indexes,])
}

filterDataCol <- function(data, filterList=c(), reservedVarList=c())
{
    indexes <- rep(FALSE, ncol(data))
    colNameList <- colnames(data)
    for (colName in colNameList)
    {
        if (is.element(colName, filterList) | 
            is.element(colName, reservedVarList))
            indexes[which(colName == colNameList)] <- TRUE
    }
    return(data[indexes])
}

# Helper function for dealing with inverse selection
inverseRangeIndex <- function(indexes, totalLength, selected)
{
    if (selected)
        return(indexes)
    else if (length(indexes) > 0)
    {
        if (length(indexes) < totalLength)
            return(seq(1, totalLength)[-indexes])
        else
            return(c())
    }
    else
        return(seq(1, totalLength))
}

# Function for getting indexes of numbers that fall into at least 
# one range specified in the range list
# data: numbers to examine
# rangeList: list of ranges (i.e. vectors of length 2 - minimum and maximum)
# selected: do normal selection or invert selection
#           if FALSE, parameter 'replicates' is forced to FALSE
# replicates: allow replicates to appear in results, i.e. numbers within
#             each range are considered independent from each other
# minNumber: minimum number of numbers to be included in each range
#            NA are added as paddings if insufficient numbers are found;
#            a positive value implies 'replicates = TRUE'; 
#            thus it is incompatible with 'selected = FALSE'
# maxNumber: maximum number of numbers to be included in each range
#            a zero value will lead to unlimited numbers to be included
rangeFilter.index <- function(data, rangeList=list(), 
                              selected = TRUE, replicates = FALSE,
                              minNumber = 0, maxNumber = 0)
{
    if (length(rangeList) < 1)
        return(inverseRangeIndex(c(), length(data), selected))
    
    indexes <- c()
    for (range in rangeList)
    {
        rangeMasks <- data >= range[1] & data <= range[2]
        
        rangeIndexes <- which(rangeMasks)
        
        if (length(rangeIndexes) < minNumber)
            rangeIndexes <- c(rangeIndexes, 
                              rep(NA, minNumber - length(rangeIndexes)))
        if (length(rangeIndexes) > maxNumber && maxNumber > 0)
            rangeIndexes <- rangeIndexes[1:maxNumber]
            
        indexes <- c(indexes, rangeIndexes)
    }
    
    if (!replicates || !selected)
        indexes <- unique(indexes)
    return(inverseRangeIndex(indexes, length(data), selected))
}

# Function for getting indexes of 2-D vectors that fall within at least 
# one range specified in the range list
# data: vectors (of length 2) to examine
# rangeList: list of ranges (i.e. vectors of length 4 - minimum 1, minimum 2, 
#            maximum 1 and maximum 2)
# selected: do normal selection or invert selection
#           if FALSE, parameter 'replicates' is forced to FALSE
# minNumber: minimum number of numbers to be included in each range
#            NA are added as paddings if insufficient numbers are found;
#            a positive value implies 'replicates = TRUE'; 
#            thus it is incompatible with 'selected = FALSE'
# maxNumber: maximum number of numbers to be included in each range
#            a zero value will lead to unlimited numbers to be included
rangeFilter2D.index <- function(data, rangeList=list(), 
                                selected = TRUE, replicates = FALSE,
                                minNumber = 0, maxNumber = 0)
{
    if (length(rangeList) < 1)
        return(inverseRangeIndex(c(), nrow(data), selected))
    
    indexes <- c()
    for (range in rangeList)
    {
        rangeMasks <- data[,1] >= range[1] & 
                      data[,1] <= range[3] &
                      data[,2] >= range[2] & 
                      data[,2] <= range[4]
        
        rangeIndexes <- which(rangeMasks)
        
        if (length(rangeIndexes) < minNumber)
            rangeIndexes <- c(rangeIndexes, 
                              rep(NA, minNumber - length(rangeIndexes)))
        if (length(rangeIndexes) > maxNumber && maxNumber > 0)
            rangeIndexes <- rangeIndexes[1:maxNumber]
        
        indexes <- c(indexes, rangeIndexes)
    }
    
    if (!replicates || !selected)
        indexes <- unique(indexes)
    return(inverseRangeIndex(indexes, nrow(data), selected))
}