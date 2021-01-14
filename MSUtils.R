# R script for find difference in m/z of two spectra


# Utility function for selecting m/z in the SOURCE list 
# that come up around the m/z in the TARGET list
# Return a list of vectors, each indicating the indexes of 
# matched m/z in the TARGET list for a given m/z in the SOURCE list
mapMZ <- function(source, target, 
                  tolerance = 1e-3, relativeTolerance = TRUE)
{
    # Calculate the m/z window for each peaks in the reference
    if (relativeTolerance)
    {
        refMin <- target - tolerance * target
        refMax <- target + tolerance * target
    }
    else
    {
        refMin <- target - tolerance
        refMax <- target + tolerance
    }
    
    matched <- sapply(source,
                      function(X, min, max)
                      {
                          return(which(
                              X > min & X < max
                          ))
                      },
                      min = refMin,
                      max = refMax)
    return(matched)
}

# Like mapMZ(), except that only the best match will be selected
mapUniqueMZ <- function(source, target, 
                        tolerance = 1e-3, relativeTolerance = TRUE,
                        noMatchAsNA = FALSE)
{
    # Find all matched target m/z for each source m/z
    if (relativeTolerance)
    {
        matched <- sapply(source,
                          target = target,
                          tolerance = tolerance,
                          function(X, target, tolerance)
                          {
                              matched <- target > X - tolerance * X &
                                         target < X + tolerance * X 
                              delta <- abs(target[matched] - X)
                              return(which(matched)[which.min(delta)])
                          })
    }
    else
    {
        matched <- sapply(source,
                          target = target,
                          tolerance = tolerance,
                          function(X, target, tolerance)
                          {
                              matched <- target > X - tolerance &
                                         target < X + tolerance
                              delta <- abs(target[matched] - X)
                              return(which(matched)[which.min(delta)])
                          })
    }
    if (noMatchAsNA)
    {
        matched <- lapply(matched,
                          function(X)
                          {
                              if (length(X) == 0)
                                  return(NA)
                              else
                                  return(X)
                          })
    }
    return(unlist(matched))
}

# Like mapMZ(), but return a mask (logical values) for the SOURCE list
selectMZ <- function(source, expected, 
                     tolerance = 1e-3, relativeTolerance = TRUE)
{
    mzIndexes <- mapMZ(source = expected,
                       target = source,
                       tolerance = tolerance, 
                       relativeTolerance = relativeTolerance)
    mzMasks <- rep(FALSE, length(source))
    mzMasks[unlist(mzIndexes)] <- TRUE
    return(mzMasks)
}

# Like selectMZ(), except that only the best match will be selected
selectUniqueMZ <- function(source, expected,
                           tolerance = 1e-3, relativeTolerance = TRUE)
{
    mzIndexes <- mapUniqueMZ(source = expected,
                             target = source,
                             tolerance = tolerance, 
                             relativeTolerance = relativeTolerance)
    mzMasks <- rep(FALSE, length(source))
    mzMasks[unlist(mzIndexes)] <- TRUE
    return(mzMasks)
}
