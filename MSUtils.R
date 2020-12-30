# R script for find difference in m/z of two spectra

# Utility function for selecting m/z in the SOURCE list 
# that come up around the m/z in the EXPECTED list
# Return a mask (logical values) for the SOURCE list
selectMZ <- function(source, expected, 
					 tolerance = 1e-3, relativeTolerance = TRUE)
{
	# Calculate the m/z window for each peaks in the reference
	if (relativeTolerance)
	{
		refMin <- expected - tolerance * expected
		refMax <- expected + tolerance * expected
	}
	else
	{
		refMin <- expected - tolerance
		refMax <- expected + tolerance
	}
	
	mzMasks <- sapply(source,
					  function(X, min, max)
					  {
					  	if (length(which(
					  		X >= min & X <= max
					  	)) > 0)
					  		return(TRUE)
					  	else
					  		return(FALSE)
					  },
					  min = refMin,
					  max = refMax)
	return(mzMasks)
}

# Utility function for selecting unique m/z in the SOURCE list 
# that come up around the m/z in the EXPECTED list
# Like selectMZ(), except that only the best match will be selected
selectUniqueMZ <- function(source, expected, 
					       tolerance = 1e-3, relativeTolerance = TRUE)
{
	# Find all matched SOURCE m/z  for each EXPECTED m/z
	if (relativeTolerance)
	{
		matched <- sapply(expected,
						  source = source,
						  tolerance = tolerance,
					      function(X, source, tolerance)
					      {
					          matched <- source > X - tolerance * X &
					 	                 source < X + tolerance * X 
					          delta <- abs(source[matched] - X)
					          return(which(matched)[which.min(delta)])
					      })
	}
	else
	{
		matched <- sapply(expected,
						  source = source,
						  tolerance = tolerance,
				    	  function(X, source, tolerance)
						  {
				    	      matched <- source > X - tolerance &
				                         source < X + tolerance
				    	 	  delta <- abs(source[matched] - X)
						      return(which(matched)[which.min(delta)])
						  })
	}
	
	# Convert indexes to masks
	mzMasks <- rep(FALSE, length(source))
	mzMasks[unlist(matched)] <- TRUE
	return(mzMasks)
}
