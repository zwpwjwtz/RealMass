# R script for reading/writing peak information from/to 
# diverse types of mass spectrum file

# Required libraries:
# mzR

# Utility functions
requirePackage <- function(name)
{
	if (requireNamespace(name))
		return(TRUE)
	stop(paste0('Library "', name, '" cannot be found. ',
				'Please install it before using.'))
}

# Handlers for the "TXT spectrum"
# This is the most commonly accept one, 
# with one m/z and its intensity per line seperated by 
# a separator (usually a whitespace)
# Using a comma (,) as separator will make it like a CSV file
readTxtSpectrum <- function(fileName, fieldSeparator = ' ')
{
    mzSpectrum <- read.csv(fileName, sep = fieldSeparator, header = FALSE)
    colnames(mzSpectrum) <- c('m/z', 'Intensity')
    return(mzSpectrum)
}
writeTxtSpectrum <- function(mzSpectrum, fileName, fieldSeparator = ' ')
{
    write.csv(fileName, sep = fieldSeparator, row.names = F, col.names = F)
}

# Handlers for the "mzXML file"
generateMZXMLHeader <- function(mzSpectrum = NULL, oldHeader = NULL)
{
	if (is.null(oldHeader))
	{
		newHeader <- 
			data.frame(stringsAsFactors = F,
					   seqNum = 0, 
					   acquisitionNum = 0, 
					   msLevel = 1, 
					   polarity = 0, 
					   peaksCount = 0, 
					   totIonCurrent = 0, 
					   retentionTime = 0, 
					   basePeakMZ = 0, 
					   basePeakIntensity = 0, 
					   collisionEnergy = 0, 
					   ionisationEnergy = 0, 
					   lowMZ = 0, 
					   highMZ = 0, 
					   precursorScanNum = 0, 
					   precursorMZ = 0, 
					   precursorCharge = 0, 
					   precursorIntensity = 0, 
					   mergedScan = 0, 
					   mergedResultScanNum = 0, 
					   mergedResultStartScanNum = 0, 
					   mergedResultEndScanNum = 0, 
					   injectionTime = 0, 
					   filterString = as.character(NA), 
					   centroided = FALSE
			)	
	} else 
		newHeader <- oldHeader
	
	if (!is.null(mzSpectrum))
	{
		newHeader$peaksCount <- nrow(mzSpectrum)
		newHeader$totIonCurrent <- sum(mzSpectrum[,2])
		newHeader$lowMZ <- min(mzSpectrum[,1])
		newHeader$lowMZ <- max(mzSpectrum[,1])
	}
	return(newHeader)
}
readMZXMLSpectrum <- function(fileName)
{
	requirePackage('mzR')
		
    mzSpectrum <- mzR::openMSfile(fileName)
    return(mzR::peaks(mzSpectrum))
}
writeMZXMLSpectrum <- function(mzSpectrum, fileName)
{
	requirePackage('mzR')
	
    # Generate a empty mzXML file for writing
    tempfileName <- tempfile()
    write('<?xml?><mzXML></mzXML>', tempfileName)
    
    mzR::copyWriteMSData(list(mzSpectrum), 
                         file = fileName,
                         original_file = tempfileName,
                         header = generateMZXMLHeader(mzSpectrum),
                         outformat = 'mzXML')
}
updateMZXMLSpectrum <- function(mzSpectrum, fileName)
{
	requirePackage('mzR')
	
    oldSpectrum <- mzR::openMSfile(fileName)
    oldHeader <- mzR::header(oldSpectrum)
    mzR::close(oldSpectrum)
    
    mzR::copyWriteMSData(list(mzSpectrum),
                         file = fileName,
                         original_file = fileName,
                         header = generateMZXMLHeader(mzSpectrum, 
                         							  oldHeader),
                         outformat = 'mzXML')
}