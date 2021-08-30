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
generateMZXMLHeader <- function(mzSpectra = NULL, oldHeader = NULL)
{
	if (is.null(oldHeader))
	{
	    template <- 
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
	    template <-  oldHeader[1,]
	
	newHeader <- data.frame()
	for (mzSpectrum in mzSpectra)
	{
	    if (!is.null(mzSpectrum))
	    {
	        newHeaderLine <- template
	        newHeaderLine$seqNum <- nrow(newHeader) + 1
	        newHeaderLine$acquisitionNum <- newHeaderLine$seqNum
	        newHeaderLine$peaksCount <- nrow(mzSpectrum)
	        newHeaderLine$totIonCurrent <- sum(mzSpectrum[,2])
	        basePeakIndex <- which.max(mzSpectrum[,2])
	        newHeaderLine$basePeakMZ = mzSpectrum[basePeakIndex,1]
	        newHeaderLine$basePeakIntensity = mzSpectrum[basePeakIndex,2]
	        newHeaderLine$lowMZ <- min(mzSpectrum[,1])
	        newHeaderLine$highMZ <- max(mzSpectrum[,1])
	        newHeader <- rbind(newHeader, newHeaderLine)
	    }
	}
	return(newHeader)
}
readMZXMLSpectra <- function(fileName)
{
	requirePackage('mzR')
		
    mzFile <- mzR::openMSfile(fileName)
    mzData <- mzR::peaks(mzFile)
    if (typeof(mzData) == 'list' && is.null(dim(mzData)))
        return(mzData)
    else
        return(list(mzData))
}
writeMZXMLSpectra <- function(spectrumList, fileName)
{
	requirePackage('mzR')
	
    # Generate a empty mzXML file for writing
    tempfileName <- tempfile()
    write('<?xml?><mzXML></mzXML>', tempfileName)
    
    if (typeof(spectrumList) != 'list')
        spectrumList <- list(spectrumList)
    mzR::copyWriteMSData(spectrumList, 
                         file = fileName,
                         original_file = tempfileName,
                         header = generateMZXMLHeader(spectrumList),
                         outformat = 'mzXML')
}
updateMZXMLSpectra <- function(spectrumList, fileName)
{
	requirePackage('mzR')
	
    oldFile <- mzR::openMSfile(fileName)
    oldHeader <- mzR::header(oldFile)
    mzR::close(oldFile)
    
    if (typeof(spectrumList) != 'list')
        spectrumList <- list(spectrumList)
    mzR::copyWriteMSData(spectrumList,
                         file = fileName,
                         original_file = fileName,
                         header = generateMZXMLHeader(spectrumList, oldHeader),
                         outformat = 'mzXML')
}


# Handlers for the "mzML file"
generateMZMLHeader <- function(mzSpectra = NULL, oldHeader = NULL)
{
    if (is.null(oldHeader))
    {
        template <- 
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
        template <- oldHeader[1,]
    
    newHeader <- data.frame()
    for (mzSpectrum in mzSpectra)
    {
        if (is.null(mzSpectrum))
            next()
        
        newHeaderLine <- template
        newHeaderLine$seqNum <- nrow(newHeader) + 1
        newHeaderLine$acquisitionNum <- newHeaderLine$seqNum
        newHeaderLine$peaksCount <- nrow(mzSpectrum)
        newHeaderLine$totIonCurrent <- sum(mzSpectrum[,2])
        basePeakIndex <- which.max(mzSpectrum[,2])
        newHeaderLine$basePeakMZ = mzSpectrum[basePeakIndex,1]
        newHeaderLine$basePeakIntensity = mzSpectrum[basePeakIndex,2]
        newHeaderLine$lowMZ <- min(mzSpectrum[,1])
        newHeaderLine$highMZ <- max(mzSpectrum[,1])
        newHeader <- rbind(newHeader, newHeaderLine)
    }
    return(newHeader)
}
readMZMLSpectra <- function(fileName)
{
    requirePackage('mzR')
    
    mzFile <- mzR::openMSfile(fileName)
    mzData <- mzR::peaks(mzFile)
    if (typeof(mzData) == 'list' && is.null(dim(mzData)))
        return(mzData)
    else
        return(list(mzData))
}
writeMZMLSpectra <- function(spectrumList, fileName)
{
    requirePackage('mzR')
    
    # Generate a empty mzML file for writing
    tempfileName <- tempfile()
    write('<?xml?>
           <mzML xmlns="http://psi.hupo.org/ms/mzml" 
               xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
               xsi:schemaLocation="http://psi.hupo.org/ms/mzml 
               http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" 
               id="" version="1.1.0">
           <cvList count="2">
              <cv id="MS" fullName="Mass spectrometry ontology" version="4.1.41" 
              URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo" />
              <cv id="UO" fullName="Unit Ontology" version="09:04:2014" 
              URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo" />
           </cvList></mzML>', tempfileName)
    
    if (typeof(spectrumList) != 'list')
        spectrumList <- list(spectrumList)
    mzR::copyWriteMSData(spectrumList, 
                         file = fileName,
                         original_file = tempfileName,
                         header = generateMZMLHeader(spectrumList),
                         outformat = 'mzML')
}
updateMZMLSpectra <- function(spectrumList, fileName)
{
    requirePackage('mzR')
    
    oldFile <- mzR::openMSfile(fileName)
    oldHeader <- mzR::header(oldFile)
    mzR::close(oldFile)
    
    if (typeof(spectrumList) != 'list')
        spectrumList <- list(spectrumList)
    mzR::copyWriteMSData(spectrumList,
                         file = fileName,
                         original_file = fileName,
                         header = generateMZMLHeader(spectrumList, oldHeader),
                         outformat = 'mzML')
}

# Read one spectrum from a data file by guessing its format
readSpectrum <- function(fileName)
{
    spectrumList <- list()
    
	# Guess file type by suffix
	if (endsWith(tolower(fileName), '.txt'))
	    spectrumList <- readTxtSpectrum(fileName)
	else if (endsWith(tolower(fileName), '.csv'))
	    spectrumList <- readTxtSpectrum(fileName, fieldSeparator = ',')
	else if (endsWith(tolower(fileName), '.mzxml'))
	    spectrumList <- readMZXMLSpectra(fileName)
    else if (endsWith(tolower(fileName), '.mzml'))
        spectrumList <- readMZMLSpectra(fileName)
	else if (endsWith(tolower(fileName), '.rds'))
	    spectrumList <- readRDS(fileName)
	else
	    spectrumList <- readTxtSpectrum(fileName)
	
	if (length(spectrumList) < 1)
	    return(NULL)
	if (typeof(spectrumList) == 'list' && is.null(dim(spectrumList)))
	    return(spectrumList[[1]])
	else
	    return(spectrumList)
}

# Read a spectra file by guessing its format
readSpectra <- function(fileName)
{
    # Guess file type by suffix
    if (endsWith(tolower(fileName), '.txt'))
        return(list(readTxtSpectrum(fileName)))
    else if (endsWith(tolower(fileName), '.csv'))
        return(list(readTxtSpectrum(fileName, fieldSeparator = ',')))
    else if (endsWith(tolower(fileName), '.mzxml'))
        return(readMZXMLSpectra(fileName))
    else if (endsWith(tolower(fileName), '.mzml'))
        return(readMZMLSpectra(fileName))
    else if (endsWith(tolower(fileName), '.rds'))
        return(readRDS(fileName))
    else
        return(NULL)
}