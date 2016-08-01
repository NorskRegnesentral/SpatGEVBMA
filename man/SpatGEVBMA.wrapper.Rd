\name{SpatGEVBMA.wrapper}
\alias{SpatGEVBMA.wrapper}

\title{Wrapper for spatial.gev.bma with netcdf files

}
\description{This is a wrapper for running the spatial.gev.bma function in the SpatGEV package and imputing the results on a spatial grid with covariate values, based on a
 netcdf-file for the covariate grid, a spreadsheet with responses at observation stations
 and their spatial location in a text-file. The output is a netcdf-file with user specified
 posterior quantiles for the user specified return period, in addition to spatial maps of 
 these and an optional Inter Quartile Range (IQR) uncertainty map. In addtion the input variables
 are written to file.

 The grid covariates values are placed in netcdf-files in the folder 'covariate.folder'.
 The spatial coordinate system can be whatever you like as long as it is decimal spaced 
 (not minutes and seconds) and the same coordinate system is used for the location of 
 the stations (see below). The coordinates must either be specified by captial 'X' and 'Y' 
 (set 'coordinate.type="XY"') or 'Lat' and 'Lon', with capital Ls (set 'coordinate.type="LatLon"').
 All netcdf-files must give covariate values at the same spatial locations (there is currently no check for this)

 The response data for stations are gathered in a spreadsheet located at 'station.annualMax.file',
 at spreadsheet named or numbered 'station.annualMax.sheet'. The first column of this spreadsheet contains
 the observation year, while first row contains the station number (starting from column 2)

 The spatial locations for the station are given in a seperate table formatted txt-file with the first row 
 containing the column names. A column named 'Stnr' must be present to denote the station number corresponding to those
 in the spreadsheet. If 'coordinate.type="XY"', then there must be columns 'X' and 'Y' specifying
 their spatial location in the same coordinate system as used in the netcdf covariate files. 
 If 'coordinate.type="LatLon"', there must be columns 'Lat' and 'Lon' corresponding to the 
 format in the netcdf-file.

 All stations needs to be within the area covered by the grid of covariate values.

 Results are imputed only in the locations where all covariates are available.
 
}
\usage{
SpatGEVBMA.wrapper(covariates.folder, station.annualMax.file,station.annualMax.sheet = 1,
	station.locations.file, output.path = getwd(), output.folder.name = "SpatGEV.res",
	return.period = 20, post.quantiles = c(0.025,0.5,0.975), show.uncertainty = TRUE, 
	coordinate.type = "XY", table.format = "html", mcmc.reps = 10^5, 
	burn.in = round(mcmc.reps*0.2), cores = 1, annualMax.name = NULL, 
    create.tempfiles = FALSE, keep.temp.files = FALSE, save.all.output = TRUE, 
	testing = FALSE)
}

\arguments{
  \item{covariates.folder}{Path to folder with covariate files in netcdf-format (see above)

}
  \item{station.annualMax.file}{File name of spreadsheet annualMax file (see above)

}
  \item{station.annualMax.sheet}{The sheet name or index containing the station annualMax to be read (exactly 1 number)

}
  \item{station.locations.file}{File name of table formatted textfile including the spatial locations of the stations 

}
  \item{output.path}{Path to the where the result folder should be stored

}
  \item{output.folder.name}{Name of result folder
}
  \item{return.period}{Return period to impute results for (single number or a vector of numbers)
}
  \item{post.quantiles}{Vector of quantiles for which the posterior should be evaluated
}
  \item{show.uncertainty}{Logical indicating whether an IQR uncertainty plot should also be provided
}
  \item{coordinate.type}{Character indicating the type/name of coordinate system being used, either "XY" or "LatLon" (see above)
}
  \item{transform.output}{Character specifying whether and how the output should be transformed. NULL corresponds to no transformation. "UTM_QQ_to_LatLon" transforms from UTM QQ (insert number) to LatLon
}
  \item{table.format}{Character indicating the format for the covariate effect summary tables. Either "html" or "latex".
}
  \item{mcmc.reps}{Number of MCMC runs for fitting the model with the station data. Should typically be at least be 10^5
}
  \item{burn.in}{The length of the initial burn-in period which is removed
}
  \item{cores}{The number of cores on the computer used for the imputation. Using detectCores()-1 is good for running locally.
}
  \item{annualMax.name}{Name of annualMax data used in output plots and netcdf files. If NULL, then the name of the specified sheet is used.
}
  \item{create.tempfiles}{Logical indicating whether temporary files should be saved in a Temp folder to perform debugging and check intermediate variables/results if the function crashes
}
  \item{keep.temp.files}{Logical indicating whether the temporary files (if written) should be kept or deleted on function completion
}
  \item{save.all.output}{Logical indicating whether all R objects should be save to file upon function completion. Allocates approx 2.5 Gb per return period for all of Norway.
}
  \item{testing}{Variable indicating whether the run is a test or not. FALSE indicates no testing, a positive number indicates the number of locations being imputed
}
}

\value{
	This returns posteriors in netcdf file format, in addition to corresponding image plots.
}

\author{
 Martin Jullum <jullum@nr.no>
}
