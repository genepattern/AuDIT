
# $Id: AuDIT.R 95 2013-03-12 21:36:04Z manidr $
print ('$Id: AuDIT.R 95 2013-03-12 21:36:04Z manidr $', quote=FALSE)



# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2009) by the
# Broad Institute. All rights are reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute cannot be responsible for its
# use, misuse, or functionality.



## suppress warnings (else GenePattern thinks an error has occurred)
options (warn = -1)


##
## Main function for Automated Detection of Inaccurate and Imprecise Transitions in
## MRM Mass Spectrometry
##


run.AuDIT <-
  function (data.file,                                        # file with data (format described below)
            libdir=NULL,                                      # GenePattern <libdir>; bypassed when NULL (for non-GP use)
            output.prefix=NULL,                               # if non-NULL, used as prefix for output files
            required.columns=c('sample', 'replicate', 'peptide', 'transition.id', 'area', 'IS.area'),
            required.columns.location=1:6,                    # column numbers for required.columns in data.file
            interference.threshold1=1e-5,                     # pvalue above which IS and/or analyte are considered interference free
            interference.threshold2=1e-5,                     # pvalue below which IS and/or analyte definitely have interference
            cv.threshold=0.2,                                 # cv threshold for combined pvalue + cv based interference decision 
            notation=list (good='good', bad='bad', ugly='?'), # notation for transition labeling
            all.pairs=FALSE,                                  # if TRUE use all possible combinations of transitions
            debug=FALSE                                       # if TRUE, extra intermediate results files will be written out
           )
{
  
  # Determine if transitions listed in the data.file have any interference.
  # Expects input in the form of a pre-assembled data table (with required.columns present in the order shown).
  # Samples and the corresponding replicates (incl. diff. concentrations) should be listed such that:
  #   sample is the actual sample id (excluding replicate notation);
  #   replicate indicates replicate number for sample.
  # (sample is usually derived from the SampleName column in Skyline, MultiQuant, etc.)
  # For a given peptide and transition, all replicates must have the same sample name
  # Note that sample must be different for different concentrations of the same sample (if any).
  # Sample replicates would generally have the sample id followed by a suffix in sample.name, with this information
  #  explicitly specified in the (sample, replicate) columns
  # Transitions are indicated in transition.id, and must be unique for each peptide. They may be indicated as 1,2,3
  #  for each peptide (as in MultiQuant) or may be derived from the fragment and charge (y7.2, y5.2, etc.) as for
  #  data from Skyline, etc. While different peptides may have the same transition.id, these must be unique for a
  #  given peptide.


  if (! is.null (libdir) ) {
    source (paste (libdir, "common.R", sep=''))
    if (libdir!='') {
      setLibPath(libdir)
      install.required.packages(libdir)
    }
  }

  # check if the file is a csv file by looking for a comma
  line <- scan(data.file, what="character", nlines=1, quiet=TRUE)
  if (length (grep(",", line)) == 0) {
    stop("Data file must be comma separated.")
  }

  data <- read.csv (data.file,
                    stringsAsFactors=FALSE, na.strings=c('NA','N/A','#N/A'))  # read so that strings are retained,
                                                                              # and numerical columns are identified
  check.columns (data, required.columns)

  if (ncol(data) < max(required.columns.location)) {
    missing.cols <- required.columns.location [required.columns.location > ncol(data)]
    missing.cols <- paste (missing.cols, collapse=", ")
    stop (paste ("Missing required columns", missing.cols, "in dataset"))
  }

  
  data <- data [, required.columns.location]
  colnames (data) <- required.columns
  data <- data.frame (data)



  # convert transition.id to transition which are identical for every peptide
  tr.old <- data [, c('sample', 'replicate', 'peptide', 'transition.id')]
  tr.new <- NULL
  temp <- by (tr.old, tr.old [,'peptide'],
              function (x) {
                trs <- unique (x[,'transition.id'])
                trs.number <- unlist (lapply (x[,'transition.id'], function (v) { which (trs %in% v) }))
                tr.new <<- rbind (tr.new, cbind (x, transition=trs.number))
              })
  data <- merge (data, tr.new)


  
  # generate list of combinations to use for calculating relative ratios
  # if all.pairs==TRUE, all possible pairs of transitions are created
  tr <- sort (unique (data [,'transition']))
  if (all.pairs) {
    library(gtools)
    tr.combo <- combinations ( length(tr), 2, tr )
  } else {
    tr.combo <- cbind (tr, c(tr[-1], tr[1]))
  }

  # calculate relative area ratios for transitions in each sample/replicate;
  # for transitions i,j in {1..n}, ratio r_ij = area(i) / area (j) for all pairs i,j;
  # if all.pairs==FALSE, ri = area(i) / area (i+1 mod n), (r1=tr1/tr2, r2=tr2/tr3, ..., rn=trn/tr1);
  # this is used to determine which transitions have interferences
  rel.ratios <- NULL
  temp <- by (data, data [, c('peptide', 'sample', 'replicate')],
              function (x) {
                ratios <- IS.ratios <- NULL
                for (r in 1:nrow(tr.combo)) {
                  ratio <- x [ x[,'transition']==tr.combo[r,1], 'area'] / x [ x[,'transition']==tr.combo[r,2], 'area']
                  is.ratio <- x [ x[,'transition']==tr.combo[r,1], 'IS.area'] / x [ x[,'transition']==tr.combo[r,2], 'IS.area']

                  if ( length (ratio) == 0 ) ratio <- NA
                  if ( length (is.ratio) == 0 ) is.ratio <- NA
                  
                  ratios <- c (ratios, ratio)
                  IS.ratios <- c (IS.ratios, is.ratio)
                }
                
                rel.ratios <<- rbind (rel.ratios, 
                                      c ( x[1,'peptide'], x[1,'sample'], x[1,'replicate'], ratios, IS.ratios ))
              })
       
  combos <- apply (tr.combo, 1, paste, collapse='.')
  ratio.colnames <- paste ('ratio', combos, sep='_')                   # transition ratios (analyte)
  ISratio.colnames <- paste ('IS.ratio', combos, sep='_')              # internal standard
  pvalue.colnames <- paste ('pvalue', combos, sep='_')
  colnames (rel.ratios) <- c ('peptide', 'sample', 'replicate',
                              ratio.colnames, ISratio.colnames)
                              

  data2 <- data.frame (rel.ratios)
             

  for (col in grep ('ratio', colnames (data2)))
    data2 [, col] <- as.numeric (unlist (lapply (data2 [, col], toString)))


  # use statistical tests to determine if the analyte and IS relative ratios indicate interference
  
  #
  # interference in each sample
  #
  # use t-test for determining (for each sample/conc) whether the analyte and IS ratios are similar
  t.result <- by.collapse (data2, data2 [,c ('peptide', 'sample')], 
                           function (x) {
                             pvals <- NULL
                             for (i in 1:length(combos)) {
                               test <- try ( t.test (x[,ratio.colnames[i]], x[,ISratio.colnames[i]]), TRUE )      # t-test
                               # test <- try ( var.test (x[,ratio.colnames[i]], x[,ISratio.colnames[i]]) )        # F-test
                               pval <- ifelse (!inherits (test, "try-error"), test$p.value , NA)
                               pvals <- c (pvals, pval)
                             }

                             return (pvals)
                           })

  t.result.data <- t.result[,3:ncol(t.result)]
  if (length (t.result.data [!is.na(t.result.data)]) == 0) {
    stop("Could not calculate any p values.")
  }
  
  colnames (t.result)[3:ncol(t.result)] <- pvalue.colnames


  # correct -test pvalues for multiple hypothesis testing (BH-FDR)
  t.pvalues <- t.result [, pvalue.colnames]
  t.corrected.pvalues <- pval2fdr ( unlist (t.pvalues), monotone=FALSE )
  t.result [, pvalue.colnames] <- matrix (t.corrected.pvalues, ncol=ncol (t.pvalues))

  
  # for each transition, combine p-values for the ratios that contain this transition
  #  (e.g., tr1 is used in ratio r1=tr1/tr2 and r3=tr3/tr1; hence combine p-values for r1 and r3 for transition 1,
  #         tr2 is used in ratio r2=tr2/tr3 and r1=tr1/tr2; combine p-values for r2 and r1,
  #         tr3 is used in ratio r3=tr3/tr1 and r2=tr2/tr3; combine p-values for r3 and r2,
  #   etc.)
  # different ratios use the same peak areas, and hence the pvalues for the ratios are DEPENDENT;
  # combination of dependent p-values is accomplished using the orginal method in:
  #  Brown (1975), A method fro combining non-independent one-sided tests of significance, Biometrics, 31:987-992
  # with modifications/improvements noted in:
  #  Kost & McDermott (2002), Combining dependent p-values, Statistics & Probablility Letters, 60:183-190
  #
  combined.p <- data.frame (matrix (rep (1, nrow(t.result) * length (tr)), ncol=length(tr)))
  colnames (combined.p) <- tr


  cov.term <- function (cov.table) {
    total.cov <- 0
    for (i in 2:nrow (cov.table))
      for (j in 1:(i-1)) {
        rho <- cov.table [i,j]
        cov.term <- 3.279 * rho + 0.711 * rho^2           # Kost & McDermott, pg. 188
        total.cov <- total.cov + cov.term
      }
    return (total.cov)
  }
  

  k <- ifelse (all.pairs, length(tr)-1, 2)   # number of tests = number of p-values
  for (i.tr in tr) {
    for (r in 1:nrow(t.result)) {
      ratio.table <- NULL
      i.combos <- which (apply (tr.combo, 1, function (x) { i.tr %in% x }))
      subdata <- data2 [ data2 [,'peptide']==t.result[r,'peptide'] & data2[,'sample']==t.result[r,'sample'], ]
      if ( nrow(subdata) == 0 ) {
        combined.p [r, i.tr] <- NA
        next
      }
      for (i.c in i.combos){
        # ratio.ic <- subdata [, ratio.colnames[i.c]] - subdata [, ISratio.colnames[i.c]]
        ratio.ic <- c (subdata [, ratio.colnames[i.c]], subdata [, ISratio.colnames[i.c]])
        ratio.table <- cbind (ratio.table, ratio.ic)
      }
      colnames (ratio.table) <- combos [i.combos]
      ratio.table [ is.na (ratio.table) ] <- 0
      
      cov.tr.r <- cor (ratio.table)
      cov.tr.r [ is.na (cov.tr.r) ] <- 0

      E.chisq <- 2*k                                      # Brown, pg. 989, 991
      Var.chisq <- 4*k + 2*cov.term (cov.tr.r)
                           
      f.dep <- 2 * E.chisq^2 / Var.chisq
      c.dep <- Var.chisq / (2 * E.chisq)

      chisq.combined <- sum (-2 * log ( t.result [r, pvalue.colnames [i.combos]] ))
      pvalue.combined <- 1 - pchisq (chisq.combined / c.dep, df=f.dep)

      combined.p [r, i.tr] <- pvalue.combined
    }
  }



  # annotate p-values to indication which transition are good, bad and ugly
  ok <- data.frame (matrix (rep (0, nrow(t.result) * length (tr)), ncol=length(tr)))
  ok [ combined.p > interference.threshold1 ] <- notation$good
  ok [ combined.p < interference.threshold2 ] <- notation$bad
  ok [ is.na (combined.p) ] <- NA
  ok [ ok==0 ] <- notation$ugly
  colnames (ok) <- paste ('transition', tr, sep='.')

  sample.result <- cbind (t.result, combined.p, ok)

  ## convert sample.result table into "long" form
  trs <- unlist (lapply (unique (data[,'transition']), toString))
  trs.status <- paste ('transition', trs, sep='.')
  res.samp <- sample.result [ , c('peptide', 'sample', trs, trs.status)]
  res.sample.level.long <- reshape (res.samp, direction='long', varying=list (trs, trs.status))
  # replace 'time' varible to reflect the correct transition
  res.sample.level.long [ ,'time'] <- trs [ res.sample.level.long [ ,'time'] ]
  colnames (res.sample.level.long) <- c('peptide', 'sample', 'transition', 'pvalue.final', 'status', 'id')

  
  ## retrieve and include the transition.id column to enable proper interpretation of results
  d.extra <- unique (data [ , c ('peptide', 'sample', 'transition', 'transition.id')])
  d.result <- merge (d.extra, res.sample.level.long)
  
  
  ## calculate cv for peak area ratio
  data.cv <- by.collapse (data [,'area'] / data [,'IS.area'],
                          data [, c ('peptide', 'sample', 'transition')], cv)
  colnames (data.cv) [ ncol(data.cv) ] <- 'cv'
  d.result <- merge (d.result, data.cv)
  cv.status <- ifelse (d.result [,'cv'] < cv.threshold, notation$good, notation$bad)
  final.call <- ifelse ((d.result[,'status']==notation$ugly & cv.status==notation$bad) |
                        d.result [,'status']==notation$bad | cv.status==notation$bad,
                        notation$bad, notation$good)    # ? are treated as bad
  d.result <- cbind (d.result, cv.status, final.call)

  # eliminate the internally created 'transition' column, and the 'id' column introduced by reshape
  d.final <- d.result [ , setdiff (colnames (d.result), c ('transition', 'id')) ]

               		
  if (!is.null (output.prefix)) {
    # save intermediate and result file
    if (debug) {
      write.table (data2, paste (output.prefix, '-intermediate1.csv', sep=''), sep=',', row.names=FALSE)
      write.table (sample.result, paste (output.prefix, '-intermediate2.csv', sep=''), sep=',', row.names=FALSE)
      write.table (d.result, paste (output.prefix, '-intermediate3.csv', sep=''), sep=',', row.names=FALSE)
    }
    write.table (d.final, paste (output.prefix, '.csv', sep=''), sep=',', row.names=FALSE)
  }

  invisible ( list (data=data2, sample.level=sample.result, result=d.final) )
}



check.columns <- function (data, reqd.cols) {
  # check that columns are exactly as required
  expected.col.names <- unlist (lapply (reqd.cols, tolower))
  col.names <- colnames(data)

  # remove any extra spaces
  col.names <- lapply (col.names, gsub, pattern=" ", replacement="")

  # make lower case
  col.names <- lapply (col.names, tolower)

  col.names <- unlist (col.names)
  if (!identical (col.names, expected.col.names)) {
    cat (paste ("Expected column: ", expected.col.names, " -- found ", col.names, ".\n", sep=''))
    stop (paste ("An error occurred which validating the input file format. \nPlease check that the columns in the input file have the correct names and are in the right order."))
  }
}



##
## Support Functions
##

cv <- function (x) {
  # CV = std. dev. / mean
  value <- try ( sd (x, na.rm=TRUE) / mean (x, na.rm=TRUE), TRUE )
  result <- ifelse ( !inherits (value, "try-error"), value, NA )
  return (result)
}



by.collapse <- function (...) {
  # an extension of the 'by' function to return a data frame
  # instead of a list/array that by usually does

  result <- by (...)
  dimensions <- expand.grid ( dimnames (result) )
  # NB: In expand.grid, the first dimension varies fastest;
  #     the same happens in the output of by, and hence the
  #     table of dimensions and results match
  
  # if the result of the applied function is a single value ... use cbind ... else use do.call (rbind, ...)
  # (look at all list elements to decide the length of the results vector
  #  any one specific list may not be representative)
  result.length <- max (unlist (lapply (result, length)), na.rm=TRUE)
  
  # not all possible levels of the indexing attributes may be present
  # fill NA for missing combination results (else, 'dimensions' will not match with 'data' below)
  for (i in 1:length(result)) 
    if ( is.null (result[[i]])) result[[i]] <- rep (NA, result.length)
 
  if (result.length==1) data <- cbind ( unlist (lapply (result, unlist)) )
  else data <- do.call (rbind, result)
  
  final.result <- cbind (dimensions, data)
  invisible (final.result)
}



# correction for multiple hypothesis testing
# converts a list of p-values to a list of fdr p-values
# courtesy of Stefano Monti

pval2fdr <- function( p, monotone=TRUE )
{
 p.ord <- order(p)
 p1 <- p[p.ord]
 fdr <- p1 * length(p) / cumineq(p1,p1)
 if (monotone)
   for ( i in (length(p)-1):1 ) {
     fdr.min <- min(fdr[i],fdr[i+1])
     if ( !is.na (fdr.min) ) fdr[i] <- fdr.min
   }
 fdr[fdr>1] <- 1
 fdr[rank(p)]
} 

# support function for pval2fdr
cumineq <- function( prm, obs, dir=1, debug=F )
{
 # INPUT:
 #  - prm    n-sized array
 #  - obs    n-sized array
 # WHAT:
 #  for each entry in obs, count how many entries in prm
 #  are <= (dir=1) or >= (dir=2) than that entry
 #
 p.ord <- order(if ( dir==1 ) prm else -prm)
 o.ord <- order(if ( dir==1 ) obs else -obs)
 o.rnk <- rank(if ( dir==1 ) obs else -obs)

 # sort entries
 #
 prm <- prm[p.ord]
 obs <- obs[o.ord]

 u.obs <- unique(obs)
 cup <- c(prm,u.obs)
 cup <- cup[order(if (dir==1) cup else -cup)]
 fp <- length(cup)+1-match(obs,rev(cup))-match(obs,u.obs)

 # return values sorted according to original order (o.rnk)
 #
 return ( if (debug)
            cbind( prm[o.rnk], obs[o.rnk], fp[o.rnk] )
          else
            fp[o.rnk] )
}




preprocess.and.run.AuDIT <- function (data.file, skyline.export=TRUE,
                                      output.prefix=NULL, row.limit=10000,
                                      required.columns=c('sample', 'replicate', 'peptide',
                                                         'transition.id', 'area', 'IS.area'),
                                      required.columns.location=1:6,       # column numbers for required.columns in data.file
                                      debug=FALSE,
                                      ...                                  # other arguments to AuDIT
                                     )
{
  # 
  # Runs AuDIT for files processed using CPTAC Study 9 (or similar) Skyline export template
  #  (containing the following columns: Sample, PeptideSequence, ReplicateName, FragmentIon,
  #                                     PrecursorCharge, ProductCharge, light.Area, heavy.Area)
  #  Sample is usually derived from SampleName, and has the same value for all replicates
  #  (see description of sample in run.AuDIT for more information).
  # Splits files into smaller fragments (to avoid memory problems in GenePattern/R).
  # Also enables processing files with peptides having different numbers of
  #  transitions per peptide
  # Finally assembles all results into a single output file specified by output.file.
  # Can also be used with non-Skyline export files (see run.AuDIT for format).
  #


  # check if the file is a csv file by looking for a comma
  line <- scan(data.file, what="character", nlines=1, quiet=TRUE)
  if (length (grep(",", line)) == 0) {
    stop("Data file must be comma separated.")
  }

  # read data
  d <- read.csv (data.file, na.strings=c('NA', '#N/A', 'N/A'))

  # check that the correct columns are present
  if (skyline.export) {
    input.columns <- c ('Sample', 'PeptideSequence', 'ReplicateName',
                        'FragmentIon', 'PrecursorCharge', 'ProductCharge',
                        'light.Area', 'heavy.Area')
    peptide.col <- 'PeptideSequence'
    if (! all (input.columns %in% colnames (d)) )
      stop (paste ("Missing required columns in dataset:\n",
                   paste ( input.columns [ which (! input.columns %in% colnames (d)) ], collapse=',')))
  } else {
    input.columns <- required.columns
    peptide.col <- 'peptide'
    check.columns (d[,required.columns.location], input.columns)
    colnames (d)[required.columns.location] <- required.columns   # to deal with case differences, etc.
  }
  

	
  # split into parts of smaller size based on row.limit (approximate)
  n <- nrow (d)
  parts <- ceiling (n / row.limit)
  peptides <- unique (d [,peptide.col])
  n.peptides <- length (peptides)
  peptides.in.part <- ceiling (n.peptides / parts)
  
	
  pid <- Sys.getpid ()
  final.result <- intermediate1 <- intermediate2 <- NULL
  for (i in 1:parts) {
    # process each part by
    # (i) creating AuDIT acceptable input file
    # (ii) createing sepatate files for peptides with specified numbers of transitions
		
    start <- peptides.in.part * (i-1) + 1
    end <- min (peptides.in.part * i, n.peptides)
    peptide.subset <- d[,peptide.col] %in% peptides [start : end]
    d.subset <- d [ peptide.subset, ]

    if (skyline.export) {
      d.new <- NULL
      tmp <- apply (d.subset, MARGIN=1,
                    function (x) {
                      transition.id <- paste (x['PrecursorCharge'], x['FragmentIon'], x['ProductCharge'], sep='.')
                      d.row <- c (toString (x['Sample']), toString (x['ReplicateName']), peptide=x[peptide.col],
                                  transition.id, area=x['light.Area'], IS.area=x['heavy.Area'])
                      d.new <<- rbind (d.new, d.row)
                    })
      d.new <- data.frame (d.new)
    } else {
      d.new <- d.subset [, required.columns.location]
    }
    colnames (d.new) <- required.columns
		
		
    # count number of transitions (K) for each peptide
    n.trs <- aggregate (d.new [,'transition.id'], list (d.new [,'peptide']),
                        function (x) { length (unique (x)) })
    colnames (n.trs) <- c ('peptide', 'n.transitions')
    data <- merge (d.new, n.trs)
    
		
    # run AuDIT for peptides with K > 2 transitions
    if ( any (data [, 'n.transitions'] <= 2) ) warning ("Peptides with <= 2 transitions removed from data set")
    data <- data [ data[,'n.transitions'] > 2, ]
    tmp <- by (data [,required.columns], list (data [,'n.transitions']),
               function (tr.data) {
                 file <- paste ('AuDITtemp.', pid, '.csv', sep='')
                 write.csv (tr.data, file, row.names=FALSE)
                 out <- run.AuDIT (file, output.prefix=NULL, required.columns=required.columns,
                                   required.columns.location=1:(length (required.columns)), debug=debug, ...)
                 final.result <<- rbind (final.result, out$result)
                 if (debug) {
                   intermediate1 <<- rbind (intermediate1, out$data)
                   intermediate2 <<- rbind (intermediate2, out$sample.level)
                 }
                 unlink (file)  # delete temp file
               })
  }
	
	
  # write final output
  if (!is.null (output.prefix)) {
    write.csv (final.result, paste (output.prefix, '.csv', sep=''), row.names=FALSE)
    if (debug) {
      write.csv (intermediate1, paste (output.prefix, '-intermediate1.csv', sep=''), row.names=FALSE)
      write.csv (intermediate2, paste (output.prefix, '-intermediate2.csv', sep=''), row.names=FALSE)
    }
  } else {
    invisible (final.result)
  }
}






##
## Command line processing to support GenePattern
##

parse.cmdline <- function (...) {
  # set up AuDIT.r for command line processing (if needed)
  # arguments are specified positionally (since there are not optional arguments) and
  # the following arguments are required
  #   libdir, dataFile, outputPrefix, pvalueThreshold, CVThreshold, allPairs, skylineExport, debug
  args <- list (...)
  if ( length (args) != 8 )
    # expected arguments not present -- error
    stop ("USAGE: R --slave --no-save --args 'libdir, dataFile, outputPrefix, pvalueThreshold, CVThreshold, allPairs, skylineExport, debug' < AuDIT.r\n")
    
  for (i in 1:8) {
    arg <- args [i]
    # remove leading and trailing blanks
    arg <- gsub ("^ *", "", arg)
    arg <- gsub (" *$", "", arg)
    # remove any embedded quotation marks
    arg <- gsub ("['\'\"]", "", arg)
    if (i==1) libdir <<- arg
    if (i==2) dataFile <<- arg
    if (i==3) outputPrefix <<- arg
    if (i==4) pvalueThreshold <<- as.numeric (arg)
    if (i==5) CVThreshold <<- as.numeric (arg)
    if (i==6) allPairs <<- ifelse (arg=='TRUE', TRUE, FALSE)
    if (i==7) skylineExport <<- ifelse (arg=='TRUE', TRUE, FALSE)
    if (i==8) debug.on <<- ifelse (arg=='TRUE', TRUE, FALSE)
  }


    
  preprocess.and.run.AuDIT (data.file = dataFile,
                            skyline.export = skylineExport,
                            libdir = libdir,
                            output.prefix = outputPrefix,
                            interference.threshold1 = pvalueThreshold,
                            interference.threshold2 = pvalueThreshold,
                            cv.threshold = CVThreshold,
                            all.pairs = allPairs,
                            debug = debug.on)
}



install.required.packages <- function (libdir) {
  if (!is.package.installed (libdir, "gtools")) {
    install.package (libdir, "gtools_2.6.1.zip", "gtools_2.6.1.tgz","gtools_2.6.1.tar.gz")
  }
}
