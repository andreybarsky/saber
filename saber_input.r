### input manipulation functions ###

# 'seqvec' (sequence vector) is the parsed character vector that drives saber.
# it is a character vector where each element is a sequence, comprising events separated by spaces.

# the following function does its best to automatically convert the input into the seqvec format.
# if given a filename (and told it is a filename), it checks whether it is csv or xls/x and converts those,
# or treats it as raw text. otherwise you can feed it a string or vector directly using the 'datatype' arg.
file_to_seqvec <- function(data=NULL, datatype='filename', seq_sep='\n', event_sep = ' ') {
  if(is.null(data))
  {
    data = file.choose()
    datatype='filename'
  }
  
  # strip whitespace from start and end:
  trim <- function( x ) {
    gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
  }
  
  
  if(datatype == 'filename' | datatype == 'file' | datatype == 'f')
  {
    if(any(grep(".csv", data))) # if the filename contains .csv it's probably a csv
    {
      seqframe = read.csv(data, stringsAsFactors=F, header=F, strip.white=TRUE)
      seqvec = gsub(" NA", "", apply(seqframe, 1, paste, collapse=event_sep)) # this is not the same as our collapse arg
      seqvec = sapply(seqvec, trim, USE.NAMES=F)
    }
    else if(any(grep(".xls",data)))  # if it contains .xls it's probably an excel file (but is it split by cell, or all in one col?)
    {
      if(!any(grep('readxl', installed.packages()[,1]))) # if the xlsx package isn't installed, we install it
      {
        warning('Please install the "readxl" package (try install.packages("readxl"))')
        #install.packages("readxl")
      }
      
      require("readxl")
      seqframe = read_excel(data, 1, col_names=F)
      if(ncol(seqframe)>1)
      { 
        seqvec = gsub(" NA", "", apply(seqframe, 1, paste, collapse=event_sep))  # collapse the df and get rid of NAs
        seqvec = sapply(seqvec, trim, USE.NAMES=F)
      }  
      else
      { 
        seqvec = seqframe[,1]
        seqvec = sapply(seqvec, trim, USE.NAMES=F)
      }
    }
    else # assume raw text, space separated
    { 
      seqvec = read.table(data, header=F, sep=seq_sep, stringsAsFactors=F)[,1]
      seqvec = sapply(seqvec, trim, USE.NAMES=F)
    }
  }
  else if(datatype == 'string' || datatype == 'str' | datatype == 's')
  {    
    seqvec = strsplit(data, seq_sep)[[1]]
    seqvec = sapply(seqvec, trim, USE.NAMES=F)
  }
  else if(datatype == 'vector' || datatype == 'v')
  {    seqvec = sapply(data, trim, USE.NAMES=F)  }
  else
  {    warning("Invalid datatype argument")  }
}




# takes seqvec and adds explicit START and END events (which should be added if they are not already implicit)
add_terminals <- function(seqvec) 
{
  add_terminal <- function(seq) {return(paste('START', seq, 'END'))}
  return(sapply(seqvec, add_terminal, USE.NAMES=FALSE))
}

# takes seqvec and vector of computed low frequency events, collapses those into a placeholder string
collapse_lfes <- function(seqvec, lfes, lfe_string="LFE")
{
  for (lfe in lfes)
  {
    seqvec = gsub(lfe, lfe_string, seqvec)
  }
  return(seqvec)
}
