### core matrix building functions:

get_pairs <- function(seqvec, lag=1) # takes sequence data vector, outputs data frame of event pairs
{
  eventpairs = data.frame(matrix(ncol=2, nrow=1)) # a growing df with 2 cols
  i = 1
  
  for(line in 1:length(seqvec))
  {
    line_events = strsplit(seqvec[line], ' ')[[1]]
    
    if(!(length(line_events) <= lag)) # ignore this line if it is too short
    {
      for (j in 1:(length(line_events)-lag))
      {
        if(line_events[j] != '' & line_events[j+lag] != '') # ignore any leftover whitespace
        {
          eventpairs[i,] = c(line_events[j], line_events[j+lag])
          i = i+1
        }
      }
    }
  }  
  return(eventpairs)
}

# # converts seqvec into table of discrete event types (count_events does this)
# gettypes <- function(seqvec, control=default_control) 
# {
#   flat_sequence = paste(seqvec, collapse = ' ')
#   eventvec = strsplit(flat_sequence, ' ')[[1]]
#   return(table(eventvec))
# }


# converts vector of sequence data to transitional frequency matrix
get_freqs <- function(seqvec, lag=1, simplify=TRUE) 
{
  event_pairs = get_pairs(seqvec, lag)
  event_table = count_events(seqvec)
  event_types = names(event_table)
  
  
  transfreqs = matrix(0, nrow = length(event_types), ncol = length(event_types), dimnames = list(event_types, event_types)) # empty trans freq matrix
  # rows are antecedents and columns are sequiturs
  # so transfreqs[a,b] is the case where a is followed by b
  
  for(i in 1:nrow(event_pairs)) # for each event pair...
  {
    #cat("counting transition #", i, "\n")
    if(is.na(event_pairs[i,1]) || is.na(event_pairs[i,2])) # if one of the pair is NA....
    {
      warning("Incomplete event pair at transition #", i, "\nRecommend setting terminals to TRUE")
    }
    else
    {
      # increment that cell
      transfreqs[event_pairs[i, 1], event_pairs[i, 2]] = transfreqs[event_pairs[i, 1], event_pairs[i, 2]] + 1
    }
  }
  
  if(simplify) # setting 'simplify' to false lets us skip the following step
  {
    # need to remove rows/cols with only zeroes, for chisquare:
    if(length(which(rowSums(transfreqs)==0)) > 0) # but only if there are those rows, else it breaks
    {
      transfreqs = transfreqs[-which(rowSums(transfreqs)==0),]
    }
    if(length(which(colSums(transfreqs)==0)) > 0)
    {
      transfreqs = transfreqs[,-which(colSums(transfreqs)==0)]  
    }
  }
  
  # else if (control$collapse > 0) # i can't remember why i need this if statement, but it breaks otherwise
  # { # this is a really hacky way to get the cluster algorithm to ignore LFEs:
  #   lfe_col = which(colnames(transfreqs)=="LFE")
  #   transfreqs = transfreqs[,-lfe_col]
  #   transfreqs = transfreqs[-lfe_col,]
  # }
  
  return(transfreqs)
}

get_exp <- function(obsfreqs) # gets expected frequency matrix from observed frequencies
{
  eventtypes = rownames(obsfreqs)
  expfreqs = matrix(0.0, nrow = nrow(obsfreqs), ncol = ncol(obsfreqs), dimnames = list(rownames(obsfreqs), colnames(obsfreqs))) # empty exp freq matrix
  
  for(i in 1:nrow(obsfreqs)) # calculating expected frequencies as ((row total * column total) / grand total)
  {
    for(j in 1:ncol(obsfreqs))
    {
      expfreqs[i,j] = (sum(obsfreqs[i,]) * sum(obsfreqs[,j])) / sum(obsfreqs)
    }
  }
  return(expfreqs)
}