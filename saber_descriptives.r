### data descriptive functions ###

# takes seqvec and outputs table of raw event frequencies, for collapsing, and for final output
count_events <- function(seqvec) 
{
  rawstring = paste(seqvec, collapse = ' ') # this 'collapse' is not the same as our collapse
  eventcounts = table(strsplit(rawstring, ' '))
  return(eventcounts)
}

# given table of eventcounts, returns default LFE cutoff value (mean - sd)
mean_minus_sd <- function(eventcounts) 
{
  return((mean(eventcounts) - (sd(eventcounts))))
}

# outputs vector of events that are considered low frequency
get_lfes <- function(seqvec, control = default_control) 
{
  lfe_min_obs = control$lfe_collapse
  
  eventcounts = count_events(seqvec)
  
  if (is.na(lfe_min_obs)) # default threshold for LFE category is mean - sd event observations
  {
    lfe_min_obs = mean_minus_sd(eventcounts)
  }
  
  lowfreqevents = c()
  
  for(e in 1:length(eventcounts))
  {
    if(eventcounts[e] < lfe_min_obs)
    {
      lowfreqevents = append(lowfreqevents, names(eventcounts)[e])
    }
  }
  return(lowfreqevents)
}
