# takes a saber object and combines particularly high frequency transitions into single events, to detect higher order dependencies
saber.recombine <- function(saber_obj, critval='default', min_obs='default', junksep='@', hftsep='-')
{
  ordered_hft = saber_obj$hft
  newdata = saber_obj$seqdata
  
  if(critval == 'default') {critval = saber_obj$critical_residual*2} # by default the recombination threshold is twice the critical value for notability
  
  if(min_obs == 'default') {min_obs = 2*mean(saber_obj$observed[saber_obj$observed > 0])} # default minimum observations is twice the mean of nonzero observed transitional frequencies. I don't know if this is good or not
  
  for(rownum in 1:nrow(ordered_hft)) # remove the spaces between high frequency events:
  {
    if((ordered_hft$Residual[rownum] > critval) & (ordered_hft$Observed[rownum] > min_obs))
    {
      original = paste(ordered_hft[rownum,1], ordered_hft[rownum,2], sep=' ')
      placeholder = paste(junksep, ordered_hft[rownum,1], hftsep, ordered_hft[rownum,2], junksep, sep='') # intermediate step to avoid over-compression
      recombined = paste(ordered_hft[rownum,1], ordered_hft[rownum,2], sep='')
      for(seqnum in 1:length(newdata))
      {
        newdata[seqnum] = gsub(original, placeholder, newdata[seqnum])
      }
      
    }
    
  }
  for(seqnum in 1:length(newdata))
  {
    newdata[seqnum] = gsub('@', '', newdata[seqnum])
  }
  return(newdata)
}
