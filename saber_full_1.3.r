# saber 1.3 - sequence analysis of behavioural events in r

# 1.2 - added ability to avoid clustering certain events, and added LFE list to output
# 1.3 - added connectivity thresholding method, and network/plotting capabilities

# check for plotting packages and install them if not available:
if(!any(grep('network', installed.packages()[,1]))) 
{  install.packages("network")}
if(!any(grep('sna', installed.packages()[,1]))) 
{  install.packages("sna")}
if(!any(grep('GGally', installed.packages()[,1]))) 
{  install.packages("GGally")}

# these packages are required for graph plotting:
library(network)
library(sna)
library(GGally)

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
      seqvec = gsub(' ', '', seqvec)
      seqvec = sapply(seqvec, trim, USE.NAMES=F)
      }
    else if(any(grep(".xls",data)))  # if it contains .xls it's probably an excel file (but is it split by cell, or all in one col?)
    {
      if(!any(grep('readxl', installed.packages()[,1]))) # if the xlsx package isn't installed, we install it
      {
        warning('Installing the "readxl" package for dealing with excel files')
        install.packages("readxl")
      }
      
      require("readxl")
      seqframe = read_excel(data, 1, col_names=F)
      if(ncol(seqframe)>1)
      { 
        seqvec = gsub(" NA", "", apply(seqframe, 1, paste, collapse=event_sep))  # collapse the df and get rid of NAs
        #seqvec = gsub(' ', '', seqvec)   # what the fuck is this for
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

### input manipulation functions ###

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


### core matrix building functions:

get_pairs <- function(seqvec, lag=1) # takes sequence data vector, outputs data frame of event pairs
{
  eventpairs = data.frame(matrix(ncol=2, nrow=1)) # a growing df with 2 cols
  i = 1
  
  for(line in 1:length(seqvec))
  {
    line_events = strsplit(seqvec[line], ' ')[[1]]

    if(!length(line_events) <= lag) # ignore this line if it is too short
    {
      for (j in 1:(length(line_events)-lag))
      {
        eventpairs[i,] = c(line_events[j], line_events[j+lag])
        i = i+1
      }
    }
  }  
  return(eventpairs)
}

# removes zero-sum rows and columns from a matrix:
simplify_matrix <- function(mat)
{
  {
    # need to remove rows/cols with only zeroes, for chisquare:
    if(length(which(rowSums(mat)==0)) > 0) # but only if there are those rows, else it breaks
    {
      mat = mat[-which(rowSums(mat)==0),]
    }
    if(length(which(colSums(mat)==0)) > 0)
    {
      mat = mat[,-which(colSums(mat)==0)]  
    }
  }
}


# converts vector of sequence data to transitional frequency matrix
get_freqs <- function(seqvec, lag=1) 
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


### statistical functions:


# calculate contingency table residuals based on difference between observed and expected frequencies
get_residuals <- function(obsfreqs, expfreqs, type='standard') 
{
  #View(obsfreqs)
  
  # eventtypes = rownames(obsfreqs)
  res_matrix = matrix(0, nrow = nrow(obsfreqs), ncol = ncol(obsfreqs), dimnames = list(rownames(obsfreqs), colnames(obsfreqs))) # empty trans freq matrix
  for(i in 1:nrow(res_matrix))
  {
    for(j in 1:ncol(res_matrix))
    {
      if (type=='sr' || type=='standard' || type=='std' || type=='s')
      {
        res = ((obsfreqs[i,j] - expfreqs[i,j]) / expfreqs[i,j]**0.5) # snr calculated as (cell.obs - cell.exp) / cell.exp**0.5 
        
      }
      else if (type=='ar' || type=='adjusted' || type=='adj' || type=='a') # why does this produce infinities?
      {
        row_proportion = obsfreqs[i,j] / sum(obsfreqs[i,])
        col_proportion = obsfreqs[i,j] / sum(obsfreqs[,j])
        if (col_proportion == 1)
        {
          col_proportion = 0 # implementation note - not sure if this is the right thing to do when AR produces infinities
        }
        
        res = ((obsfreqs[i,j] - expfreqs[i,j]) / 
                 sqrt(expfreqs[i,j] * (1+row_proportion) * (1-col_proportion))) # adj res calculated as (cell.obs - cell.exp) / cell.exp * (1 + row_prop) * (1 - col prop) ** 0.5
        #cat(colnames(obsfreqs)[i], colnames(obsfreqs)[j], ':', obsfreqs[i,j], 'minus', expfreqs[i,j], 'divided by sqrt of (', expfreqs[i,j], 'times (1 +', row_proportion, ') * (1 -', col_proportion, ')) is', res, '\n')
      }
      else
      {
        warning("Invalid residual type (try residual='standard' or 'adjusted')")
      }
      
      if(!is.nan(res)) # catch divide by zero errors
      {
        res_matrix[i,j] = res
      }
      else
      {
        res_matrix[i,j] = 0
      }
      
    }
  }
  return(res_matrix)
}

#### -------------------------------------------------------------------------------


### network stuff:

make_network <- function(full_res_matrix, res_crit)
{
  # from a non-simplified res matrix and critical cutoff value, output a network object
  # and strip non-connected events
  adj_matrix = full_res_matrix
  adj_matrix[adj_matrix < res_crit] = 0
  adj_matrix = round(adj_matrix, 1) # round for display purposes
  
  # clean up adj_matrix to remove non-connected events:
  nce = numeric()
  for (i in 1:ncol(adj_matrix))
  {
    if(sum(adj_matrix[,i]) == 0)
    {
      if(sum(adj_matrix[i,]) == 0)
      {
        nce = c(nce, i)
      }
    }
  }
  if(length(nce) != 0){ # the adj_matrix gets nulled if this is 0
    adj_matrix = adj_matrix[-nce, -nce]
  }
  
  net = network(adj_matrix, directed=T, ignore.eval=F, names.eval='res')
  return(list(adj_matrix, net))
}



### Calculating residual critical value:

res_crit <- function(obs, type='acceleration', plots=T)
{ # the 'plots' local variable is inherited by the subfunctions because R doesn't like recursive default argument references
    # 'obs' should be an unsimplified observed freq matrix
  
  # first, some quick geometric functions:
  
  get_distance <- function(x1, y1, x2, y2) # given two points, get euclidean distance using pythagorean theorem
  { sqrt((x1 - x2)^2 + (y1-y2)^2) }
  
  get_line <- function(x1, y1, x2, y2)  # given two points, get slope and intercept of the line that connects them
  { slope = (y2-y1)/(x2-x1) ; intercept = y1-(x1*slope) ; list(slope=slope, intercept=intercept) }
  
  get_intersection <- function(m1, c1, m2, c2) # given two lines, get the point at which they intersect
  {  x = (c2-c1)/(m1-m2) ; y = (m1 * x) + c1 ; list(x=x, y=y) }
  
  # inner line method:
  res_crit_innerline <- function(obs)  # adjusted so as to be scale invariant, because it treats the graph as a square; but in doing so, is sensitive to the tail
  {
    obs = simplify_matrix(obs)
    res_matrix = get_residuals(obs, get_exp(obs))
    
    res_vector = sort(as.vector(res_matrix), decreasing=T)
    res_vector = res_vector[res_vector>0]
    scaled_res_vector = res_vector*(length(res_vector)/max(res_vector)) # we have to scale the y-axis - to basically make the graph a square (trust me)
    if(plots)
    {
      plot(res_vector, type="S", xlab="Transitions", ylab="Residual", main="Inner Line")
    }
      inner_lines = numeric()
      for (i in 1:length(scaled_res_vector))
      {
        if(plots)
          {segments(0,0,i,res_vector[i], col="light gray")}
        inner_lines[i] = get_distance(0,0,i,scaled_res_vector[i])
      }
      critpoint = which.min(inner_lines)
      if(plots)
        {segments(0,0,critpoint, res_vector[critpoint], col="blue")
        abline(v=critpoint, col="red")}
    
    return(res_vector[critpoint])
  }
  
  # vector projection method:
  res_crit_projection <- function(obs)
  {
    obs = simplify_matrix(obs)
    res_matrix = get_residuals(obs, get_exp(obs))
    
    res_vector = sort(as.vector(res_matrix), decreasing=T)
    res_vector = res_vector[res_vector>0]
    if(plots)
      {plot(res_vector, type="S", xlab="Transitions", ylab="Residual", main="Vector Projection")}
    inv_slope = 1
    baseline = get_line(0,max(res_vector), length(res_vector),0)
    if(plots)
      {abline(baseline$intercept, baseline$slope, col="dark green")}
    outer_lines = numeric()
    for (i in 1:length(res_vector))
    {
      inv_intercept = res_vector[i] - (i*inv_slope)
      base_intersect = get_intersection(inv_slope, inv_intercept, baseline$slope, baseline$intercept)
      this_distance = get_distance(base_intersect$x, base_intersect$y, i, res_vector[i])
      if(plots)
        {segments(base_intersect$x, base_intersect$y, i, res_vector[i], col="light gray")}
      outer_lines = append(outer_lines, this_distance)
    }
    critpoint = which.max(outer_lines)
    inv_intercept = res_vector[critpoint] - (critpoint*inv_slope)
    base_intersect = get_intersection(inv_slope, inv_intercept, baseline$slope, baseline$intercept)
    if(plots)
      {segments(base_intersect$x, base_intersect$y, critpoint, res_vector[critpoint], col="blue")
      abline(v=critpoint, col="red")}
    
    
    return(res_vector[critpoint])  }
  
  # curve acceleration method:
  res_crit_acceleration <- function(obs, smoothing=0.5) # looks for the sharpest turn in the smoothed scree curve. smoothing is the loess smooth parameter
  {
    obs = simplify_matrix(obs)
    res_matrix = get_residuals(obs, get_exp(obs))
    
    res_vector = sort(as.vector(res_matrix), decreasing=T)
    res_vector = res_vector[res_vector>0]
    second_derivative <- function(x,i) {return(x[i+1] + x[i-1] - 2 * x[i])}
    res_len = 1:length(res_vector)
    res_smoothed = loess(res_vector ~ res_len, span=smoothing)$fitted
    second_derivatives = diff(diff(c(NA, res_smoothed, NA))) # second derivative of smoothed loess curve
    scale_factor = max(res_vector) / max(second_derivatives, na.rm=T) # allows us to plot the second derivative on the same scale
    
    critpoint = which.max(second_derivatives)
    # critpoint = which(second_derivatives<0)[1] # is this better?
    
    if(plots)
    {
      plot(res_vector, type="S", ylab="Residual", xlab="Transitions", main='Curve Acceleration')
      lines(res_smoothed, col="blue", lty=2)
      lines(second_derivatives*scale_factor, lty=3, col="dark green")
      abline(v=critpoint, col="red")
      
      legend('topright', c('Smoothed Curve', 'Second Derivative'), lty=c(2,3), col=c('blue', 'dark green'))
    }
    return(res_vector[critpoint])  }
  
  # chi square method:
  res_crit_chisq <- function(obs, chi2, matsize) # sqrt ( chisquare / number of cells)
  {
    obs = simplify_matrix(obs)
    res_matrix = get_residuals(obs, get_exp(obs))
    
    critval = (chi2 / matsize)**0.5
    res_vector = sort(as.vector(res_matrix), decreasing=T)
    res_vector = res_vector[res_vector>0]
    # plot(res_vector, type="S", ylab="Residual", xlab="Transitions", main="Chi Square")
    critpoint = which(res_vector<critval)[1]
    # abline(v=critpoint, col="red")
    return(critpoint)
  }
  
  chisq <- function(obsfreqs, expfreqs) # calculate chi-squared statistic of transitional frequencies
  {
    chi2 = 0
    for(i in 1:nrow(obsfreqs))
    {
      for(j in 1:ncol(obsfreqs))
      {
        chi2 = chi2 + (((obsfreqs[i,j] - expfreqs[i,j])**2) / expfreqs[i,j])
      }
    }
    return(chi2)
  }
  
  # connectivity method:
  res_crit_connectivity <- function(obs, min_transitions = 5)
  {
    # we require a symmetrical adjacency matrix, so don't simplify:
    res_matrix = get_residuals(obs, get_exp(obs))
    
    #res_matrix = fullres
    res_vector = sort(as.vector(res_matrix), decreasing=T)
    
    # if min_transitions is 1 then you just get a two-event plot, so we start a little way in:
    #min_transitions = 5
    
    for (rval in res_vector[min_transitions:length(res_vector)])
    {
      net = make_network(res_matrix, rval)[[2]]
      if (is.connected(net, 'weak'))
      {bestval = (rval); break}
    }
    
    if(plots)
    {
      if(any(grep('GGally', installed.packages()[,1]))) # check for ggnet2 plotting capability
      {
      print(ggnet2(net, arrow.size=8, arrow.gap = 0.04, label=T)) #, edge.label='res')
      }
      else
      {
        warning('Unable to display network plot - requires package GGally')
      }
    }
    #critpoint = which(res_vector<bestval)[1]
    
    return(bestval)
  }

          if(type=='acceleration')
  {    return(res_crit_acceleration(obs))
  }  else if(type=='innerline')
  {    return(res_crit_innerline(obs))
  }  else if(type=='projection')
  {    return(res_crit_projection(obs))
  }  else if(type=='connectivity')
  {    return(res_crit_connectivity(obs))
  }  else if(type=='chisq')
  {    return(res_crit_chisq(obs))
  }  else
  {
    warning("Invalid method for automatically calculating critical residual value (try 'acceleration','projection', or 'connectivity', or a numeric value)")
    warning("Using 1 as default.")
    return(1)
  }
}

get_high_freqs <- function(obs_matrix, exp_matrix, res_matrix, critical, min_obs, type='high')
{
  pairnum = 1 # iterator
  notable_transitions = data.frame(matrix(ncol=7, nrow=1), stringsAsFactors=F)
  colnames(notable_transitions) = c('Antecedent', 'Sequitur', 'Observed', 'Expected', 'Residual', 'Proportion Out', 'Proportion In')
  notable_matrix=obs_matrix
  
  for(i in 1:nrow(res_matrix))
  {
    for(j in 1:ncol(res_matrix))
    {
      if(type=='high')
      {
        if((res_matrix[i,j] >= critical) & (obs_matrix[i,j] >= min_obs))
        {
          notable_transitions[pairnum,] = c(rownames(res_matrix)[i], colnames(res_matrix)[j], obs_matrix[i,j], round(exp_matrix[i,j],3), round(res_matrix[i,j],3), round((obs_matrix[i,j] / sum(obs_matrix[i,])), 3), round((obs_matrix[i,j] / sum(obs_matrix[,j])), 3))
          pairnum = pairnum + 1
        }
        else
        {
          notable_matrix[i,j] = 0
        }
      }
      else if(type=='low') # can also be used to fetch transitions that occur significantly less often than expected
      {
        if(res_matrix[i,j] <= -critical)
        {
          notable_transitions[pairnum,] = c(rownames(res_matrix)[i], colnames(res_matrix)[j], obs_matrix[i,j], round(exp_matrix[i,j],3), round(res_matrix[i,j],3), round((obs_matrix[i,j] / sum(obs_matrix[i,])), 3), round((obs_matrix[i,j] / sum(obs_matrix[,j])), 3))
          pairnum = pairnum + 1
        }
        else
        {
          notable_matrix[i,j] = 0
        }
      }
    }
    
  }
  #force data types:
  notable_transitions$Observed = as.integer(notable_transitions$Observed)
  notable_transitions$Expected = as.double(notable_transitions$Expected)
  notable_transitions$Residual = as.double(notable_transitions$Residual)
  return(list(notable_transitions[order(-notable_transitions$Residual),], notable_matrix))
  
}

# performs hierarchical clustering on event types based on their correlations in the transition matrix:
saber_clust <- function(seqvec, control, events_not_to_cluster)
{
  raw_obs = get_freqs(seqvec, control$lag)
  # raw_obs = simplify_matrix(raw_obs)
  
  # eventrows is rows of event antecedents and sequiturs, so that they can be hierarchically clustered (treating events as cases)
  
  eventrows = as.data.frame(matrix(ncol=2*ncol(raw_obs), nrow=ncol(raw_obs)))
  rownames(eventrows) = colnames(raw_obs)
  
  for (k in colnames(raw_obs))
  {
    eventrows[k,] = c(raw_obs[k,], raw_obs[,k])
  }
  
  eventcols = t(eventrows)
  
  indices_not_to_cluster = numeric()
  
  for (i in events_not_to_cluster)
  {
    indices_not_to_cluster = c(indices_not_to_cluster, which(colnames(eventcols) == i))
  }
  
  newcols = eventcols[,-indices_not_to_cluster]
  
  rowdist = as.dist(1 - cor(newcols)) # this distance metric is somewhat arbitrary
  eventclust = hclust(rowdist, method='average')
  
  ####
  # cut_height is the correlation cutoff
  
  cut_height = control$cluster_height
  ####
  if(!control$quiet){
  plot(eventclust, sub=NA, xlab = NA, main='Event Cluster Dendrogram') ; abline(h=cut_height, lty=2, col="red")
  }
    
  if(cut_height > 0)
  {
    event_tree = cutree(eventclust, h=cut_height)
    
    x_clusters = which(table(event_tree)>1) # which hclust IDs have more than one element?
    
    clust_id = 1
    events_to_cluster = list()
    for (x in x_clusters)
    {
      events_to_cluster[[clust_id]] = labels(which(event_tree==x))
      cluster_name = paste("X", clust_id, sep='')
      names(events_to_cluster)[[clust_id]] = cluster_name
      clust_id = clust_id + 1
    }
    
    
    ## now we actually do the string replacement:
    
    new_input = seqvec
    
    
    clust_num = 1
    for (x in events_to_cluster)
    {
      for (event in x)
      {
        #cat("replacing", event, "with", names(events_to_cluster)[clust_num], "\n")
        new_input = gsub(event,names(events_to_cluster)[clust_num],new_input, fixed=TRUE)
      }
      clust_num = clust_num + 1
    }
    
    
    
    return(list(new_input = new_input, clusters = events_to_cluster, event_cor = cor(eventcols)))
  }
  else
  {
    return(list(new_input = saber_input, clusters = NULL, event_cor = NULL))
  }
}




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


########## main function:
saber <- function(data=NULL, # input data
                  datatype='filename', # is the data variable referring to a filename, or is it itself a sequence string, or a character vector
                  terminals=TRUE, # whether or not to add explicit START and END events to sequence data
                  lag=1, # time offset for event pairs
                  event_sep=' ', # separator character between events
                  seq_sep = '\n', # separator character between sequences
                  residual='standard', # which type of residual to calculate (standard or adjusted)
                  critical='connectivity', # how to calculate critical residual cutoff value for notability
                  min_obs=5, # minimum number of observations for transition notability
                  lfe_collapse='default', # default value is (mean - sd)
                  cluster_height=0, # height at which to cut the hierarchical cluster tree for correlated events (equivalent to 1 - pearson's R)
                  structural_events=c('START', 'END'), # do not cluster these
                  n_recombine=0, # number of times to recombine high-frequency transitions
                  recombine_crit='default', # default is twice the critical residual value
                  recombine_obs='default', # default is twice the mean of nonzero observations
                  quiet=FALSE # suppress plots and output
                  )  
{
  control = list(terminals=terminals, 
                 lag=lag,
                 event_sep=event_sep,
                 seq_sep=seq_sep,
                 residual=residual,
                 critical=critical,
                 min_obs=min_obs, 
                 lfe_collapse=lfe_collapse, 
                 cluster_height=cluster_height,
                 structural_events=structural_events,
                 n_recombine=n_recombine,
                 recombine_crit=recombine_crit, 
                 recombine_obs=recombine_obs, 
                 quiet=quiet)
  
  # process the input data:

  raw_seqvec = file_to_seqvec(data, datatype, control$seq_sep, control$event_sep)
  seqvec = raw_seqvec
  
  if (control$terminals) # if you want explicit START and END events, we add those here
  {
    seqvec = add_terminals(seqvec)
  }
  
  # perform hierarchical clustering if asked for:
  if (control$cluster_height > 0)
  {
    #par(mfrow=c(1,2))
    cluster_obj = saber_clust(seqvec, control, structural_events)
    seqvec = cluster_obj$new_input
    saber_clusters = cluster_obj$clusters
    event_cor = cluster_obj$event_cor
  }
  else
  {
    #par(mfrow=c(1,1))
    saber_clusters = NULL
    event_cor = NULL
  }
  
  # collapse low-frequency-events if we are doing that:
  if (control$lfe_collapse == 'default')
  {
    control$lfe_collapse = mean_minus_sd(count_events(seqvec))
  }
  if (control$lfe_collapse > 0)
  {
    lfes = get_lfes(seqvec, control)
    saber_clusters$LFE = lfes
    seqvec = collapse_lfes(seqvec, lfes)
  }
  
  
  # calculate transition matrices:
  full_obs_matrix = get_freqs(seqvec, lag)
  obs_matrix = simplify_matrix(full_obs_matrix)
  exp_matrix = get_exp(obs_matrix)
  res_matrix = get_residuals(obs_matrix, exp_matrix, control$residual)
  
  if (!is.numeric(control$critical)) # the user can pass a numeric value to 'critical', or a string that acts as a method for extracting the value automatically
  {
    crit_residual = res_crit(full_obs_matrix, control$critical, plots=T)
  } else
  {
    crit_residual = control$critical
  }
  
  full_res_matrix = get_residuals(full_obs_matrix, get_exp(full_obs_matrix), control$residual)
  adj_matrix_and_net = make_network(full_res_matrix, crit_residual)
  adj_matrix = adj_matrix_and_net[[1]]
  net = adj_matrix_and_net[[2]]
  
  hft = get_high_freqs(obs_matrix, exp_matrix, res_matrix, crit_residual, control$min_obs)[[1]]
  
  saber_obj = list(input = raw_seqvec,
                   seqdata = seqvec,
                   alphabet = count_events(seqvec),
                   observed = obs_matrix,
                   expected = exp_matrix,
                   residuals = res_matrix,
                   adjacency = adj_matrix,
                   critical_residual = crit_residual,
                   event_cor = event_cor,
                   clusters = saber_clusters,
                   hft = hft
  )
  
  class(saber_obj) = 'saber'
  
  if(!control$quiet & n_recombine == 0) # show plots and output
  {
    Observed_Transitions = obs_matrix
    Residual_Matrix = round(res_matrix,2)
    View(Observed_Transitions)
    View(Residual_Matrix)
    print(hft)
  }
  
  #par(mfrow=c(1,1))
  
  # recursive higher order sequence analysis:
  if(n_recombine > 0)
  {
    new_input = saber.recombine(saber_obj, critval=recombine_crit, min_obs=recombine_obs)
    saber_obj = saber(new_input, 'v', terminals=FALSE, control$lag, control$event_sep, control$seq_sep, control$residual, control$critical, lfe_collapse=0, cluster_height=0, n_recombine = n_recombine-1, quiet=quiet)

  }
  
  if(!quiet)
  {
    if(!any(grep('GGally', installed.packages()[,1]))) # check for ggnet2 plotting capability
    {plot(saber_obj)}
  }
  return(saber_obj)
}

# plot a saber object (really just the adjacency matrix)
plot.saber <- function(saber_obj, node_labels=TRUE, weights = 'labels', edge_label_size=3)
{ # weights should be one of 'width', 'labels', or NA
  
  if(!any(grep('GGally', installed.packages()[,1]))) # check for ggnet2 plotting capability
  {
    warning('Installing package GGally for network plotting')
    install.packages("GGally")
  }
  require("GGally")
  
  if (weights == 'width')
  {
    # rescale residuals to arrow size:
    weight_max = max(saber_obj$adjacency)
    weight_min = min(saber_obj$adjacency[saber_obj$adjacency != 0])
    a = 0.1 # min arrow width
    b = 1.1 # max arrow width
    new_adj = round(((((b-a) * (saber_obj$adjacency - weight_min)) / (weight_max - weight_min)) + a),2)
    # get rid of rescaled zeroes:
    zeroval = min(new_adj)
    new_adj[new_adj == zeroval] = 0
    
    net = network(new_adj, directed=T, ignore.eval=F, names.eval='res')
    print(ggnet2(net, arrow.size=8, arrow.gap = 0.04, label=node_labels, edge.size = 'res'))
  }
  else if (weights == 'labels')
  {
  
    net = network(saber_obj$adjacency, directed=T, ignore.eval=F, names.eval='res')
    print(ggnet2(net, arrow.size=8, arrow.gap = 0.04, label=node_labels, edge.label='res', edge.label.size = edge_label_size))
  }
  else
  {
    net = network(saber_obj$adjacency, directed=T, ignore.eval=F, names.eval='res')
    print(ggnet2(net, arrow.size=8, arrow.gap = 0.04, label=node_labels))
  }
  
}

# the control variable is a list of parameters passed to the main
# function that specifies how to parse the data, as well as how
# to treat and collapse it etc

# only required if you want to use the internal functions that make up the main saber loop

# default_control = list(datatype='filename', 
#                        terminals=T, 
#                        event_sep=' ',
#                        seq_sep='\n', 
#                        critical='projection', 
#                        residual='std', 
#                        min_obs=2, 
#                        lfe_collapse='default', # default value (mean - sd)
#                        cluster_height=0,
#                        quiet = F)


######### end of program

