########## main function:

### to do: combine LFE and cluster membership variable

source("./saber_input.r")
source("./saber_descriptives.r")
source("./saber_frequencies.r")
source("./saber_stats.r")
source("./saber_highfreqs.r")
source("./saber_clust.r")
source("./saber_recombine.r")

saber <- function(data=NULL, # input data
                  datatype='filename', # is the data variable referring to a filename, or is it itself a sequence string, or a character vector
                  terminals=TRUE, # whether or not to add explicit START and END events to sequence data
                  lag=1, # time offset for event pairs
                  event_sep=' ', # separator character between events
                  seq_sep = '\n', # separator character between sequences
                  residual='standard', # which type of residual to calculate (standard or adjusted)
                  critical='acceleration', # how to calculate critical residual cutoff value for notability
                  min_obs=5, # minimum number of observations for transition notability
                  lfe_collapse='default', # default value is (mean - sd)
                  cluster_height=0, # height at which to cut the hierarchical cluster tree for correlated events (equivalent to 1 - pearson's R)
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
                 n_recombine=n_recombine,
                 recombine_crit=recombine_crit, 
                 recombine_obs=recombine_obs, 
                 quiet=quiet)
  
  # process the input data:
  seqvec = file_to_seqvec(data, datatype, control$seq_sep, control$event_sep)
  
  if (control$terminals) # if you want explicit START and END events, we add those here
  {
    seqvec = add_terminals(seqvec)
  }
  
  # collapse low-frequency-events if we are doing that:
  if (control$lfe_collapse == 'default')
  {
    control$lfe_collapse = mean_minus_sd(count_events(seqvec))
  }
  if (control$lfe_collapse > 0)
  {
    lfes = get_lfes(seqvec, control)
    seqvec = collapse_lfes(seqvec, lfes)
  }
  
  # perform hierarchical clustering if asked for:
  if (control$cluster_height > 0)
  {
    par(mfrow=c(1,2))
    cluster_obj = saber_clust(seqvec, control)
    seqvec = cluster_obj$new_input
    saber_clusters = cluster_obj$clusters
    event_cor = cluster_obj$event_cor
  }
  else
  {
    par(mfrow=c(1,1))
    saber_clusters = NULL
    event_cor = NULL
  }
  
  
  # calculate transition matrices:
  obs_matrix = get_freqs(seqvec, lag)
  exp_matrix = get_exp(obs_matrix)
  res_matrix = get_residuals(obs_matrix, exp_matrix, control$residual)
  
  if (!is.numeric(control$critical)) # the user can pass a numeric value to 'critical', or a string that acts as a method for extracting the value automatically
  {
    crit_residual = res_crit(res_matrix, control$critical, plots=T)
  } else
  {
    crit_residual = control$critical
  }
  
  hft = get_high_freqs(obs_matrix, exp_matrix, res_matrix, crit_residual, control$min_obs)[[1]]
  
  
  saber_obj = list(input = data,
                   seqdata = seqvec,
                   alphabet = count_events(seqvec),
                   observed = obs_matrix,
                   expected = exp_matrix,
                   residuals = res_matrix,
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
  
  par(mfrow=c(1,1))
  
  # recursive higher order sequence analysis:
  if(n_recombine > 0)
  {
    new_input = saber.recombine(saber_obj, critval=recombine_crit, min_obs=recombine_obs)
    saber_obj = saber(new_input, 'v', terminals=FALSE, control$lag, control$event_sep, control$seq_sep, control$residual, control$critical, lfe_collapse=0, cluster_height=0, n_recombine = n_recombine-1, quiet=quiet)
    
  }
  
  return(saber_obj)
}