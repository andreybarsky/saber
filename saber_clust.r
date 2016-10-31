# performs hierarchical clustering on event types based on their correlations in the transition matrix:
saber_clust <- function(seqvec, control)
{
  raw_obs = get_freqs(seqvec, control$lag, simplify=FALSE)
  
  # eventrows is rows of event antecedents and sequiturs, so that they can be hierarchically clustered (treating events as cases)
  
  eventrows = as.data.frame(matrix(ncol=2*ncol(raw_obs), nrow=ncol(raw_obs)))
  rownames(eventrows) = colnames(raw_obs)
  
  for (k in colnames(raw_obs))
  {
    eventrows[k,] = c(raw_obs[k,], raw_obs[,k])
  }
  
  eventcols = t(eventrows)
  
  rowdist = as.dist(1 - cor(eventcols)) # this distance metric is somewhat arbitrary
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
