
get_high_freqs <- function(obs_matrix, exp_matrix, res_matrix, critical, min_obs, type='high')
{
  pairnum = 1 # iterator
  notable_transitions = data.frame(matrix(ncol=7, nrow=1), stringsAsFactors=F)
  colnames(notable_transitions) = c('Antecedent', 'Sequitur', 'Observed', 'Expected', 'Residual', 'Proportion Out', 'Proportion In')
  notable_matrix=obs_matrix
  
  # if (!is.numeric(control$critical)) # the user can pass a numeric value to 'critical', or a string that acts as a method for extracting the value automatically
  # {
  #   control$critical = res_crit(res_matrix, control$critical, plots=!control$quiet)
  # }
  
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