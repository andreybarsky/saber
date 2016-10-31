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


### Calculating residual critical value:

res_crit <- function(res_matrix, type='acceleration', plots=T)
{ # the 'plots' local variable is inherited by the subfunctions because R doesn't like recursive default argument references
  
  # first, some quick geometric functions:
  
  get_distance <- function(x1, y1, x2, y2) # given two points, get euclidean distance using pythagorean theorem
  { sqrt((x1 - x2)^2 + (y1-y2)^2) }
  
  get_line <- function(x1, y1, x2, y2)  # given two points, get slope and intercept of the line that connects them
  { slope = (y2-y1)/(x2-x1) ; intercept = y1-(x1*slope) ; list(slope=slope, intercept=intercept) }
  
  get_intersection <- function(m1, c1, m2, c2) # given two lines, get the point at which they intersect
  {  x = (c2-c1)/(m1-m2) ; y = (m1 * x) + c1 ; list(x=x, y=y) }
  
  
  # vector projection method:
  res_crit_projection <- function(res_matrix)
  {
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
  res_crit_acceleration <- function(res_matrix, smoothing=0.5) # looks for the sharpest turn in the smoothed scree curve. smoothing is the loess smooth parameter
  {
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
  

  
  if(type=='acceleration')
  {
    return(res_crit_acceleration(res_matrix))
  }  else if(type=='projection')
  {
    return(res_crit_projection(res_matrix))
  }  else
  {
    warning("Invalid method for automatically calculating critical residual value (try 'acceleration' or 'projection', or a numeric value)")
  }
}
