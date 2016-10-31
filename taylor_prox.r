min_dist <- function(seq, a, b) # for each A, how close is the next B
{
  splitseq = strsplit(seq, ' ')[[1]]
  a_pos = which(splitseq == a)
  b_pos = which(splitseq == b)
  distances = numeric()
  
  if(a == b) # special case where A and B are the same
  {
    return(diff(a_pos))
  }
  
  for (a in 1:length(a_pos)) # for every occurrence of A
  {
    if (!(length(which(b_pos > a_pos[a])) == 0)) # not if no B occurs after this A...
    {
      next_b = b_pos[which(b_pos > a_pos[a])[1]]
      
      difference = next_b - a_pos[a]
      distances = append(distances, difference)
    }
  }
return (distances)
}

proximity <- function(seq, a, b, w=1)
{

  splitseq = strsplit(seq, ' ')[[1]]
  min_dists = min_dist(seq, a, b)
  
  if (length(min_dists) == 0)
  {
    return(NA)
  }
  
  proximity = 1 - (
    sum(w * (min_dists-1)) / ( # optional weighting parameter
      length(min_dists) * (
        length(splitseq)-2
      )
    )
  )
  
  return(proximity)
}


prox_matrix <- function(seq, weight=1) # calculate proximity matrix for a single sequence
{
  splitseq = strsplit(seq, ' ')[[1]]
  eventtypes = names(table(strsplit(seq, ' ')[[1]]))
  
  pmatrix = matrix(0, nrow = length(eventtypes), ncol = length(eventtypes), dimnames = list(eventtypes, eventtypes)) # empty trans freq matrix
  
  for (r in 1:length(eventtypes))
  {
    for (c in 1:length(eventtypes))
    {
      pmatrix[r,c] = proximity(seq, eventtypes[r], eventtypes[c], weight)
    }
  }
  return(pmatrix)
}

exp_prox <- function(seq, n) # calculate expected proximity matrix by resampling sequences many times
{
  splitseq = strsplit(seq, ' ')[[1]]
  eventtypes = names(table(strsplit(seq, ' ')[[1]]))
  
  # empty matrix:
  exp_matrix = matrix(0, 
                      nrow = length(eventtypes), 
                      ncol = length(eventtypes), 
                      dimnames = list(eventtypes, eventtypes)) 
  
  for (i in 1:n)
  {
    new_matrix = prox_matrix(paste(sample(splitseq), collapse=' '))
    new_matrix[which(is.na(new_matrix))] = 0 # replace NAs with 0 
    
    exp_matrix = exp_matrix + new_matrix
  }
  
  return (exp_matrix/n)
}