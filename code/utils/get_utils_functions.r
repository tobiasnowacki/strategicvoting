return_max_name <- function(x){
     x
     pos <- which(x == max(x))
     if (length(pos) > 1){pos <- 99}
     return(pos)
}

return_sec_name <- function(y){
     out <- which(y == as.numeric(sort(y)[2]))
     if (length(out) > 1){out <- 99}
     return(out)
}

return_third_name <- function(y){
     out <- which(y == as.numeric(sort(y)[1]))
     if (length(out) > 1){out <- 99}
     return(out)
}

first_second_difference <- function(x){
  # Function that yields the difference between the first and the second preference;
  # in terms of score (NVM utilities)
  max(x) - max(x[-which(x == max(x))])
}
