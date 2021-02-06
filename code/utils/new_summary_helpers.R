# Function to get summary statistics
get_sum_stats = function(obj){
  w = obj$rcv[[i]]$weights
  map_dfr(c("rcv", "plur"), function(y){
    names(obj[[y]]) = 1:length(obj[[y]])
    map(obj[[y]], ~
          tibble(
            prev = weighted.mean(.x$tau > 0, w),
            mag  = weighted.mean(.x$tau[.x$tau > 0], w[.x$tau > 0]),
            eb = prev * mag)) %>%
      bind_rows(.id = "iter") %>%
      mutate(system = y)
  })
}