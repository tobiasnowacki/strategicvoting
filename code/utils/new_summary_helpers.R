# Function to get summary statistics
get_sum_stats = function(obj){
  w = obj$rcv[[1]]$weights
  map_dfr(c("rcv", "plur"), function(y){
    names(obj[[y]]) = 1:length(obj[[y]])
    map(obj[[y]], ~
          tibble(
            prev = weighted.mean(tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere, w),
            mag  = weighted.mean(.x$tau[tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere], w[tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere]),
            eb = prev * mag)) %>%
      bind_rows(.id = "iter") %>%
      mutate(system = y)
  })
}

# # Function to get summary statistics
# get_sum_stats = function(obj){
#   w = obj$rcv[[i]]$weights
#   map_dfr(c("rcv", "plur"), function(y){
#     names(obj[[y]]) = 1:length(obj[[y]])
#     map(obj[[y]], ~
#           tibble(
#             prev = weighted.mean(.x$tau > 0, w),
#             mag  = weighted.mean(.x$tau[.x$tau > 0], w[.x$tau > 0]),
#             eb = prev * mag)) %>%
#       bind_rows(.id = "iter") %>%
#       mutate(system = y)
#   })
# }