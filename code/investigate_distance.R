library(tidyverse)

# Load vvec data
load("output/files/1/85_vvec.Rdata")
load("output/files/1/85_sum.Rdata")


plot_v_vec_distance <- function(obj, filepath, 
                                custom_red = 0.0001, 
                                n_lag = NULL,
                                avg_span){
  
  ################################################
  # Takes the full object and creates the following
  # six distance plots:
  #   "/dist_irv" -- v-vec(t - 1) to v-vec(t) under IRV
  #   "/dist_plur" -- v-vec(t - 1) to v-vec(t) under Plur
  #   "/dist_irv_lagged" -- v-vec(lag) to v-vec(t) under IRV
  #   "/dist_plur_lagged" -- v-vec(lag) to v-vec(t) under Plur
  #   "/dist_br" -- v-vec(t) to BR(t) under both
  # "/dist_br_lagged" -- v-vec(t) to BR(t) under both
  ################################################
  
  cat("V.VEC TO V.VEC DISTANCE \n")
  cat("Calculating distances. \n")
  
  conv_dist <- list()
  
  for (j in 1:length(obj)){
    # print(j)
    df <- obj[[j]][[1]] 
    df_dist <- euclid_together(df) %>%
      select(diff:diff_first) %>%
      mutate(system = 'IRV',
             case = nn[j],
             iter = row_number())
    if(is.numeric(n_lag) == TRUE){
      df_dist$avg_dist <- calc_lag_dist(df, n_lag, avg_span)
    }
    conv_dist[[(2 * j) - 1]] <- df_dist
    
    df <- obj[[j]][[2]] 
    df_dist <- euclid_together(df) %>%
      select(diff:diff_first) %>%
      mutate(system = 'Plurality',
             case = nn[j],
             iter = row_number())
    if(is.numeric(n_lag) == TRUE){
      df_dist$avg_dist <- calc_lag_dist(df, n_lag, avg_span)
    }
    conv_dist[[2 * j]] <- df_dist
  }
  
  # return(conv_dist)
  conv_dist_df <- do.call(rbind, conv_dist)
  
  cat("Distances done. Plotting now. \n")
  
  # Plot IRV convergence
  p1 <- ggplot(conv_dist_df %>% filter(system == "IRV"), aes(iter, log(diff))) +
    geom_line(aes(group = interaction(case, system)), 
              alpha = 0.1,
              color = "blue") +
    labs(x = "Iteration (Strategic Responsiveness)",
         y = "ln(Distance)")
  ggsave(paste0(filepath, "/dist_irv.pdf"), 
         p1, 
         width = 4, 
         height = 4)
  
  # Plot distance to lagged average
  p1_lag <- ggplot(conv_dist_df %>% filter(system == "IRV"), 
                   aes(iter, log(avg_dist))) +
    geom_line(aes(group = interaction(case, system)), 
              alpha = 0.1,
              color = "blue") +
    labs(x = "Iteration (Strategic Responsiveness)",
         y = "ln(Distance to lagged average)") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_irv_lagged.pdf"), 
         p1_lag, 
         width = 4,
         height = 4)
  
  # Plot plurality convergence
  p2 <- ggplot(conv_dist_df %>% filter(system == "Plurality"), aes(iter, log(diff))) +
    geom_line(aes(group = interaction(case, system)), 
              alpha = 0.1, 
              color = "blue") +
    labs(x = "Iteration (Strategic Responsiveness)",
         y = "ln(Distance)") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_plur.pdf"), 
         p2, 
         width = 4,
         height = 4)
  
  ################################################
  
  cat("V.VEC TO BEST RESPONSE DISTANCE \n")
  
  # compare best response vector to polling v_vec
  br_dist_list <- list()
  for(c in 1:length(obj)){
    case_name <- names(big_list_na_omit)[c]
    br_dist_list[[c]] <- br_distance(obj[[c]],
                                     n_lag = 20,
                                     avg_span = 10,
                                     case = case_name)
  }
  br_dist <- do.call(rbind, br_dist_list)
  names(br_dist)[1] <- 'd'
  return(br_dist)
  # plot results
  cat("Distances done. Plotting now. \n")
  
  # t to t-1
  plot_br_1 <- ggplot(br_dist %>% filter(iter > 0, case %in% c("FRA_2002", "SWE_2014")), aes(iter, log(d))) +
    geom_line(aes(group = case), alpha = 0.1) +
    facet_wrap(. ~ system) +
    labs(x = "Iteration", y = "ln(Distance)") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_br.pdf"),
         plot_br_1,
         height = 4,
         width = 8)
  
  # t to avg(t-10 : t-20)
  plot_br_2 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, log(l_avg))) +
    geom_line(aes(group = case), alpha = 0.1) +
    facet_wrap(. ~ system) +
    labs(x = "Iteration", y = "ln(Lagged Distance)") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_br_lag.pdf"),
         plot_br_2,
         height = 4,
         width = 8)
  
  # linear, rather than logged
  plot_br_3 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, d)) +
    geom_line(aes(group = case), alpha = 0.1) +
    facet_wrap(. ~ system) +
    labs(x = "Iteration", y = "Distance") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_br_lin.pdf"),
         plot_br_3,
         height = 4,
         width = 8)
  
  plot_br_4 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, l_avg)) +
    geom_line(aes(group = case), alpha = 0.1) +
    facet_wrap(. ~ system) +
    labs(x = "Iteration", y = "Lagged Distance") +
    theme_tn()
  ggsave(paste0(filepath, "/dist_br_lag_lin.pdf"),
         plot_br_4,
         height = 4,
         width = 8)
  
  # save(br_dist, file = here(paste0(filepath, "v_vecs_br.Rdata")))
  
  # if(is.numeric(n_lag) == TRUE ){
  #   return(conv_dist_df)
  # }
  
}

# produce a vector of distance elements for each row of a given v_vec DF
get_distance_by_row <- function(df){
  out <- df^2 %>% rowSums %>% sqrt
  return(out)
}

# produce DF with averaged values
get_lagged_avg <- function(df, avg_span, n_lag){
  new_mat <- matrix(NA, nrow = nrow(df), ncol = ncol(df))
  for(i in 1:(nrow(df) - n_lag)){
    new_mat[i + n_lag, ] <- colMeans(df[i:(i + avg_span), ])
  }
  return(new_mat)
}

# Neat function to calculate distances: vvec and best response
# could modify to accommodate vvec to vvec too...
br_distance <- function(obj, n_lag, avg_span, case = "NA"){
  br_rcv <- obj[[3]] %>% as.matrix
  d_rcv <- (obj[[1]] - br_rcv) %>%
    get_distance_by_row(.) %>%
    as.data.frame %>%
    mutate(iter = 0:(length(.) - 1),
           system = "IRV",
           case = case)
  avg_rcv <- get_lagged_avg(obj[[1]], 
                            avg_span = avg_span, 
                            n_lag = n_lag)
  d_rcv_avg <- (avg_rcv - br_rcv) %>% 
    get_distance_by_row(.)
  d_rcv$l_avg <- d_rcv_avg
  
  br_plur <- obj[[4]] %>% as.matrix
  # br_plur = br_plur[, c(1, 3, 5)] + br_plur[, c(2, 4, 6)]
  d_plur <- (obj[[2]] - br_plur) %>%
    get_distance_by_row(.) %>%
    as.data.frame %>%
    mutate(iter = 0:(length(.) - 1),
           system = "Plurality",
           case = case)     
  avg_plur <- get_lagged_avg(obj[[2]],
                             avg_span = avg_span,
                             n_lag = n_lag)
  d_plur$l_avg <- (avg_plur - br_plur) %>%
    get_distance_by_row(.)
  return(rbind(d_rcv, d_plur))
}

euclid <- function(x) {
  # Euclidean distance between one iteration and next
  if (isdf <- is.data.frame(x)) {
    x <- data.matrix(x)
  }
  dij <- c(NA, sqrt(rowSums((tail(x, -1) - head(x, -1))^2)))
  x <- cbind(x, diff = dij)
  if (isdf) {
    x <- as.data.frame(x)
  }
  x
}


euclid_first <- function(df){
  # Euclidean distance between one iteration and first
  df_follow <- df[2:nrow(df), ]
  dist <- apply(df_follow, 1, function(x) sqrt(sum((x - df[1, ])^2)))
  return(dist)
}

euclid_together <- function(df){
  a <- euclid(df)
  diff_first <- c(NA, euclid_first(df))
  return(cbind(a, diff_first))
}

calc_lag_dist <- function(df, n_lag, avg_span){
  # calculates
  sqsums <- rowSums(df^2)
  max_dist <- nrow(df) - n_lag
  d <- rep(NA, n_lag)
  for(i in 1:max_dist){
    avg_d <- sum(colMeans(df[i:(i + avg_span), ])^2) 
    avg_d <- sqsums[i + n_lag - 1]
    d[i + n_lag] <- sqrt(abs(avg_d - sqsums[i + n_lag]))
  }
  return(d)
}

joint_v_vec_plot <- function(obj, filepath){
  rcv_df <- list()
  plur_df <- list()
  # Pull cases together into one DF
  for (i in 1:length(obj)){
    rcv_df[[i]] <- obj[[i]][[1]] %>% mutate(iter = 1:nrow(.),
                                            case = names(obj)[i],
                                            system = "IRV",
                                            state = c("first", rep("middle", nrow(.) - 2), "last"))
    plur_df[[i]] <- obj[[i]][[2]] %>% mutate(iter = 1:nrow(.),
                                             case = names(obj)[i],
                                             system = "Plurality",
                                             state = c("first", rep("middle", nrow(.) - 2), "last"))
  }
  rcv_df <- do.call(rbind, rcv_df) %>%
    mutate(A = abc + acb,
           B = bac + bca,
           C = cba + cab) %>%
    dplyr::select(c(iter, case, system, state, A, B, C))
  plur_df <- do.call(rbind, plur_df) %>%
    rename(A = a, B = b, C = c) %>%
    dplyr::select(c(iter, case, system, state, A, B, C))
  v_vec_df <- rbind(rcv_df, plur_df)
  
  
  library(ggtern)
  # Create plot
  p1 <- ggtern(v_vec_df %>% filter(case == "SWE_2014" | case == "FRA_2002"), 
               aes(A, B, C)) +
    geom_line(aes(group = interaction(system, case)),
              alpha = 0.1) +
    geom_point(data = v_vec_df %>% filter(state %in% c("first", "last")) %>% filter(case == "SWE_2014" | case == "FRA_2002"),
               aes(colour = state),
               size = 1.2,
               alpha = 0.5) +
    facet_wrap(~ system) +
    scale_colour_manual(values = c("red", "blue")) +
    # geom_ternary_labels(vertex_labels = c("A", "B", "C")) +
    theme(panel.spacing = unit(0, "cm"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
  return(p1)
}

brobj = plot_v_vec_distance(vvecdf, 
                    paste0("output/figures/", lambda, "/", s), 
                    n_lag = 20, avg_span = 10)

joint_v_vec_plot(vvecdf ,  paste0("output/figures/", lambda, "/", s))

ggplot(sum_df %>% filter(case %in% c("SWE_2014", "FRA_2002"), system == "rcv"), aes(as.numeric(iter), prev)) +
  geom_point(aes(colour = case)) +
  geom_line(aes(colour = case))
  
