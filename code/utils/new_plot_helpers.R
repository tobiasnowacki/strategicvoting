source("code/utils/ternary_functions.R")

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
           C = cba + cab,
           C = C + .5 * B,
           B = sqrt(3/4)*B) %>%
    dplyr::select(c(iter, case, system, state, A, B, C))
  plur_df <- do.call(rbind, plur_df) %>%
    rename(A = a, B = b, C = c) %>%
    mutate(C = C + .5 * B,
           B = sqrt(3/4)*B) %>%
    dplyr::select(c(iter, case, system, state, A, B, C))
  v_vec_df <- rbind(rcv_df, plur_df)

  

  # Create plot
  p1 <- ggplot(v_vec_df, aes(C, B)) +
    geom_line(aes(group = interaction(system, case)),
              alpha = 0.2) +
    geom_point(data = v_vec_df %>% filter(state %in% c("first", "last")),
               aes(colour = state),
               size = 0.5,
               alpha = 0.6) +
    facet_wrap(~ system) +
    theme_void() +
    geom_ternary_boundary() +
    geom_ternary_gridlines(alpha = 0.1, size = 0.3) +
    # geom_ternary_labels(vertex_labels = c("A", "B", "C")) +
    theme(legend.position = "bottom") 
   ggsave(paste0(filepath, "/v_vec_path.pdf"), 
         p1, 
         device = cairo_pdf,
         height = 4,
         width = 6)
}

non_conv_v_vec_plot <- function(obj, filepath, max_iter){
  which_cases <- which(sapply(obj, function(x) x[[7]] == 0))
  rcv_df <- list()

  if(length(which_cases) < 1){
    return(cat("No non-convergent cases"))
  }

  for(i in which_cases){
    rcv_df[[i]] <- obj[[i]][[2]] %>% mutate(iter = 1:nrow(.),
                                       case = names(obj)[i],
                                       system = "IRV",
                                       state = c("first", rep("middle", nrow(.) - 2), "last"))
  }
  rcv_df <- do.call(rbind, rcv_df)
  ggtern(rcv_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
    geom_line(alpha = .2) +
    geom_point(data = rcv_df %>% filter(state %in% c("first", "last")), aes(colour = state)) +
    facet_wrap(~ case) +
    theme_sv() +
    labs(x = "A", y = "B", z = "C",
         phase = "Iteration")

}