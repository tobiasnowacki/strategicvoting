transform_ternary = function(df){
    names(df) = c("A", "B", "C")
    df %>% 
        mutate(C = C + .5*B,
               B = sqrt(3/4)*B)
}

geom_ternary_boundary <- function(alpha = 1, col = "black"){
  vertex_df <- data.frame(x = c(1, 0, 0), y = c(0, 1, 0), z = c(0, 0, 1)) %>% 
    mutate(x = x + 0.5 * y, y = sqrt(3/4) * y)

  list(geom_polygon(data = vertex_df, ggplot2::aes(x = x, y = y), show.legend = F, alpha = alpha, col = col, fill = NA))
}

one_line <- function(val, overhang = .05){
  out <- data.frame(
    rbind(c(val, 1 - val + overhang, -overhang),
          c(val, -overhang, 1 - val + overhang)))
  colnames(out) <- c("x", "y", "z")
  out
}

one_df <- function(val, overhang = .05){
  data.frame(value = val, one_line(val, overhang))
}  

gridlines_df_for_one <- function(vals, overhang = .01){
  lapply(vals, one_df, overhang = overhang) %>% dplyr::bind_rows()
}

gridlines_df <- function(vals = c(.25, .5, .75), overhang = .01){
  g1 <- gridlines_df_for_one(vals, overhang = overhang)
  dplyr::bind_rows(
    g1 %>% dplyr::mutate(vertex = "x"),
    g1 %>% dplyr::mutate(vertex = "y") %>% 
      dplyr::select(vertex, value, x = y, y = x, z),
    g1 %>% dplyr::mutate(vertex = "z") %>% 
      dplyr::select(vertex, value, x = z, y, z = x)
  ) %>% dplyr::select(vertex, dplyr::everything())
}

geom_ternary_gridlines <- function(at = (1:3)/4, line_overhang = 0.05, text_overhang = 0.1, alpha = 0.75, col = "darkgray", size = 2){
    line_df <- gridlines_df(vals = at, overhang = line_overhang) %>% 
        dplyr::mutate(x = x + 0.5 * y, y = sqrt(3/4) * y)
    text_df <- gridlines_df(vals = at, overhang = text_overhang) %>% 
        dplyr::mutate(x = x + 0.5 * y, y = sqrt(3/4) * y) %>% 
        dplyr::distinct(vertex, value, .keep_all = TRUE)
    list(ggplot2::geom_line(data = line_df, ggplot2::aes(x = x, 
        y = y, linetype = "dotted", group = interaction(vertex, 
            value)), show.legend = FALSE, alpha = alpha, col = col), 
        ggplot2::geom_text(data = text_df, ggplot2::aes(x = x, 
            y = y, label = value), size = size))
}

geom_ternary_labels <- function(vertex_labels = c("A", "B", "C"), padding = .1, label_offset = .05){
  # vertex_labels are in "lower-left, top, right" order
  list(ggplot2::expand_limits(x = c(-padding, 1 + padding), y = sqrt(3/4)*c(-padding, 1 + padding)),
    ggplot2::annotate(geom = "text", x = c(0,.5,1), y = c(0,sqrt(3/4),0) + label_offset*c(-1,1,-1), label = vertex_labels))
}


# New function w/o ggtern for replication
joint_v_vec_plot <- function(obj, filepath){
  rcv_df <- list()
  plur_df <- list()
  # Pull cases together into one DF
  for (i in 1:length(obj)){
    rcv_df[[i]] <- obj[[i]][[2]] %>% mutate(iter = 1:nrow(.),
                                       case = names(obj)[i],
                                       system = "IRV",
                                       state = c("first", rep("middle", nrow(.) - 2), "last"))
    plur_df[[i]] <- obj[[i]][[4]] %>% mutate(iter = 1:nrow(.),
                                       case = names(obj)[i],
                                       system = "Plurality",
                                       state = c("first", rep("middle", nrow(.) - 2), "last"))
  }
  rcv_df <- do.call(rbind, rcv_df)
  plur_df <- do.call(rbind, plur_df)
  v_vec_df <- rbind(rcv_df, plur_df)

  # Create plot

  v_vec_df_tern = v_vec_df %>%
    mutate(A = V1 + V2,
           B = V3 + V4,
           C = V5 + V6,
           # ternary transformation
           C = C + .5 *B,
           B = sqrt(3/4)*B)
  return(v_vec_df_tern)
  p1 = ggplot(v_vec_df_tern, aes(x = C, y = B)) +
    geom_line(aes(group = interaction(system, case)),
              alpha = 0.2) +
    geom_point(data = v_vec_df %>% filter(state %in% c("first", "last")),
               aes(colour = state),
               size = 0.5,
               alpha = 0.6) +
    coord_fixed() +
    geom_ternary_boundary() +
    facet_wrap(~ system) +
    theme_tn()
  return(p1)
  # p1 <- ggtern(v_vec_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
  #   geom_line(aes(group = interaction(system, case)),
  #             alpha = 0.2) +
  #   geom_point(data = v_vec_df %>% filter(state %in% c("first", "last")),
  #              aes(colour = state),
  #              size = 0.5,
  #              alpha = 0.6) +
  #   facet_wrap(~ system) +
  #   theme_sv() +
  #   theme(legend.position = "bottom") +
  #   labs(x = "A", y = "B", z = "C",
  #        colour = "Iteration")
  #  ggsave(here(paste0(filepath, "/v_vec_path.pdf")), 
  #        p1, 
  #        device = cairo_pdf,
  #        height = 4,
  #        width = 6)
}
