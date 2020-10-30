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

