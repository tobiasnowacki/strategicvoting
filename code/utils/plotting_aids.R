# plotting_aids.r
# Two functions that help create plots

combine <- function(share, point1, point2){
  # Point between two points, at share weighted distance
  pointnew <- (1 - share) * point1 + share * point2
  return(pointnew)
}

dial_points <- function(m_AB, m_BA, m_CB){
  # Given ratios between preferences for three FP, return co-ords for 'dials'
  centre <- c(1/3, 1/3, 1/3)
  extreme_right <- c(0.25, 0.25, 0.5)
  extreme_left <- c(0.5, 0.25, 0.25)
  extreme_top <- c(0.25, 0.5, 0.25)
  
  dial_origin_right <- c(0, 0.5, 0.5)
  dial_origin_left <- c(0.5, 0.5, 0)
  dial_origin_bottom <- c(0.5, 0, 0.5)
  
  ## A/C dials for single-peaked cases:
  if (m_CB > 0.5){
    ## Right-hand side
    d_right <- abs(m_CB - 0.5)/0.5
    bend_right <- (1 - d_right) * centre + d_right * extreme_right
    
    
    ## Left-hand side
    d_left <- abs(m_AB - 0.5)/0.5
    bend_left <- (1 - d_left) * centre + d_left * extreme_left
    
    dials1 <- data.frame(rbind(dial_origin_left, bend_left, centre, bend_right, dial_origin_right),
                         type = as.character("side"), stringsAsFactors = FALSE)
  }
  
  ## A/C dials for uniform or divided cases:
  if (m_CB <= 0.5){
    ## Right-hand side
    d_right <- abs(m_CB - 0.5)/0.5
    bend_right <- (1 - d_right) * centre + d_right * extreme_top
    
    ## Left-hand side
    d_left <- abs(m_AB - 0.5)/0.5
    bend_left <- (1 - d_left) * centre + d_left * extreme_top
    
    dials1 <- data.frame(rbind(dial_origin_left, bend_left, centre, bend_right, dial_origin_right),
                         type = as.character("side"), stringsAsFactors = FALSE)
  }
  
  ## Bottom side
  if(m_BA < 0.5){
    # Bottom dial veers to the left
    d_bottom <- abs(m_BA - 0.5) / 0.5
    bend_bottom <- (1 - d_bottom) * centre + d_bottom * extreme_left
  }
  if(m_BA >= 0.5){
    # Bottom dial veers to the right (or is straight)
    d_bottom <- abs(m_BA - 0.5) / 0.5
    bend_bottom <- (1 - d_bottom) * centre + d_bottom * extreme_right
  }
  dials1[6:8, 1:3] <- data.frame(rbind(dial_origin_bottom, bend_bottom, centre))
  dials1[6:8, 4] <- as.character("bottom")
  return(dials1)
}