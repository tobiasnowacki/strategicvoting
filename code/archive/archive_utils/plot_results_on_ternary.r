### plotting results in different systems 

simplex.x <- function(x, omitted = 0){
  if(sum(x) == 0){return(.5)}
  return( (1 - omitted)*(x[2] + 0.5 * x[3]) / sum(x) + omitted/2)
}
simplex.y <- function(x, omitted = 0){
  if(sum(x) == 0){return(sqrt(.75)*(1/3))}  
  return( (1 - omitted)*(sqrt(0.75) *  x[3]) / sum(x) + omitted/2)
} 

add.ternary.point = function(point, col = "black", pch = 19, cex = 1){
  points(x = simplex.x(point), y = simplex.y(point), pch = pch, col = col, cex = cex)
}

add.ternary.text = function(point, labels, col = "black", cex = 1, x.offset = 0, y.offset = 0){
  text(x = simplex.x(point) + x.offset, y = simplex.y(point) + y.offset, labels = labels, col = col, cex = cex)
}

add.ternary.lines = function(point.1, point.2, col = "black", lwd = 1, lty = 1){
  lines(x = c(simplex.x(point.1), simplex.x(point.2)), y = c(simplex.y(point.1), simplex.y(point.2)), col = col, lwd = lwd, lty = lty)
}


add.ternary.polygon = function(point.mat, border = NULL, col = NA){
  xs = apply(point.mat, 1, simplex.x)
  ys = apply(point.mat, 1, simplex.y)
  polygon(xs, ys, border = border, col = col)
}

add.ternary.boundary = function(k = 0){
  bca.vertex = c(1 - k, 0, 0)
  cba.vertex = c(0, 1 - k, 0)
  bac.vertex = c(0, 0, 1-k)
  bac.v.x = simplex.x(bac.vertex, k); bac.v.y = simplex.y(bac.vertex, k)
  bca.v.x = simplex.x(bca.vertex, k); bca.v.y = simplex.y(bca.vertex, k)
  cba.v.x = simplex.x(cba.vertex, k); cba.v.y = simplex.y(cba.vertex, k) 
  
  lines(c(bac.v.x, cba.v.x), c(bac.v.y, cba.v.y))
  lines(c(bca.v.x, cba.v.x), c(bca.v.y, cba.v.y))
  lines(c(bac.v.x, bca.v.x), c(bac.v.y, bca.v.y))
}


library(colorspace)


plot.plurality.result = function(v.vec, from.v.vec = NULL, vertex.labels = c("A", "B", "C"), shading.cols = rainbow_hcl(3, alpha = .4), main = NULL, new = T, border = "black", fp.result.col = "black", fp.result.cex = 1, space = .1, add.fp.result = F){
  # the basic plot
  if(new){
    xs = c(0, 1) + c(-space, space); ys = c(0, sqrt(3/4)) + sqrt(3/4)*c(-space, space)
    plot(xs, ys, type = "n", xlab = "", ylab = "", axes = F, main = main)
    add.ternary.boundary()
  }
  
  # the regions
  A.vertex = c(1,0,0)
  B.vertex = c(0,0,1)
  C.vertex = c(0,1,0)
  AB.tie.C.0 = c(1/2, 0, 1/2)
  AC.tie.B.0 = c(1/2, 1/2, 0)
  BC.tie.A.0 = c(0, 1/2, 1/2)
  midpoint = c(1/3, 1/3, 1/3)
  A.point.mat = rbind(A.vertex, 
                    AB.tie.C.0,
                    midpoint,
                    AC.tie.B.0,
                    A.vertex
  )		
  B.point.mat = rbind(B.vertex, 
                      AB.tie.C.0,
                      midpoint,
                      BC.tie.A.0,
                      B.vertex
  )		
  C.point.mat = rbind(C.vertex, 
                      AC.tie.B.0,
                      midpoint,
                      BC.tie.A.0,
                      C.vertex
  )		
  add.ternary.polygon(A.point.mat, col = shading.cols[1], border = border)
  add.ternary.polygon(B.point.mat, col = shading.cols[2], border = border)
  add.ternary.polygon(C.point.mat, col = shading.cols[3], border = border)
  
  # the result
  fp.vec = c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
  
  # FP result
  if(add.fp.result){
    add.ternary.point(fp.vec[c(1,3,2)], pch = 19, cex = fp.result.cex, col = fp.result.col)
  }
  
  if(new){
    label.offset = .05
    add.ternary.text(c(1,0,0), vertex.labels[1], x.offset = -label.offset, y.offset = -sqrt(3/4)*label.offset)
    add.ternary.text(c(0, 0,1), vertex.labels[2], x.offset = 0, y.offset = sqrt(3/4)*label.offset)
    add.ternary.text(c(0,1,0), vertex.labels[3], x.offset = label.offset, y.offset = -sqrt(3/4)*label.offset)
  }
  
}

# was called plot.v.vec in original version
plot.av.result = function(v.vec, from.v.vec = NULL, secondary.line.col = "gray", secondary.line.lwd = 2, vertex.labels = c("A", "B", "C"), shading.cols = rainbow_hcl(3, alpha = .4), main = NULL, new = T, border = "black", fp.result.col = "black", fp.result.cex = 1, space = .1, clipped = F, add.fp.result = F){
  # v.vec is in AB, AC, BA, BC, CA, CB, AX, BX, CX order 
  if(length(v.vec == 6)){
    v.vec = c(v.vec, 0, 0, 0)
  }
  # B is the top vertex
  
  # the basic plot
  if(new){
    xs = c(0, 1) + c(-space, space); ys = c(0, sqrt(3/4)) + sqrt(3/4)*c(-space, space)
    if(clipped){
      xs = c(1/4, 3/4) + c(-space/2, space/2)
      ys = c(0, sqrt(3/4)/2) + sqrt(3/4)*c(-space/2, space/2)
    }
    plot(xs, ys, type = "n", xlab = "", ylab = "", axes = F, main = main)
    add.ternary.boundary()
    
    # majority thresholds
    add.ternary.lines(c(1/2, 1/2, 0), c(0, 1/2, 1/2), col = secondary.line.col, lwd = secondary.line.lwd, lty = 3)
    add.ternary.lines(c(0, 1/2, 1/2), c(1/2, 0, 1/2), col = secondary.line.col, lwd = secondary.line.lwd, lty = 3)
    add.ternary.lines(c(1/2, 0, 1/2), c(1/2, 1/2, 0), col = secondary.line.col, lwd = secondary.line.lwd, lty = 3)
    
    # first-round pivotal events
    add.ternary.lines(c(1/2, 1/4, 1/4), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)
    add.ternary.lines(c(1/4, 1/2, 1/4), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)
    add.ternary.lines(c(1/4, 1/4, 1/2), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)		
    
  }
  
  
  # the first preference shares
  fp.vec = c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
  
  if(!is.null(from.v.vec)){
    from.fp.vec = c(sum(from.v.vec[c(1,2,7)]), sum(from.v.vec[c(3,4,8)]), sum(from.v.vec[c(5,6,9)]))
    plot.av.result(v.vec = from.v.vec, new = F, shading.cols = c(NULL, NULL, NULL), border = "gray", fp.result.col = "black", fp.result.cex = .5)
    add.ternary.lines(fp.vec[c(1,3,2)], from.fp.vec[c(1,3,2)], col = "black")		
  }
  
  
  # second-round pivotal events 
  mAB = (v.vec[1] - v.vec[2])/fp.vec[1]
  mBA = (v.vec[3] - v.vec[4])/fp.vec[2]
  mCA = (v.vec[5] - v.vec[6])/fp.vec[3]
  
  # Determining key points
  AB.tie.C.0 = c(1/2, 0, 1/2)
  AC.tie.B.0 = c(1/2, 1/2, 0)
  BC.tie.A.0 = c(0, 1/2, 1/2)
  
  if(mAB > 0){
    denom.pos = 3 + mAB
    BC.tie.A.max = c(1/denom.pos, 1 - 2/denom.pos, 1/denom.pos)
  }else{
    denom.neg = 3 - mAB
    BC.tie.A.max = c(1/denom.neg, 1/denom.neg, 1 - 2/denom.neg)	
  }
  
  if(mBA > 0){
    denom.pos = 3 + mBA
    AC.tie.B.max = c(1/denom.pos, 1 - 2/denom.pos, 1/denom.pos)
  }else{
    denom.neg = 3 - mBA
    AC.tie.B.max = c(1 - 2/denom.neg, 1/denom.neg, 1/denom.neg)	
  }
  
  if(mCA > 0){
    denom.pos = 3 + mCA
    AB.tie.C.max = c(1/denom.pos, 1/denom.pos, 1 - 2/denom.pos)
  }else{
    denom.neg = 3 - mCA
    AB.tie.C.max = c(1 - 2/denom.neg, 1/denom.neg, 1/denom.neg)	
  }
  
  # now to draw this. 
  
  # A winning area 
  point.mat = rbind(
    AC.tie.B.0,
    c(1,0,0), 
    AB.tie.C.0,
    AB.tie.C.max) # this will always be the order
  if(mCA > 0 | mBA > 0){
    point.mat = rbind(point.mat, c(1/3, 1/3, 1/3))
  }
  point.mat = rbind(point.mat, 
                    AC.tie.B.max,
                    AC.tie.B.0
  )		
  add.ternary.polygon(point.mat, col = shading.cols[1], border = border)
  
  # B winning area 
  point.mat = rbind(
    AB.tie.C.0,
    c(0,0,1), 
    BC.tie.A.0,
    BC.tie.A.max) # this will always be the order
  if(mAB > 0 | mCA < 0){
    point.mat = rbind(point.mat, c(1/3, 1/3, 1/3))
  }
  point.mat = rbind(point.mat, 
                    AB.tie.C.max,
                    AB.tie.C.0
  )		
  add.ternary.polygon(point.mat, col = shading.cols[2], border = border)
  
  # C winning area 
  point.mat = rbind(
    BC.tie.A.0,
    c(0,1,0), 
    AC.tie.B.0,
    AC.tie.B.max) # this will always be the order
  if(mBA < 0 | mAB < 0){
    point.mat = rbind(point.mat, c(1/3, 1/3, 1/3))
  }
  point.mat = rbind(point.mat, 
                    BC.tie.A.max,
                    BC.tie.A.0
  )		
  add.ternary.polygon(point.mat, col = shading.cols[3], border = border)
  
  # FP result
  if(add.fp.result){
    add.ternary.point(fp.vec[c(1,3,2)], pch = 19, cex = fp.result.cex, col = fp.result.col)
  }
  
  if(new){
    label.offset = .05
    if(clipped){
      add.ternary.text(c(3/4,1/4,0), vertex.labels[1], x.offset = -label.offset, y.offset = -sqrt(3/4)*label.offset)
      # masking box
      y.inc = label.offset/2
      polygon(c(-100, 100, 100, -100), sqrt(3/4)*c(1/2 + y.inc, 1/2 + y.inc, 100, 100), col = "white", border = F)
      add.ternary.text(c(1/4, 1/4,1/2), vertex.labels[2], x.offset = 0, y.offset = sqrt(3/4)*label.offset)
      add.ternary.text(c(1/4,3/4,0), vertex.labels[3], x.offset = label.offset, y.offset = -sqrt(3/4)*label.offset)
    }else{
      add.ternary.text(c(1,0,0), vertex.labels[1], x.offset = -label.offset, y.offset = -sqrt(3/4)*label.offset)
      add.ternary.text(c(0, 0,1), vertex.labels[2], x.offset = 0, y.offset = sqrt(3/4)*label.offset)
      add.ternary.text(c(0,1,0), vertex.labels[3], x.offset = label.offset, y.offset = -sqrt(3/4)*label.offset)
    }
    
  }
  
}

plot.condorcet.result = function(v.vec, from.v.vec = NULL, secondary.line.col = "gray", secondary.line.lwd = 2, vertex.labels = c("A", "B", "C"), shading.cols = rainbow_hcl(3, alpha = .4), main = NULL, new = T, border = "black", fp.result.col = "black", fp.result.cex = 1, space = .1, clipped = F, add.fp.result = F){
  
  # basic plot
  if(new){
    xs = c(0, 1) + c(-space, space); ys = c(0, sqrt(3/4)) + sqrt(3/4)*c(-space, space)
    plot(xs, ys, type = "n", xlab = "", ylab = "", axes = F, main = main)
    add.ternary.boundary()
  }
  
  # FP vote shares
  fp.vec = c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
  
  # margins 
  mBC = (v.vec[1] - v.vec[2])/fp.vec[1]  # to what extend does A's second pref favor B over C?
  mAC = (v.vec[3] - v.vec[4])/fp.vec[2]  # to what extend does B's second pref favor A over C?
  mAB = (v.vec[5] - v.vec[6])/fp.vec[3]  # to what extend does C's second pref favor A over B?
  
  # majority two-way tie points
  AB.tie.C.0 = c(1/2, 0, 1/2)
  AC.tie.B.0 = c(1/2, 1/2, 0)
  BC.tie.A.0 = c(0, 1/2, 1/2)
  
  # vertices
  A.vertex = c(1,0,0)
  B.vertex = c(0,0,1)
  C.vertex = c(0,1,0)
  
  # AB majority tie line 
  ab.point = c(0, 1 - mAB/(1 + mAB), mAB/(1 + mAB))
  if(mAB < 0){
    ab.point = c(1 - 1/(1-mAB), 1/(1-mAB), 0)
  }
  add.ternary.lines(AB.tie.C.0, ab.point, col = secondary.line.col)
  
  # AC majority tie line 
  ac.point = c(0, mAC/(1 + mAC), 1 - mAC/(1 + mAC))
  if(mAC < 0){
    ac.point = c(1 - 1/(1-mAC), 0, 1/(1-mAC))
  }
  add.ternary.lines(AC.tie.B.0, ac.point, col = secondary.line.col)
  
  # BC majority tie line 
  bc.point = c(1/(1+mBC), 1 - 1/(1+mBC), 0)
  if(mBC < 0){
    bc.point = c(1/(1 - mBC), 0, 1 - 1/(1-mBC))
  }
  add.ternary.lines(BC.tie.A.0, bc.point, col = secondary.line.col)
  
  # intersections 
  AB.AC.int.B = ((1 + mAB)/2)/(1 + mAB + (1+mAC)*(1-mAB)/2)
  AB.AC.int.A = (AB.AC.int.B*(1 + mAB) - mAB)/(1 - mAB)
  AB.AC.int = c(AB.AC.int.A, 1 - AB.AC.int.A - AB.AC.int.B, AB.AC.int.B)
  
  AB.BC.int.A = (1 - (2*mAB/(1 + mAB)))/(2*(1-mAB)/(1 + mAB) + 1 + mBC)
  AB.BC.int.B = 1/2*(1 - AB.BC.int.A*(1 + mBC))
  AB.BC.int = c(AB.BC.int.A, 1 - AB.BC.int.A - AB.BC.int.B, AB.BC.int.B)
  
  AC.BC.int.A = ((1 - mAC)/2)/(2 - (1 + mBC)*(1 + mAC)/2)
  AC.BC.int.B = (1/2)*(1 - AC.BC.int.A*(1 + mBC))
  AC.BC.int = c(AC.BC.int.A, 1 - AC.BC.int.A - AC.BC.int.B, AC.BC.int.B)
    
  # winning areas 
  add.ternary.polygon(rbind(A.vertex, AC.tie.B.0, AB.AC.int, AB.tie.C.0, A.vertex), col = shading.cols[1], border = border)
  add.ternary.polygon(rbind(B.vertex, AB.tie.C.0, AB.BC.int, BC.tie.A.0, B.vertex), col = shading.cols[2], border = border)
  add.ternary.polygon(rbind(C.vertex, AC.tie.B.0, AC.BC.int, BC.tie.A.0, C.vertex), col = shading.cols[3], border = border)
  
  # FP result
  if(add.fp.result){
    add.ternary.point(fp.vec[c(1,3,2)], pch = 19, cex = fp.result.cex, col = fp.result.col)
  }
  
  if(new){
    label.offset = .05
    add.ternary.text(c(1,0,0), vertex.labels[1], x.offset = -label.offset, y.offset = -sqrt(3/4)*label.offset)
    add.ternary.text(c(0, 0,1), vertex.labels[2], x.offset = 0, y.offset = sqrt(3/4)*label.offset)
    add.ternary.text(c(0,1,0), vertex.labels[3], x.offset = label.offset, y.offset = -sqrt(3/4)*label.offset)
  }
  
}

# side-by-side comparison 
#basic.v.vec = c(20, 6, 13,10, 4, 21, 2,4,3)
#v.vec = basic.v.vec/sum(basic.v.vec)
#compare_three_systems(v.vec)

compare_three_systems = function(v.vec){
  par(mfrow = c(1,3))
  plot.plurality.result(v.vec, main = "plurality")
  plot.av.result(v.vec, main = "AV")
  plot.condorcet.result(v.vec, main = "Condorcet")
}


## Next steps: 
# shiny app (for these three or more -- text entries for counts, or arrows?)
# Borda count (or positional methods in general)
# overlays to show where they disagree

# writeup: 
