
make.plot.of.saturated.regressions = function(reg.list, 
                                              x.lab = "Relative tactical incentive (deciles)", 
                                              y.lab = "Proportion voting tactically", 
                                              y.lims = c(0,1), 
                                              y.line.loc = 7.5, 
                                              pchs = rep(19, length(reg.list)), 
                                              offsets = seq(-.15, .15, length = length(reg.list)), 
                                              pt.cex = .5, 
                                              point.cols = rep("black", length(reg.list)), 
                                              confint.cols = point.cols, 
                                              line.cols = point.cols, 
                                              confint.ltys = rep(2, length(reg.list)), 
                                              ltys = 1:(length(reg.list)), 
                                              x.vals = 1:length(coef(reg.list[[1]])),  # equal spacing of the points
                                              plot.confints = F,
                                              error.bar.confint.plot = T,
                                              default.axis = F,
                                              axis.with.no.labels = F,
                                              x.labels = paste0("D", x.vals),
                                              tones = c(.3, .4, .5, .6, .7),
                                              lwds = rep(1, length(reg.list))){
  
  plot(x = range(x.vals), y = c(0,1), type = "n", axes = F, xlab = x.lab, ylab = y.lab, ylim = y.lims)
  if(default.axis){
    axis(1)
  }else if(axis.with.no.labels){
    axis(1, at = c(-100000000, 10000000), labels = c("Should not", "see these"))
  }else{
    axis(1, at = x.vals, labels = x.labels) 		
  }
  axis(2)
  abline(v = y.line.loc, lty = 2, col = "gray") # identifies where tau turns positive
  
  for(j in 1:length(reg.list)){
    # plot point estimates 
    points(x.vals + offsets[j], coef(reg.list[[j]]), pch = pchs[j], col = point.cols[j], cex = pt.cex)
    # connect them
    lines(x = x.vals + offsets[j], y = coef(reg.list[[j]]), col = line.cols[j], lty = ltys[j], lwd = lwds[j])  # 1:10
    # plot the confidence intervals
    if(plot.confints){
      # I want to do this with polygon instead. 
      if(error.bar.confint.plot){
        for(i in 1:length(x.vals)){
          lines(x = c(x.vals[i], x.vals[i]) + offsets[j], y = confint(reg.list[[j]])[i,], lty = confint.ltys[j], col = confint.cols[j])
        }							
      }else{
        ci.mat = confint(reg.list[[j]])
        the.tone = tones[j] # .25 + .12*j
        polygon(x = c(x.vals, rev(x.vals)) + offsets[j], y = c(ci.mat[,1], rev(ci.mat[,2])), col = rgb(the.tone, the.tone, the.tone, alpha = .3), border = F)  # 1:nrow(ci.mat), nrow(ci.mat):1) + offsets[j]
      }
    }				
  }
  
}