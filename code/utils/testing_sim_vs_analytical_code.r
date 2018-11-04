compute.sim.based.and.analytical.piv.probs.for.J.scenarios = function(hyper.v.vec = c(6,4,5,5,4,6,1,1,1), J = 50, tol = .005, hyper.s = 20, lower.s = 15, higher.s = 60, types = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.AB", "AC.BC", "BC.BC", "BC.BA", "BC.AC"), Ms = c(100, 500, 1000, 5000, 10000, 50000, 100000), increments = 50, increments.at.end = 10, increments.for.dir.mat = 20, different.second.pref.precision = F, skip.sims = F, expected.result.mat = NULL, s.mat = NULL){
	stopifnot(length(hyper.v.vec) == 9)
	stopifnot(length(hyper.s) == 1)
	if(is.null(expected.result.mat)){expected.result.mat = rdirichlet(n = J, alpha = hyper.s*hyper.v.vec/sum(hyper.v.vec))}else{cat("Using expected result mat passed in arguments.\n")}
	if(is.null(s.mat)){s.mat = matrix(runif(n = J*4, lower.s, higher.s), nrow = J, ncol = 4)}else{cat("Using s.mat from arguments.\n")}
	# types.sim = tolower(types)
	stub = matrix(NA, nrow = J, ncol = length(types))
	colnames(stub) = types
	sim.results = analytical.results = stub
	ar.list = sim.list = list()
	for(M in Ms){
		cat("=============", M, "============\n")
		for(j in 1:J){
			if(j%%10 == 0){cat(j)}else{cat(".")}
			v.vec = expected.result.mat[j,]
			if(different.second.pref.precision){
				s.vec = s.mat[j,]
			}else{
				s.vec = rep(s.mat[j,1], 4)
			}
			if(M == Ms[1]){
				analytical.results[j,] = unlist(av.pivotal.event.probs.general(v.vec, s.vec, increments = increments, increments.at.end = increments.at.end, increments.for.dir.mat = increments.for.dir.mat))[types]				
			}
			if(!skip.sims){
				sim.results[j, ] = unlist(av.pivotal.probs.from.sims(v.vec = v.vec, s.vec = s.vec, M = M, tol = tol))[types.sim] 				
			}
		}
		key = paste0("M", M)
		ar.list[[key]] = analytical.results # we paste this for all Ms even though it's only calculated for the first one
		sim.list[[key]] = sim.results
		cat("Done.\n\n")
	}
	list(ar.list = ar.list, sim.list = sim.list, expected.result.mat = expected.result.mat, s.mat = s.mat)	
}


## 

plot.rmse.by.pivotal.event.and.M = function(out, types = NULL, Ms = NULL){
	if(is.null(Ms)){
		Ms = as.integer(gsub("M", "", names(out$sim.list)))
	}
	if(is.null(types)){
		types = colnames(out$sim.list[[paste0("M", Ms[1])]])
	}
	par(mfrow = c(4,3), mar = c(2,2,1,1))
	for(j in 1:length(types)){
		plot(log(range(Ms)), c(0, 1), type = "n", xlab = "Number of draws (M)", main = types[j], ylab = "RMSE (sim vs analytical)", axes = F)
		rmses = c()
		axis(1, at = log(Ms), labels = paste0(Ms/1000, "k")) # c("1K", "10k", "100k", "1M"))
		axis(2)
		for(M in Ms){
			key = paste0("M", M)
			sim.result = out$sim.list[[key]][,j]
			analytical.result = out$ar.list[[key]][,j]
			rmses = c(rmses, sqrt(mean((analytical.result - sim.result)^2)))
		}
		abline(h = 0, col = "gray")
		lines(log(Ms), rmses, pch = 19, cex = .3)
	}
}


plot.sim.results.vs.analytical.results = function(sim.results, analytical.results, types = NULL, plot.rows = 4, plot.cols = 3, log = T, lims = NULL, add.line.with.slope = NULL, text.labels = NULL, new = T, mar = c(2,2,1,1), report.cor.coef = T, frp.lims = NULL, srp.lims = NULL){
	if(is.null(types)){
		types = colnames(sim.results)
	}
	if(length(types) > plot.rows*plot.cols){cat("you won't have enough panels in your plot.\n")}
	par(mfrow = c(plot.rows, plot.cols), mar = mar)
	if(log){
		for(j in 1:ncol(analytical.results)){
			use = which(sim.results[,j] > 0 & !is.na(analytical.results[,j]))
			if(is.null(lims)){lims = c(-6, 3)}
			the.lm = lm(analytical.results[use,j] ~ sim.results[use,j])
			the.ln.lm = lm(I(log(analytical.results[use,j])) ~ I(log(sim.results[use,j])))
			if(new){
				plot(log(sim.results[,j]), log(analytical.results[,j]), type = "n", main = paste0(types[j], ifelse(report.cor.coef, paste0(": cor = ", round(cor(sim.results[,j], analytical.results[,j]), 2), ", coef = ", round(coef(the.ln.lm)[2], 2), " (sim on x)"), "")), xlim = lims, ylim = lims)
				abline(a = 0, b = 1, lty = 3)
			}
			if(!is.null(text.labels)){
				text(log(sim.results[,j]), log(analytical.results[,j]), labels = text.labels, cex = .7)
			}else{
				points(log(sim.results[,j]), log(analytical.results[,j]), pch = 19, col = rgb(.6, .6, .6, alpha = .5))
			}
		 	abline(the.ln.lm, col = "red")
			if(!is.null(add.line.with.slope)){
				abline(a = 0, b = add.line.with.slope, lty = 2, col = "red")
			}
		}			
	}else{
		for(j in 1:ncol(analytical.results)){
			use = !is.na(analytical.results[,j])
			the.lm = lm(analytical.results[use,j] ~ sim.results[use,j])
			if(!is.null(frp.lims) & j > 3){
				lims = frp.lims
			}else if(!is.null(srp.lims) & j <= 3){
				lims = srp.lims
			}else{
				lims = range(c(analytical.results[use,j], sim.results[use,j]), na.rm = T)				
			}
			if(new){
				plot(sim.results[,j], analytical.results[,j], type = "n", main = paste0(types[j], ifelse(report.cor.coef, paste0(": cor = ", round(cor(sim.results[,j], analytical.results[,j]), 2), ", coef = ", round(coef(the.lm)[2], 2), " (sim on x)"), "")), xlim = lims, ylim = lims)
				abline(a = 0, b = 1, lty = 3)
			}
			if(!is.null(text.labels)){
				text(sim.results[,j], analytical.results[,j], labels = text.labels, cex = .7)
			}else{
				points(sim.results[,j], analytical.results[,j], pch = 19, col = rgb(.6, .6, .6, alpha = .5)) 
			}
			abline(the.lm, col = "red")
			if(!is.null(add.line.with.slope)){
				abline(a = 0, b = add.line.with.slope, lty = 2, col = "red")
			}
		}		
	}
}

plot.rmse.boxplots = function(out, Ms = NULL, main = NULL){
	if(is.null(Ms)){
		Ms = as.integer(gsub("M", "", names(out$sim.list)))
	}
	all.rmses = matrix(NA, nrow = 0, ncol = 2)
	colnames(all.rmses) = c("M", "rmse")
	for(M in Ms){
		key = paste0("M", M)
		sim.mat = out$sim.list[[key]]
		analytical.mat = out$ar.list[[key]]
		rmses = c()
		for(j in 1:nrow(sim.mat)){
			rmses = c(rmses, sqrt(mean((sim.mat[j,] - analytical.mat[j,])^2)))
		}
		all.rmses = rbind(all.rmses, cbind(M, rmses))
	}
	plot(c(0.5, length(Ms) + .5), c(0, 2), type = "n", xlab = "Number of draws (M)", main = main, ylab = "Distribution of RMSE (sim vs analytical)", axes = F)
	axis(1, at = 1:length(Ms), labels = paste0(Ms/1000, "k"))
	axis(2)
	boxplot(all.rmses[,2] ~ all.rmses[,1], add = T, at = 1:length(Ms), names = rep("", length(Ms)))
}

## How about a single RMSE plot for all pivotal events

plot.rmse.by.pivotal.event.and.M.one.plot = function(out, types = NULL, Ms = NULL, main = NULL){
	if(is.null(Ms)){
		Ms = as.integer(gsub("M", "", names(out$sim.list)))
	}
	if(is.null(types)){
		types = colnames(out$sim.list[[paste0("M", Ms[1])]])
	}
	plot(log(range(Ms)), c(0, 1), type = "n", xlab = "Number of draws (M)", main = main, ylab = "RMSE (sim vs analytical)", axes = F)
	axis(1, at = log(Ms), labels = paste0(Ms/1000, "k")) # c("1K", "10k", "100k", "1M"))
	axis(2)
	abline(h = 0, col = "gray")
	for(j in 1:length(types)){
		rmses = c()
		for(M in Ms){
			key = paste0("M", M)
			sim.result = out$sim.list[[key]][,j]
			analytical.result = out$ar.list[[key]][,j]
			rmses = c(rmses, sqrt(mean((analytical.result - sim.result)^2)))
		}
		lines(log(Ms), rmses, lty = j, col = rgb(.3, .3, .3, alpha = .5))
	}
}
