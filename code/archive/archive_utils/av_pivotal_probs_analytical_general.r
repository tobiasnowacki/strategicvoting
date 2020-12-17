
get.increment.midpoints = function(start, end, increments, increments.at.end = NULL){
	inc.length = (end - start)/increments
	out = seq(start + inc.length/2, end - inc.length/2, length = increments)
	if(!is.null(increments.at.end)){
		to.add.at.end = get.increment.midpoints(mean(out[c(length(out) - 1, length(out))]), end, increments.at.end) 
		out = c(out[-length(out)], to.add.at.end)
	}
	out
}


#####################################
#### second-round pivotal events ####
#####################################


## calculate the probability that A and B are tied in the second-round of a three-candidate AV contest. 
# the "general" approach takes a length-6 vector of expected results (v.vec) and a length-4 vector of precision parameters (s.vec)
# -- FP precision and precision on each second-preference ratio.   
pr.second.round.pivotal.cands.1.and.2.general = function(v.vec, s.vec, increments = 50){
	stopifnot(length(v.vec) == 6)
	stopifnot(length(s.vec) == 4)
	fp.alpha = s.vec[1]*c(sum(v.vec[1:2]), sum(v.vec[3:4]), sum(v.vec[5:6]))
	midpoints.y = get.increment.midpoints(1/4, 1/2, increments)  # the y values we will be evaluating  
	y.inc = midpoints.y[2] - midpoints.y[1]
	sum = 0
	for(y in midpoints.y){ 
		this.sum = 0
		the.from = 1/2 - y
		the.to = ifelse(y < 1/3, y, (1-y)/2)
		how.many.x = round((the.to - the.from)/y.inc) # roughly getting squares
		if(how.many.x == 0){next} # skip if it's really small
		midpoints.x = get.increment.midpoints(the.from, the.to, how.many.x)
		for(x in midpoints.x){
			dir.part = ddirichlet(c(y, 1-y-x, x), alpha = fp.alpha)
			alpha_31 = s.vec[4]*v.vec[5]
			alpha_32 = s.vec[4]*v.vec[6]
			beta.part = dbeta((1/2-y)/x, alpha_31, alpha_32)*(1/x)
			this = dir.part*beta.part
			this.sum = this.sum + this
		}
		sum = sum + this.sum*(the.to - the.from)/how.many.x
	}
	sum*y.inc
}

# permutations of the above so we get second-round pivotal probabilities AB, AC, BC
second.round.pivotal.event.probs.general = function(v.vec, s.vec, increments = 50){
	list("AB" = pr.second.round.pivotal.cands.1.and.2.general(v.vec = v.vec, s.vec = s.vec, increments = increments), 
	# switching C and B
	"AC" = pr.second.round.pivotal.cands.1.and.2.general(v.vec = v.vec[c(2,1,5,6,3,4)], s.vec = s.vec[c(1,2,4,3)], increments = increments),
	# moving A last, moving B and C up
	"BC" = pr.second.round.pivotal.cands.1.and.2.general(v.vec = v.vec[c(4,3,6,5,1,2)], s.vec = s.vec[c(1,3,4,2)], increments = increments))
}




# a faster approach for the Dirichlet special case where there is just one precision parameter 
pr.second.round.pivotal.cands.2.and.3 = function(alpha.vec, increments = 200){
	# now we reverse the alpha.vec so that we can calculate as a tie between 1 and 2, as in the appendix calculations
	alpha.vec = alpha.vec[c(6,5,4,3,2,1)] 
	midpoints.y = get.increment.midpoints(0, 1/2, increments)  # the y values we will be evaluating  
	sum = 0
	for(y in midpoints.y){ 
		first.part = ddirichlet(c(y, 1/2 - y, 1/2), alpha = c(sum(alpha.vec[1:2]), alpha.vec[5], sum(alpha.vec[c(3:4, 6)])))
		if(y < 1/3){int.limit = 4*y -1}else{int.limit = y}
		second.part = pbeta(int.limit, alpha.vec[6], sum(alpha.vec[3:4]))
		sum = sum + first.part*second.part
	}
	(midpoints.y[2] - midpoints.y[1])*sum
}

# permutations of the Dirichlet special case 
second.round.pivotal.event.probs = function(v.vec, s.vec, increments = 100){
	stopifnot(length(unique(s.vec)) == 1) # this is only safe when s.vec is rep(x, 4)
	# first we constitute the alpha.vec from v.vec and s.vec
	alpha.vec = v.vec*s.vec[1]
	list("AB" = pr.second.round.pivotal.cands.2.and.3(alpha.vec[c(5,6,2,1,4,3)], increments = increments), 
	"AC" = pr.second.round.pivotal.cands.2.and.3(alpha.vec[c(3,4,1,2,6,5)], increments = increments),
	"BC" = pr.second.round.pivotal.cands.2.and.3(alpha.vec, increments = increments))
}


## the truncated case where we have ballots that just list 1 candidate. 
## AB, AC, BA, BC, CA, CB, A, B, C  

pr.second.round.pivotal.cands.1.and.2.truncated = function(v.vec, s.vec, increments = 50){
	# this function allows truncated ballots. slower because another dimension of integration. 
	stopifnot(length(unique(s.vec)) == 1)
	stopifnot(length(v.vec) == 9)
	alpha.vec = s.vec[1]*v.vec
	midpoints.y = get.increment.midpoints(1/4, 1/2, increments)  # the y values we will be evaluating
	y.inc = midpoints.y[2] - midpoints.y[1]
	value.at.y = 0
	values.at.y = c()
	for(y in midpoints.y){
		to.x = 4*y - 1 
		if(y < 1/3){to.x = 1 - 2*y}
		how.many.x = round(to.x/y.inc)
		if(how.many.x == 0){next}
		midpoints.x = get.increment.midpoints(0, to.x, how.many.x)
		for(x in midpoints.x){
			dir.part = ddirichlet(c(y, (1-x)/2 - y, x, (1-x)/2), alpha = c(sum(alpha.vec[c(1,2,7)]), alpha.vec[5], alpha.vec[9], sum(alpha.vec[c(3,4,8,6)])))
			beta.limit = (4*y-1-x)/(1-x)
			if(y > 1/3){beta.limit = (y - x)/(1 - x)}
			beta.part = pbeta(beta.limit, alpha.vec[6], sum(alpha.vec[c(3,4,8)]))
			value.at.y = value.at.y + y.inc*dir.part*beta.part
		}
	values.at.y = c(values.at.y, value.at.y)
	}
	y.inc*sum(values.at.y)
}


pr.second.round.pivotal.cands.1.and.2.truncated.2 = function(v.vec, s.vec, increments = 50, noisy = F){
	stopifnot(length(v.vec) == 9)
	stopifnot(length(unique(s.vec)) == 1)
	fp.alpha = s.vec[1]*c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
	alpha.vec = s.vec[1]*v.vec
	midpoints.y = get.increment.midpoints(1/4, 1/2, increments)  # the y values we will be evaluating  
	y.inc = midpoints.y[2] - midpoints.y[1]
# 	if(noisy){cat("y.inc: ", y.inc, " ")}
	sum = 0
# 	if(noisy){cat(length(midpoints.y))}
	dir.parts = dir.sums = c()
	for(y in midpoints.y){ 
		this.sum = 0
		the.from = 1/2 - y
		the.to = ifelse(y < 1/3, y, (1-y)/2)
		how.many.x = round((the.to - the.from)/y.inc) # roughly getting squares
		if(how.many.x == 0){next} # skip if it's really small
		midpoints.x = get.increment.midpoints(the.from, the.to, how.many.x)
		# if(noisy & y == midpoints.y[5]){cat(" x1 ", length(midpoints.x))}
		# x.inc = midpoints.x[2] - midpoints.x[1]  # should be roughly y.inc
		# if(noisy){cat("x.inc: ", x.inc, " ")}
		for(x in midpoints.x){
			dir.parts = c(dir.parts, y.inc*y.inc*ddirichlet(c(y, 1-y-x, x), alpha = fp.alpha)) # it's a tile. the area of the tile is y.inc^2
			# now to get the inner integral
			the.from = 0
			the.to = ifelse(x < 1 - 2*y, 2*x + 2*y - 1, 1-2*y) # to keep things within the right area -- shouldn't change anything but will keep us out of trouble and be faster
			how.many.z = round((the.to - the.from)/y.inc)
			if(how.many.z == 0){next}  # skippable
			midpoints.z = get.increment.midpoints(the.from, the.to, how.many.z)
			# if(noisy & x == midpoints.x[1]){cat(" x2 ", length(midpoints.z))}
			# z.inc = midpoints.z[2] - midpoints.z[1]  # should be roughly y.inc
			# if(noisy){cat("z.inc: ", z.inc, " ")}
			dir.sum.parts = c()
			for(z in midpoints.z){
				dir.sum.parts = c(dir.sum.parts, y.inc*ddirichlet(c(((1 - z)/2 - y)/x, (x - (1 + z)/2 + y)/x, z/x), alpha = alpha.vec[c(5,6,9)]))
			}
			dir.sums = c(dir.sums, sum(dir.sum.parts, na.rm = T))
		}
	}
	list(sum = sum(dir.parts*dir.sums, na.rm = T), 
		dir.parts = dir.parts, dir.sums = dir.sums)
}


# this is a more elegant way of doing it, but it's giving me the same answers I think. to figure out tomorrow. 


pr.second.round.pivotal.cands.1.and.2.truncated.3 = function(v.vec, s.vec, increments.grid = 50, increments.at.end = 20, noisy = F){
	stopifnot(length(v.vec) == 9)
	stopifnot(length(unique(s.vec)) == 1)
	alpha.vec = s.vec[1]*v.vec
	midpoints.xy = get.increment.midpoints(1/4, 1/2, increments.grid)  # the x values we will be evaluating in the FP plot
	xy.inc = midpoints.xy[2] - midpoints.xy[1]
	fp.grid = sp.grid = matrix(0, length(midpoints.xy), length(midpoints.xy))
	# now we cycle through the grid
	for(i in 1:length(midpoints.xy)){
		x = midpoints.xy[i]
		for(j in 1:length(midpoints.xy)){
			y = midpoints.xy[j]
			if(y < .5 - x/2 | x < .5 - y/2){next} # skipping FP vecs where C is not last -- easier than specifying limits of integration; only slightly less efficient I believe. 
			fp.grid[i,j] = xy.inc^2*ddirichlet(c(x, y, 1-x-y), alpha = c(sum(alpha.vec[c(1,2,7)]), sum(alpha.vec[c(3,4,8)]), sum(alpha.vec[c(5,6,9)])))
			sp.grid[i,j] = prob.being.along.line.analytical(alpha.vec = alpha.vec[c(5,6,9)], x = x, y = y, pica.increment.length = xy.inc, increments.at.end = increments.at.end) 
		}
	}
	list(sum = sum(fp.grid*sp.grid), fp.grid = fp.grid, sp.grid = sp.grid)
}

# analytical version 
prob.being.along.line.analytical = function(alpha.vec, x, y, pica.increment.length = .005, increments.at.end = 100, picx.cutoff = 1){
	if(alpha.vec[3] < picx.cutoff){
		return((1/(2*(1-x-y)))*dbeta((1-2*y)/(2*(1-x-y)), alpha.vec[2], alpha.vec[1]))  # not sure why there is a 2 in the denom -- TODO: worry about this. 
	}
	# first locate the line, so we can integrate along it
	to.pica = 1/2 - x # \pi_{CA} at end point
	from.pica = ifelse(x > y, 0, y-x)
	how.many.pica = round((to.pica - from.pica)/pica.increment.length)
	if(how.many.pica <= 1){return(0)}
	midpoints.pica = get.increment.midpoints(from.pica, to.pica, how.many.pica, increments.at.end = increments.at.end) # oversampling at end
	basic.inc = midpoints.pica[2] - midpoints.pica[1]
	increments.pica = c(rep(basic.inc, how.many.pica - 1), rep(basic.inc/increments.at.end, increments.at.end))
	stopifnot(length(midpoints.pica) == length(increments.pica))
	pica.dirs = c()
	for(i in 1:length(midpoints.pica)){
		pica = midpoints.pica[i]; this.inc = increments.pica[i]
		this.vec = c(pica, x - y + pica, 1 - 2*(x + pica)) # the raw vote share vec we are considering, e.g. .1, .05, .01 
		this.vec = this.vec/sum(this.vec) # normalize (divide by 1 - x - y)
		pica.dirs = c(pica.dirs, this.inc*ddirichlet(this.vec, alpha.vec)) # calculate and store the dirichlet density at this point.  
	}
	(1/((1-x-y)^2))*sum(pica.dirs, na.rm = T) # the normalization: sqrt(2) was because our increments are horizontal, but the lengths are along the 45 degree line; pica.inc because this is the size of the increments; divide by (1-x-y)^2 because the area of the full dirichlet is 1/2 but the area of the smaller one is (1 - x - y)^2/2. 
}

plot.pica = function(x, y, alpha.vec = NULL, plot.draws = NULL, pica.increment.length = .005){
	plot(c(0, 1), c(0,1), type = "n", xlim = c(0, 1-x-y), ylim = c(0,1-x-y), axes = F, xlab = "", ylab = "")
	axis(1); axis(2)
	mat = get.pica.and.picb.vectors.given.x.and.y(x,y,pica.increment.length)
	lines(mat[,1], mat[,2])
	abline(a = 1-x-y, b = -1, lty = 2)
	if(!is.null(alpha.vec)){
		exp.vec = (1 - x - y)*alpha.vec/sum(alpha.vec)
		cat("plotting at x = ", round(exp.vec[1], 2), " and y = ", round(exp.vec[2], 2), "..\n")
		points(exp.vec[1], exp.vec[2], pch = 19, col = "red")
	}
	if(!is.null(plot.draws)){
		draws = rdirichlet(n = plot.draws, alpha = alpha.vec)*(1 - x - y)
		points(draws[,1], draws[,2], pch = 19, cex = .25, col = rgb(.5, .5, .5, alpha = .5))
	}
}




# just some checking -- this looks fine. 
get.pica.and.picb.vectors.given.x.and.y = function(x,y,inc){
	to.pica = (1 - 2*x)/2
	from.pica = ifelse(x > y, 0, y-x)
	how.many.pica = round((to.pica - from.pica)/inc)
	if(how.many.pica == 0){return(NA)}
	picas = get.increment.midpoints(from.pica, to.pica, how.many.pica)
	picbs = x - y + picas
	cbind(picas, picbs)
}

get.segment.length.given.x.and.y = function(x,y){
	to.pica = (1 - 2*x)/2
	from.pica = ifelse(x > y, 0, y-x)
	sqrt(2)*(to.pica - from.pica)
}




# permutations of the Dirichlet special case, when we have truncated ballots 
second.round.pivotal.event.probs.truncated = function(v.vec, s.vec, increments.grid = 40, increments.at.end = 10){
	stopifnot(length(unique(s.vec)) == 1) # this is only safe when s.vec is rep(x, 4)
	ab.thing = pr.second.round.pivotal.cands.1.and.2.truncated.3(v.vec, s.vec, increments.grid = increments.grid, increments.at.end = increments.at.end)
	ac.thing = pr.second.round.pivotal.cands.1.and.2.truncated.3(v.vec[c(2,1,5,6,3,4,7,9,8)], s.vec, increments.grid = increments.grid, increments.at.end = increments.at.end)
	bc.thing = pr.second.round.pivotal.cands.1.and.2.truncated.3(v.vec[c(4,3,6,5,1,2,8,9,7)], s.vec, increments.grid = increments.grid, increments.at.end = increments.at.end)
	list("AB" = ab.thing$sum, 
	"AC" = ac.thing$sum,
	"BC" = bc.thing$sum,
	"A.last" = sum(bc.thing$fp.grid, na.rm = T),
	"B.last" = sum(ac.thing$fp.grid, na.rm = T),
	"C.last" = sum(ab.thing$fp.grid, na.rm = T))
}


# permutations of the Dirichlet special case, when we have truncated ballots 
second.round.pivotal.event.probs.truncated.2 = function(v.vec, s.vec, increments = 100){
	stopifnot(length(unique(s.vec)) == 1) # this is only safe when s.vec is rep(x, 4)
	ab.thing = pr.second.round.pivotal.cands.1.and.2.truncated.2(v.vec, s.vec, increments = increments)
	ac.thing = pr.second.round.pivotal.cands.1.and.2.truncated.2(v.vec[c(2,1,5,6,3,4,7,9,8)], s.vec, increments = increments)
	bc.thing = pr.second.round.pivotal.cands.1.and.2.truncated.2(v.vec[c(4,3,6,5,1,2,8,9,7)], s.vec, increments = increments)
	list("AB" = ab.thing$sum, 
	"AC" = ac.thing$sum,
	"BC" = bc.thing$sum,
	"A.last" = sum(bc.thing$dir.parts, na.rm = T),
	"B.last" = sum(ac.thing$dir.parts, na.rm = T),
	"C.last" = sum(ab.thing$dir.parts, na.rm = T))
}




####################################
#### first-round pivotal events ####
####################################

# these functions are used below in first.round.pivotal.events.cands.1.and.2.truncated()
# we want to compute the dirichlet once over the whole possible area, and then add up the sums as appropriate. 
dirichlet.matrix.to.beat.3 = function(alpha.segment, increments = 20){
	paxy = get.increment.midpoints(0, 1, increments) # \pi_{AX}/y, thus paxy
	pacy = get.increment.midpoints(0, 1, increments) # \pi_{AC}/y, thus pacy
	dir.mat = matrix(0, ncol = length(paxy), nrow = length(pacy)) # this matrix covers the simplex of \pi_{AX}/y, \pi_{AC}/y.  
	colnames(dir.mat) = round(paxy, 4); rownames(dir.mat) = round(pacy, 4)
	# now calculate the dirichlet at each point in the grid 
	for(i in 1:length(pacy)){
		for(j in 1:length(paxy)){
			eval.at = c(1 - paxy[j] - pacy[i], pacy[i], paxy[j])
			if(min(eval.at) <= 0 | max(eval.at) >= 1){next} # don't need to calculate when something is below or above 1.
			dir.mat[i,j] = ddirichlet(eval.at, alpha = alpha.segment)
		}
	}
	# normalize -- this should add up to 1. previously multiplied by the tile size
	dir.mat/sum(dir.mat) # *(paxy[2] - paxy[1])*(pacy[2] - pacy[1])
}

use.mat.to.beat.3 = function(y, increments = 20){
	# a matrix of integers: "For this y, should this combination of parameters be included in the Dirichlet sum?"
	paxy = get.increment.midpoints(0, 1, increments) # \pi_{AX}/y, thus paxy
	pacy = get.increment.midpoints(0, 1, increments) # \pi_{AC}/y, thus pacy
	use.mat = matrix(NA, ncol = 0, nrow = length(pacy))
	for(j in 1:length(paxy)){
		use.mat = cbind(use.mat, as.integer(pacy <= 2 - 1/(2*y) - paxy[j]/2))
	}
	use.mat	
}


first.round.pivotal.events.cands.1.and.2.truncated = function(v.vec, s.vec, increments = 50, increments.2 = 20){
	# this function allows truncated ballots. slower because another dimension of integration. 
	stopifnot(length(unique(s.vec)) == 1)
	stopifnot(length(v.vec) == 9)
	alpha.vec = s.vec[1]*v.vec
	midpoints.y = get.increment.midpoints(1/4, 1/3, increments)  # the y values we will be evaluating
	y.inc = midpoints.y[2] - midpoints.y[1]
	# for each y, we will calculate the probability of 1 beating 3 and 2 beating 3; we will then sum up with weights based on the probability of a tie for second between 1 and 2 at each y
	dir.mat.13 = dirichlet.matrix.to.beat.3(alpha.vec[c(3,4,8)], increments.2) # this matrix shows the probability of getting entries in a simplex with BX along the horizontal and BC on the vertical. 
	dir.mat.23 = dirichlet.matrix.to.beat.3(alpha.vec[c(1,2,7)], increments.2)
	pr.tie.for.second.at.y = pr.1.beats.3.given.y = pr.2.beats.3.given.y = c()
	for(y in midpoints.y){
		pr.tie.for.second.at.y = c(pr.tie.for.second.at.y, ddirichlet(c(y, y, 1-2*y), alpha = c(sum(alpha.vec[c(1,2,7)]), sum(alpha.vec[c(3,4,8)]), sum(alpha.vec[c(5,6,9)]))))
		this.use.mat = use.mat.to.beat.3(y, increments.2)
		pr.1.beats.3.given.y = c(pr.1.beats.3.given.y, sum(this.use.mat*dir.mat.13))
		pr.2.beats.3.given.y = c(pr.2.beats.3.given.y, sum(this.use.mat*dir.mat.23))
	} 
	list(
		"12.12" = y.inc*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*pr.2.beats.3.given.y),
		"12.32" = y.inc*sum(pr.tie.for.second.at.y*(1 - pr.1.beats.3.given.y)*pr.2.beats.3.given.y),
		"12.13" = y.inc*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*(1 - pr.2.beats.3.given.y)) #,
#		"12" = y.inc*sum(pr.tie.for.second.at.y),
#		"tie.vec" = pr.tie.for.second.at.y,
#		"1.beats.3.vec" = pr.1.beats.3.given.y,
#		"2.beats.3.vec" = pr.2.beats.3.given.y,
#		"y.vec" = midpoints.y
		)
}


first.round.pivotal.event.probs.truncated = function(v.vec, s.vec, increments = 100, increments.2 = 20){
	frpe.12 = first.round.pivotal.events.cands.1.and.2.truncated(v.vec, s.vec, increments = increments, increments.2 = increments.2)
	frpe.13 = first.round.pivotal.events.cands.1.and.2.truncated(v.vec = v.vec[c(2,1,5,6,3,4,7,9,8)], s.vec = s.vec, increments = increments, increments.2 = increments.2)
	frpe.23 = first.round.pivotal.events.cands.1.and.2.truncated(v.vec = v.vec[c(4,3,6,5,1,2,8,9,7)], s.vec = s.vec, increments = increments, increments.2 = increments.2)
	list("AB.AB" = frpe.12[["12.12"]], "AB.AC" = frpe.12[["12.13"]], "AB.CB" = frpe.12[["12.32"]], 
	"AC.AC" = frpe.13[["12.12"]], "AC.AB" = frpe.13[["12.13"]], "AC.BC" = frpe.13[["12.32"]], 
	"BC.BC" = frpe.23[["12.12"]], "BC.BA" = frpe.23[["12.13"]], "BC.AC" = frpe.23[["12.32"]] #,
	#"AB.1" =  frpe.12[["12"]], "AC.1" =  frpe.13[["12"]], "BC.1" =  frpe.23[["12"]]
	)
}


first.round.pivotal.events.cands.1.and.2.general = function(v.vec, s.vec, increments = 100){
	stopifnot(length(v.vec) == 6)
	stopifnot(length(s.vec) == 4)
	fp.alpha = s.vec[1]*c(sum(v.vec[1:2]), sum(v.vec[3:4]), sum(v.vec[5:6]))
	midpoints.y = get.increment.midpoints(1/4, 1/3, increments)  # the y values we will be evaluating
	pr.tie.for.second.at.y = pr.1.beats.3.given.y = pr.2.beats.3.given.y = c()
	for(y in midpoints.y){
		pr.tie.for.second.at.y = c(pr.tie.for.second.at.y, ddirichlet(c(y, y, 1-2*y), alpha = fp.alpha))
		pr.1.beats.3.given.y = c(pr.1.beats.3.given.y, pbeta(2 - 1/(2*y), s.vec[3]*v.vec[4], s.vec[3]*v.vec[3])) 
		pr.2.beats.3.given.y = c(pr.2.beats.3.given.y, pbeta(2 - 1/(2*y), s.vec[2]*v.vec[2], s.vec[2]*v.vec[1])) 
	} 
	inc.size = (midpoints.y[2] - midpoints.y[1])
	list(
		"12.12" = inc.size*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*pr.2.beats.3.given.y),
		"12.32" = inc.size*sum(pr.tie.for.second.at.y*(1 - pr.1.beats.3.given.y)*pr.2.beats.3.given.y),
		"12.13" = inc.size*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*(1 - pr.2.beats.3.given.y))
		)
}



first.round.pivotal.event.probs.general = function(v.vec, s.vec, increments = 100){
	frpe.12 = first.round.pivotal.events.cands.1.and.2.general(v.vec, s.vec, increments = increments)
	frpe.13 = first.round.pivotal.events.cands.1.and.2.general(v.vec = v.vec[c(2,1,5,6,3,4)], s.vec = s.vec[c(1,2,4,3)], increments = increments)
	frpe.23 = first.round.pivotal.events.cands.1.and.2.general(v.vec = v.vec[c(4,3,6,5,1,2)], s.vec = s.vec[c(1,3,4,2)], increments = increments)
	list("AB.AB" = frpe.12[["12.12"]], "AB.AC" = frpe.12[["12.13"]], "AB.CB" = frpe.12[["12.32"]], 
	"AC.AC" = frpe.13[["12.12"]], "AC.AB" = frpe.13[["12.13"]], "AC.BC" = frpe.13[["12.32"]], 
	"BC.BC" = frpe.23[["12.12"]], "BC.BA" = frpe.23[["12.13"]], "BC.AC" = frpe.23[["12.32"]] 
	)
}



####################################
######### combining them ###########
####################################

av.pivotal.event.probs.general = function(v.vec, s.vec, increments = 100){
	if(length(unique(s.vec)) == 1){
		c(second.round.pivotal.event.probs(v.vec, s.vec, increments), first.round.pivotal.event.probs.general(v.vec, s.vec, increments))
	}else{
		c(second.round.pivotal.event.probs.general(v.vec, s.vec, increments), first.round.pivotal.event.probs.general(v.vec, s.vec, increments))		
	}
}

av.pivotal.event.probs.truncated = function(v.vec, s.vec, increments = 100, increments.2 = 20, increments.grid = 20, increments.at.end = 10){
	c(second.round.pivotal.event.probs.truncated(v.vec, s.vec, increments.grid, increments.at.end), first.round.pivotal.event.probs.truncated(v.vec, s.vec, increments, increments.2))
}







