#### plurality pivotal probabilities 

# three-party version, paralleling the approach for AV 

## utility copied from AV...v2 used below 

library(gtools)

get.increment.midpoints = function(start, end, increments, increments.at.end = NULL){
	inc.length = (end - start)/increments
	out = seq(start + inc.length/2, end - inc.length/2, length = increments)
	if(!is.null(increments.at.end)){
		to.add.at.end = get.increment.midpoints(mean(out[c(length(out) - 1, length(out))]), end, increments.at.end) 
		out = c(out[-length(out)], to.add.at.end)
	}
	out
}

pr.tie.for.first.1.and.2 = function(v.vec, s, increments = 50){
	ys = get.increment.midpoints(1/3, 1/2, increments)
	incr.length = ys[2] - ys[1]
	incr.length*sum(apply(cbind(ys, ys, 1-2*ys), 1, ddirichlet, s*v.vec))
}

plurality.pivotal.probabilities = function(v.vec, s, increments = 50){
	list("AB" = pr.tie.for.first.1.and.2(v.vec, s, increments), 
	"AC" = pr.tie.for.first.1.and.2(v.vec[c(1,3,2)], s, increments), 
	"BC" = pr.tie.for.first.1.and.2(v.vec[c(2,3,1)], s, increments) 
	)	
}