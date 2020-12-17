## choose weight to put on beta for correction of numerical integration near edges of dirichlet 

# when the dirichlet gets very steep near edge (which happens when one of the alphas is very small), the numerical estimation is biased. 

# I tried getting a numerical solution, but did not make much progress. 

# the idea is that at the limit this is a beta. so, let's blend them optimally to address this problem. 

source("utils/plurality_pivotal_probabilities_analytical.r")
source("utils/plurality_pivotal_probabilities_from_sims.r")

# first illustrate the problem 

pdf("./../output/figs/20190410/diagnosing_problem_role_of_increments_plus_2x2_50_100_triangular_kernel.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
xs = seq(3, .1, -.05)
incs = c(50, 100)
for(s in c(20, 85)){
  for(big.x in c(20, 22)){
    xs.scale = (big.x + 18)/(big.x + 18 + xs)
    sto = stoa = stob = stoc = stod = stoe = stof = rep(NA, length(xs))
    for(i in 1:length(xs)){
      x = xs[i]
      this.v.vec = normalize.to.1(c(big.x, 18, x))
      sto[i] = plurality.pivotal.probs.from.sims(v.vec = this.v.vec, s = s)$AB
      stoa[i] = plurality.pivotal.probabilities(v.vec = this.v.vec, s = s, increments = incs[1], mix.with.beta = F)$AB
      stob[i] = plurality.pivotal.probabilities(v.vec = this.v.vec, s = s, increments = incs[2], mix.with.beta = F)$AB
      stoc[i] = plurality.pivotal.probabilities(v.vec = this.v.vec, s = s, increments = incs[2], mix.with.beta = T)$AB
      stod[i] = plurality.pivotal.probabilities(v.vec = this.v.vec, s = s, increments = incs[2], mix.with.beta = T, mix.start = 3)$AB
      stoe[i] = plurality.pivotal.probabilities(v.vec = this.v.vec, s = s, increments = incs[2], mix.with.beta = T, mix.start = 5)$AB
      # stof[i] = simpson_rule_tie_for_first_1_and_2(v.vec = this.v.vec, s = s, increments = incs[2])
    }
    plot(xs.scale, sto, type = "l", ylim = c(1, 4), xlab = "Sum of top two vote shares", ylab = "")
    lines(xs.scale, stoa, lty = 2, col = "red")
    lines(xs.scale, stob, lty = 2, col = "blue")
    lines(xs.scale, stoc, lty = 2, col = "orange", lwd = 2)
    lines(xs.scale, stod, lty = 2, col = "pink", lwd = 2)
    lines(xs.scale, stoe, lty = 2, col = "purple", lwd = 2)
    # lines(xs.scale, stof, lty = 1, col = rgb(.5, .5, .5, alpha = .5), lwd = 2)
    # it should converge to a beta, right? 
    v.vec.2 = c(big.x, 18)/sum(c(big.x, 18))
    abline(h = .5*dbeta(x = .5, shape1 = s*v.vec.2[1], shape2 = s*v.vec.2[2]))
  }
}
dev.off()

# let's say switch to beta when alpha_3 - 1 = 0

# interesting -- so it goes off the rails. 
# I think the blend with beta may also be best.

## that's a reasonable compromise. could test it more comprehensively. 

## another approach is to adapt the intervals in the numerical approach to deal with the curvature. 

## looking at Simpson's rule and other things for numerical integration. 

# let's take a case close to the end 
x = xs[i - 4]
this.v.vec = normalize.to.1(c(20, 18, x))
ys = seq(.33, .5, by = .001)
where.to.evaluate = cbind(ys, ys, 1-2*ys) 
densities.along.locus = ddirichlet(where.to.evaluate, this.v.vec*50)
plot(ys, densities.along.locus, type = "l")
use = !is.infinite(densities.along.locus)
summary(lm(densities.along.locus[use] ~ poly(ys[use], 8)))

# some solutions: switch to beta. use different ddirichlet function
# mixtools, gtools, MCMCpack all do the same. 
# gtools and MCMCpack code is identical.
## MCMCpack says "Code is taken from Greg's Miscellaneous Functions (gregmisc)"
## they probably all use that. 
## well the density must be correct. 
## increasing the increments does reduce the problem. so it seems to be about using midpoints when the function is very convex. 
## alternatives: 
## mix in the beta. in plurality case, this is pretty simple: looking for tie between A and B? this is a mixture of the analytical dirichlet answer and .5*dbeta(.5, alpha_A, alpha_B) where the weight on the beta is (1 - alpha_C/(alpha_A + alpha_B + alpha_C))^k, for maybe k = 2?  
## change the way we do increments -- I already messed with that in the analytical probs for the AV case.
## solve analytically 


## OK check a bit more, try to implement for AV too. 




v.vec = normalize.to.1(c(50, 40, .1))
s = 80
pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 20)
pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 50)
pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 500)
pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 5000)
pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 50000)
pr.tie.for.first.1.and.2.analytical.beta(v.vec, s)
plurality.pivotal.probs.from.sims(v.vec, s, M = 1000000)$AB

m = 100
M = 1000000
v.vecs = rdirichlet(m, alpha = c(.48, .46, .06)*20)
ss = runif(m, 40, 40) # for now, keeping stable
mat = matrix(NA, nrow = nrow(v.vecs), ncol = 3)
colnames(mat) = c("sim", "analytical", "beta")
for(i in 1:nrow(mat)){
  v.vec = v.vecs[i,]
  s = ss[i]
  mat[i,1] = plurality.pivotal.probs.from.sims(v.vec, s, M = 100000)$AB
  mat[i,2] = pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments = 100)
  mat[i,3] = pr.tie.for.first.1.and.2.analytical.beta(v.vec, s)
}
cbind(v.vecs[,3], mat)

# so we can at least specify a class of functions and try out different values. 
wob = (1 - v.vecs[,3])^10  # weight on beta 
try.combo = (1-wob)*mat[,2] + wob*mat[,3]

plot(mat[,"sim"], mat[,"analytical"], pch = 19, col = rgb(.6, .6, .8, alpha = .5))
points(mat[,"sim"], mat[,"beta"], pch = 20, col = rgb(.6, .6, .4, alpha = .5))
abline(lm(mat[,"sim"] ~ try.combo), col = "red")
abline(lm(mat[,"sim"] ~ mat[,"analytical"]), col = "blue")
abline(lm(mat[,"sim"] ~ mat[,"beta"]), col = "green")
abline(a = 0, b = 1, lty = 2)

# maybe use optim to choose among a set of parameters

