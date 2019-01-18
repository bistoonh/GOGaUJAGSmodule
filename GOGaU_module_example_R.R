

library(rjags)
# load module
load.module("GOGaU")

dGOGaU = function(x,alpha, beta, a = 0, b = 1, log = FALSE)
	{
	par = c(alpha, beta, a, b)
	G = punif(x,par[3],par[4])
	g = dunif(x,par[3],par[4])
	Gb = G^par[2]
	d = par[2]*g*G^(par[1]*par[2]-1)*exp(-Gb/(1-Gb))/(gamma(par[1])*(1-Gb)^(par[1]+1))
	d[!is.finite(d)] = NA
	if(log == TRUE) d = log(d)
	return(d)
	}

rGOGaU = function(n, alpha, beta, a = 0, b = 1)
	{
	par = c(alpha, beta, a, b)
	GI=rgamma(n,par[1],1)
	q = qunif((GI/(1+GI))^(1/par[2]),par[3],par[4])
	return(q)
	}
	
model.string =  
		"
		model
			{
			alpha ~ dgamma(0.001, 0.001)
			beta ~ dgamma(0.001, 0.001)
			a ~ dunif(-40,minx)
			b ~ dunif(maxx,50)

			for (i in 1:N)
				{
				x[i] ~ dGOGaU(alpha, beta, a, b)
				}
			}
		"

	# simulation data
coli = c('#5B9BD5', '#ED7D31', '#9B9B9B', '#70AD47', '#FFC000', '#A532A5', '#4076C4', '#255E91')
N = 300
set.seed(100)
par = c(5, 0.28, 0, 10)
x <- rGOGaU(N, par[1], par[2], par[3], par[4])
dat <- list(x=x, N=N, minx = min(x), maxx = max(x))
# inits
inits1 <- list(alpha = par[1] , beta = par[2], a=min(x)-1, b = max(x)+1)
inits <- list(inits1)
# sample
model.spec <- textConnection(model.string)
j.model = NA
j.model <- jags.model(model.spec, data = dat, inits=inits, n.chains=1, n.adapt=1000)
j.samples <- coda.samples(j.model, c("alpha", "beta", "a", "b"), n.iter=6000, thin=1)


# plot
s = summary(j.samples)
jagsparh = s$statistics[c("alpha", "beta", "a", "b"), 1]
d1 = dGOGaU(sort(x), par[1], par[2], par[3], par[4])
d2 = dGOGaU(sort(x), jagsparh[1], jagsparh[2], jagsparh[3], jagsparh[4])
hist(x, prob = T, ylim = c(0, max(d1,d2)))
	lines(sort(x), d1, lwd = 2, col = '#5B9BD5')
	lines(sort(x), d2, lwd = 2, col = '#ED7D31')
	legend(0.5, max(d1,d2),
			c('Real density', 'Bayes estimator'),
			col = c('#5B9BD5', '#ED7D31'),
			lwd = c(2, 2)
			) 


