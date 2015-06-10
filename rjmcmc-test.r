
### parametres de RJMCMC

kmax <- 30        ### nombre maximum de nucleosomes par regions
lambda <- 3       ### parametre de la loi de poisson pour k. On peut essayer lambda = 1,2,3,...,
ecartmin <- 146   ### ecart minimal entre deux nucleosomes
ecartmax <- 292   ### ecart maximal entre deux nucleosomes.  
                  ### A choisir experimentalement en prenant la distance maximale entre les pics des donnees
minReads <- 3     ### nombre ninimum de brins dans une région candidate
niter <- 25000    ### nombre d'iterations. On peut aller jusqu'a 100 000


### En utilisant la segmentation de PING
### seg@List[[i]]@yF est le vecteur contenant les brins +
### seg@List[[i]]@yR est le vecteur contenant les brins -

### Appel de la fonction RJMCMC pour une région candidate i

source("rjmcmc.r")

ptm <- proc.time()
liste  <- RJMCMC(seg@List[[i]]@yF, seg@List[[i]]@yR, niter, kmax, lambda, ecartmin, ecartmax, minReads)
proc.time() - ptm

liste$mu    ### retourne les positions de nucléosomes dans la region candidate i
