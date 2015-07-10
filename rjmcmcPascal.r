
#library("PING")
library("MCMCpack")
library("pscl")
#library("GenomicRanges")
#library(parallel)


###############################################################################
###### Fonctions auxiliaires  #################################################
#######################################################################


### Probabilit? des mouvements Death et Birth

Dk = function(k,lambda,kmax) {ifelse((k==1||k>kmax),0,0.5*min(1,dpois(k-1,lambda)/dpois(k,lambda)))}
Bk = function(k,lambda,kmax) {ifelse((k>=kmax),0,0.5*min(1,dpois(k+1,lambda)/dpois(k,lambda)))}

### Fonction qui g?n?re d'une loi normale tronqu?e sur [a,b]

tnormale = function(mu,sigma,a,b)
{
    repeat {
        y <- rnorm(1,mu,sd=sqrt(sigma))
        if (y > a & y < b) break() }
    return(y)
}


student.mixture = function(i,k,w,mu,sigma,dfr)
{
    v <- c(0,w)
    u <- runif(1,0,1)
    for (j in 1:k) {
        if(sum(v[1:j])<u & u<=sum(v[1:(j+1)]))     {
            mixte <- mu[j]+sqrt(sigma[j])*rt(1,dfr[j]) }}
    return(mixte)
}


normal.mixture = function(i,k,w,mu,sigma)
{
    v <- c(0,w)
    u <- runif(1,0,1)
    for (j in 1:k) {
        if (sum(v[1:j])<u & u<=sum(v[1:(j+1)]))
        {
            mixte <- rnorm(1,mu[j],sd=sqrt(sigma[j])) }}
    return(mixte)
}


priormu = function(mu,y)
{
    k <- length(mu)
    V <- matrix(nrow=k,ncol=k)
    for (i in 1:k) {
        for (j in 1:k) {
            if (j == i) {V[i,j] <- 1}
            else if (j == i-1) {V[i,j] <- -1}
            else {V[i,j] <- 0 }}}
    omega <- t(V)%*%V
    R <- max(y) - min(y)
    E <- (max(y) + min(y))/2
    tau <- 1/R^2
    M <- rep(E,k)
    const <- (pi/(2*tau))^{-k/2}
    prior <- const * exp(-(tau/2) * (t(mu - M) %*% omega %*% (mu - M)) )
    return(prior)
}


### Fonction qui calcule le mode d'une s?rie de donn?es


mode = function(sample) {
    tabsample <- tabulate(sample)
    samplemode <- which(tabsample== max(tabsample))[1]
    if(sum(tabsample == max(tabsample))>1) samplemode <- NA
    samplemode
}


#################################################################
### Fonction qui sert ? fusionner deux r?gions nucl?osomales ####
#################################################################

merge = function(yf, yr, y, liste, ecartmin, ecartmax, minReads)
{
    k <- length(liste$mu)
    if (k>1) {
#         ecart.min <- min(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]}))

        ## Find the smallest distance between 2 nucleosomes
        ecart     <- diff(liste$mu)
        ecart.min <- min(ecart)

        if (ecart.min < ecartmin)
        {
            repeat
            {
#                 p <- which(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]})==ecart.min)[1]

                ## Find the first position with minimum gap
                p <- which.min(ecart)

                classes <- y[y>=liste$mu[p] & y<liste$mu[p+1]]
                classesf <- yf[yf>=liste$mu[p] & yf<liste$mu[p+1]]
                classesr <- yr[yr>=liste$mu[p] & yr<liste$mu[p+1]]

                if (length(classes) > minReads){
                    mu <- mean(round(classes))
                } else {
                    mu <- mean(c(liste$mu[p], liste$mu[p+1]))
                }

                liste$mu        <- sort(liste$mu[-p])
                liste$mu[p]     <- mu
                liste$mu        <- sort(liste$mu)
                liste$sigmaf    <- liste$sigmaf[-p]
                liste$sigmar    <- liste$sigmar[-p]
                liste$delta     <- liste$delta[-p]
                liste$dl        <- liste$dl[-p]
                liste$w         <- liste$w[-p]/sum(liste$w[-p])

                k <- k - 1

                # Update the smallest distance between 2 nucleosomes
                if (k > 1) {
#                     ecart.min <- min(sapply(1:(k-1),function(i){liste$mu[i+1]-liste$mu[i]}))
                    ecart <- diff(liste$mu)
                    ecart.min <- min(ecart)
                }
                if ( k== 1 || ecart.min > ecartmin) break()
            } ### end of boucle repeat
            liste <- list(k=k, mu=liste$mu, sigmaf=liste$sigmaf, sigmar=liste$sigmar, delta=liste$delta, dl=liste$dl, w=liste$w)
        } ### end of condition if (ecart.min < ecartmin)
#         else {liste=liste}
        liste <- split(yf, yr, y, liste, ecartmin, ecartmax, minReads)
    } ### condition else if (k > 1)
    return(liste)
}



#################################################################
### Fonction qui sert ? subdiviser une r?gion nucl?osomale ######
#################################################################

split = function(yf, yr, y, liste, ecartmin, ecartmax, minReads)
{
    k <- length(liste$mu)
    if (k>1) {
#         ecart.max <- max(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]}))

        ## Find the largest distance between 2 nucleosomes
        # Its position and its value
        ecart       <- diff(liste$mu)
        p           <- order(ecart, decreasing = T)[1]
        ecart.max   <- ecart[p]

        if (ecart.max > ecartmax) {
            j <- 1
            repeat {
#                 p <- which(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]})==ecart.max)

                classes <- y[y>=liste$mu[p] & y<liste$mu[p+1]]
                classesf <- yf[yf>=liste$mu[p] & yf<liste$mu[p+1]]
                classesr <- yr[yr>=liste$mu[p] & yr<liste$mu[p+1]]

#                 j <- 1

                if (length(classes) > minReads)
                {
                    new.mu <- sort(c(liste$mu[1:k],mean(round(classes))))
                    new.sigmaf <- c(liste$sigmaf[1:k], (liste$sigmaf[p]+liste$sigmaf[p+1])/2)
                    new.sigmaf[p+1] <- (liste$sigmaf[p]+liste$sigmaf[p+1])/2
                    new.sigmaf[k+1] <- liste$sigmaf[k]
                    new.sigmar <- c(liste$sigmar[1:k], (liste$sigmar[p]+liste$sigmar[p+1])/2)
                    new.sigmar[p+1] <- (liste$sigmar[p]+liste$sigmar[p+1])/2
                    new.sigmar[k+1] <- liste$sigmar[k]
                    new.delta <- c(liste$delta[1:k], (liste$delta[p]+liste$delta[p+1])/2)
                    new.delta[p+1] <- (liste$delta[p]+liste$delta[p+1])/2
                    new.delta[k+1] <- liste$delta[k]
                    new.dl <- round(c(liste$dl[1:k], (liste$dl[p]+liste$dl[p+1])/2))
                    new.dl[p+1] <- (liste$dl[p]+liste$dl[p+1])/2
                    new.dl[k+1] <- liste$dl[k]
                    new.w <- c(liste$w[1:k], (liste$w[p]+liste$w[p+1])/2)
                    new.w[p+1] <- (liste$w[p]+liste$w[p+1])/2
                    new.w[k+1] <- liste$w[k]
                    k <- length(new.mu)
                    liste <- list(k=k, mu=new.mu, sigmaf=new.sigmaf, sigmar=new.sigmar, delta=new.delta, dl=new.dl, w=new.w/sum(new.w))

#                     ecart.max <- max(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]}))

                    ## Update the largest distance between 2 nucleosomes
                    j           <- 1
                    ecart       <- diff(liste$mu)
                    p           <- order(ecart, decreasing = T)[j]
                    ecart.max   <- ecart[p]
                }
                else
                {
#                     liste <- liste
#                     ecart.max <- sort(sapply(1:(k-1),function(j){liste$mu[j+1]-liste$mu[j]}))[k-1-j]

                    ## Select the next maximum distance between 2 nucleosomes
                    j           <- j + 1
                    p           <- order(ecart, decreasing = T)[j]
                    ecart.max   <- ecart[p]
                }
                if ( j >= (k-1) || ecart.max <= ecartmax) break()
            } ### end of boucle repeat
        } ### end of condition if (ecart.max > ecartmax)
    } ### end of condition if (k>1)
    return(liste)
}


###############################################################################
###### Fonction principale  #########
#######################################
############################################


RJMCMC = function(yf, yr, niter, kmax, lambda, ecartmin, ecartmax, minReads)
{

y <- sort(c(yf,yr))
n <- length(y)
size <- n
nf <- length(yf)
nr <- length(yr)
d <- sapply(1:size,function(m){ifelse(min(abs(yr-y[m]))==0,-1,1)})
zeta <- 147
deltamin <- 142
deltamax <- 152

  ##############################################################
  #### Initialisation des param?tres############################
  ##############################################################

  k <- rep(0,niter)
  ktilde <- rep(0,niter)
  mu <- matrix(0,nrow=niter,ncol=kmax)
  mutilde <- matrix(0,nrow=niter,ncol=kmax)
  sigmaftilde <- matrix(0,nrow=niter,ncol=kmax)
  sigmaf <- matrix(0,nrow=niter,ncol=kmax)
  sigmartilde <- matrix(0,nrow=niter,ncol=kmax)
  sigmar <- matrix(0,nrow=niter,ncol=kmax)
  deltatilde <- matrix(0,nrow=niter,ncol=kmax)
  delta <- matrix(0,nrow=niter,ncol=kmax)
  wtilde <- matrix(0,nrow=niter,ncol=kmax)
  w <- matrix(0,nrow=niter,ncol=kmax)
  a <- matrix(0,nrow=niter,ncol=kmax+1)
  atilde <- matrix(0,nrow=niter,ncol=kmax+1)
  dimtilde <- matrix(0,nrow=niter,ncol=kmax)
  dim <- matrix(0,nrow=niter,ncol=kmax)
  dl <- matrix(0,nrow=niter,ncol=kmax)
  dltilde <- matrix(3,nrow=niter,ncol=kmax)

  k[1] <- 1

  mu[1,1] <- runif(1,min(y),max(y))
  sigmaf[1,1] <- 1
  sigmar[1,1] <- 1
  delta[1,1] <- runif(1,0, 2*(mu[1,1]-min(y)))
  w[1,1] <- 1
  dl[1,1] <- 3

  a[1,1] <- min(y)
  a[1,k[1]+1] <- max(y)

  dim[1,1] <- length(y[a[1,1]<=y & y<=max(y)])

  rhob <- rep(0,niter)
  rhod <- rep(0,niter)
  rhomh <- rep(0,niter)
  Kn1 <- rep(0, niter)
  Kn2 <- rep(0, niter)
  Kn <-  rep(0, niter)
  Ln1 <- rep(0, niter)
  Ln2 <- rep(0, niter)
  Ln <-  rep(0, niter)

  Kaf <- matrix(0,nrow=nf,ncol=kmax)
  Kbf <- matrix(0,nrow=nf,ncol=kmax)
  Kar <- matrix(0,nrow=nr,ncol=kmax)
  Kbr <- matrix(0,nrow=nr,ncol=kmax)

  Y1f <- rep(0,nf)
  Y2f <- rep(0,nf)
  Y1r <- rep(0,nr)
  Y2r <- rep(0,nr)

 niter <-  ifelse((nf+nr)<=10, 1000, niter)

  for (i in 2:niter) {

    if (k[i-1]==1) {

      u<-runif(1)

         if (u<=0.5) {

          ktilde[i] <- k[i-1] + 1
          compteur <- 1
          repeat {
            j <- sample(1:k[i-1],1)
            mutilde[i,j] <- runif(1,min(y),mu[i-1,j])
            mutilde[i,1:ktilde[i]] <- sort(c(mu[i-1,1:k[i-1]],mutilde[i,j]))

            atilde[i,j+1] <- runif(1,mutilde[i,j],mutilde[i,j+1])
            atilde[i,1:(ktilde[i]+1)] <- sort(c(a[i-1,1:ktilde[i]],atilde[i,j+1]))
            atilde[i,1] <- min(y)
            atilde[i,(ktilde[i]+1)] <- max(y)

            dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
            dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=max(y)])
            if (ktilde[i]>2) {
            for (m in 2: (ktilde[i]-1)) {
            dimtilde[i,m] <- length(y[(atilde[i,m]<=y & y<atilde[i,m+1])]) }}
            Pr <- min(dimtilde[i,1:ktilde[i]])

            ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
            classesf <- y[atilde[i,j]<=y & y<=ybar]
            classesr <- y[ybar<=y & y<=atilde[i,j+1]]

            Lf <- length(classesf[!duplicated(classesf)])
            Lr <- length(classesr[!duplicated(classesr)])
            compteur <- compteur + 1

           if ( (Pr>1 & Lf>1 & Lr>1)  || compteur==1000) break()

          }

         if (compteur==1000) {rhob[i] <- 0}

         else  {

          dltilde[i,j] <- sample(3:30,1)

          sigmaftilde[i,j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
          sigmartilde[i,j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

          sigmaftilde[i,1:ktilde[i]] <- c(sigmaf[i-1,1:k[i-1]],sigmaftilde[i,j])
          sigmartilde[i,1:ktilde[i]] <- c(sigmar[i-1,1:k[i-1]],sigmartilde[i,j])

          deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)
          deltatilde[i,1:ktilde[i]]<- c(deltatilde[i,j],delta[i-1,1:k[i-1]])

          alpha <- rep(1,k[i-1])
          alphatilde <- rep(1,ktilde[i])
          alphaproptilde <- rep(1,ktilde[i])
          alphaprop <- rep(1,k[i-1])
          wtilde[i,1:ktilde[i]] <- rdirichlet(1,alphaproptilde)
          ennetilde <- dimtilde[i,1:ktilde[i]]
          enne <- rmultinom(1,n,w[i-1,1:k[i-1]])

          #Rapport de vraisemblance

          for (m in 1:k[i-1]) {
            Kaf[,m] <- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt((yf-mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
            Kbf[,m] <- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt((yf-mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
          }
          Kaf[,ktilde[i]] <- (wtilde[i,ktilde[i]]*(1/sqrt(sigmaftilde[i,ktilde[i]]))*dt(( yf -mutilde[i,ktilde[i]]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,ktilde[i]]),dltilde[i,m]))
          for (s in 1:nf) {
            Y1f[s] <- log(sum(Kaf[s,1:ktilde[i]]))
            Y2f[s] <- log(sum(Kbf[s,1:k[i-1]]))
          }

          for (m in 1:k[i-1]) {
            Kar[,m] <- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt((yr-mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
            Kbr[,m] <- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt((yr-mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
          }
          Kar[,ktilde[i]] <- (wtilde[i,ktilde[i]]*(1/sqrt(sigmartilde[i,ktilde[i]]))*dt(( yr -mutilde[i,ktilde[i]]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,ktilde[i]]),dltilde[i,m]))
          for (s in 1:nr) {
            Y1r[s] <- log(sum(Kar[s,1:ktilde[i]]))
            Y2r[s] <- log(sum(Kbr[s,1:k[i-1]]))
          }

          Kn <- sum(Y1f)+sum(Y1r)
          Kd <- sum(Y2f)+sum(Y2r)

          q<-Kn-Kd
          rap.q<-exp(q)

          rap.vrais=rap.q

          if (j==1) {qalloc <- 1/(mu[i-1,j]-min(y))}       #densit? de mutilde[i,j]
          else {qalloc <- 1/(mu[i-1,j]-mu[i-1,j-1])}

          rap.priormu = (priormu(mutilde[i,1:ktilde[i]],y)/priormu(mu[i-1,1:k[i-1]],y))
          rap.priorw = (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
          rap.priorenne = dmultinom(ennetilde,n,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]],n,w[i-1,1:k[i-1]])
          rap.priork = (dpois(ktilde[i],lambda)/dpois(k[i-1],lambda))
          rap.propmu = (1/(qalloc))
          rap.propw = (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))

          rap.prior = rap.priormu * rap.priorw * rap.priorenne * rap.priork
          rap.prop =  rap.propmu  * rap.propw

          rhob[i] = min(1,(rap.vrais) * (rap.prior)  *  (rap.prop ) * (Dk(ktilde[i],lambda,kmax)/Bk(k[i-1],lambda,kmax)))

        }

         rhob[i] = ifelse ( is.na(rhob[i])==FALSE, rhob[i], 0)

        v=runif(1)      #Acceptation/rejet du Birth move

        if (rhob[i]>=v & ktilde[i]<=kmax) {
          k[i]<-ktilde[i]
          mu[i,1:k[i]]<- mutilde[i,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaftilde[i,1:k[i]]
          sigmar[i,1:k[i]]<- sigmartilde[i,1:k[i]]
          delta[i,1:k[i]]<- deltatilde[i,1:k[i]]
          dl[i,1:(k[i])] <- dltilde[i,1:k[i]]
          w[i,1:k[i]]<- wtilde[i,1:k[i]]
          dim[i,1:k[i]]<- dimtilde[i,1:k[i]]
          a[i,1:(k[i]+1)]<- atilde[i,1:(k[i]+1)]

        }

        else {
          k[i]<-k[i-1]
          mu[i,1:k[i]]<- mu[i-1,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaf[i-1,1:k[i]]
          sigmar[i,1:k[i]]<- sigmar[i-1,1:k[i]]
          delta[i,1:k[i]]<- delta[i-1,1:k[i]]
          dl[i,1:(k[i])] <- dl[i-1,1:k[i]]
          w[i,1:k[i]]<- w[i-1,1:k[i]]
          dim[i,1:k[i]]<- dim[i-1,1:k[i]]
          a[i,1:(k[i]+1)]<- a[i-1,1:(k[i]+1)]

         }

      } ### end of B move in case k=1

      else {

        ###Metropolis-Hastings move

         ktilde[i] <- k[i-1]
         compteur<-1
          repeat {
           j <- sample(1:k[i-1],1)
           mutilde[i,j] <- runif(1,mu[i-1,j],max(y))
           mutilde[i,1:ktilde[i]] <- sort(c(mutilde[i,1:ktilde[i]]))

            atilde[i,j] <-  min(y)
            atilde[i,j+1] <- max(y)

            dimtilde[i,1] <- length(y[atilde[i,1]<=y & y< atilde[i,2]])
            dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=max(y)])
            if (ktilde[i]>2) {
            for (m in 2: (ktilde[i]-1)) {
            dimtilde[i,m] <-length(y[(atilde[i,m]<=y & y<atilde[i,m+1])])}}
            Pr <- min(dimtilde[i,1:ktilde[i]])

            ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
            classesf <- y[atilde[i,j]<=y & y<=ybar]
            classesr <- y[ybar<=y & y<=atilde[i,j+1]]

            Lf <- length(classesf[!duplicated(classesf)])
            Lr <- length(classesr[!duplicated(classesr)])
            compteur <- compteur + 1

            if ( (Pr>1 & Lf>1 & Lr>1)  || compteur==1000) break()

           }

          if (compteur==1000) {rhomh[i] <- 0}

          else {

          sigmaftilde[i,1:ktilde[i]] <- sigmaf[i-1,1:k[i-1]]
          sigmartilde[i,1:ktilde[i]] <- sigmar[i-1,1:k[i-1]]

          dltilde[i,j] <- sample(3:30,1)

          sigmaftilde[i,j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
          sigmartilde[i,j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

          deltatilde[i,1:ktilde[i]] <- delta[i-1,1:k[i-1]]

          deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)

          alpha <- rep(1,k[i-1])
          alphatilde <- rep(1,ktilde[i])
          ennetilde <- dimtilde[i,1:ktilde[i]]
          alphaproptilde <- rep(1,ktilde[i])
          alphaprop <- rep(1,k[i-1])
          wtilde[i,1:ktilde[i]] <- rdirichlet(1,alphaproptilde)

          ### calcul du rapport de vraisemblance de M-H move

          for (m in 1:ktilde[i]) {
            Kaf[,m]<- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( yf -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
            Kbf[,m]<- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( yf -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
          }
          for (s in 1:nf) {
            Y1f[s] <-log(sum(Kaf[s,1:ktilde[i]]))
            Y2f[s] <-log(sum(Kbf[s,1:k[i-1]]))
          }

          for (m in 1:ktilde[i]) {
            Kar[,m]<- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt((yr -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
            Kbr[,m]<- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt((yr -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
          }
          for (s in 1:nr) {
            Y1r[s] <-log(sum(Kar[s,1:ktilde[i]]))
            Y2r[s] <-log(sum(Kbr[s,1:k[i-1]]))
          }

          Kn <- sum(Y1f)+sum(Y1r)
          Kd <- sum(Y2f)+sum(Y2r)

          q <- Kn-Kd
          rap.q <- exp(q)

          rap.vrais = rap.q

          rap.priormu= (priormu(mutilde[i,1:ktilde[i]],y)/priormu(mu[i-1,1:k[i-1]],y))
          rap.priorw= (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
          rap.priorenne=dmultinom(ennetilde,n,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]],n,w[i-1,1:k[i-1]])
          rap.priork= 1
          rap.propmu=1
          rap.propw= (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))

          rap.prior= rap.priormu * rap.priorw * rap.priorenne * rap.priork
          rap.prop=  rap.propmu  * rap.propw
          rhomh[i]=min(1, rap.vrais * (rap.prior)  *  (rap.prop ))

        }

         rhomh[i]= ifelse ( is.na(rhomh[i])==FALSE, rhomh[i], 0)

        v = runif(1)

        if (rhomh[i]>=v & ktilde[i]<=kmax) {
          k[i]<-ktilde[i]
          mu[i,1:k[i]]<- mutilde[i,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaftilde[i,1:k[i]]
          sigmar[i,1:k[i]]<- sigmartilde[i,1:k[i]]
          delta[i,1:k[i]]<- deltatilde[i,1:k[i]]
          dl[i,1:(k[i])] <- dltilde[i,1:k[i]]
          w[i,1:k[i]]<- wtilde[i,1:k[i]]
          dim[i,1:k[i]]<- dimtilde[i,1:k[i]]
           a[i,1:(k[i]+1)]<- atilde[i,1:(k[i]+1)]

         }

        else {
          k[i]<-k[i-1]
          mu[i,1:k[i]]<- mu[i-1,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaf[i-1,1:k[i]]
          sigmar[i,1:k[i]]<- sigmar[i-1,1:k[i]]
          delta[i,1:k[i]]<- delta[i-1,1:k[i]]
          dl[i,1:(k[i])] <- dl[i-1,1:k[i]]
          w[i,1:k[i]]<- w[i-1,1:k[i]]
          dim[i,1:k[i]]<- dim[i-1,1:k[i]]
         a[i,1:(k[i]+1)]<- a[i-1,1:(k[i]+1)]


        }


      }    ### end of M-H move in case k=1

    }  ### end of test of k=1

    else {

      u<-runif(1)

      if (u<= Dk(k[i-1],lambda,kmax)) {

        ### Death move

        ktilde[i]<-k[i-1]-1
         compteur=1
          repeat {

            j <- sample(1:k[i-1],1)
            X <- mu[i-1,1:k[i-1]]
            mutilde[i,1:ktilde[i]] <- sort(c(X[-j]))

            Yf <- sigmaf[i-1,1:k[i-1]]
            sigmaftilde[i,1:ktilde[i]] <- Yf[-j]

            Yr <- sigmar[i-1,1:k[i-1]]
            sigmartilde[i,1:ktilde[i]] <- Yr[-j]

            Delta <- delta[i-1,1:k[i-1]]
            deltatilde[i,1:ktilde[i]] <- Delta[-j]

            Dl <-dl[i-1,1:k[i-1]]
            dltilde[i,1:ktilde[i]] <- Dl[-j]

            Z <- a[i-1,1:(k[i-1]+1)]
            if (j==k[i-1]) {atilde[i,1:(ktilde[i]+1)] <- Z[-j]}
            else { atilde[i,1:(ktilde[i]+1)] <- Z[-(j+1)] }
            atilde[i,1] <- min(y)
            atilde[i,ktilde[i]+1] <- max(y)

            dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
            dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=max(y)])
            if (ktilde[i]>2) {
            for (m in 2: (ktilde[i]-1)) {
            dimtilde[i,m] <- length(y[atilde[i,m]<=y & y<atilde[i,m+1]]) }}
            Pr <- min(dimtilde[i,1:ktilde[i]])

            ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
            classesf <- y[atilde[i,j]<=y & y<=ybar]
            classesr <- y[ybar<=y & y<=atilde[i,j+1]]

            Lf <- length(classesf[!duplicated(classesf)])
            Lr <- length(classesr[!duplicated(classesr)])
            compteur = compteur + 1

          if ( (Pr>1 & Lf>1 & Lr>1)  || compteur==1000) break()

          }

          if (compteur==1000) {rhod[i] <- 0}

          else {

          alpha <- rep(1,k[i-1])
          alphatilde <- rep(1,ktilde[i])
          alphaproptilde <- rep(1,ktilde[i])
          alphaprop <- rep(1,k[i-1])
          ennetilde <- dimtilde[i,1:ktilde[i]]
          enne <- rmultinom(1,n,w[i-1,1:k[i-1]])
          wtilde[i,1:ktilde[i]]<-rdirichlet(1,alphaproptilde)

          ### Rapport de vraisemblance ###

          for (m in 1:ktilde[i]) {
            Kaf[,m]<- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( yf -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
            Kbf[,m]<- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( yf -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
          }
          Kbf[,k[i-1]]<- (w[i-1,k[i-1]]*(1/sqrt(sigmaf[i-1,k[i-1]]))*dt(( yf -mu[i-1,k[i-1]]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,k[i-1]]),dltilde[i,m]))

          for (s in 1:nf) {
            Y1f[s] <-log(sum(Kaf[s,1:ktilde[i]]))
            Y2f[s] <-log(sum(Kbf[s,1:k[i-1]]))
          }

          for (m in 1:ktilde[i]) {
            Kar[,m]<- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( yr -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
            Kbr[,m]<- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( yr -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
          }
          Kbr[,k[i-1]]<- (w[i-1,k[i-1]]*(1/sqrt(sigmar[i-1,k[i-1]]))*dt(( yr -mu[i-1,k[i-1]]-delta[i-1,m]/2)/sqrt(sigmar[i-1,k[i-1]]),dl[i-1,m]))

          for (s in 1:nr) {
            Y1r[s] <-log(sum(Kar[s,1:ktilde[i]]))
            Y2r[s] <-log(sum(Kbr[s,1:k[i-1]]))
          }

          Kn<-sum(Y1f)+sum(Y1r)
          Kd<-sum(Y2f)+sum(Y2r)

          q=Kn-Kd
          rap.q<-exp(q)

          rap.vrais=rap.q

          if (j==1) {qalloc <- 1/(mu[i-1,j+1]-min(y))}  #densit? de mutilde[i,j]
          else { if (j==k[i-1]) { qalloc <- 1/(max(y)-mu[i-1,j-1])}
          else { qalloc <- 1/(mu[i-1,j+1]-mu[i-1,j-1]) }}

          rap.priormu = (priormu(mutilde[i,1:ktilde[i]],y)/priormu(mu[i-1,1:k[i-1]],y))
          rap.priorw = (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
          rap.priorenne=dmultinom(ennetilde,n,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]],n,w[i-1,1:k[i-1]])
          rap.priork= (dpois(ktilde[i],lambda)/dpois(k[i-1],lambda))
          rap.propmu= (qalloc)
          rap.propw= (ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde)/ddirichlet(w[i-1,1:k[i-1]],alphaprop))

          rap.prior = rap.priormu * rap.priorw * rap.priorenne * rap.priork
          rap.prop =  rap.propmu  * rap.propw

          rhod[i]= min(1,(rap.vrais) * (rap.prior)  *  (rap.prop) * (Bk(ktilde[i],lambda,kmax)/Dk(k[i-1],lambda,kmax)))

        }

        rhod[i] =  ifelse ( is.na(rhod[i])==FALSE, rhod[i], 0)

        v = runif(1)      #Acceptation/rejet du Death move

        if (rhod[i]>=v & ktilde[i]<=kmax ) {
          k[i]<-ktilde[i]
          mu[i,1:k[i]]<- mutilde[i,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaftilde[i,1:k[i]]
          sigmar[i,1:k[i]]<- sigmartilde[i,1:k[i]]
          delta[i,1:k[i]]<- deltatilde[i,1:k[i]]
          dl[i,1:k[i]] <- dltilde[i,1:k[i]]
          w[i,1:k[i]]<- wtilde[i,1:k[i]]
          dim[i,1:k[i]]<- dimtilde[i,1:k[i]]
           a[i,1:(k[i]+1)]<- atilde[i,1:(k[i]+1)]

          }

        else {
          k[i]<-k[i-1]
          mu[i,1:k[i]]<- mu[i-1,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaf[i-1,1:k[i]]
          sigmar[i,1:k[i]]<- sigmar[i-1,1:k[i]]
          delta[i,1:k[i]]<- delta[i-1,1:k[i]]
          dl[i,1:k[i]] <- dl[i-1,1:k[i]]
          w[i,1:k[i]]<- w[i-1,1:k[i]]
          dim[i,1:k[i]]<- dim[i-1,1:k[i]]
          a[i,1:(k[i]+1)]<- a[i-1,1:(k[i]+1)]

        }

      }

      else {  if (u<= (Dk(k[i-1],lambda,kmax)+Bk(k[i-1],lambda,kmax))) {

        #Birth move

         ktilde[i] <- k[i-1] + 1
         compteur <- 1
          repeat {
          j <- sample(1:k[i-1],1)
           if (j==1) {mutilde[i,j] <- runif(1,min(y),mu[i-1,j])}
           else {
           mutilde[i,j]<-runif(1,mu[i-1,j-1],mu[i-1,j])}
           mutilde[i,1:ktilde[i]]=sort(c(mu[i-1,1:k[i-1]],mutilde[i,j]))

            atilde[i,j+1]<- ifelse(j<ktilde[i-1], runif(1,mutilde[i,j],mutilde[i,j+1]), runif(1,mutilde[i,j],max(y)))
            atilde[i,1:(ktilde[i]+1)]<- sort(c(a[i-1,1:ktilde[i]],atilde[i,j+1]))
            if (j==1){
            atilde[i,j]<-min(y)}
            else {
            atilde[i,j]<-runif(1,mutilde[i,j-1],mutilde[i,j])}
            atilde[i,1]<-min(y)
            atilde[i,ktilde[i]+1]<-max(y)

            dimtilde[i,1]=length(y[atilde[i,1]<=y & y<atilde[i,2]])
            dimtilde[i,ktilde[i]]=length(y[atilde[i,ktilde[i]]<=y & y<=max(y)])
            if (ktilde[i]>2) {
            for (m in 2: (ktilde[i]-1)) {
            dimtilde[i,m]<-length(y[atilde[i,m]<=y & y<atilde[i,m+1]]) }}
            Pr <- min(dimtilde[i,1:ktilde[i]])

            ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
            classesf <- y[atilde[i,j]<=y & y<=ybar]
            classesr <- y[ybar<=y & y<=atilde[i,j+1]]

            Lf <- length(classesf[!duplicated(classesf)])
            Lr <- length(classesr[!duplicated(classesr)])
            compteur = compteur+1
           if ( (Pr>1 & Lf>1 & Lr>1)  || compteur==1000) break()

          }

          if (compteur==1000) {rhob[i] <- 0}

         else {

          dltilde[i,j]=sample(3:30,1)
          sigmaftilde[i,j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
          sigmartilde[i,j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

          if (j==1) {
            sigmaftilde[i,1:ktilde[i]]<- c(sigmaf[i-1,1:k[i-1]],sigmaftilde[i,j])
            sigmartilde[i,1:ktilde[i]]<- c(sigmar[i-1,1:k[i-1]],sigmartilde[i,j] ) }
            else {
            sigmaftilde[i,1:ktilde[i]] <- c(sigmaf[i-1,1:(j-1)],sigmaftilde[i,j],sigmaf[i-1,j:k[i-1]])
            sigmartilde[i,1:ktilde[i]] <- c(sigmar[i-1,1:(j-1)],sigmartilde[i,j],sigmar[i-1,j:k[i-1]]) }

           deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)
          if (j==1) {
            deltatilde[i,1:ktilde[i]] <- c(deltatilde[i,j],delta[i-1,1:k[i-1]]) }
          else if (j==k[i-1]) {
            deltatilde[i,1:ktilde[i]] <- c(delta[i-1,1:k[i-1]],deltatilde[i,j]) }
          else {
            deltatilde[i,1:ktilde[i]] <- c(delta[i-1,1:(j-1)],deltatilde[i,j],delta[i-1,j:k[i-1]]) }

          alpha <- rep(1,k[i-1])
          alphatilde <- rep(1,ktilde[i])
          alphaproptilde <- rep(1,ktilde[i])
          alphaprop <- rep(1,k[i-1])
          ennetilde <- dimtilde[i,1:ktilde[i]]
          enne <- rmultinom(1,n,w[i-1,1:k[i-1]])
          wtilde[i,1:ktilde[i]] <- rdirichlet(1,alphaproptilde)

          ### Rapport de vraisemblance ###

          for (m in 1:k[i-1]) {
            Kaf[,m]<- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( yf -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
            Kbf[,m]<- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( yf -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
          }
          Kaf[,ktilde[i]]<- (wtilde[i,ktilde[i]]*(1/sqrt(sigmaftilde[i,ktilde[i]]))*dt(( yf -mutilde[i,ktilde[i]]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,ktilde[i]]),dltilde[i,m]))
          for (s in 1:nf) {
            Y1f[s] <-log(sum(Kaf[s,1:ktilde[i]]))
            Y2f[s] <-log(sum(Kbf[s,1:k[i-1]]))
          }

          for (m in 1:k[i-1]) {
            Kar[,m]<- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( yr -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
            Kbr[,m]<- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( yr -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
          }
          Kar[,ktilde[i]]<- (wtilde[i,ktilde[i]]*(1/sqrt(sigmartilde[i,ktilde[i]]))*dt(( yr -mutilde[i,ktilde[i]]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,ktilde[i]]),dltilde[i,m]))
          for (s in 1:nr) {
            Y1r[s] <-log(sum(Kar[s,1:ktilde[i]]))
            Y2r[s] <-log(sum(Kbr[s,1:k[i-1]]))
          }

          Kn <- sum(Y1f)+sum(Y1r)
          Kd <- sum(Y2f)+sum(Y2r)

          q=Kn-Kd
          rap.q<-exp(q)

          rap.vrais=rap.q

          if (j==1) { qalloc=1/(mu[i-1,j]-min(y)) }       #densit? de mutilde[i,j]
          else { qalloc<-1/(mu[i-1,j]-mu[i-1,j-1]) }

          rap.priormu= (priormu(mutilde[i,1:ktilde[i]],y)/priormu(mu[i-1,1:k[i-1]],y))
          rap.priorw= (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
          rap.priorenne=dmultinom(ennetilde,n,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]],n,w[i-1,1:k[i-1]])
          rap.priork= (dpois(ktilde[i],lambda)/dpois(k[i-1],lambda))
          rap.propmu= (1/(qalloc))
          rap.propw= (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))

          rap.prior = rap.priormu * rap.priorw * rap.priorenne * rap.priork
          rap.prop =  rap.propmu  * rap.propw
          rhob[i]= min(1, (rap.vrais) * (rap.prior)  *  (rap.prop ) * (Dk(ktilde[i],lambda,kmax)/Bk(k[i-1],lambda,kmax)))
        }

        v = runif(1)      # Acceptation/rejet du Birth move
        rhob[i] = ifelse ( is.na(rhob[i])==FALSE, rhob[i], 0)

          if (rhob[i]>=v & ktilde[i]<=kmax) {
          k[i]<-ktilde[i]
          mu[i,1:k[i]]<- mutilde[i,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaftilde[i,1:k[i]]
          sigmar[i,1:k[i]]<- sigmartilde[i,1:k[i]]
          delta[i,1:k[i]]<- deltatilde[i,1:k[i]]
          dl[i,1:(k[i])] <- dltilde[i,1:k[i]]
          w[i,1:k[i]]<- wtilde[i,1:k[i]]
          dim[i,1:k[i]]<- dimtilde[i,1:k[i]]
          a[i,1:(k[i]+1)]<- atilde[i,1:(k[i]+1)]

        }

        else {
          k[i]<-k[i-1]
          mu[i,1:k[i]]<- mu[i-1,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaf[i-1,1:k[i]]
          sigmar[i,1:k[i]]<- sigmar[i-1,1:k[i]]
          delta[i,1:k[i]]<- delta[i-1,1:k[i]]
          dl[i,1:(k[i])] <- dl[i-1,1:k[i]]
          w[i,1:k[i]]<- w[i-1,1:k[i]]
          dim[i,1:k[i]]<- dim[i-1,1:k[i]]
          a[i,1:(k[i]+1)]<- a[i-1,1:(k[i]+1)]

        }

      }

      else {

        ### Metropolis-Hastings move

           ktilde[i] <- k[i-1]
            compteur <- 1
            repeat {
            j <- sample(2:k[i-1],1)
            mutilde[i,1:ktilde[i]] <- mu[i-1,1:k[i-1]]
            if (j==1) {
            mutilde[i,j] <- runif(1,min(y),mu[i-1,j+1])}
            else {if (j==ktilde[i]) {
            mutilde[i,j] <- runif(1,mu[i-1,j-1],max(y))}
            else { mutilde[i,j] <- runif(1,mu[i-1,j-1],mu[i-1,j+1]) }}
            mutilde[i,1:ktilde[i]] <- sort(c(mutilde[i,1:ktilde[i]]))

            atilde[i,1:(ktilde[i]+1)] <- sort(c(a[i-1,1:(k[i-1]+1)]))
            if (j==ktilde[i]) {
            atilde[i,j] <- runif(1,mutilde[i,j],max(y))
            atilde[i,j+1] <- max(y) }
            else { if (j==1) {
            atilde[i,j] <- min(y)
            atilde[i,j+1] <- runif(1,mutilde[i,j],mutilde[i,j+1])}
            else {
            atilde[i,j] <- runif(1,mutilde[i,j-1],mutilde[i,j])
            atilde[i,j+1] <- runif(1,mutilde[i,j],mutilde[i,j+1])  }}
            atilde[i,1] <- min(y)
            atilde[i,ktilde[i]+1] <- max(y)

            dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
            dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=max(y)])
            if (ktilde[i]>2) {
            for (m in 2: (ktilde[i]-1)) {
            dimtilde[i,m]<-length(y[atilde[i,m]<=y & y<atilde[i,m+1]]) }}
            Pr <- min(dimtilde[i,1:ktilde[i]])

            ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
            classesf <- y[atilde[i,j]<=y & y<=ybar]
            classesr <- y[ybar<=y & y<=atilde[i,j+1]]

            Lf <- length(classesf[!duplicated(classesf)])
            Lr <- length(classesr[!duplicated(classesr)])
            compteur <- compteur + 1
           if ( (Pr>1 & Lf>1 & Lr>1)  || compteur==1000) break()
          }

          if (compteur==1000) {rhomh[i] <- 0}

          else {

          sigmaftilde[i,1:ktilde[i]] <- sigmaf[i-1,1:k[i-1]]
          sigmartilde[i,1:ktilde[i]] <- sigmar[i-1,1:k[i-1]]

          dltilde[i,j]=sample(3:30,1)

          sigmaftilde[i,j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
          sigmartilde[i,j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

          deltatilde[i,1:ktilde[i]] <- delta[i-1,1:k[i-1]]
           deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)

          alpha <- rep(1,k[i-1])
          alphatilde <- rep(1,ktilde[i])
          ennetilde <- dimtilde[i,1:ktilde[i]]
          alphaproptilde <- rep(1,ktilde[i])
          alphaprop <- rep(1,k[i-1])
          wtilde[i,1:ktilde[i]] <- rdirichlet(1,alphaproptilde)

          ###calcul du rapport de vraisemblance de M-H move

          for (m in 1:ktilde[i]) {
            Kaf[,m]<- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( yf -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
            Kbf[,m]<- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( yf -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
          }
          for (s in 1:nf) {
            Y1f[s] <-log(sum(Kaf[s,1:ktilde[i]]))
            Y2f[s] <-log(sum(Kbf[s,1:k[i-1]]))
          }

          for (m in 1:ktilde[i]) {
            Kar[,m]<- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( yr -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
            Kbr[,m]<- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( yr -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
          }
          for (s in 1:nr) {
            Y1r[s] <-log(sum(Kar[s,1:ktilde[i]]))
            Y2r[s] <-log(sum(Kbr[s,1:k[i-1]]))
          }

          Kn<-sum(Y1f)+sum(Y1r)
          Kd<-sum(Y2f)+sum(Y2r)

          q=Kn-Kd
          rap.q<-exp(q)

          rap.vrais = rap.q

          rap.priormu= (priormu(mutilde[i,1:ktilde[i]],y)/priormu(mu[i-1,1:k[i-1]],y))
          rap.priorw= (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha))
          rap.priorenne=dmultinom(ennetilde,n,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]],n,w[i-1,1:k[i-1]])
          rap.priork=(1)
          rap.propmu=(1)
          rap.propw= (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))
          rap.prior= rap.priormu * rap.priorw * rap.priorenne * rap.priork
          rap.prop=  rap.propmu  * rap.propw
          rhomh[i]=min(1, rap.vrais * (rap.prior)  *  (rap.prop ))

        }

        rhomh[i] <- ifelse ( is.na(rhomh[i])==FALSE, rhomh[i], 0)

        v=runif(1)

        if (rhomh[i]>=v & ktilde[i]<=kmax ) {
          k[i]<-ktilde[i]
          mu[i,1:k[i]]<- mutilde[i,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaftilde[i,1:k[i]]
          sigmar[i,1:k[i]]<- sigmartilde[i,1:k[i]]
          delta[i,1:k[i]]<- deltatilde[i,1:k[i]]
          dl[i,1:(k[i])] <- dltilde[i,1:k[i]]
          w[i,1:k[i]]<- wtilde[i,1:k[i]]
          dim[i,1:k[i]]<- dimtilde[i,1:k[i]]
         a[i,1:(k[i]+1)]<- atilde[i,1:(k[i]+1)]

         }

          else {
          k[i]<-k[i-1]
          mu[i,1:k[i]]<- mu[i-1,1:k[i]]
          sigmaf[i,1:k[i]]<- sigmaf[i-1,1:k[i]]
          sigmar[i,1:k[i]]<- sigmar[i-1,1:k[i]]
          delta[i,1:k[i]]<- delta[i-1,1:k[i]]
          dl[i,1:(k[i])] <- dl[i-1,1:k[i]]
          w[i,1:k[i]]<- w[i-1,1:k[i]]
          dim[i,1:k[i]]<- dim[i-1,1:k[i]]
          a[i,1:(k[i]+1)]<- a[i-1,1:(k[i]+1)]

         }
       }
       } #end of else of B move and M-H move
      } ###end of moves in case k>=2
      } ###end of boucle RJMCMC

liste = rep(list(NULL),niter)
 for (i in 1:niter)
 {
   new.list <- list(k=k[i],
                    mu=mu[i,1:k[i]],
                    sigmaf=sigmaf[i,1:k[i]],
                    sigmar=sigmar[i,1:k[i]],
                    delta=delta[i,1:k[i]],
                    dl=dl[i,1:k[i]],
                    w=w[i,1:k[i]]
                 )

   liste[[i]] <- merge(yf, yr, y, new.list, ecartmin, ecartmax, minReads)
 }

  kmax = max(kmax,sapply(1:niter,function(i){liste[[i]]$k}))
  mu <- matrix(0,nrow=niter,ncol=kmax)
  sigmaf <- matrix(0,nrow=niter,ncol=kmax)
  sigmar <- matrix(0,nrow=niter,ncol=kmax)
  delta <- matrix(0,nrow=niter,ncol=kmax)
  w <- matrix(0,nrow=niter,ncol=kmax)
  dl <- matrix(0,nrow=niter,ncol=kmax)

 for (i in 1:niter)
 {
   k[i]             <- liste[[i]]$k
   mu[i,1:k[i]]     <- liste[[i]]$mu
   sigmaf[i,1:k[i]] <- liste[[i]]$sigmaf
   sigmar[i,1:k[i]] <- liste[[i]]$sigmar
   delta[i,1:k[i]]  <- liste[[i]]$delta
   dl[i,1:k[i]]     <- liste[[i]]$dl
}

  km <- mode(k)
  K <- which(k==km)

  mu_hat     <- sapply(1:km,function(j){mean(mu[K,j])})
  sigmaf_hat <- sapply(1:km,function(j){mean(sigmaf[K,j])})
  sigmar_hat <- sapply(1:km,function(j){mean(sigmar[K,j])})
  w_hat      <- sapply(1:km,function(j){mean(w[K,j])})
  delta_hat  <- sapply(1:km,function(j){mean(delta[K,j])})
  dl_hat     <- round(sapply(1:km,function(j){mean(dl[K,j])}))

  qmu <- matrix(0,nrow=km,ncol=2)
  colnames(qmu) <- c("2.5%", "97.5%")
  qsigmaf <- matrix(0,nrow=km,ncol=2)
  colnames(qsigmaf) <- c("2.5%", "97.5%")
  qsigmar <- matrix(0,nrow=km,ncol=2)
  colnames(qsigmar) <- c("2.5%", "97.5%")
  qdelta <- matrix(0,nrow=km,ncol=2)
  colnames(qdelta) <- c("2.5%", "97.5%")
  qdl <- matrix(0,nrow=km,ncol=2)
  colnames(qdl) <- c("2.5%", "97.5%")
  qw <- matrix(0,nrow=km,ncol=2)
  colnames(qw) <- c("2.5%", "97.5%")

  for (j in 1:km)
  {
  qmu[j,] <- quantile(mu[,j], probs=c(0.025,0.975), names=FALSE)
  qsigmaf[j,] <- quantile(sigmaf[,j], probs=c(0.025,0.975), names=FALSE)
  qsigmar[j,] <- quantile(sigmar[,j], probs=c(0.025,0.975), names=FALSE)
  qdelta[j,] <- quantile(delta[,j], probs=c(0.025,0.975), names=FALSE)
  qdl[j,] <- quantile(dl[,j], probs=c(0.025,0.975), names=FALSE)
  qw[j,] <- quantile(w[,j], probs=c(0.025,0.975), names=FALSE)
  }

  liste <- list(
                K=k,
                k=km,
                mu=mu_hat,
                sigmaf=sigmaf_hat,
                sigmar=sigmar_hat,
                delta=delta_hat,
                dl=dl_hat,
                w=w_hat,
                qmu=qmu,
                qsigmaf=qsigmaf,
                qsigmar=qsigmar,
                qdelta=qdelta,
                qdl=qdl,
                qw=qw
               )
 return(liste)

}


