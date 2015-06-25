n<-100000
set.seed(101)
system.time({
    for (i in 1:n) {
        rjmcmc:::tnormale(10, 2, 9, 11)
    }
})

