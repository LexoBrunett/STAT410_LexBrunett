   w <- 48
    n <- 64
    alpha <- 0.05
    (pi.hat <- w/n)
    (pi.tilde <- (w + qnorm(p = 1-alpha/2)^2 /2) / (n + qnorm(p = 1-alpha/2)^2))
    wilson <- pi.tilde + qnorm(p=c(alpha/2, 1-alpha/2))*sqrt(n)/
              (n+qnorm(p=1-alpha/2)^2)*sqrt(pi.hat*(1-pi.hat) +
                                      qnorm(p = 1-alpha/2)^2/(4*n))
    round(wilson, digits=4)
    library(package=binom)
    binom.confint(x = w, n = n, conf.level = 1-alpha, methods = "wilson")

#ConfLevelWaldOnly.R



  # Initial settings and calculations
    alpha <- 0.05
    n <- 40
    w <- 0:n
    pi.hat <- w/n
    pi.seq <- seq(from = 0.001, to = 0.999, by = 0.0005)
  # Wald
    var.wald   <- pi.hat*(1-pi.hat)/n
    lower.wald <- pi.hat - qnorm(p = 1-alpha/2) * sqrt(var.wald)
    upper.wald <- pi.hat + qnorm(p = 1-alpha/2) * sqrt(var.wald)
  # Save true confidence levels in a matrix
    save.true.conf <- matrix(data = NA, nrow = length(pi.seq), ncol = 2)
  # Create counter for the loop
    counter <- 1
  # Loop over each pi that the true confidence level is calculated on
    for(pi in pi.seq) {
      pmf <- dbinom(x=w, size=n, prob=pi)
      save.wald <- ifelse(test=pi>lower.wald,
                          yes=ifelse(test = pi<upper.wald, yes=1, no=0),
                          no = 0)
      wald <- sum(save.wald*pmf)
      save.true.conf[counter,] <- c(pi, wald)
      # print(save.true.conf[counter,])  # Uncomment to view results one-by-one
      counter <- counter+1
    }
    head(save.true.conf, 150)
    tail(save.true.conf)
    # Plot
    #x11(width = 7, height = 6, pointsize = 12)
    plot(x=save.true.conf[,1], y=save.true.conf[,2],
         main="Wald", xlab=expression(pi),
         ylab="True confidence level", type="l", ylim=c(0.85,1))
    abline(h = 1-alpha, lty = "dotted", lwd=3)

  # Alternatively
    library(binom)
    binom.coverage(p = 0.16, n = n, conf.level = 0.95, method = "asymptotic")
  
    binom.plot(n = 40, method = binom.asymp, np = 500, conf.level = 0.95)

#ConfLevelTwoProb.R

 # Initial settings
    alpha <- 0.05
    pi1   <- 0.2
    pi2   <- 0.4
    n1    <- 10
    n2    <- 10
    numb.bin.samples <- 1000

  # Estimated true confidence level
  # Simulate w1 and w2
    set.seed(2349)
    w1 <- rbinom(n=numb.bin.samples, size = n1, prob = pi1)
    w2 <- rbinom(n=numb.bin.samples, size = n2, prob = pi2)
    pi.hat1 <- w1/n1
    pi.hat2 <- w2/n2
  # Wald
    var.wald <- pi.hat1*(1-pi.hat1) / n1 + pi.hat2*(1-pi.hat2) / n2
    lower <- pi.hat1 - pi.hat2 - qnorm(p = 1-alpha/2) * sqrt(var.wald)
    upper <- pi.hat1 - pi.hat2 + qnorm(p = 1-alpha/2) * sqrt(var.wald)
  # Intervals 1-5
    data.frame(w1, w2, lower, upper)[1:5,]
    dim(data.frame(w1, w2, lower, upper))
  # Calculate estimated true confidence level
    save <- ifelse(test=pi1-pi2 > lower,
                   yes=ifelse(test = pi1-pi2 < upper, yes = 1, no = 0),
                   no =0)
    save[1:5]
    true.conf <- mean(save)
    round(true.conf, 4)
  # Agresti-Caffo
    pi.tilde1 <- (w1+1)/(n1+2)
    pi.tilde2 <- (w2+1)/(n2+2)
    var.AC <- pi.tilde1*(1-pi.tilde1) / (n1+2) +
              pi.tilde2*(1-pi.tilde2) / (n2+2)
    lower.AC <- pi.tilde1 - pi.tilde2 - qnorm(p = 1-alpha/2) * sqrt(var.AC)
    upper.AC <- pi.tilde1 - pi.tilde2 + qnorm(p = 1-alpha/2) * sqrt(var.AC)
    save.AC <- ifelse(test = pi1-pi2 > lower.AC,
                    yes=ifelse(test = pi1-pi2 < upper.AC, yes = 1, no = 0),
                    no = 0)
    save.AC[1:10]
    true.conf.AC <- mean(save.AC)
    round(true.conf.AC, 4)

  # Estimated true confidence level holding pi2 fixed at 0.3
    numb.bin.samples<-10000
    pi1seq <- seq(from = 0.001, to = 0.999, by = 0.0005)
  # pi1seq<-0.2  # Testing
  # pi1seq<-seq(from = 0.1, to = 0.9, by = 0.1)  # Testing
  # Save true confidence levels in a matrix
    save.true.conf <- matrix(data = NA, nrow = length(pi1seq), ncol = 3)
  # Create counter for the loop
    counter <- 1
    set.seed(2114)
  # Loop over each pi1 that the true confidence level is calculated on
    for(pi1 in pi1seq) {
      w1 <- rbinom(n = numb.bin.samples, size = n1, prob = pi1)
      w2 <- rbinom(n = numb.bin.samples, size = n2, prob = pi2)
      pi.hat1 <- w1/n1
      pi.hat2 <- w2/n2
    # Wald
      lower <- pi.hat1 - pi.hat2 - qnorm(p = 1-alpha/2) *
        sqrt(pi.hat1*(1-pi.hat1) / n1 + pi.hat2*(1-pi.hat2) / n2)
      upper <- pi.hat1 - pi.hat2 + qnorm(p = 1-alpha/2) *
        sqrt(pi.hat1*(1-pi.hat1) / n1 + pi.hat2*(1-pi.hat2) / n2)
      save <- ifelse(test = pi1-pi2 > lower,
                   yes=ifelse(test=pi1-pi2 < upper, yes=1, no=0),
                   no =0)
      wald <- mean(save)
    # Agresti-Caffo
      pi.tilde1 <- (w1+1)/(n1+2)
      pi.tilde2 <- (w2+1)/(n2+2)
      lower.AC <- pi.tilde1 - pi.tilde2 - qnorm(p = 1-alpha/2) *
        sqrt(pi.tilde1*(1-pi.tilde1) / (n1+2) +
             pi.tilde2*(1-pi.tilde2) / (n2+2))
      upper.AC <- pi.tilde1 - pi.tilde2 + qnorm(p = 1-alpha/2) *
        sqrt(pi.tilde1*(1-pi.tilde1) / (n1+2) +
             pi.tilde2*(1-pi.tilde2) / (n2+2))
      save.AC <- ifelse(test = pi1-pi2 > lower.AC,
                        yes=ifelse(test=pi1-pi2 < upper.AC, yes=1, no=0),
                        no = 0)
      AC <- mean(save.AC)
      save.true.conf[counter,] <- c(pi1, wald, AC)
      counter <- counter+1
    }
  # Plot
    plot(x=save.true.conf[,1], y=save.true.conf[,2],
         xlab = expression(pi[1]),
         ylab = "Estimated true confidence level",
         type = "l", ylim = c(0.85,1), lty = "solid", col = "blue")
    lines(x = save.true.conf[,1], y = save.true.conf[,3],
          lty = "dashed", col = "red")
    abline(h = 1-alpha, lty = "dotted")
    legend(x = 0.1, y = 0.88, legend=c("Wald", "Agresti-Caffo"),
           lty = c("solid", "dashed"), bty = "n", col = c("blue", "red"))



###########################################################################
# True confidence level

  # All possible combinations of w1 and w2
  w.all<-expand.grid(w1 = 0:n1, w2 = 0:n2)
  
  # All possible combinations of pi^_1 and pi^_2
  pi.hat1<-(0:n1)/n1
  pi.hat2<-(0:n2)/n2
  pi.hat.all<-expand.grid(pi.hat1 = pi.hat1, pi.hat2 = pi.hat2)
 
  # Find joint probability for w1 and w2
  prob.w1<-dbinom(x = 0:n1, size = n1, prob = pi1)
  prob.w2<-dbinom(x = 0:n2, size = n2, prob = pi2)
  prob.all<-expand.grid(prob.w1 = prob.w1, prob.w2 = prob.w2)
  pmf<-prob.all$prob.w1*prob.all$prob.w2
  
  # Joint probability of observing w1 and w2 (i.e., P(W1 = w1, W2 = w2))
  head(data.frame(w.all, pmf = round(pmf,4)))
  tail(data.frame(w.all, pmf = round(pmf,4)))
  
  # Wald
  var.wald<-pi.hat.all[,1]*(1-pi.hat.all[,1]) / n1 + pi.hat.all[,2]*(1-pi.hat.all[,2]) / n2
  lower<-pi.hat.all[,1] - pi.hat.all[,2] - qnorm(p = 1-alpha/2) * sqrt(var.wald)
  upper<-pi.hat.all[,1] - pi.hat.all[,2] + qnorm(p = 1-alpha/2) * sqrt(var.wald)
  save<-ifelse(test = pi1-pi2 > lower,
               yes = ifelse(test = pi1-pi2 < upper, yes = 1, no = 0), no = 0)
  sum(save*pmf)
  data.frame(w.all, round(data.frame(pmf, lower, upper),4), save)[1:15,] #Example
  
  # Agresti-Caffo
  pi1tilde<-(0:n1+1)/(n1+2)
  pi2tilde<-(0:n2+1)/(n2+2)
  pi.all.tilde<-expand.grid(pi1tilde = pi1tilde, pi2tilde = pi2tilde)
  var.ac<-pi.all.tilde[,1]*(1-pi.all.tilde[,1]) / (n1+2) +
          pi.all.tilde[,2]*(1-pi.all.tilde[,2]) / (n2+2)
  lower.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] - qnorm(p = 1-alpha/2) * sqrt(var.ac)
  upper.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] + qnorm(p = 1-alpha/2) * sqrt(var.ac)
  save.AC<-ifelse(test = pi1-pi2 > lower.AC,
                  yes = ifelse(test = pi1-pi2 < upper.AC, yes = 1, no = 0), no = 0)
  sum(save.AC*pmf)
  data.frame(w.all, round(data.frame(pmf, lower, upper),4), save)[1:15,]  #Example

  
  
###########################################################################
# True confidence level holding pi2 fixed
 
  pi1seq<-seq(from = 0.001, to = 0.999, by = 0.0005)
  # pi1seq<-0.2  # Testing
  # pi1seq<-seq(from = 0.1, to = 0.9, by = 0.1)  # Testing

  # Save true confidence levels in a matrix
  save.true.conf<-matrix(data = NA, nrow = length(pi1seq), ncol = 3)

  # Create counter for the loop
  counter<-1

  # All possible combinations of w1 and w2
  w.all<-expand.grid(w1 = 0:n1, w2 = 0:n2)

  # All possible combinations of pi^_1 and pi^_2
  pi.hat1<-0:n1/n1
  pi.hat2<-0:n2/n2
  pi.hat.all<-expand.grid(pi.hat1 = pi.hat1, pi.hat2 = pi.hat2)
  
  # Wald
  lower<-pi.hat.all[,1] - pi.hat.all[,2] - qnorm(p = 1-alpha/2) * 
         sqrt(pi.hat.all[,1]*(1-pi.hat.all[,1]) / n1 + pi.hat.all[,2]*(1-pi.hat.all[,2]) / n2)
  upper<-pi.hat.all[,1] - pi.hat.all[,2] + qnorm(p = 1-alpha/2) * 
         sqrt(pi.hat.all[,1]*(1-pi.hat.all[,1]) / n1 + pi.hat.all[,2]*(1-pi.hat.all[,2]) / n2)

  # Agresti-Caffo
  pi1tilde<-(0:n1+1)/(n1+2)
  pi2tilde<-(0:n2+1)/(n2+2)
  pi.all.tilde<-expand.grid(pi1tilde = pi1tilde, pi2tilde = pi2tilde)
  lower.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] - qnorm(p = 1-alpha/2) *
            sqrt(pi.all.tilde[,1]*(1-pi.all.tilde[,1]) / (n1+2) +
              pi.all.tilde[,2]*(1-pi.all.tilde[,2]) / (n2+2))
  upper.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] + qnorm(p = 1-alpha/2) *
            sqrt(pi.all.tilde[,1]*(1-pi.all.tilde[,1]) / (n1+2) +
              pi.all.tilde[,2]*(1-pi.all.tilde[,2]) / (n2+2))


  # Loop over each pi1 that the true confidence level is calculated on
  for(pi1 in pi1seq) {

    # Find joint probability for w1 and w2
    prob.w1<-dbinom(x = 0:n1, size = n1, prob = pi1)
    prob.w2<-dbinom(x = 0:n2, size = n2, prob = pi2)
    prob.all<-expand.grid(prob.w1 = prob.w1, prob.w2 = prob.w2)
    pmf<-prob.all$prob.w1*prob.all$prob.w2
   
    # Wald
    save<-ifelse(test = pi1-pi2 > lower,
                 yes = ifelse(test = pi1-pi2 < upper, yes = 1, no = 0), no = 0)
    wald<-sum(save*pmf)

    # Agresti-Caffo
    save.AC<-ifelse(test = pi1-pi2 > lower.AC,
                    yes = ifelse(test = pi1-pi2 < upper.AC, yes = 1, no = 0), no = 0)
    AC<-sum(save.AC*pmf)
  
    save.true.conf[counter,]<-c(pi1, wald, AC)
    counter<-counter+1
  }
  
  
  # Plot
  x11(width = 7, height = 6, pointsize = 12)
  # pdf(file = "c:\\figures\\Figure1.4color.pdf", width = 7, height = 6, colormodel = "cmyk")   # Create plot for book
  plot(x = save.true.conf[,1], y = save.true.conf[,2], xlab = expression(pi[1]),
    ylab = "True confidence level", type = "l", ylim = c(0.85,1), lty = "solid", col = "blue")
  lines(x = save.true.conf[,1], y = save.true.conf[,3], lty = "dashed", col = "red")
  abline(h = 1-alpha, lty = "dotted")
  legend(x = 0.1, y = 0.88, legend = c("Wald", "Agresti-Caffo"), lty = c("solid", "dashed"),
    bty = "n", col = c("blue", "red"))
  # dev.off()  # Create plot for book

  # Black-and-white version of plot
  # pdf(file = "c:\\figures\\Figure1.4BW.pdf", width = 7, height = 6, colormodel = "cmyk")   # Create plot for book
  plot(x = save.true.conf[,1], y = save.true.conf[,2], xlab = expression(pi[1]),
    ylab = "True confidence level", type = "l", ylim = c(0.85,1), lty = "solid", col = "black")
  lines(x = save.true.conf[,1], y = save.true.conf[,3], lty = "dashed", col = "black")
  abline(h = 1-alpha, lty = "dotted")
  legend(x = 0.1, y = 0.88, legend = c("Wald", "Agresti-Caffo"), lty = c("solid", "dashed"),
    bty = "n", col = c("black", "black"))
  # dev.off()  # Create plot for book



###########################################################################
# 3D plot of the true confidence level
#   NOTE: This code can take a significant amount of time to run. We recommend
#         using the test cases of 0.1 to 0.9 by 0.1 to obtain an estimate of how
#         long it should run for the 0.001 to 0.999 by 0.0005 case. 

  # Find start time
  start.time<-proc.time()

  pi1seq<-seq(from = 0.001, to = 0.999, by = 0.0025)  # using a smaller "by" argument value can lead to slow rendering when rotating plots
  # pi1seq<-seq(from = 0.1, to = 0.9, by = 0.1)  # Testing

  pi2seq<-seq(from = 0.001, to = 0.999, by = 0.0025)
  # pi2seq<-seq(from = 0.1, to = 0.9, by = 0.1)  # Testing


  # Save true confidence levels in a matrix
  save.true.conf<-matrix(data = NA, nrow = length(pi1seq)*length(pi2seq), ncol = 4)

  # Create counter for the loop
  counter<-1

  # All possible combinations of w1 and w2
  w.all<-expand.grid(w1 = 0:n1, w2 = 0:n2)


  # All possible combinations of pi^_1 and pi^_2
  pi.hat1<-0:n1/n1
  pi.hat2<-0:n2/n2
  pi.hat.all<-expand.grid(pi.hat1 = pi.hat1, pi.hat2 = pi.hat2)
  
  # Wald
  lower<-pi.hat.all[,1] - pi.hat.all[,2] - qnorm(p = 1-alpha/2) * 
         sqrt(pi.hat.all[,1]*(1-pi.hat.all[,1]) / n1 + pi.hat.all[,2]*(1-pi.hat.all[,2]) / n2)
  upper<-pi.hat.all[,1] - pi.hat.all[,2] + qnorm(p = 1-alpha/2) * 
         sqrt(pi.hat.all[,1]*(1-pi.hat.all[,1]) / n1 + pi.hat.all[,2]*(1-pi.hat.all[,2]) / n2)

  # Agresti-Caffo
  pi1tilde<-(0:n1+1)/(n1+2)
  pi2tilde<-(0:n2+1)/(n2+2)
  pi.all.tilde<-expand.grid(pi1tilde = pi1tilde, pi2tilde = pi2tilde)
  lower.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] - qnorm(p = 1-alpha/2) *
            sqrt(pi.all.tilde[,1]*(1-pi.all.tilde[,1]) / (n1+2) +
              pi.all.tilde[,2]*(1-pi.all.tilde[,2]) / (n2+2))
  upper.AC<-pi.all.tilde[,1] - pi.all.tilde[,2] + qnorm(p = 1-alpha/2) *
            sqrt(pi.all.tilde[,1]*(1-pi.all.tilde[,1]) / (n1+2) +
              pi.all.tilde[,2]*(1-pi.all.tilde[,2]) / (n2+2))


  # Loop over each pi1 and pi2 that the true confidence level is calculated on
  for(pi1 in pi1seq) {
    for(pi2 in pi2seq) {

      # Find joint probability for w1 and w2
      prob.w1<-dbinom(x = 0:n1, size = n1, prob = pi1)
      prob.w2<-dbinom(x = 0:n2, size = n2, prob = pi2)
      prob.all<-expand.grid(prob.w1 = prob.w1, prob.w2 = prob.w2)
      pmf<-prob.all$prob.w1*prob.all$prob.w2
   
      # Wald
      save<-ifelse(test = pi1-pi2 > lower,
                   yes = ifelse(test = pi1-pi2 < upper, yes = 1, no = 0), no = 0)
      wald<-sum(save*pmf)

      # Agresti-Caffo
      save.AC<-ifelse(test = pi1-pi2 > lower.AC,
                      yes = ifelse(test = pi1-pi2 < upper.AC, yes = 1, no = 0), no = 0)
      AC<-sum(save.AC*pmf)
  
      save.true.conf[counter,]<-c(pi1, pi2, wald, AC)
      counter<-counter+1
    }
    print(pi1)
  }
  
  # Find end time and total time elapsed
  end.time<-proc.time()
  save.time<-end.time-start.time
  cat("\n Number of minutes running:", save.time[3]/60, "\n \n")

  # Write file out with results to save for later (if needed)
  # write.table(x = save.true.conf, file = "c:\\chris\\save.true.conf.txt", quote = FALSE, row.names = FALSE)
  # save.true.conf<-read.table(file = "c:\\chris\\save.true.conf.txt", header = TRUE)
  
  # 3D plot package
  library(rgl) 
  
  # Wald plot with plane at 0.95
  open3d()
  persp3d(x = pi1seq, y = pi2seq, z = save.true.conf[,3], xlim = c(0,1), ylim =
      c(0,1), zlim = c(0.85, 1), xlab = "pi1", ylab = "pi2", zlab = "True confidence level", ticktype = "detailed", col="red")
  # grid3d(side = c("x-", "y-", "z"), col = "lightgray")
  true.conf<-data.frame(x = c(0,0,0.1,0.1), y = c(0,0.1,0,0.1), z = c(0.95, 
      0.95, 0.95, 0.95))
  persp3d(x = c(0,1), y = c(0,1), z = matrix(data = c(0.95,0.95, 0.95, 
      0.95), nrow = 2, ncol = 2), add = TRUE, col = "green")

  # AC plot with plane at 0.95
  open3d()
  persp3d(x = pi1seq, y = pi2seq, z = save.true.conf[,4], xlim = c(0,1), ylim =
      c(0,1), zlim = c(0.85, 1), xlab = "pi1", ylab = "pi2", zlab = "True confidence level",
      aspects = 10, ticktype = "detailed", col="red")
  # grid3d(side = c("x-", "y-", "z"), col = "lightgray")
  true.conf<-data.frame(x = c(0,0,0.1,0.1), y = c(0,0.1,0,0.1), z = c(0.95, 
      0.95, 0.95, 0.95))
  persp3d(x = c(0,1), y = c(0,1), z = matrix(data = c(0.95,0.95, 0.95, 
      0.95), nrow = 2, ncol = 2), add = TRUE, col = "green")


  # The zlim option in persp3d does not fix the axis limits like it should. Below is a fix
  #  to the problem in order to get both the Wald and AC plots on the same scale.
  test.true.conf.wald<-ifelse(test = save.true.conf[,3]<0.85, yes = NA, no = 1)
  save.true.conf.wald2<-test.true.conf.wald*save.true.conf[,3]  # Put NA's in vector if true confidence level is < 0.85
  open3d()
  persp3d(x = pi1seq, y = pi2seq, z = save.true.conf.wald2, xlim = c(0,1), ylim =
      c(0,1), xlab = "pi1", ylab = "pi2",  zlim = c(0.85, 1),
      zlab = "True confidence level", ticktype = "detailed", col="red")
  true.conf<-data.frame(x = c(0,0,0.1,0.1), y = c(0,0.1,0,0.1), z = c(0.95, 
      0.95, 0.95, 0.95))
  persp3d(x = c(0,1), y = c(0,1), z = matrix(data = c(0.95,0.95, 0.95, 
      0.95), nrow = 2, ncol = 2), add = TRUE, col = "green")
 
  # Note that there is not a par(mfrow = c(1,2)) type of option yet in the rgl package, so I used a
  #  graphics editor to get the two plots side-by-side. Also, I used a print screen to obtain
  #  the plot graphics to put into the graphics editor.  

#ConfLevel.R

# Initial settings
    alpha <- 0.05
    pi    <- 0.0735  # 0.0730 better than 0.0735, 0.156 better than 0.157
    n     <- 40
    w     <- 0:n
    pi.hat <- w/n
    pmf    <- dbinom(x=w, size=n, prob=pi)
    var.wald <- pi.hat*(1-pi.hat)/n
    lower <- pi.hat - qnorm(p = 1-alpha/2) * sqrt(var.wald)
    upper <- pi.hat + qnorm(p = 1-alpha/2) * sqrt(var.wald)
    save <- ifelse(test = pi>lower,
                   yes=ifelse(test = pi<upper, yes=1, no=0), no=0)
    sum(save*pmf)
    data.frame(w, pi.hat, round(data.frame(pmf, lower, upper), 4), save)[1:13,]
  # Following is only for pi = 0.157
    sum(dbinom(x=4:11, size=n, prob=pi))

  # Estimated true confidence level
    numb.bin.samples <- 1000
    set.seed(4516)
    w <- rbinom(n = numb.bin.samples, size = n, prob = pi)
    (counts <- table(x = w))
    sum(counts[4:11])/numb.bin.samples
    pi.hat <- w/n
    pi.hat[1:10]
    var.wald <- pi.hat*(1-pi.hat)/n
    lower <- pi.hat - qnorm(p = 1-alpha/2) * sqrt(var.wald)
    upper <- pi.hat + qnorm(p = 1-alpha/2) * sqrt(var.wald)
    data.frame(w, pi.hat, lower, upper)[1:10,]
    save <- ifelse(test=pi>lower,
                   yes = ifelse(test = pi<upper, yes = 1, no = 0), no = 0)
    save[1:10]
    mean(save)
    true.conf <- mean(save)
    cat("An estimate of the true confidence level is:",
        round(true.conf,4), "\n")

    library(package = binom)
    binom.confint(x=sum(save), n=numb.bin.samples,
                  conf.level = 1-alpha, methods = "wilson")

  # Compare the two ways
  # Simulate same samples again
    set.seed(4516)
    w <- rbinom(n = numb.bin.samples, size = n, prob = pi)
    table(w)
    prop.w <- table(w)/numb.bin.samples
    obs.w <- as.integer(names(table(w)))
    binom.prob <- round(dbinom(x = obs.w, size = n, prob = pi),4)
    data.frame(w=obs.w, obs.prop=prop.w, binom.prob=binom.prob)
    sum(prop.w)
    sum(binom.prob)
  # Note: not equal to 1 because some possible values of w were not observed


#ConfLevel4Intervals.R

  # Initial settings and calculations
    alpha <- 0.05
    n <- 40
    w <- 0:n
    pi.hat <- w/n
    p.tilde <- (w + qnorm(p = 1-alpha/2)^2 /2) / (n+qnorm(1-alpha/2)^2)

  # Wald
    var.wald   <- pi.hat*(1-pi.hat)/n
    lower.wald <- pi.hat - qnorm(p = 1-alpha/2) * sqrt(var.wald)
    upper.wald <- pi.hat + qnorm(p = 1-alpha/2) * sqrt(var.wald)

  # Agresti-Coull
    lower.AC<-p.tilde - qnorm(p=1-alpha/2)*sqrt(p.tilde*(1-p.tilde) /
                                                  (n+qnorm(1-alpha/2)^2))
    upper.AC<-p.tilde + qnorm(p=1-alpha/2)*sqrt(p.tilde*(1-p.tilde) /
                                                  (n+qnorm(1-alpha/2)^2))

  # Wilson
    lower.wilson <- p.tilde -
      qnorm(p=1-alpha/2)*sqrt(n) / (n+qnorm(1-alpha/2)^2) *
      sqrt(pi.hat*(1-pi.hat) + qnorm(1-alpha/2)^2/(4*n))
    upper.wilson <- p.tilde +
      qnorm(p=1-alpha/2)*sqrt(n) / (n+qnorm(1-alpha/2)^2) *
      sqrt(pi.hat*(1-pi.hat) + qnorm(1-alpha/2)^2/(4*n))

  # Clopper-Pearson - Little more complicated due to the y = 0 and n cases
    lower.CP <- numeric(n+1) # Initializes vector to save the lower bounds
    upper.CP <- numeric(n+1) # Initializes vector to save the upper bounds
    # y = 0
    w0 <- 0  # Set here for emphasis
    lower.CP[1] <- 0
    upper.CP[1] <- qbeta(p = 1-alpha/2, shape1 = w0+1, shape2 = n-w0)
    # y = n
    wn <- n  # Set here for emphasis
    lower.CP[n+1] <- qbeta(p = alpha/2, shape1 = wn, shape2 = n-wn+1)
    upper.CP[n+1] <- 1
    # y = 1, ..., n-1
    w.new<-1:(n-1)
    lower.CP[2:n]<-qbeta(p = alpha/2,   shape1=w.new,   shape2=n-w.new+1)
    upper.CP[2:n]<-qbeta(p = 1-alpha/2, shape1=w.new+1, shape2=n-w.new)

  # All pi's
    pi.seq <- seq(from=0.001, to=0.999, by=0.0005)
  # pi.seq<-0.16 #Testing
  # pi.seq<-seq(from = 0.1, to = 0.9, by = 0.1) #Testing
  # Save true confidence levels in a matrix
    save.true.conf <- matrix(data=NA, nrow=length(pi.seq), ncol=5)
  # Create counter for the loop
    counter <- 1
  # Loop over each pi that the true confidence level is calculated on
    for(pi in pi.seq) {
      pmf <- dbinom(x=w, size=n, prob=pi)
      save.wald <- ifelse(test=pi > lower.wald,
                          yes= ifelse(test=pi<upper.wald, yes=1, no=0),
                          no=0)
      wald <- sum(save.wald*pmf)
      save.AC <- ifelse(test=pi > lower.AC,
                        yes=ifelse(test=pi<upper.AC, yes=1, no=0),
                        no=0)
      AC <- sum(save.AC*pmf)
      save.wilson <- ifelse(test=pi > lower.wilson,
                            yes=ifelse(test=pi<upper.wilson, yes=1, no=0),
                            no=0)
      wilson <- sum(save.wilson*pmf)
      save.CP <- ifelse(test=pi > lower.CP,
                        yes=ifelse(test=pi<upper.CP, yes=1, no=0),
                        no=0)
      CP <- sum(save.CP*pmf)
      save.true.conf[counter,] <- c(pi, wald, AC, wilson, CP)
      counter <- counter+1
    }
  # Plots
    #x11(width = 7, height = 6, pointsize = 12)
    par(mfrow = c(2,2))
    plot(x=save.true.conf[,1], y=save.true.conf[,2],
         main="Wald", xlab=expression(pi),
         ylab="True confidence level", type="l", ylim=c(0.85,1))
    abline(h = 1-alpha, lty="dotted", lwd=3)
    plot(x=save.true.conf[,1], y=save.true.conf[,3],
         main="Agresti-Coull", xlab=expression(pi),
         ylab="True confidence level", type="l", ylim = c(0.85,1))
    abline(h = 1-alpha, lty="dotted", lwd=3)
    plot(x=save.true.conf[,1], y=save.true.conf[,4],
         main="Wilson", xlab=expression(pi),
         ylab="True confidence level", type="l", ylim = c(0.85,1))
    abline(h = 1-alpha, lty="dotted", lwd=3)
    plot(x=save.true.conf[,1], y=save.true.conf[,5],
         main="Clopper-Pearson", xlab=expression(pi),
         ylab="True confidence level", type="l", ylim = c(0.85,1))
    abline(h = 1-alpha, lty="dotted", lwd=3)

  # pi = 0.157
    save.true.conf[save.true.conf[,1]==0.157,]
  # AC and Wilson same true conf. levels at pi=0.157: doesn't always happen
    sum(save.true.conf[,3] != save.true.conf[,4])  # Number of differences
    length(pi.seq)  # Number of true confidence levels

  # Alternatively, these true confidence levels can be found using the
  # binom.coverage() function in the binom package
    library(binom)
    binom.coverage(p=0.16, n=n, conf.level=0.95,
                   method=c("asymptotic","agresti-coull","wilson","exact"))
  # The above produces a warning message, but still works

  # Alternatively, plots of the true confidence levels can be found
  # from the binom.plot() function in the binom package. Note that
  # the method argument uses a function name that calculates the
  # individual confidence interval. Also, the y-axis is the true
  # confidence level and the x-axis is pi. The lattice package is
  # used for plotting so the usual arguments (like xlab = ) for
  # basic R plots do not always work.
    binom.plot(n=40, method=binom.asymp,         np=500, conf.level=0.95)
    binom.plot(n=40, method=binom.agresti.coull, np=500, conf.level=0.95)
    binom.plot(n=40, method=binom.wilson,        np=500, conf.level=0.95)
    binom.plot(n=40, method=binom.exact,         np=500, conf.level=0.95)

#Bird.R

# Create contingency table - notice the data is entered by columns
    (c.table <- array(data=c(251, 48, 34, 5), dim=c(2,2),
                      dimnames=list(First =c("made", "missed"),
                                   Second=c("made", "missed"))))
    list(First = c("made", "missed"), Second = c("made", "missed"))
    c.table[1,1]      # w1
    c.table[1,]       # w1 and n1-w1
    sum(c.table[1,])  # n1
    rowSums(c.table)  # n1 and n2
  # Find the estimated pi^j
    (pi.hat.table <- c.table/rowSums(c.table))
    sum(pi.hat.table[1,])
    prop.table(c.table, 1)
    addmargins(prop.table(c.table, 1), 2)

  # Another way to create a contingency table
    (c.table2 <- array(data=c(251, 48, 34, 5), dim=c(2,2),
                       dimnames=list(c("first made", "first missed"),
                                     c("second made", "second missed"))))

  # What if the data did not already come in a contingency table format?
  # Create "raw" data
    miss.miss <- matrix(rep(c("missed", "missed"), 5),   5, 2, byrow=T)
    miss.make <- matrix(rep(c("missed", "made"),  48),  48, 2, byrow=T)
    make.miss <- matrix(rep(c("made", "missed"),  34),  34, 2, byrow=T)
    make.make <- matrix(rep(c("made", "made"),   251), 251, 2, byrow=T)
  # Put "raw" data into one data.frame
    all.data <- rbind(miss.miss, miss.make, make.miss, make.make)
    all.data2 <- data.frame(all.data)
  # Gives new names to columns
    names(all.data2) <- c("first", "second")
  # Rearrange rows to "simulate" how the data may have been observed as
    set.seed(9212)
    all.data2 <- all.data2[sample(x=1:nrow(all.data2), replace=FALSE),]
    row.names(all.data2) <- NULL  # Remove original row numbers
    head(all.data2)
    tail(all.data2)

  # Find contingency table two different ways
    (bird.table1 <- table(all.data2$first, all.data2$second))
    bird.table1[1,1]  # w1
    (bird.table2 <- xtabs(formula = ~ first + second, data=all.data2) )
    bird.table2[1,1]  # w1
    bird.table2/rowSums(bird.table2)
  # summary(bird.table2)  # Gives Pearson chi-square test for independence

  # Confidence interval for difference of two probabilities
    alpha   <- 0.05
    pi.hat1 <- pi.hat.table[1,1]
    pi.hat2 <- pi.hat.table[2,1]
  # Wald
    var.wald <- pi.hat1*(1-pi.hat1)/sum(c.table[1,]) +
                pi.hat2*(1-pi.hat2)/sum(c.table[2,])
    pi.hat1-pi.hat2 + qnorm(p=c(alpha/2, 1-alpha/2))*sqrt(var.wald)
  # Agresti-Caffo
    pi.tilde1 <- (c.table[1,1]+1)/(sum(c.table[1,])+2)
    pi.tilde2 <- (c.table[2,1]+1)/(sum(c.table[2,])+2)
    var.AC <- pi.tilde1*(1-pi.tilde1)/(sum(c.table[1,])+2) +
              pi.tilde2*(1-pi.tilde2)/(sum(c.table[2,])+2)
    pi.tilde1 - pi.tilde2 + qnorm(p = c(alpha/2, 1-alpha/2))*sqrt(var.AC)

  # Note: Each interval limit could be calculated one at a time as well.
  # For example, the Wald interval is
    lower <- pi.hat1 - pi.hat2 - qnorm(p = 1-alpha/2) *
             sqrt(pi.hat1*(1-pi.hat1) / sum(c.table[1,]) +
                  pi.hat2*(1-pi.hat2) / sum(c.table[2,]))
    upper <- pi.hat1 - pi.hat2 + qnorm(p = 1-alpha/2) *
             sqrt(pi.hat1*(1-pi.hat1) / sum(c.table[1,]) +
                  pi.hat2*(1-pi.hat2) / sum(c.table[2,]))
    data.frame(lower, upper)

  # Other ways to get the C.I.s
  # We could also avoid using tables and obtain the Wald interval
    w1 <- 251
    n1 <- 285
    w2 <-  48
    n2 <-  53
    alpha <- 0.05
    pi.hat1 <- w1/n1
    pi.hat2 <- w2/n2
    var.wald <- pi.hat1*(1-pi.hat1) / n1 +  pi.hat2*(1-pi.hat2) / n2
    pi.hat1 - pi.hat2 + qnorm(p = c(alpha/2, 1-alpha/2)) *  sqrt(var.wald)
  # C.I. and also hypothesis test for Ho: pi_1|1 - pi_1|2
    prop.test(x=c.table[,1], n=rowSums(c.table),
              conf.level=0.95, correct = FALSE)
    prop.test(x=c.table, conf.level=0.95, correct=FALSE)
  # Above two prop.test's give the same results.

  # Wald statistic
    (Z.W <- (pi.hat1-pi.hat2)/sqrt( pi.hat1*(1-pi.hat1)/sum(c.table[1,]) +
                                    pi.hat2*(1-pi.hat2)/sum(c.table[2,]) ))
    Z.W^2
  # Incorporate null hypothesis into variance in the denominator
    pi.hat.Ho <- sum(c.table[,1])/sum(c.table)
    Z.0 <- (pi.hat1-pi.hat2)/sqrt( pi.hat.Ho*(1-pi.hat.Ho)*
                              (1/sum(c.table[1,])+1/sum(c.table[2,])) )
    Z.0^2

  # Calculations using the PropCIs package
    library(package = PropCIs)
  # Wald
    wald2ci(x1=c.table[1,1], n1=sum(c.table[1,]),
            x2=c.table[2,1], n2=sum(c.table[2,]),
            conf.level=0.95, adjust="Wald")
  # Agresti-Caffo
    wald2ci(x1 = c.table[1,1], n1 = sum(c.table[1,]),
            x2 = c.table[2,1], n2 = sum(c.table[2,]),
            conf.level = 0.95, adjust = "AC")

  # Hypothesis test for difference of two probabilities
    prop.test(x = c.table, conf.level = 0.95, correct = FALSE)

  # LRT
    pi.bar <- colSums(c.table)[1]/sum(c.table)
    log.Lambda <- c.table[1,1]*log(pi.bar/pi.hat.table[1,1]) +
                  c.table[1,2]*log((1-pi.bar)/(1-pi.hat.table[1,1])) +
                  c.table[2,1]*log(pi.bar/pi.hat.table[2,1]) +
                  c.table[2,2]*log((1-pi.bar)/(1-pi.hat.table[2,1]))
    test.stat <- -2*log.Lambda
    crit.val  <- qchisq(p=0.95, df=1)
    p.val     <- 1-pchisq(q=test.stat, df=1)
    round(data.frame(pi.bar, test.stat, crit.val, p.val, row.names=NULL), 4)
    library(package = vcd)
    assocstats(x = c.table)

  # Other ways to do the hypothesis tests
    chisq.test(x = c.table, correct = FALSE)
    summary(bird.table2)
    summary(as.table(c.table))
    class(bird.table2)
    class(as.table(c.table))
    summary.table(bird.table2)

  # Relative risk
    cat("The sample relative risk is", round(pi.hat1/pi.hat2, 4), "\n \n")
    alpha <- 0.05
    n1 <- sum(c.table[1,])
    n2 <- sum(c.table[2,])
  # Wald confidence interval
    ci <- exp(log(pi.hat1/pi.hat2) + qnorm(p=c(alpha/2, 1-alpha/2)) *
                sqrt((1-pi.hat1)/(n1*pi.hat1) + (1-pi.hat2)/(n2*pi.hat2)))
    round(ci, 4)
    rev(round(1/ci, 4))  # inverted

  # Change contingency table so that we are looking at the "missed" column
    (1-pi.hat1)/(1-pi.hat2)
    exp(log((1-pi.hat1)/(1-pi.hat2)) + qnorm(p = c(alpha/2, 1-alpha/2)) *
    sqrt((pi.hat1)/(n1*(1-pi.hat1)) + (pi.hat2)/(n2*(1-pi.hat2))))

  # OR
    OR.hat <- c.table[1,1]*c.table[2,2] / (c.table[2,1]*c.table[1,2])
    round(OR.hat, 6)
    round(1/OR.hat, 6)
    alpha <- 0.05
    var.log.or <- 1/c.table[1,1] + 1/c.table[1,2] +
                  1/c.table[2,1] + 1/c.table[2,2]
    OR.CI <- exp(log(OR.hat) + qnorm(p = c(alpha/2, 1-alpha/2)) *
             sqrt(var.log.or))
    round(OR.CI, 6)
    rev(round(1/OR.CI, 6))

  # Another way to get the OR
  # Note that the function below automatically adds 0.5 to each cell
  # for variance but not for the odds ratio itself unless there are 0
  # counts (see code in function).
    library(vcd)
    save.OR <- oddsratio(x = c.table, log = TRUE)
    attributes(save.OR)  # names( ) does not work
    summary(save.OR) 
    confint(save.OR, level = 0.95)
    OR.tilde <- (c.table[1,1]+0.5)*(c.table[2,2]+0.5)/
                ((c.table[1,2]+0.5)*(c.table[2,1]+0.5))
    log(OR.tilde)  # Does not match vcd's log(OR^)
    sqrt(1/(c.table[1,1]+0.5) + 1/(c.table[2,2]+0.5) +
         1/(c.table[1,2]+0.5) + 1/(c.table[2,1]+0.5))
    # Above matches vcd's standard error (sqrt(var^(log(OR^)))


#Binomial.R

 # P(W = 1) with n = 5 and pi = 0.6
    dbinom(x=1, size=5, prob=0.6) 

  # The entire pmf
    (pmf <- dbinom(x=0:5, size=5, prob=0.6) )

  # Make the printing look a little nicer
    pmf <- dbinom(x=0:5, size=5, prob=0.6) 
    (save <- data.frame(w=0:5, prob=round(x=pmf, digits=4)))

  # Plot PMF
    # x11(width=6, height=6, pointsize=12)
    plot(x=save$w, y=save$prob, type="h", xlab="w", ylab="P(W=w)",
         main="Plot of a binomial PMF for n=5, pi=0.6", lwd=3)
    #     panel.first=grid(col="gray", lty="dotted"), lwd=3)
    abline(h=0)

  # Alternative plot using the expression() function
    plot(x=save$w, y=save$prob, type="h", xlab="w", ylab="P(W=w)",
         main=expression(paste("Plot of a binomial PMF for ",
                         italic(n)==5, " and ", italic(pi)==0.6)), 
         panel.first=grid(col="gray", lty="dotted"), lwd=3)
    abline(h=0)

  # Simulate observations from a Binomial distribution
    set.seed(4848)
    bin5 <- rbinom(n=1000, size=5, prob=0.6)
    bin5[1:10]
    mean(bin5)
    var(bin5)

  # Frequency distribution
    table(x=bin5) 
    # x11(width=6, height=6, pointsize=12)
  # Relative frequency histogram
    hist(x=bin5, main="Binomial with n=5, pi=0.6, 1000 observations",
         probability=TRUE, ylab="Relative frequency")
  # pdf(file="c:\\figures\\Figure1.2.pdf", width=6, height=6, colormodel="cmyk")
    hist(x=bin5, main="Binomial with n=5, pi=0.6, 1000 observations",
         probability=TRUE, breaks=c(-0.5:5.5), ylab="Relative frequency")
  # dev.off()  # Create plot for book
  # Better plot
    (save.count <- table(bin5))

    barplot(height=save.count,
            names=c("0", "1", "2", "3", "4", "5"),
            main="Binomial with n=5, pi=0.6, 1000 bin. observations",
            xlab="x")



#Beta.R

  # Basic computations
    #x11(width = 6, height = 6, pointsize = 12)
    par(mfrow = c(2,2))
    curve(expr=dbeta(x=x, shape1=1, shape2=1), ylab="f(v)", 
          xlab="v", from=0, to=1, ylim=c(0,5), main="beta(1,1)")
    curve(expr=dbeta(x=x, shape1=1, shape2=5), ylab="f(w)", 
          xlab="v", from=0, to=1, ylim=c(0,5), main="beta(1,5)")
    curve(expr=dbeta(x=x, shape1=5, shape2=1), ylab="f(w)", 
          xlab="v", from=0, to=1, ylim=c(0,5), main="beta(5,1)")
    curve(expr=dbeta(x=x, shape1=5, shape2=5), ylab="f(w)", 
          xlab="v", from=0, to=1, ylim=c(0,5), main="beta(5,5)")
    par(mfrow=c(1,1))



    


