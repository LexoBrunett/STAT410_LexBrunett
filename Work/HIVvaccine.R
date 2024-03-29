

  # This program is not referenced in the book, but it is referenced in
  # workshops that the authors give corresponding to the book.

    (c.table <- array(data = c(51, 74, 8146, 8124), dim=c(2,2),
                      dimnames=list(Trt=c("vaccine", "placebo"),
                                    Response=c("HIV", "No HIV"))))

  # Relative risk
    n1 <- sum(c.table[1,])
    n2 <- sum(c.table[2,])
    pi.hat1 <- c.table[1,1] / sum(c.table[1,])
    pi.hat2 <- c.table[2,1] / sum(c.table[2,])
    round(pi.hat1/pi.hat2, 4)
    alpha <- 0.05
  # Wald confidence interval
    var.log.RR <- 1/c.table[1,1] - 1/n1 + 1/c.table[2,1] - 1/n2
    RR.ci <- exp(log(pi.hat1/pi.hat2) + qnorm(p = c(alpha/2, 1-alpha/2)) *
                   sqrt(var.log.RR))
    round(RR.ci, 4)
    rev(round(1/RR.ci, 4))

  # Odds ratio
    OR.hat <- c.table[1,1]*c.table[2,2] / (c.table[2,1]*c.table[1,2])
    round(OR.hat, 4)
    round(1/OR.hat, 4)
    alpha<-0.05
    var.log.or <- 1/c.table[1,1] + 1/c.table[1,2] +
                  1/c.table[2,1] + 1/c.table[2,2]
    OR.CI <- exp(log(OR.hat) +
                   qnorm(p = c(alpha/2, 1-alpha/2)) *sqrt(var.log.or))
    round(OR.CI, 4)
    rev(round(1/OR.CI, 4))

  # Using Epi package
    library(package = Epi)
    twoby2(c.table, alpha=0.05)
    library(package = epitools)
    oddsratio(c.table, conf.level = 0.95, method = "wald")
