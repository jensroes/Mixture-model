# Create fake data
mog <- function (n, theta, mu1, mu2, sig1, sig2) {
  y0 <- rlnorm(n, mean=mu1, sd = sig1)
  y1 <- rlnorm(n, mean=mu2, sd = sig2)
  flag <- rbinom(n, size=1, prob=theta)
  y <- y0*(1 - flag) + y1*flag 
}

set.seed(123)
Nsubj = 100 # number of subjects
K = 10 # number of observations per subject per condition
# assumed data were transformed to proportions
# Population parameters

beta_mean = 1000; log(beta_mean)
beta_sd = 2
theta <- c(.25, .35) # mixting proportion for condition 1 and 2
delta_mean = 20; log(delta_mean)
delta_sd = .5
sigma = log(c(2, 3)) # trial-by-trial error

SubjBeta = log(rnorm(Nsubj, beta_mean, beta_sd))
hist(SubjBeta)
SubjDelta = log(rnorm(Nsubj, delta_mean, delta_sd))
hist(SubjDelta)

# iterate over subject to generate data for each one
data = NULL

for(subnum in 1:Nsubj){
  data = rbind(
    data
    , data.frame(
      id = subnum
      , condition = 1
      , value = mog(n = K, 
                    theta = theta[1],
                    mu1 = SubjBeta[subnum],
                    mu2 = SubjBeta[subnum] + SubjDelta[subnum],
                    sig1 = sigma[1],
                    sig2 = sigma[2])
    )
    , data.frame(
      id = subnum
      , condition = 2
      , value = mog(n = K, 
                    theta = theta[2],
                    mu1 = SubjBeta[subnum],
                    mu2 = SubjBeta[subnum] + SubjDelta[subnum],
                    sig1 = sigma[1],
                    sig2 = sigma[2])
    )
  )
}

data <- as_tibble(data)
#print(theta)