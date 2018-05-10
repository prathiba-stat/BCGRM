#set working directory on line 2 between quotes
setwd("")
require(runjags)
#dataset 1 - only one teacher/student
#dataset 2 - more than 1 teacher/student
#dataset 3 - all students
#model 1 is GRM for dataset 1 - standard priors
#model 2 is GRM for dataset 2 - standard priors
#model 3 is GRM for dataset 2 - priors from posteriors of model 1 
#model 4 is GRM for dataset 3 - standard priors
#model 5 is modified GRM to include teacher information
#for dataset 3 - standard priors
model_1 <- "model {
  # model for observations
  for (i in 1:N) {
    for (j in 1:L) {
      for (k in 1:(K-1)) {
        logit(p.star[i, j, k]) <- a[j] * (theta[i] - b[j, k])
      }
      p[i, j, 1] <- 1 - p.star[i, j, 1]
      for (k in 2:(K - 1)) {
        p[i, j, k] <- p.star[i, j, k - 1] - p.star[i, j, k]
      }
      p[i, j, K] <- p.star[i, j, K - 1]
      y[i, j] ~ dcat(p[i, j, 1:K])
    }
  }
  
  # theta parameters
  for (i in 1:N) {
    theta[i] ~ dnorm(0, 1)
  }
  
  # item parameters
  for (j in 1:L) {
    a[j] ~ dnorm(0, 1)I(0,)
    
    for (k in 1:(K - 1)) {
      b.star[j, k] ~ dnorm(0, 1)  
    }
    b[j, 1:(K - 1)] <- sort(b.star[j, 1:(K - 1)])
  }
}"


data <- read.csv("all-students-not-reversed.csv")[,2:9]
data <- data[,c(1, 2, 6, 7, 8, 3, 4, 5)]

#obtain only the first instance of the ratings. The rest are duplicates
newdata <- data[row.names(unique(data[c("teacherID", "studentID")])),]

#index students who have more than 2 teacher ratings
double.students <- which(colSums(table(newdata$teacherID, newdata$studentID))>1)
#extract students indexed in the previous step
mult.data <- newdata[newdata[,2] %in% double.students, ]

#index students with only one rating
single.students <- which(colSums(table(newdata$teacherID, newdata$studentID))==1)
#extract students indexed in the previous step
single.data <- newdata[newdata[,2] %in% single.students, ]

model <- model_1
y <- single.data
teacher1 <- y[ , 1]
student1 <- y[ , 2]
n.students1 <- length(unique(student1))
n.teachers1 <- length(unique(teacher1))
y <- y[ , -c(1, 2)]  # drop teachers, students

N1 <- nrow(y)
L1 <- ncol(y)
K <- 5
# recode negatively worded items
for (j in 1:3) {
  y[ , j] <- K - y[ , j] + 1
}
y <- as.matrix(y)

data1 <- list(
  y = y, N = N1, L = L1, K = K
)
results1 <- run.jags(
  model = model_1,
  data = data1,
  monitor = c("a", "b", "theta", "DIC"),
  n.chains = 2,
  burnin = 1000,
  adapt = 1000, 
  sample = 10000,
  #method = "rjparallel",
  summarise = TRUE,
  inits = function() {
    list(
      "a" = exp(rlnorm(data$L, 0, 0.1)),
      "b.star" = matrix(rnorm(data1$L * (data1$K - 1), 0, 1), 
                        nrow=data1$L, ncol=data1$K-1),
      "theta" = rnorm(N1, 0, 1)
    )
  }
)

results1$draws <- combine.mcmc(results1$mcmc)
s1 <- summary(results1)
DIC <- c(results1$dic, rep(NA, 3))
s1 <- rbind(s1, DIC)
write.csv(s1, "data1-estimates.csv")
##################END OF MODEL 1

#obtain priors for model 3
prior.a <- s1[1:6,4]
prior.b <- matrix(s1[7:30, 4], 6, 4)

################BEGIN MODEL 2
#select students with more than one teacher rating
y <- mult.data
n.teachers <- length(unique(y[,1]))
new.tids <- seq(1:n.teachers)
j <- 0
for (i in unique(y[,1])){
  j <- j + 1
  y[y[,1]==i,1] <- new.tids[j]
}

n.students <- length(unique(y[,2]))
new.sids <- seq(1:n.students)
k <- 0
for (i in unique(y[,2])){
  k <- k + 1
  y[y[,2]==i, 2] <- new.sids[k]
}
teacher <- y[ , 1]
student <- y[ , 2]

y <- y[ , -c(1, 2)]  # drop teachers, students

N <- nrow(y)
L <- ncol(y)
K <- 5
# recode negatively worded items
for (j in 1:3) {
  y[ , j] <- K - y[ , j] + 1
}
y <- as.matrix(y)

model_2 <- "model {
  # model for observations
  for (i in 1:N) {
    for (j in 1:L) {
      for (k in 1:(K-1)) {
        logit(p.star[i, j, k]) <- a[j] * (theta[i] - b[j, k])
      }
      p[i, j, 1] <- 1 - p.star[i, j, 1]
      for (k in 2:(K - 1)) {
        p[i, j, k] <- p.star[i, j, k - 1] - p.star[i, j, k]
      }
      p[i, j, K] <- p.star[i, j, K - 1]
      y[i, j] ~ dcat(p[i, j, 1:K])
    }
  }
  
  # student & teacher parameters
  for (i in 1:N) {
    theta[i] <- theta.s[student[i]] + theta.t[teacher[i]]
  }
  
  for (s in 1:n.students) {
    theta.s[s] ~ dnorm(0, 1)
  }
  
  for (t in 1:n.teachers) {
    theta.t[t] ~ dnorm(0, 1)
  }
  
  # item parameters
  for (j in 1:L) {
    a[j] ~ dnorm(0, 1)I(0,)
    
    for (k in 1:(K - 1)) {
      b.star[j, k] ~ dnorm(0, 1)  
    }
    b[j, 1:(K - 1)] <- sort(b.star[j, 1:(K - 1)])
  }
}"

data2 <- list(
  y = y, N = N, L = L, K = K,
  student = student, n.students = n.students,
  teacher = teacher, n.teachers = n.teachers
)
results2 <- run.jags(
  model = model_2,
  data = data2,
  monitor = c("a", "b", "theta.s", "theta.t", "dic"),
  n.chains = 2,
  burnin = 1000,
  adapt = 1000, 
  sample = 10000,
  # method = "rjparallel",
  summarise = TRUE,
  inits = function() {
    list(
      "a" = exp(rlnorm(data$L, 0, 0.1)),
      "b.star" = matrix(rnorm(data$L * (data$K - 1), 0, 1), nrow=data$L, ncol=data$K-1),
      "theta.s" = rnorm(data$n.students, 0, 1),
      "theta.t" = rnorm(data$n.teachers, 0, 1)
    )
  }
)

results2$draws <- combine.mcmc(results2$mcmc)
s2 <- summary(results2)
DIC2 <- c(results2$dic, rep(NA, 3))
s2 <- rbind(s2, DIC2)
write.csv(s2, "model2-estimates.csv")

##################MODEL 3
model_3 <- "model {
  # model for observations
  for (i in 1:N) {
    for (j in 1:L) {
      for (k in 1:(K-1)) {
        logit(p.star[i, j, k]) <- a[j] * (theta[i] - b[j, k])
      }
      p[i, j, 1] <- 1 - p.star[i, j, 1]
      for (k in 2:(K - 1)) {
        p[i, j, k] <- p.star[i, j, k - 1] - p.star[i, j, k]
      }
      p[i, j, K] <- p.star[i, j, K - 1]
      y[i, j] ~ dcat(p[i, j, 1:K])
    }
  }
  
  # student & teacher parameters
  for (i in 1:N) {
    theta[i] <- theta.s[student[i]] + theta.t[teacher[i]]
  }
  
  for (s in 1:n.students) {
    theta.s[s] ~ dnorm(0, 1)
  }
  
  for (t in 1:n.teachers) {
    theta.t[t] ~ dnorm(0, 1)
  }
  
  # item parameters
  for (j in 1:L) {
    a[j] ~ dnorm(prior.a[j], 1)I(0,)
    
    for (k in 1:(K - 1)) {
      b.star[j, k] ~ dnorm(prior.b[j, k], 1)  
    }
    b[j, 1:(K - 1)] <- sort(b.star[j, 1:(K - 1)])
  }
}"
data3 <- list(
  y = y, N = N, L = L, K = K,
  student = student, n.students = n.students,
  teacher = teacher, n.teachers = n.teachers,
  prior.a = prior.a, prior.b = prior.b
)
results3 <- run.jags(
  model = model_3,
  data = data3,
  monitor = c("a", "b", "theta.s", "theta.t", "dic"),
  n.chains = 2,
  burnin = 1000,
  adapt = 1000, 
  sample = 10000,
  #method = "rjparallel",
  summarise = TRUE,
  inits = function() {
    list(
      "a" = exp(rlnorm(data$L, 0, 0.1)),
      "b.star" = matrix(rnorm(data$L * (data$K - 1), 0, 1), 
                        nrow=data$L, ncol=data$K-1),
      "theta.s" = rnorm(data$n.students, 0, 1),
      "theta.t" = rnorm(data$n.teachers, 0, 1)
    )
  }
)

results3$draws <- combine.mcmc(results3$mcmc)
s3 <- summary(results3)
DIC3 <- c(results3$dic, rep(NA, 3))
s3 <- rbind(s3, DIC3)
# student / teacher means
theta.s <- s3[which(startsWith(rownames(s3), "theta.s")), "Mean"]
theta.t <- s3[which(startsWith(rownames(s3), "theta.t")), "Mean"]

hist(theta.s)
hist(theta.t)

summary(theta.s)
summary(theta.t)

write.csv(s3, "model3-estimates.csv")
###############END of MODEL 3

#Gather  all data for model 4 which is basically model 1
#but estimated for all data, students with 1 or more teachers

y <- newdata
n.teachers <- length(unique(y[,1]))
new.tids <- seq(1:n.teachers)
j <- 0
for (i in unique(y[,1])){
  j <- j + 1
  y[y[,1]==i,1] <- new.tids[j]
}

n.students <- length(unique(y[,2]))
new.sids <- seq(1:n.students)
k <- 0
for (i in unique(y[,2])){
  k <- k + 1
  y[y[,2]==i, 2] <- new.sids[k]
}
teacher <- y[ , 1]
student <- y[ , 2]

y <- y[ , -c(1, 2)]  # drop teachers, students

N <- nrow(y)
L <- ncol(y)
K <- 5
# recode negatively worded items
for (j in 1:3) {
  y[ , j] <- K - y[ , j] + 1
}
y <- as.matrix(y)

data4 <- list(
  y = y, N = N, L = L, K = K
)
results4 <- run.jags(
  model = model_1,
  data = data4,
  monitor = c("a", "b", "theta", "dic"),
  n.chains = 2,
  burnin = 1000,
  adapt = 1000, 
  sample = 10000,
  # method = "rjparallel",
  summarise = TRUE,
  inits = function() {
    list(
      "a" = exp(rlnorm(data$L, 0, 0.1)),
      "b.star" = matrix(rnorm(data$L * (data$K - 1), 0, 1), nrow=data$L, ncol=data$K-1),
      "theta" = rnorm(N, 0, 1)
    )
  }
)

results4$draws <- combine.mcmc(results4$mcmc)
s4 <- summary(results4)
DIC4 <- c(results4$dic, rep(NA, 3))
s4 <- rbind(s4, DIC4)
write.csv(s4, "model4-estimates.csv")
#########END OF MODEL 4
#BEGIN MODEL 5 WHICH MODELS TEACHER INFORMATION
#Model 5
model_5 <- "model {
  # model for observations
  for (i in 1:N) {
    for (j in 1:L) {
      for (k in 1:(K-1)) {
        logit(p.star[i, j, k]) <- a[j] * (theta[i] - b[j, k])
      }
      p[i, j, 1] <- 1 - p.star[i, j, 1]
      for (k in 2:(K - 1)) {
        p[i, j, k] <- p.star[i, j, k - 1] - p.star[i, j, k]
      }
      p[i, j, K] <- p.star[i, j, K - 1]
      y[i, j] ~ dcat(p[i, j, 1:K])
    }
  }
  
  # student & teacher parameters
  for (i in 1:N) {
    theta[i] ~ dnorm((theta.s[student[i]] + theta.t[teacher[i]]), 1)
  }
  
  for (s in 1:n.students) {
    theta.s[s] ~ dnorm(0, 1)
  }
  
  for (t in 1:n.teachers) {
    theta.t[t] ~ dnorm(0, 1)
  }
  
  # item parameters
  for (j in 1:L) {
    a[j] ~ dnorm(0, 1)I(0,)
    
    for (k in 1:(K - 1)) {
      b.star[j, k] ~ dnorm(0, 1)  
    }
    b[j, 1:(K - 1)] <- sort(b.star[j, 1:(K - 1)])
  }
}"

data5 <- list(
  y = y, N = N, L = L, K = K,
  student = student, n.students = n.students,
  teacher = teacher, n.teachers = n.teachers
)
results5 <- run.jags(
  model = model_5,
  data = data5,
  monitor = c("a", "b", "theta.s", "theta.t", "dic"),
  n.chains = 2,
  burnin = 1000,
  adapt = 1000, 
  sample = 10000,
  # method = "rjparallel",
  summarise = TRUE,
  inits = function() {
    list(
      "a" = exp(rlnorm(data$L, 0, 0.1)),
      "b.star" = matrix(rnorm(data$L * (data$K - 1), 0, 1), nrow=data$L, ncol=data$K-1),
      "theta.s" = rnorm(data$n.students, 0, 1),
      "theta.t" = rnorm(data$n.teachers, 0, 1)
    )
  }
)

results5$draws <- combine.mcmc(results5$mcmc)
s5 <- summary(results5)
DIC5 <- c(results5$dic, rep(NA, 3))
s5 <- rbind(s5, DIC5)
write.csv(s5, "model5-estimates.csv")
