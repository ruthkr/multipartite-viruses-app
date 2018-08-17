# # Third Fixed Point
# R2.star3 <- tail(data$V3,1)
# s.star3 <- tail(data$V5,1)
#
# a.for.p3 <- sigma1*kappa1/omega + kappa1*(1-R2.star3)
# b.for.p3 <- -(epsilon1*s.star3 + gamma1 + kappa1*(1-s.star3)*(1-R2.star3))
# c.for.p3 <- (epsilon1*s.star3 + gamma1)*(1-s.star3)
#
# p.star3.plus <- (-b.for.p3 + sqrt(b.for.p3^2 - 4*a.for.p3*c.for.p3)) / (2*a.for.p3)
# p.star3.minus <- (-b.for.p3 - sqrt(b.for.p3^2 - 4*a.for.p3*c.for.p3)) / (2*a.for.p3)
#
# R1.star3.plus <- 1 - R2.star3 - epsilon1*s.star3 / (kappa1*p.star3.plus) - gamma1 / (kappa1*p.star3.plus)
# R1.star3.minus <- 1 - R2.star3 - epsilon1*s.star3 / (kappa1*p.star3.minus) - gamma1 / (kappa1*p.star3.minus)
#
# v1.star3.plus <- epsilon1*R1.star3.plus*s.star3 / delta1
# v1.star3.minus <- epsilon1*R1.star3.minus*s.star3 / delta1
#
# v2.star3 <- epsilon2*R2.star3*s.star3 / delta2

# ################################################################################################################
# Third Fixed Point
R1.star3 <- tail(data$V2,1)

x <- sigma1*(kappa2*omega*R1.star3*(1 - R1.star3) - gamma2*(sigma1 + omega*R1.star3))
y <- sigma2*(sigma1 + omega*R1.star3)*kappa2*R1.star3 + epsilon2*sigma1*(sigma1 + omega*R1.star3) + sigma1*kappa2*omega*R1.star3*(1 - R1.star3)
s.star3 <- x / y
p.star3 <- omega*R1.star3*(1-s.star3) / (sigma1 + omega*R1.star3)
R2.star3 <- 1 - R1.star3 - epsilon2*s.star3 / (kappa2*p.star3) - gamma2 / (kappa2*p.star3)
v1.star3 <- epsilon1*R1.star3*s.star3 / delta1
v2.star3 <- epsilon2*R2.star3*s.star3 / delta2
################################################################################################################

# Fixed point 3 to try
R2.star3 <- tail(data$V3,1)

a.p3 <- sigma1*kappa1/(sigma2*omega) * (sigma2 + omega*R2.star3) + kappa1 * (1-R2.star3) + epsilon1*omega*R2.star3 / (sigma2+omega*R2.star3)
b.p3 <- -(kappa1*(1-R2.star3) + 2*epsilon1*omega*R2.star3 / (sigma2+omega*R2.star3) + gamma1 )
c.p3 <- gamma1 + epsilon1*omega*R2.star3 / (sigma2+omega*R2.star3)

p.star3.plus <- (-b.p3 + sqrt(b.p3^2 - 4*a.p3*c.p3)) / (2*a.p3)
p.star3.minus <- (-b.p3 - sqrt(b.p3^2 - 4*a.p3*c.p3)) / (2*a.p3)

s.star3.plus <- omega*R2.star3*(1-p.star3.plus) / (sigma2 + omega*R2.star3)
s.star3.minus <- omega*R2.star3*(1-p.star3.minus) / (sigma2 + omega*R2.star3)

R1.star3.plus <- 1 - R2.star3 - epsilon1*s.star3.plus / (kappa1*p.star3.plus) - gamma1 / (kappa1*p.star3.plus)
R1.star3.minus <- 1 - R2.star3 - epsilon1*s.star3.minus / (kappa1*p.star3.minus) - gamma1 / (kappa1*p.star3.minus)

v1.star3.plus <- epsilon1*R1.star3.plus*s.star3.plus / delta1
v1.star3.minus <- epsilon1*R1.star3.minus*s.star3.minus / delta1

v2.star3.plus <- epsilon2*R2.star3*s.star3.plus / delta2
v2.star3.minus <- epsilon2*R2.star3*s.star3.minus / delta2

############################################################

s.star3 <- tail(data$V5,1)

a3 <- (sigma1 / omega + 1 ) * kappa1
b3 <- (sigma1 / omega) * kappa1 * s.star3 - kappa1 * (1 - s.star3) - epsilon1 * s.star3 - gamma1
c3 <- (epsilon1*s.star3 + gamma1) * (1 - s.star3)

p.star3.p <- (-b3 + sqrt(b3^2 - 4*a3*c3)) / (2*a3)
R2.star3.p <- sigma2 * kappa2 * p.star3.p * s.star3 / (omega * (1 - p.star3.p - s.star3))
R1.star3.p <- 1 - R2.star3.p - epsilon1*s.star3 / (kappa1*p.star3.p) - gamma1 / (kappa1*p.star3.p)
v1.star3.p <- epsilon1*R1.star3.p*s.star3 / delta1
v2.star3.p <- epsilon2*R2.star3.p*s.star3 / delta2

p.star3.m <- (-b3 - sqrt(b3^2 - 4*a3*c3)) / (2*a3)
R2.star3.m <- sigma2 * kappa2 * p.star3.m * s.star3 / (omega * (1 - p.star3.m - s.star3))
R1.star3.m <- 1 - R2.star3.m - epsilon1*s.star3 / (kappa1*p.star3.m) - gamma1 / (kappa1*p.star3.m)
v1.star3.m <- epsilon1*R1.star3.m*s.star3 / delta1
v2.star3.m <- epsilon2*R2.star3.m*s.star3 / delta2



