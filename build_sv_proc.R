# y_t = 

phi = 0.95
eta =  0.3
ups = 0.4

num = 1000

h = vector("numeric", num)
h[1] = rnorm(1, 0, sqrt(eta^2 / (1 - phi^2)))

for(i in 2:num) {
	h[i] = ups + phi * h[i - 1] + rnorm(1, 0, eta^2)
}
plot(h, type="l")

eps = rnorm(num)
y = eps * exp(h / 2)

plot(cumsum(y), type="l")
