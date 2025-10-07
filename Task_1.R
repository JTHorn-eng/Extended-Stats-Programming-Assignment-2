

n <- 1000; h <- rep(NA,100); hmax <- 5
h


rep(1:(n%/%sum(1:hmax)+1), rep(1:5,6), lenght.out=)

# Suppose n = 90 and hmax = 5
# sum(1:hmax) = 15
# so every cycle we would assign 15 people 5 numbers
# eg 122333444455555
# This happens a total of n/(sum(hmax)) = 6 times
# so we would need 5 numbers * 6 cycles = 30 numbers.
rep(1:30, times = sample(1:hmax, 30, replace = TRUE))

# however, what if n is not a multiple of hmax?
sample(rep(1:(n%/%sum(hmax))), 1000)

rep(1:hmax, rep(1:hmax))

# If n is a multiple of sum(hmax) it is easy. for example if hmax=5 and n=15:
rep(1:5,1:5)

# Or hmax=5 and n = 30
rep(1:10,rep(1:5,2))
n = 45; hmax = 6
# more generally, if n is not a multiple of hmax it is tricky because the rep function 
# only allows the same number of repeats

n = 100000; hmax = 5
# is it best to cut off the repeats to n and knowingly be left with more smaller households than I should,
# or is it better to have an equal number of households of all the sizes and just sample 1000 values from them?
h <- sample(rep(rep(1:(hmax*((n%/%sum(hmax))+1)), rep(1:hmax, n%/%sum(hmax)+1)),length.out=n),n, replace=FALSE)
h
h <- sort(h)
h
rep((rep(1:20,rep(1:5,4))),length.out = 50)
