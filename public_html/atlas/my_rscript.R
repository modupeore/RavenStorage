args <- commandArgs(TRUE)
 
L <- as.numeric(args[1])
S <- args[2]
cnames <- c("Ross", "Illinois")
N <- args[3]
mymatrix <- matrix(L, ncol=2, byrow=TRUE, dimnames=list(S,cnames))

L
typeof(L)
#x <- rnorm(N,0,1)
#png(filename="OUTPUT/temp.png", width=500, height=500)
#barplot(mymatrix, beside=TRUE)
#dev.off()

