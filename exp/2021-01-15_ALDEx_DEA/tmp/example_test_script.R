library(ALDEx2)
data(selex) #subset only the last 400 features for efficiency
selex.sub <- selex[1:400,]
conds <- c(rep("NS", 7), rep("S", 7))
x.all <- aldex(selex.sub, conds, mc.samples=16, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",ylab="Difference")