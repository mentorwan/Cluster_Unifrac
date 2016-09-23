load("./otu.RData")

alpha <- c(0,0.5,1)
n <- nrow(otu.tab.rff)

dimname3 <- c(paste("d", alpha, sep="_"), "d_UW", "d_VAW")
unifracs <- array(NA, c(n, n, length(alpha)+2),
                  dimnames=list(rownames(otu.tab), rownames(otu.tab), dimname3))


for (i in 2:n){
  
con <- file(paste("/data/wany/R/unifrac",i,".RData",sep=""),"rb")

load(con)
unifracs[i,,] <- temp

}

for (i in 1:(length(alpha)+2)){
  for (j in 1:n){
    unifracs[j, j, i] <- 0
  }
} 


for (i in 1:(n-1)) {
  for (j in (i+1):n)
  {
    for (k in 1:(length(alpha)+2)){
      unifracs[i,j,k] <- unifracs[j,i,k]
    }
    
  }
}
  

