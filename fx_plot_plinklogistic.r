#=================================
# This function relies on R/ggplot2 and will generate Manhattan plot and QQ plot 
# from plink association analysis output
# Original source: https://github.com/freeseek/gwaspipeline; free software by G. Genovese
# Used and modified under the terms of the GNU General Public License as published by
# the Free Software Foundation,
#====================================

plot_plinklogistic <- function(fname){
#e.g. fname="xyz.assoc.logistic"

library(ggplot2)

df = read.table(paste(fname, '.tsv',sep=''), as.is=T, head=T)
colnames(df) = c('CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'OR', 'STAT', 'P') 

#=======MA plot
#hg19
chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067,
            159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
            115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
            59128983, 63026520, 48129895, 51305566, 155270560, 59373566)
ticks <- c(0, cumsum(chrlen))[1:23] + chrlen[1:23]/2

df$CHR[df$CHR==25] <- 23
df <- df[df$CHR %in% 1:23 & !is.na(df$P) & df$P>0 & df$P<=1, ]
df$POS <- c(0, cumsum(chrlen))[df$CHR] + df$BP

p <- ggplot(df, aes(x = POS, y = -log10(P), color = as.factor(CHR %% 2))) + geom_point() +
     scale_x_continuous(name = 'Chromosome', breaks = ticks, labels = c(1:22, 'X'), limits = c(0, sum(chrlen[1:23]))) +
     scale_y_continuous(expression(-log[10](italic(p))), breaks = 0:8, limits = c(0, .5 + -log10(min(df$P)))) +
     scale_colour_manual(guide = FALSE, values = c('gray10', 'gray60')) +
     geom_hline(yintercept = -log10(5e-8), color = 'red') +
     theme_bw()

png(paste(fname, '_MA.png',sep=''), width = 960, height = 480, type = 'cairo')
print(p)
dev.off()

#======QQ plot
p=df$P
N <- length(p) # number of p-values
# compute the null distribution 
MAX <- -log10(min(min(p), 1/N))
# compute the confidence intervals
if (N>100) {
  ind <- c(1:100, seq(101, N, (N-101)/round(sqrt(N))))
} else {
  ind <- c(1:N)
}
c95 <- rep(0, length(ind))
c05 <- rep(0, length(ind))
# the jth order statistic from a 
# uniform(0,1) sample 
# has a beta(j,n-j+1) distribution 
# (Casella & Berger, 2002, 
# 2nd edition, pg 230, Duxbury)
for(i in 1:length(ind)) {
  c95[i] <- qbeta(0.95, ind[i], N-ind[i]+1)
  c05[i] <- qbeta(0.05, ind[i], N-ind[i]+1)
}
#==start the plot
png(paste(fname, '_QQ.png',sep=''), width = 480, height = 480, type = 'cairo')
# plot the two confidence lines
par(mar = c(5,4,4,2) + 0.1 + c(0, 2, 0, 0))
plot(-log10(ind/N), -log10(c95), ylim=c(0,MAX), xlim=c(0,MAX), type='l', axes=FALSE, xlab='', ylab='')
par(new=T)
plot(-log10(ind/N), -log10(c05), ylim=c(0,MAX), xlim=c(0,MAX), type='l', axes=FALSE, xlab='', ylab='')
# add the diagonal
abline(0, 1, col='red')
par(new=T)
p <- sort(p)
idx1 <- rev(!duplicated(p))
p <- sort(p, decreasing=T)
idx2 <- !duplicated(p)
x <- -log10(1-(which(idx1) + which(idx2) - 2)/2/N)
y <- -log10(p[idx2])
plot(x, y, ylim = c(0, MAX), xlim = c(0, MAX), xlab = expression('Expected -log'[10]*'(p)'), ylab = expression('Observed -log'[10]*'(p)'), cex.axis = 1.5, cex.lab = 1.5)
dev.off()

}