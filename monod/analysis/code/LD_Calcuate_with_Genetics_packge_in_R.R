library("genetics")

r1<-genotype(c("C/C","C/C","T/T","C/T"))
LD(r1)
LD

g1 <- genotype( c('T/A', NA, 'T/T', NA, 'T/A', NA, 'T/T', 'T/A','T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',NA, 'T/A', 'T/A',NA) )
g2 <- genotype( c('C/A','C/A','C/C','C/A','C/C','C/A','C/A','C/A','C/A','C/C','C/A','A/A','C/A','A/A','C/A','C/C','C/A', 'C/A', 'C/A', 'A/A') )

# Compute LD on a single pair
g1 <- genotype( c('C/C','T/T','T/C',"C/T") )
g1 <- genotype( c('C/C','T/T','T/C',"C/T") )

LD(g1,g2)

x<-c(0,0,0,0)
y<-c(1,1,1,1)

x<-c(0,1,0,1)
y<-c(0,1,0,1)

summary(lm(x~y))

cor(x,y)
