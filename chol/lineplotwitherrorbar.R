#install.packages("ggplot2")
library(ggplot2)

df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ToothGrowth, varname="len",groupnames=c("supp", "dose"))

# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)
head(df2)

p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9)) 
print(p)
# Finished bar plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+ 
  theme_classic() + scale_fill_manual(values=c('#999999','#E69F00'))

# Keep only upper error bars
ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=len, ymax=len+sd), width=.2,
                position=position_dodge(.9)) 

# Default line plot
p<- ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))
print(p)
# Finished line plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'))

# Use geom_pointrange
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))
# Use geom_line()+geom_pointrange()
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line()+
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))

