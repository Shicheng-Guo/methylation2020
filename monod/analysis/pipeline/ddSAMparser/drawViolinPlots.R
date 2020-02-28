require(ggplot2)
require(gridExtra)
pdf("TestMixturesOut.pdf", width=15, height=10)
data <- read.table("TestMixturesOut.txt", header=T, sep="\t");

p1 <- ggplot(data, aes(factor(Mix), Frac.Minor), size=0.5) + 
	geom_violin(trim=FALSE, fill = "grey80", colour = "#3366FF", adjust=0.5)+geom_jitter(height=0) + 
	theme(
		axis.title.x = element_text(face="bold", colour="#990000", size=20), 
		axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
		axis.title.y = element_text(face="bold", colour="#990000", size=20),
		axis.text.y  = element_text(size=16)
		)

p2 <- ggplot(data, aes(factor(Mix), ScoresMajor), size=0.5) +
	geom_violin(trim=FALSE, fill = "grey80", colour = "#3366FF", adjust=0.5)+geom_jitter(height=0) +
        theme(
                axis.title.x = element_text(face="bold", colour="#990000", size=20),
                axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
                axis.title.y = element_text(face="bold", colour="#990000", size=20),
                axis.text.y  = element_text(size=16)
                )

p3 <- ggplot(data, aes(factor(Mix), ScoresMinor), size=0.5) + 
	geom_violin(trim=FALSE, fill = "grey80", colour = "#3366FF", adjust=0.5)+geom_jitter(height=0) +
        theme(
                axis.title.x = element_text(face="bold", colour="#990000", size=20),
                axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
                axis.title.y = element_text(face="bold", colour="#990000", size=20),
                axis.text.y  = element_text(size=16)
                )
grid.arrange(p1, p2, p3, ncol=3, nrow=2)

dev.off()
