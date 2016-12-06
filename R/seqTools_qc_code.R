### R code from vignette source 'seqTools_qual_report.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: seqTools_qual_report.Rnw:88-95
###################################################
library(seqTools)
fqdir<-"C:/Users/jerem/Dropbox/Bioinformatics/Reports/macpipe_trans/data"
fqq<-fastqq(file.path(fqdir,
                      c("Mcla_R1_subsample.fastq.gz", "Mcla_R2_subsample.fastq.gz")),
            k=4,probeLabel=c("g4","g5"))
# Set this to some other location
basedir<-getwd()


###################################################
### code chunk number 2: seqTools_qual_report.Rnw:99-100
###################################################
fqq


###################################################
### code chunk number 3: seqTools_qual_report.Rnw:105-109
###################################################
dfr<-data.frame(file=basename(fileNames(fqq)),
                sample=probeLabel(fqq),
                reads=format(nReads(fqq), big.mark=Sys.localeconv()[7]))
print(dfr)


###################################################
### code chunk number 4: seqTools_qual_report.Rnw:117-118
###################################################
plotNucCount(fqq)


###################################################
### code chunk number 5: seqTools_qual_report.Rnw:123-124
###################################################
plotGCcontent(fqq)


###################################################
### code chunk number 6: seqTools_qual_report.Rnw:130-142
###################################################
for(i in 1:nFiles(fqq))
{
  file<-file.path(basedir, paste("nucFreq_", i, ".pdf", sep=""))
  pdf(file ,width=6, height=6)
  plotNucFreq(fqq, i)
  invisible(dev.off())
  cat("\\begin{figure}[H]\n")
  cat("\\begin{center}\n")
  cat("\\includegraphics{", file, "}\n",sep="")
  cat("\\end{center}\n")
  cat("\\end{figure}\n\n")
}


###################################################
### code chunk number 7: seqTools_qual_report.Rnw:152-164
###################################################
for(i in 1:nFiles(fqq))
{
  file<-file.path(basedir, paste("phredQuant_", i, ".pdf",sep=""))
  pdf(file)
  plotPhredQuant(fqq, i)
  dev.off()
  cat("\\begin{figure}[H]\n")
  cat("\\begin{center}\n")
  cat("\\includegraphics{", file, "}\n", sep="")
  cat("\\end{center}\n")
  cat("\\end{figure}\n\n")
}


###################################################
### code chunk number 8: seqTools_qual_report.Rnw:172-198
###################################################
fqi<-fqq
probeLabel(fqi)<-paste(1:nFiles(fqi), probeLabel(fqi), sep="_")
lbl<-probeLabel(fqi)
# May set another palette
cols<-terrain.colors(4)

col_label<-function(n)
{
  if(is.leaf(n))
  {
    a<-attributes(n)
    i<-which(a$label==lbl)
    cat(a$label, "\t", i, "\n")
    attr(n,"nodePar") <- c(a$nodePar, list(lab.col=cols[i%%4+1],
                                           pch="", lab.cex=1.2))
  }
  return(n)
}

cbm<-cbDistMatrix(fqi)
hc<-as.dendrogram(hclust(as.dist(cbm)))
hcd<-dendrapply(hc, col_label)
#
op <- par(oma=c(1, 1, 1, 1), mar=c(1, 1, 1, 12) + 0.1)
plot(hcd, horiz=TRUE, edgePar=list(lwd=2, lty=1))
par(op)


