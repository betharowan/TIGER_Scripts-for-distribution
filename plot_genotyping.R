
args <- commandArgs(trailingOnly = TRUE)
sample<-args[1]
outfile<-args[2]
finalfile<-args[3]
finalbreaksfile<-args[4]
finalbreaksrefinedfile<-args[5]
preBaseCalls<-args[6]
sliding_window_breaksfile<-args[7]
print(sample)
print(outfile)
print(finalfile)
print(finalbreaksfile)

pdf(file=outfile, width=18, height=10)


d=read.table(finalfile)
#print(finalfile)
#print(finalbreaksfile)
#print ("done_reading")
b=read.table(finalbreaksfile)
print(finalbreaksrefinedfile)
r=read.table(finalbreaksrefinedfile)
print(preBaseCalls)
pre<-read.table(preBaseCalls)
print(sliding_window_breaksfile)
s<-read.table(sliding_window_breaksfile)

lm=matrix(data=c(1,2,3,4), nrow=4, ncol=1)
layout(lm)
chr_maxi<-max(b$V2)

for (chr in c(1:chr_maxi)) {
	
	#plot genotyping with HMM
	
	plot(1,1,type="n",xlim=c(1,max(b$V4)), ylim=c(0,1), axes=F, xlab="", ylab="", main=paste("Sample ", sample, " -- chr", chr, sep=""))
	rect(b$V3[b$V2==chr], 0.35, b$V4[b$V2==chr], 0.65, col=ifelse(b$V5[b$V2==chr]=="CC", "red", ifelse(b$V5[b$V2==chr]=="LL", "blue", "purple")))
	
	rect(r$V3[r$V2==chr], 0, r$V4[r$V2==chr], 0.3, col=ifelse(r$V5[r$V2==chr]=="CC", "red", ifelse(r$V5[r$V2==chr]=="LL", "blue", "purple")))
	
	
	#plot genotyping data
	
	plot(d$V3[d$V2==chr & d$V8 != 0], d$V8[d$V2==chr & d$V8!=0], col="red", type="h", ylim=c(-10,10), xlim=c(0, max(d$V3)), axes=F, xlab="Position", ylab="Allele count")
	axis(2)
	points(d$V3[d$V2==chr & d$V10!=0], d$V10[d$V2==chr & d$V10 != 0]*-1, col="blue", type="h")
	
	plot(1,1,type="n",xlim=c(1,max(b$V4)), ylim=c(0,1), axes=F, xlab="", ylab="", main="")
	rect(s$V3[s$V2==chr], 0.35, s$V4[s$V2==chr], 0.65, col=ifelse(s$V5[s$V2==chr]=="CC", "red", ifelse(s$V5[s$V2==chr]=="LL", "blue", "purple")))
	
	plot(pre[pre[,1]==chr,2],pre[pre[,1]==chr,3],xlab="Position", axes=F,xlim=c(0, max(d$V3)),ylim=c(-100,100),ylab="allele freq in %",type="l",lty=1)
	lines(c(0,max(d$V3[d$V2==chr])),c(0,0),col="purple",lty=2)
	lines(c(0,max(d$V3[d$V2==chr])),c(100,100),col="red",lty=2)
	lines(c(0,max(d$V3[d$V2==chr])),c(-100,-100),col="blue",lty=2)
	axis(1, at=c(1, max(d$V3[d$V2==chr])/2,max(d$V3[d$V2==chr])))
	axis(2)
}


dev.off()

