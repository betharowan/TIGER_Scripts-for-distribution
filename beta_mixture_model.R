# TODO: Add comment
# 
# Author: patel
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
sliding_window_outcome<-args[1]
outfile<-args[2]

cutoffs<-function(p1,p2,alpha,beta)
{
	return(beta(alpha+p1,beta+p2)/beta(alpha,beta));
}

beta_prob<-function(y,alpha_v,beta_v)
{
	## if(alpha_v <0 | beta_v<0)
	## {
	##     return (-10000000000000);
	## }
	
	return ((y^(alpha_v-1)*(1-y)^(beta_v-1))/beta(alpha_v, beta_v))
}

beta_prob_2<-function(y,alpha_v,beta_v,p)
{
	## if(alpha_v <0 | beta_v<0)
	## {
	##     return (-10000000000000);
	## }
	loglike<-(log(y^(alpha_v-1)*(1-y)^(beta_v-1)/beta(alpha_v, beta_v)))+log(p)
	return (loglike)
}
mix_beta_model_log_function<-function(x,data,z,p_)
{
	alpha<-x[1:3];
	beta<-x[4:6];
	#total_sum<-c();
	local_sum<-sapply(data,beta_prob_2,alpha=x[1:3],beta=x[4:6],p=p_[])	
	ll<-sum(t(local_sum)*z)*-1
	return(ll)
}


mix_beta_model<-function(x,data,z,p_)
{
	#cat("now running it with ",x," ",p,"\n",file="data1.txt", sep="\t", append=TRUE)
	alpha<-x[1:3];
	beta<-x[4:6];
	
	midterm<-c();
	for(j in 1:3)
	{			
		local_sum<-0;
		local_sum<-log(beta_prob(data,alpha[j],beta[j]))	
		local_sum<-(local_sum+log(p_[j]));
		midterm<-append(midterm,local_sum);
	}
	ll<-sum(z*midterm)*(-1);
	#ll<-sumi*(-1);
	return (ll);	
}


dominator_nominator_func<-function(y,p,alpha,beta)
{
	return((p*beta_prob(y,alpha,beta)))
}


e_step<-function(x,data,p)
{
	alpha<-x[1:3];
	beta<-x[4:6];
	new_z<-c()
	dominator<-sapply(data,dominator_nominator_func,p=p[],alpha=x[1:3],beta=x[4:6])
	new_z<-apply(dominator, 1, function(x) x / colSums(dominator))
	return (new_z)
}

prob<-function(z)
{
	return(colSums(z)/dim(z)[1])
}

m_step<-function(x,data,z,p)
{
	v<-optim(par=x,fn=mix_beta_model_log_function,data=data,p=p,z=z)
	
	return (v$par);
}
#####

main_f<-function(filename)
{
	
	d<-read.table(filename);
	
	#sample down to 50%
	k<-as.data.frame(table(d[,3]))
	k[,2]<-k[,2]*0.5
	data<-rep(as.numeric(levels(k[,1]))[k[,1]],k[,2])
	data <- data + 100
	data <- data/200
	data[data==0] = 0.00001
	data[data==1] = 1-0.00001
	#d[,3]<-d[,3]+100
	#d[,3]<-d[,3]/200
	#d[d[,3]==0,3]<-0.00001 #beta function is not defined at 0 and 1 
	#d[d[,3]==1,3]<-1-0.00001
	
	#data<-d[,3];
	z<-c();
	
	
	z<-matrix(0,length(data),3)
	for(i in 1:length(data))
	{
		if(data[i]>0.75)
		{
			z[i,]<-c(0.001,0.001,0.998);
		}else if(data[i]<0.25)
		{
			z[i,]<-c(0.998,0.001,0.001);
		}else{
			z[i,]<-c(0.001,0.998,0.001);
		}
	}
	x<-c(1,1,1,1,1,1);
	
	
	
	counter<-0;
	v<-c();
	diff_log=1000000;
	while(TRUE & counter < 100)
	{
		#cat("counter: ",counter," \n");
		#z<-model.matrix(~0 + as.factor(z));
		p<-prob(z);
		#cat("p ",p,"\n");
		#m-step
		#v<-nlm(f=mix_beta_model,p=x,data=data,p_=p,z=z)
		v<-nlm(f=mix_beta_model,p=x,data=data,p_=p,z=z,print.level=0)
		if(length(v$estimate[v$estimate>50]))
		{
			break;
		}
		
		#v<-nlm(f=mix_beta_model_log_function,p=x,data=data,p_=p,z=z)
		#e-step
		#z<-e_step(v$par,data,p)
		z<-e_step(v$estimate,data,p)
		x<-v$estimate;
		if(abs(v$minimum-diff_log)<0.2)
		{
			break;
		}
		diff_log<-v$minimum;
		counter<-counter+1;
		#cat(v)
		#produce the histogram 
	}
	#return (v)
	#while
	#output

	#rownames(z)<-d[,3];
	#t<-z[z[,1]<z[,2] & z[,2]>z[,3],]
	s<-seq(0.001,0.99,by=0.001)
	k<-matrix(0,4,length(s))
	k[1,]<-s
	k[2,]<-dbeta(s,x[1],x[1+3])
	k[3,]<-dbeta(s,x[2],x[2+3])
	k[4,]<-dbeta(s,x[3],x[3+3])
	k<-t(k);
	t<-k[k[,2]<k[,3] & k[,3]>k[,4],]

	upper_bound<-max(t[,1])
	lower_bound<-min(t[,1])
	
	upper_bound<-upper_bound*200-100
	lower_bound<-lower_bound*200-100
	#now we have to take care, if values are to extreme.......
	#based on expierence
	if(length(lower_bound)==0 || lower_bound >0)
	{
		lower_bound<--25;
		
	}else if(lower_bound < -70)
	{
		lower_bound<- -50;
	}
	
	if(length(upper_bound)==0 || upper_bound> 90)
	{
		upper_bound<-50;
	}else if(upper_bound < 10)
	{
		upper_bound<- 25;
	}

	border_matrix<-matrix(0,2,1)
	border_matrix[1,]<-lower_bound
	border_matrix[2,]<-upper_bound
#	write.table( border_matrix, sep="\n", file="allele_count_population_corrected.base_call_pre_sw_allele_border", row.names=F,col.names=F)
	write.table( border_matrix, sep="\n", file=outfile, row.names=F,col.names=F)


}
main_f(sliding_window_outcome)





## expensiveBigLibraryFunction <- function(x,warning=function(w) {print(paste('warning:',w));browser()},error=function(e) {print(paste('e:',e));browser()}                                         ) 
## {
##    print(paste("big expensive step we don't want to repeat for x:",x))
##    z <- x  # the "expensive operation" 
##            # (not really, just standing in for computation)
##    repeat 
##     withRestarts(
##            withRestarts(
##                  tryCatch(   # you could call withCallingHandlers  # with identical arguments here, too
##                  {
##                       print(paste("attempt cheap operation for z:",z))
##                       return(log(z))
##                  },warning = warning,error = error),
##                 flipArg = function() {z <<- -z} ),
##             zapOutArg = function() {z <<- 1} ) 
## }




## s<-seq(0.001,0.99,by=0.001)
## xx<-v$estimate;
## plot(s,dbeta(s,xx[3],xx[3+3]),col="blue")
## points(s,dbeta(s,xx[1],xx[1+3]),col="red")
## points(s,dbeta(s,xx[2],xx[2+3]),col="green")
## points(s,dbeta(s,xx[3],xx[3+3]),col="blue")
