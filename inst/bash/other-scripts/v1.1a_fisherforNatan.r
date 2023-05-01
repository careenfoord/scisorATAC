
fT<-function(X,direc){
  #cat(X[2],X[3], X[4], X[5],X[6], X[7],"\n", file=stderr()) 
  m=matrix(c(as.numeric(X[2]), as.numeric(X[3]),as.numeric(X[4]), as.numeric(X[5])), nrow=2,ncol=2,byrow=TRUE);
  #print(X)
  ##rownames(m)=c("ex1In","ex1Out")
  ##colnames(m)=c("ex2In","ex2Out")
  #print(m)  
  mp=fisher.test(m,alternative=direc);
  return(mp$p.value);  
}

LOR<-function(X){
  #cat(X[2],X[3], X[4], X[5],X[6], X[7],"\n", file=stderr()) 
  m=matrix(c(as.numeric(X[2]), as.numeric(X[3]),as.numeric(X[4]), as.numeric(X[5])), nrow=2,ncol=2,byrow=TRUE);
  if(m[1,1] ==0){m[1,1]=0.5;}  
  if(m[2,1] ==0){m[2,1]=0.5;}  
  if(m[1,2] ==0){m[1,2]=0.5;}  
  if(m[2,2] ==0){m[2,2]=0.5;}  
  #rownames(m)=c("ex1In","ex1Out")
  #colnames(m)=c("ex2In","ex2Out")
  #mp=fisher.test(m,alternative=direc);
  return(log(m[1,1],base=2) + log(m[2,2],base=2) - log(m[1,2],base=2) - log(m[2,1],base=2));  # log((a*d)/(b*c)) = log(a) + log(d) - log(b) - log(c)
}

# written by Hagen Tilgner,forgot to record date!
args<-commandArgs(trailingOnly=TRUE);



cat("# 1. options\n", file=stderr())

inputFile=args[1];
cat("inputFile=",inputFile,"\n", file=stderr())
direction=args[2];
cat("direction=",direction,"\n", file=stderr())
correction=args[3];
cat("correction=",correction,"\n", file=stderr())

if(direction!="two.sided" && direction!="greater" && direction!="less"){
	stop(paste("ERROR:",direction," is not a legal direction for the fishe test"));
}

if(correction!="BH" && correction!="bonferroni" && correction!="BY" && correction!="holm" && correction!="hochberg" && correction!="hommel" && correction!="fdr"){
	stop(paste("ERROR:",correction, " is not a legal correction method"));
}

cat("# 2. reading data\n", file=stderr()) 
data2=read.table(inputFile);
#data2=data[2:length(data[,1]),]

cat("# 3. data structures\n", file=stderr())
pVals=apply(data2,MARGIN=1,fT,direction);
#print(pVals);

cat("# 4. p-value correction\n", file=stderr())
pvalsCorrected=p.adjust(pVals, method = correction);
#if(correction=="BH"){
#  pvalsCorrected=p.adjust(pVals, method = "BH");
#}
#if(correction=="bonferroni"){
#  pvalsCorrected=p.adjust(pVals, method = "bonferroni");
#}

cat("# 5. LOG-ODDS ratios\n", file=stderr())
lor=apply(data2,MARGIN=1,LOR);


cat("# 6. printing\n", file=stderr()) 
for(i in 1:length(data2[,1])){    
  cat(as.character(data2[i,1]),"\t", as.character(data2[i,2]),"\t",  as.character(data2[i,3]),"\t", as.character(data2[i,4]),"\t",  as.character(data2[i,5]), "\t", pVals[i],"\t", pvalsCorrected[i],"\t", lor[i],"\n");	      
}
