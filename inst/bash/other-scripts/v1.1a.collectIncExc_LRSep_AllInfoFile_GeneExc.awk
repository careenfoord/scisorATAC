function abs(x){
    if(x>=0){return(x);}
    x=(-1)*x;
    return(x);
}
function min(a,b){    
    if(a<=b){return(a);}
    return(b);
}
function max(a,b){   
    if(a>=b){return(a);}
    return(b);
}
function allNeccesaryOptionsThere(){
          
    if(!chr){
	print "-v chr=??? must be provided; exiting" > "/dev/stderr";
	exit;
    }
    else{
	print "OPTION chr="chr  > "/dev/stderr";
    }
    if(!strand){
	print "-v strand=??? must be provided; exiting" > "/dev/stderr";
	exit;
    }
    else{
	print "OPTION strand="strand  > "/dev/stderr";
    }
    if(!fileGZ){
	print "-v fileGZ=??? must be provided; exiting" > "/dev/stderr";
	exit;
    }
    else{
	print "OPTION fileGZ="fileGZ  > "/dev/stderr";
    }
    if(!exons2BeConsidered){
	print "-v exons2BeConsidered=??? must be provided; exiting" > "/dev/stderr";
	exit;
    }
    else{
	print "OPTION exons2BeConsidered="exons2BeConsidered  > "/dev/stderr";
    }
    if(!unzipCommand){
	print "-v unzipCommand=??? must be provided; exiting" > "/dev/stderr";
	exit;
    }
    else{
	print "OPTION unzipCommand="unzipCommand  > "/dev/stderr";
    }


}


BEGIN{

    allNeccesaryOptionsThere();

    OFS="\t";
    totExons=0;
    totTranscripts=0;

    print "# 0. all the exons we are interested in" > "/dev/stderr";
    comm=unzipCommand" "exons2BeConsidered;
    while(comm | getline){	
	split($1,a,"_");
	if(toupper(a[1])!=toupper(chr) || toupper(a[4])!=toupper(strand)){continue;}
	EXON_IN[a[2]"_"a[3]"_"$2]=0;
	EXON_OUT[a[2]"_"a[3]"_"$2]=0;
	EXON_OTHER[a[2]"_"a[3]"_"$2]=0;
	if(!(a[2] in pos2IntExon)){	    
	    pos2IntExon[a[2]]=a[2]"_"a[3]"_"$2;
	}
	else{
	    pos2IntExon[a[2]]=pos2IntExon[a[2]]";.;"a[2]"_"a[3]"_"$2;
	}	
    }

    
    print "# A. treating input file "file" to read exons and introns" > "/dev/stderr";
    comm=unzipCommand" "fileGZ;
    cou=0;
    while(comm | getline){	
	
	#if(cou % 1000000 ==0){print cou > "/dev/stderr";}
	#cou++;

	if($2 ~/@/ || $2 =="none"){continue;}
	moleculeKey=$2"\t"$4"\t"$5;
	if(moleculeKey in seenMolecules){continue;}
	seenMolecules[moleculeKey]=1;
	
	nExonsPlus1=split($9,a,";%;");
	if(nExonsPlus1<3){
	    print "ERROR:n="n" for read "$1 > "/dev/stderr";
	    exit(0);
	}
   	
	split(a[2],b,"_");
	
	if(toupper(b[1])!=toupper(chr)){continue;}
	if(toupper(b[4])!=toupper(strand)){continue;}
	#print "bb:"a[2]
	#if(cou % 10 ==0){print cou > "/dev/stderr";}
	#cou++;

	# if we get here this is a read on the right chromosome and strand and that has been assigned to one gene only
	trID=$1;
	tr2Gene[$1]=$2;
	for(i=2;i<=nExonsPlus1;i++){
	    split(a[i],b,"_");
	    exIDwithChromAndGene=b[2]"_"b[3]"_"$2;

	    if(trID in tr2Exons){
		tr2Exons[trID]=tr2Exons[trID]";.;"exIDwithChromAndGene;
	    }
	    else{
		tr2Exons[trID]=exIDwithChromAndGene;
	    }	    
	}
    }
       
    print "# B. saving all internal exons" > "/dev/stderr";
    cou=0;
    for(t in tr2Exons){	
	n=split(tr2Exons[t],a,";.;");
	if(verbose){print "n="n > "/dev/stderr";}
	split(a[1],b,"_");
	trStart[t]=b[1];
	if(b[2] in siteToStartReads){
	    siteToStartReads[b[2]]=siteToStartReads[b[2]]";.;"t;
	    siteToStartCorrespondingPos[b[2]]=siteToStartCorrespondingPos[b[2]]";.;"b[1];
	}
	else{
	    siteToStartReads[b[2]]=t;
	    siteToStartCorrespondingPos[b[2]]=b[1];
	}
	
	
	#EXTERNAL_EXON_IN[a[1]]++;
	split(a[n],b,"_");
	trEnd[t]=b[2];
	if(b[1] in siteToEndReads){
	    siteToEndReads[b[1]]=siteToEndReads[b[1]]";.;"t;
	    siteToEndCorrespondingPos[b[1]]=siteToEndCorrespondingPos[b[1]]";.;"b[2];
	}
	else{
	    siteToEndReads[b[1]]=t;
	    siteToEndCorrespondingPos[b[1]]=b[2];
	}
	

	for(i=2;i<=n;i++){
	    # treaing internal exons
	    if(i<=n-1){
		if(!(a[i] in EXON_IN)){
		    print "ERROR: on "chr" "strand" we have an unknown exon "a[i]" at "i" and n="n" and t="t > "/dev/stderr";
		    exit(0);
		}
				
		EXON_IN[a[i]]++;
		
		# remebering the read for this exon
		if(a[i] in exon2Reads){exon2Reads[a[i]]=exon2Reads[a[i]]";.;"t;}
		else{
		    exon2Reads[a[i]]=t;
		    exon2ExclusionReads[a[i]]="none";
		}
	    }
	    # finding introns
	    split(a[i-1],b,"_");
	    split(a[i],c,"_");	    
	    intStart=b[2]+1;
	    intEnd=c[1]-1;
	    if(intStart>intEnd){
		print "ERROR:"intStart">"intEnd";a[i-1]="a[i-1]";a[i]="a[i] > "/dev/stderr";
		exit(0);
	    }
	    intron=intStart"_"intEnd"_"tr2Gene[t];
	    # remebering the read for this intron
	    if(intron in intron2Reads){intron2Reads[intron]=intron2Reads[intron]";.;"t;}
	    else{intron2Reads[intron]=t;}
	    INTRON_COUNT[intron]++;
	    #if(intron=="1010895_1011365"){
	    #	print "INTRON_COUNT[1010895_1011365]="INTRON_COUNT[intron]";t="t > "/dev/stderr";
	    #}
	    if(intStart in pos2Intron){
		if(!(pos2Intron[intStart] ~ intron)){
		    pos2Intron[intStart]=pos2Intron[intStart]";.;"intron;
		}
	    }
	    else{
		pos2Intron[intStart]=intron;	    
	    }
	    #if(intron=="1010895_1011365"){
	    #	print pos2Intron[intStart]";t="t > "/dev/stderr";
	    #}
	}    
    }
    
    print "# C. going over all introns and checking whether they contain an exon (+50 adjacent bases)" > "/dev/stderr";
    for(ipos in pos2Intron){
	if(verbose)print "ipos="ipos";pos2Intron[ipos]="pos2Intron[ipos] > "/dev/stderr";
	n=split(pos2Intron[ipos],introns,";.;")
	for(i=1;i<=n;i++){
	    if(verbose)print "i="i";introns[i]="introns[i] > "/dev/stderr";
	    split(introns[i],posInt,"_");
	    #if(introns[i]=="1010895_1011365"){
	    #	print "checking region:",posInt[1]+50,posInt[2]-50 > "/dev/stderr";
	    #}
	    for(j=posInt[1]+50;j<=posInt[2]-50;j++){		
		if(verbose)print "j="j > "/dev/stderr";
		if(j in pos2IntExon){
		    #if(introns[i]=="1010895_1011365" && j==1011095){
		    #	print "got here and found:",pos2IntExon[j],";t="t > "/dev/stderr";
		    #}
		    if(verbose)print j, pos2IntExon[j] > "/dev/stderr";
		    m=split(pos2IntExon[j],exons,";.;");
		    for(k=1;k<=m;k++){
			split(exons[k],posExon,"_");
			if(posExon[1]!=j){
			    print "ERROR:posExon[1]="posExon[1]";j="j > "/dev/stderr";
			    exit(0);
			}
			# this intron does skip the exon in question
			if(posExon[2]<=posInt[2]-50 && posInt[3]==posExon[3]){
			    if(!(exons[k] in EXON_OUT)){
				print "ERROR: unknown exon "a[i]" for EXON_OUT" > "/dev/stderr";
				exit(0);
			    }
			    EXON_OUT[exons[k]]+=INTRON_COUNT[introns[i]];
			    #exon2Reads[exons[k]]=exon2Reads[exons[k]]";.;"intron2Reads[introns[i]];
			    if(exon2ExclusionReads[exons[k]]=="none")
				exon2ExclusionReads[exons[k]]=intron2Reads[introns[i]]
			    else
				exon2ExclusionReads[exons[k]]=exon2ExclusionReads[exons[k]]";.;"intron2Reads[introns[i]]
			}
		    }
		}
	    }
	}
    } 
    
    print "# D. adding half inclusion reads and printing inclusion and exclusion reads" > "/dev/stderr";
    OFS="\t";
    for(ex in EXON_IN){

	split(ex,b,"_");
	OL = 0;
	for(t in trStart){
	    if(trStart[t]>b[2] || trEnd[t]<b[1]){continue;}
	    OL++;
	}
	
	
	halfIncReadLeft=0;
	halfIncReadRight=0;

	
	if(b[2] in siteToStartReads){
	    n=split(siteToStartReads[b[2]],locReads,";.;");
	    m=split(siteToStartCorrespondingPos[b[2]],locPos,";.;");
	    if(n!=m){
		print "ERROR:n="n";m="m > "/dev/stderr";
		exit(0);
	    }

	    for(i=1;i<=m;i++){
		if(locPos[i]>=b[1]){
		    exon2Reads[ex]=exon2Reads[ex]";.;"locReads[i];
		    halfIncReadRight++;
		}		
	    }	
	}
	if(b[1] in siteToEndReads){
	    n=split(siteToEndReads[b[1]],locReads,";.;");
	    m=split(siteToEndCorrespondingPos[b[1]],locPos,";.;");
	    if(n!=m){
		print "ERROR:n="n";m="m > "/dev/stderr";
		exit(0);
	    }
	    for(i=1;i<=m;i++){
		if(locPos[i]<=b[2]){
		    exon2Reads[ex]=exon2Reads[ex]";.;"locReads[i];
		    halfIncReadLeft++;
		}		
	    }	
	}
	if(substr(exon2ExclusionReads[ex],1,3)==";.;"){
	    exon2ExclusionReads[ex]=substr(exon2ExclusionReads[ex],4,length(exon2ExclusionReads[ex])-3);
	}
	print chr"_"ex"_"strand, EXON_IN[ex],halfIncReadLeft,halfIncReadRight,EXON_OUT[ex],OL,exon2Reads[ex],exon2ExclusionReads[ex];
    }
}
