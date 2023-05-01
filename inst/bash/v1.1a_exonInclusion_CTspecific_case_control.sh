#!/bin/bash
# by Hagen Tilgner for mapping of 454/pacBio.MOLECULO reads : 5/2014

function checkForFile {
    file=$1;
    if [ ! -f $file ]
    then
	echo "ERROR:"$file" is not a regular file ... exiting "
	exit;
    fi
}

function checkForDirectory {
    dir=$1;
    if [ ! -d $dir ]
    then
	echo "ERROR:"$dir" is not a regular directory ... exiting "
	exit;
    fi
}


function getReadNumbersInParallel2 {

    echo "executing getReadNumbersInParallel2"

    infileGZ=$1;
    echo "infileGZ="$infileGZ
    if [ ! -f $infileGZ ]
    then
	echo "ERROR:"$infileGZ" is not a regular file"
	exit;
    fi
    numberLines=`$unzipCommand $infileGZ | wc -l`
    echo $infileGZ" has "$numberLines" lines"

    outfile=$2;
    echo "outfile="$outfile
    if [ -e $outfile ]
    then
	echo "ERROR:"$outfile" exists in parallel"
	exit;
    fi

    n=$3;
    echo "n="$n
    tmp=$4;
    echo "tmp="$tmp
    chrs=$5;
    echo "chrs="$chrs
    collectCommand="cat "
    rmCommand="rm "

    exons2BeConsidered=$6
    echo "exons2BeConsidered="$exons2BeConsidered

    for c in `cat $chrs`; do
	for s in - +; do
	    echo -e $c"\t"$s;
	    str="awk -f "${scriptDir}"/v1.1a.collectIncExc_LRSep_AllInfoFile_GeneExc.awk -v chr="$c" -v strand="$s" -v fileGZ="$infileGZ" -v exons2BeConsidered="$exons2BeConsidered" -v unzipCommand="$unzipCommand" > "$tmp"/"$$.$c.$s
	    echo $str >> $tmp"/"$$"parallel.comp.guide";
	    collectCommand=$collectCommand" "$tmp"/"$$.$c.$s
	    rmCommand=$rmCommand" "$tmp"/"$$.$c.$s
	done
    done

    time python2 ${scriptDir}/v0.2.executeInParallel.py --commandFile $tmp"/"$$"parallel.comp.guide" --n $n
    echo $collectCommand;
    echo $rmCommand;
    `$collectCommand > $outfile`;
    gzip $outfile
    `$rmCommand`;
    rm $tmp"/"$$"parallel.comp.guide"

}



# when and where did we do this ?
#############################
# when and where did we do this ?
{ echo "#################################";
  echo "# RUNNING [$0]";
#  echo "##### on params [$1]";
  echo "# Current date:`date`";
  echo "# Current dir: `pwd`"; } 1>&2;

###############
# 0. the arguments
echo "+++++++++++++++++++++++++ 1. arguments";

echo "++++++++++++++++++ 1a. flexible";

ALLINFOcaseList=$1
echo "ALLINFOcaseList="$ALLINFOcaseList
checkForFile $ALLINFOcaseList

ALLINFOcontrolList=$2
echo "ALLINFOcontrolList="$ALLINFOcontrolList
checkForFile $ALLINFOcontrolList

tmpDirFile=$3
echo "tmpDirFile="$tmpDirFile
checkForFile $tmpDirFile

chroms=${4};
echo "chroms="$chroms
checkForFile $chroms

numCPUsHighMemory=$5
echo "numCPUsHighMemory="$numCPUsHighMemory

annotationGZ=$6
echo "annotationGZ="$annotationGZ

scriptDir=${7}
echo "scriptDir="$scriptDir

minPSI=${8}
echo "minPSI="$minPSI

maxPSI=${9}
echo "maxPSI="$maxPSI

minReadNumber=${10}
echo "minReadNumber="$minReadNumber

unzipCommand=${11}
echo "unzipCommand="$unzipCommand

minOLfrac=${12}
echo "minOLfrac="$minOLfrac

cellTypeFile=${13}
echo "cellTypeFile="$cellTypeFile


echo "++++++++++++++++++ 1b. deduced from arguments";
tmpdir1=`head -1 $tmpDirFile`
echo "tmpdir1="$tmpdir1
tmpdir1=$tmpdir1"tmp."$$
echo "tmpdir1="$tmpdir1
mkdir $tmpdir1
echo "tmpdir1="$tmpdir1



echo "+++++++++++++++++++++++++ 2. data organization";

echo "++++++++++++++++++ 2.a. organizing data";

echo "+++++++++++++++ 2.a.1 case reads into one file";
for file in `cat $ALLINFOcaseList`; do $unzipCommand $file; done | gzip -c > $tmpdir1"/cases.allInfo.gz"

echo "+++++++++++++++ 2.a.2 control reads into one file";
for file in `cat $ALLINFOcontrolList`; do $unzipCommand $file; done | gzip -c > $tmpdir1"/controls.allInfo.gz"

echo "+++++++++++++++ 2.a.3 combined";
$unzipCommand $tmpdir1"/cases.allInfo.gz" $tmpdir1"/controls.allInfo.gz" | gzip -c > $tmpdir1"/combined.allInfo.gz"


echo "+++++++++++++++ 2.a.1 read stats";

n=`$unzipCommand $tmpdir1"/cases.allInfo.gz" | wc -l`;
echo "# total number of reads for cases: "$n;

n=`$unzipCommand $tmpdir1"/controls.allInfo.gz" | wc -l`;
echo "# total number of reads for controls: "$n;




echo "++++++++++++++++++ 2.b. getting inclusion and exclusion reads";

echo "+++++++++++++++ 2.b_1 getting all exons for all reads";

$unzipCommand $tmpdir1"/combined.allInfo.gz" | awk '{if($2!~/@/ && $2!="none"){print $1"\t"$2"\t"$9;}}' | awk '{n=split($3,a,";%;"); if(n<3){print "ERROR:n="n" for read "$1 > "/dev/stderr";} for(i=3;i<=n-1;i++){print a[i]"\t"$2;}}' | sort | uniq -c | awk '{print $2"\t"$3"\t"$1;}' | gzip -c > all.internalExons.tab.gz

echo "+++++++++++++++ 2.b_2 getting inclusion and exclusion reads and remving intron retention";
getReadNumbersInParallel2 $tmpdir1"/combined.allInfo.gz" sampleBoth.inc.exc.tab $numCPUsHighMemory $tmpdir1 $chroms all.internalExons.tab.gz

$unzipCommand sampleBoth.inc.exc.tab.gz | awk -v anno=$annotationGZ -v unzipCommand=$unzipCommand 'BEGIN{comm=unzipCommand" "anno; while(comm|getline){if($3=="exon"){exon[$1"_"$4"_"$5"_"$7]=1;  for(i=$4;i<=$5;i++){exonic[$1"_"i"_"$7]=1;}     st[$1"_"$4"_"$7]=1; en[$1"_"$5"_"$7]=1;}}}{split($1,a,"_"); x=0; for(i=a[2];i<=a[3];i++){ s=a[1]"_"i"_"a[5]; if(!(s in exonic)){x++;}  }   if(!(a[1]"_"a[2]"_"a[3]"_"a[5] in exon) && x>=70 && (a[1]"_"a[2]"_"a[5] in st || a[1]"_"a[3]"_"a[5] in en)){next;} print $0;}' | gzip -c > sampleBoth.inc.exc.tab.NoRetention.gz

echo "+++++++++++++++ 2.b_3 getting genes that have >= 2 exons with >=minPSI and <=maxPSI exon inclusion";
echo "min PSI " $minPSI
echo "max PSI " $maxPSI
echo "min OL frac " $minOLfrac
echo "min read number "  $minReadNumber
$unzipCommand sampleBoth.inc.exc.tab.NoRetention.gz | awk -v minPSI=$minPSI -v maxPSI=$maxPSI -v minOLfrac=$minOLfrac -v minReadNumber=$minReadNumber '{if($2+$3+$5==0 || $2+$4+$5==0){next;} if($2+$3+$4+$5<minReadNumber){next;} psi=($2+$3+$4)/($2+$3+$4+$5); leftPSI=($2+$3)/($2+$3+$5); rightPSI=($2+$4)/($2+$4+$5); split($1,a,"_"); if(($2+$3+$4+$5)/$6<minOLfrac){next;} if(rightPSI>=minPSI && rightPSI<=maxPSI && leftPSI>=minPSI && leftPSI<=maxPSI && psi>=minPSI && psi<=maxPSI){print $0; }}' | gzip -c > all.altExons.tab.gz

echo "++++++++++++++++++ 2.c. getting inclusion and exclusion reads for each cell type";
echo "+++++++++++++++ 2.c_1 in separate output files";

for CT in `cat $cellTypeFile`; do
    $unzipCommand $tmpdir1"/combined.allInfo.gz" | awk -v CT=$CT '$3==CT' | gzip -c > tmp.gz
    getReadNumbersInParallel2 tmp.gz sampleBoth.inc.exc.tab.$CT $numCPUsHighMemory $tmpdir1 $chroms all.internalExons.tab.gz
    #gzip sampleBoth.inc.exc.tab.$CT
    #$unzipCommand sampleBoth.inc.exc.tab.$CT.gz | awk -v exonsToUse=all.altExons.tab.gz -v unzipCommand=$unzipCommand 'BEGIN{comm=unzipCommand" "exonsToUse; while(comm|getline){use[$1]=1;}}{if($1 in use){print $0;}}' | gzip -c > all.altExons.inc.exc.tab.$CT.gz
done
rm tmp.gz



echo "+++++++++++++++ 2.c_2 in combing the data";
awk -v cellTypeFile=$cellTypeFile -v unzipCommand=$unzipCommand -v ASexons=all.altExons.tab.gz -v CTspecificTrunk=sampleBoth.inc.exc.tab -v minReadNumber=$minReadNumber 'BEGIN{comm="cat "cellTypeFile; ctNumber=0; while(comm|getline){ctNumber++; CT[ctNumber]=$1;}  comm=unzipCommand" "ASexons; while(comm|getline){exonsToUse[$1]=1; for(i=1;i<=ctNumber;i++){inc[$1"\t"CT[i]]=0; exc[$1"\t"CT[i]]=0; skip[$1"\t"CT[i]]=0; psi[$1"\t"CT[i]]="NA";}} for(i=1;i<=ctNumber;i++){comm=unzipCommand" "CTspecificTrunk"."CT[i]".gz"; while(comm|getline){if(!($1 in exonsToUse)){continue;} if($2+$3+$4+$5<minReadNumber){continue;} psi[$1"\t"CT[i]]=($2+$3+$4)/($2+$3+$4+$5);} close(comm);}  commentLine="exonID"; for(i=1;i<=ctNumber;i++){commentLine=commentLine"\t"CT[i];} print "#"commentLine; for(ex in exonsToUse){str=ex; for(i=1;i<=ctNumber;i++){str=str" "psi[ex"\t"CT[i]]; }      print str;}}' | gzip -c > all.altExons.matrix.tab.gz

echo "+++++++++++++++++++++++++ 3. splitting cases and controls: all.altExons.matrix.tab.gz is our universe of exons for the cell types called"

echo "++++++++++++++++++++++ 3.0 reformatting the universe";
$unzipCommand all.altExons.matrix.tab.gz | awk '{split($1,a,"_"); if($1~/chr/){print a[1]"_"a[2]"_"a[3]"_"a[5]"\t"a[4];}}' | gzip -c > universe.tab.gz


echo "++++++++++++++++++++++ 3.1 cases";
getReadNumbersInParallel2 $tmpdir1"/cases.allInfo.gz" CASES.sampleBoth.inc.exc.tab $numCPUsHighMemory $tmpdir1 $chroms all.internalExons.tab.gz


echo "++++++++++++++++++++++ 3.2 controls";
getReadNumbersInParallel2 $tmpdir1"/controls.allInfo.gz" CONTROLS.sampleBoth.inc.exc.tab $numCPUsHighMemory $tmpdir1 $chroms  all.internalExons.tab.gz

echo "++++++++++++++++++++++ 3.3 collect data from both files to get inclusion and exclusion reads from both sources";
awk -v uni=universe.tab.gz -v unzipCommand=zcat -v cases=CASES.sampleBoth.inc.exc.tab.gz -v controls=CONTROLS.sampleBoth.inc.exc.tab.gz 'BEGIN{comm=unzipCommand" "uni; while(comm|getline){split($1,a,"_"); exonGene=a[1]"_"a[2]"_"a[3]"_"$2"_"a[4]; caseInc[exonGene]=0; caseExc[exonGene]=0; conInc[exonGene]=0; conExc[exonGene]=0;}  comm=unzipCommand" "cases; while(comm|getline){if($1 in caseInc){caseInc[$1]=$2+$3+$4; caseExc[$1]=$5;}} comm=unzipCommand" "controls; while(comm|getline){if($1 in conInc){conInc[$1]=$2+$3+$4; conExc[$1]=$5;}}     for(k in caseInc){print k,caseInc[k],caseExc[k], conInc[k],conExc[k];}     }' | gzip -c > cases_vs_controls.counts.tab.gz

echo "++++++++++++++++++++++ 3.4 check chi-squared criterion and only output those that pass";
zcat cases_vs_controls.counts.tab.gz | awk '{x=0; n=$2+$3+$4+$5; if(n==0){next;} f1=($2+$3)/n;  f2=($2+$4)/n; if(f1*f2*n>=5){x++;} if(f1*(1-f2)*n>=5){x++;} if((1-f1)*f2*n>=5){x++;} if((1-f1)*(1-f2)*n>=5){x++;} if(x>=3){print $0;}}' | gzip -c > cases_vs_controls.counts.passed-chi-sq-crit.tab.gz

echo "++++++++++++++++++++++ 3.5 testing 2 x 2 tables and correcting for multiple testing";
R --vanilla --slave --args cases_vs_controls.counts.passed-chi-sq-crit.tab.gz two.sided BY < ${scriptDir}/v1.1a_fisherforNatan.r | gzip -c > cases_vs_controls.counts.pVal.FDR.LOR.tab.gz


echo "+++++++++++++++++++++++++ 4. cleaning up ";
rm $tmpdir1"/cases.allInfo.gz"
rm $tmpdir1"/controls.allInfo.gz"
rm $tmpdir1"/combined.allInfo.gz"


rmdir $tmpdir1;
exit;
