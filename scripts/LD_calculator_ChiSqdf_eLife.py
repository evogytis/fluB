#####################################################################################
##
##  This script calculates linkage disequilibrium (LD) between sites in different alignments (nucleotide or amino acid).
##  It is meant to accompany the manuscript "Reassortment between influenza B lineages and the emergence of a co-adapted PB1-PB2-HA gene complex"
##  and has only been tested with alignments produced specifically for it.
##  You must have python, numpy and biopython installed on your machine to use this script.
##
##  You also require either protein or nucleotide alignments (in nexus format) with the same number of taxa that have the same names.
##  The sequences used for this project are listed in the acknowledgment table accompanying the manuscript:
##  https://github.com/evogytis/fluB/tree/master/acknowledgement%20tables
##
##  and can be downloaded from GISAID available at:
##  http://platform.gisaid.org
##  
##  open command line, go to folder where the script is located and type in:
##  "python LD_calculatorChiSqdf_eLife.py" (ignore ")
##
##  To switch between nucleotide and amino acid modes comment 
##
##  The following figures in the manuscript were made using this script:
##  Figures 10, S7-S8
##
##
##  Alternatively, instead of running the script the output is available here:
##  https://github.com/evogytis/fluB/tree/master/data
##  
##
##  Gytis Dudas
##  Institute of Evolutionary Biology
##  University of Edinburgh
##  EH9 3JT
##  Edinburgh, UK
##
#####################################################################################


import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
from datetime import datetime as dt
import time
import calendar
from Bio import SeqIO
from Bio import AlignIO

def unique(o, idfun=repr):
    """
    Returns a list of unique values in a list.
    """
    seen = {}
    return [seen.setdefault(idfun(e),e) for e in o if idfun(e) not in seen]

def column(data,col):
    """
    Returns a list that is a column in a list of lists.
    """
    return [row [col] for row in data]

def index(data,item):
    """
    Returns a list of indices where an item can be found within a given list.
    """
    return [i for i,x in enumerate(data) if x == item]

def frequency(data):
    """
    Returns a list of counts of each value in a list.
    """
    uvals=unique(data)
    out = dict(zip(uvals, list(0 for u in range(len(uvals)))))
    for i in range(len(uvals)):
        out[uvals[i]]=data.count(uvals[i])
    return out

def decimalDate(sequenceName):
    """
    Returns a yyyy-mm-dd date format as a decimal date.
    """
    dateRegex='[A-Z]\/[A-Za-z\/\-\_\.0-9]+\/[0-9]+\_([0-9\-]+)'
    C=re.search('%s'%(dateRegex),sequenceName)
    if C is not None:
        try:
            yrs = int(C.group(1).split('-')[0])
            mons = int(C.group(1).split('-')[1])
            dds = int(C.group(1).split('-')[2])
            return (float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))
        except IndexError:
            try:
                u = int(calendar.monthrange(int(C.group(1).split('-')[0]),int(C.group(1).split('-')[1]))[1])
                l = 1
                yrs = int(C.group(1).split('-')[0])
                mons = int(C.group(1).split('-')[1])
                dds=int((u+l)/float(2))
                return (float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))
            except IndexError:
                yrs = int(C.group(1))
                return float(yrs)+0.5

def removeItem(someList,itemList):
    """
    Returns a list with an item excluded.
    """
    return [x for x in someList if x not in itemList]

def removeItem2(someList,itemList):
    """
    Returns a list if an item is present.
    """
    return [x for x in someList if x in itemList]

def flatten(someList):
    """
    Reduces a list of lists into a single list.
    """
    return [item for sublist in someList for item in sublist]

import collections
def overlap(a,b):
    """
    Finds the overlap between two lists.
    """
    a_multiset = collections.Counter(a)
    b_multiset = collections.Counter(b)

    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    return overlap, a_remainder, b_remainder


## Makes a list of dates if partitioning of LD is desired
# Default value makes a single bin for all sequences isolated between 1900 and 2014
step=3
timeline=np.arange(1984,2014,step)

## Segments to analyze                                                                       
segments=['PB1','PB2','PA','HA','NP','NA','M1','NS1']

## Taxa names to exclude from analysis
# LD analysis in manuscript did not exclude any taxa
#mixedLineage=['B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994']
mixedLineage=[]

## Select mode - nt for nucleotide alignments, aa for amino acid alignments
#mode='nt'
mode='aa'


#segments=['NA','M1','NS1']

## Iterate over all pairs of segments
for mem1 in range(len(segments)):
    seen={}
    for mem2 in range(mem1,len(segments)):
        print 'Comparing %s and %s'%(segments[mem1],segments[mem2])

        dictA=list({} for u in range(len(timeline)))
        dictB=list({} for u in range(len(timeline)))
        alignmentA=list([] for u in range(len(timeline)))
        alignmentB=list([] for u in range(len(timeline)))
        strains=list([] for u in range(len(timeline)))

        segment1=segments[mem1]
        segment2=segments[mem2]

        if mode=='aa':
            ## Point to folder containing alignments
            path='/Users/admin/Documents/Viral sequences/InfB reassortment/Alignments/Translated/'
            ## Define alignment name
            seg1='InfB_%s translation alignment.nex'%(segment1)
            seg2='InfB_%s translation alignment.nex'%(segment2)
            f = open('/Users/admin/Documents/Viral sequences/InfB reassortment/Alignments/Translated/%s_%s_rDChicomparisonNoSelfTimed3.csv'%(segment1,segment2),'w')
        elif mode=='nt':
            path='/Users/admin/Documents/Viral sequences/InfB reassortment/Alignments/Without ancient seqs and pruned/'
            seg1='InfB_%s.nex'%(segment1)
            seg2='InfB_%s.nex'%(segment2)
            f = open('/Users/admin/Documents/Viral sequences/InfB reassortment/Alignments/Without ancient seqs and pruned/%s_%s_rDChiNtcomparisonNoSelfRemake.csv'%(segment1,segment2),'w')
        
        handle1 = open(path+seg1, "rU")
        handle2 = open(path+seg2, "rU")

        seqCheck=[]

        ## Iterate over all sequences, assign to bins, exclude some if necessary
        for record1,record2 in zip(AlignIO.read(handle1, "nexus"),AlignIO.read(handle2, "nexus")):
            seqCheck.append([record1,record2])
            for i in range(len(timeline)):
                #if (record1.id not in strains[i] and record1.id not in mixedLineage):
                    
                #if (record2.id not in strains[i] and record2.id not in mixedLineage):
                    
                
                if timeline[i]<=decimalDate(record1.id)<timeline[i]+step and record1.id not in mixedLineage:
                    dictA[i][record1.id]=record1.seq
                    alignmentA[i].append(record1.seq)
                    if record1.id not in mixedLineage and record1.id not in strains[i]:
                        strains[i].append(record1.id)
                else:
                    pass
                    
                if timeline[i]<=decimalDate(record2.id)<timeline[i]+step and record2.id not in mixedLineage:
                    dictB[i][record2.id]=record2.seq
                    alignmentB[i].append(record2.seq)
                    if record2.id not in mixedLineage and record2.id not in strains[i]:
                        strains[i].append(record2.id)
                else:
                    pass

        print 'Sequence number check:',len(seqCheck)
        for tp in range(len(timeline)):
            print 'Numbers of sequences left:',timeline[tp],len(alignmentA[tp]),len(alignmentB[tp]),len(dictA[tp]),len(dictB[tp]),len(strains[tp])
            
            ## Iterate over all pairs of sites in both alignments
            for i in range(len(alignmentA[tp][0])):
                #print timeline[tp],'site',i+1
                for j in range(len(alignmentB[tp][0])):
                    if mode=='aa':
                        totalA=len(removeItem(column(alignmentA[tp],i),['-','?']))
                        totalB=len(removeItem(column(alignmentB[tp],j),['-','?']))
                    elif mode=='nt':
                        totalA=len(removeItem2(column(alignmentA[tp],i),['A','C','T','G']))
                        totalB=len(removeItem2(column(alignmentB[tp],j),['A','C','T','G']))
                    

                    total=min([totalA,totalB])
                    if mode=='aa':
                        poly1=unique(removeItem(column(alignmentA[tp],i),['-','?']))
                        poly2=unique(removeItem(column(alignmentB[tp],j),['-','?']))
                    elif mode=='nt':
                        poly1=unique(removeItem2(column(alignmentA[tp],i),['A','C','T','G']))
                        poly2=unique(removeItem2(column(alignmentB[tp],j),['A','C','T','G']))
                    
                    
                    chiSq=None
                    chi=0

                    ## Count polymorphisms
                    polies=[]
                    for p1 in poly1:
                        if mode=='aa':
                            polies.append(removeItem(column(alignmentA[tp],i),['-','?']).count(p1))
                        elif mode=='nt':
                            polies.append(removeItem2(column(alignmentA[tp],i),['A','C','T','G']).count(p1))
                    for p2 in poly2:
                        if mode=='aa':
                            polies.append(removeItem(column(alignmentB[tp],j),['-','?']).count(p2))
                        elif mode=='nt':
                            polies.append(removeItem2(column(alignmentB[tp],j),['A','C','T','G']).count(p2))

                    ## check the count of the rarest polymorphism, proceed if minimum number met
                    if len(polies)==0:
                        allHap=False
                    ## this is the cutoff for how many times a polymorphism has to be seen to be included in the analysis
                    elif min(polies)>=1:
                        allHap=True
                    else:
                        allHap=False


                    hapCount=0
                    for x in strains[tp]:
                        ## Exclude haplotypes that have gaps
                        if mode=='aa':
                            if dictA[tp].has_key(x) and dictB[tp].has_key(x) and ('-' in [dictA[tp][x][i],dictB[tp][x][j]] or '?' in [dictA[tp][x][i],dictB[tp][x][j]]):
                                allHap=False
                            elif dictA[tp].has_key(x)==False or dictB[tp].has_key(x)==False:
                                print 'missing key:',x,dictA[tp].has_key(x),dictB[tp].has_key(x),decimalDate(x)
                        elif mode=='nt':
                            if dictA[tp].has_key(x) and dictB[tp].has_key(x) and (('A' or 'C' or 'T' or 'G') not in [dictA[tp][x][i],dictB[tp][x][j]]):
                                allHap=False
                            elif dictA[tp].has_key(x)==False or dictB[tp].has_key(x)==False:
                                print 'missing key:',x,dictA[tp].has_key(x),dictB[tp].has_key(x),decimalDate(x)
                                
                    ## Ignore comparisons where a site is compared to itself
                    if seg1==seg2 and i==j:
                        allHap=False

                    ## Construct contingency tables if there are at least 2 isolates with a haplotype
                    if len(poly1)>1 and len(poly2)>1 and allHap==True:
                        chiSq=0
                        contingencyTable=list(list(0 for t in range(len(poly2))) for p in range(len(poly1)))
                        observedTable=list(list(0 for t in range(len(poly2))) for p in range(len(poly1)))
                        expectedTable=list(list(0 for t in range(len(poly2))) for r in range(len(poly1)))

                        ## Count haplotypes
                        for x in strains[tp]:
                            if mode=='aa':
                                if dictA[tp].has_key(x) and dictB[tp].has_key(x) and '-' not in [dictA[tp][x][i],dictB[tp][x][j]] and '?' not in [dictA[tp][x][i],dictB[tp][x][j]]:
                                    contingencyTable[index(poly1,dictA[tp][x][i])[0]][index(poly2,dictB[tp][x][j])[0]]+=1
                                else:
                                    pass
                            elif mode=='nt':
                                if dictA[tp].has_key(x) and dictB[tp].has_key(x) and len(removeItem2([dictA[tp][x][i],dictB[tp][x][j]],['A','C','T','G']))==2:
                                    contingencyTable[index(poly1,dictA[tp][x][i])[0]][index(poly2,dictB[tp][x][j])[0]]+=1
                                else:
                                    pass
                                
                        ## Output stored chiSqdf value if contingency table already seen
                        if seen.has_key(str(contingencyTable))==True:
                            print>>f,'%s,%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],len(strains[tp]),i+1,j+1,len(poly1),len(poly2),seen[str(contingencyTable)]*(len(poly1)-1)*(len(poly2)-1)*total,seen[str(contingencyTable)])
                            #print '%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],i+1,j+1,len(poly1),len(poly2),seen[str(contingencyTable)]*(len(poly1)-1)*(len(poly2)-1)*total,seen[str(contingencyTable)])
                            print '%s\tseqs %s\t%s site %s\t%s site %s\t%s'%(timeline[tp],len(strains[tp]),segments[mem1],i+1,segments[mem2],j+1,seen[str(contingencyTable)])
                        else:
                            total1=0
                            total2=0

                            ## Calculate expected numbers of haplotypes under a null model with given polymorphism frequencies
                            for w in range(len(poly1)):
                                freqs1=0
                                for level1 in range(len(poly2)):
                                    freqs1+=contingencyTable[w][level1]
                                freqs1=freqs1/float(total)
                                total1+=freqs1
                                for q in range(len(poly2)):
                                    freqs2=0
                                    for level1 in range(len(poly1)):
                                        freqs2+=contingencyTable[level1][q]
                                    freqs2=freqs2/float(total)
                                    total2+=freqs2
                                    expectedTable[w][q]=freqs1*freqs2*total
                                    chi+=np.absolute(contingencyTable[w][q]-expectedTable[w][q])
                                    if expectedTable[w][q]==0:
                                        print poly1[w],poly2[q],contingencyTable,expectedTable
                                        print column(alignmentA[tp],i)
                                        print column(alignmentB[tp],j)
                                        print i,j,poly1,poly2
                                        haps=[]
                                        for x in strains[tp]:
                                            print x,dictA[tp].has_key(x),dictB[tp].has_key(x)
                                            if dictA[tp].has_key(x)==True and dictB[tp].has_key(x)==True:
                                                haps.append([dictA[tp][x][i],dictB[tp][x][j]])
                                            else:
                                                print 'dicts are missing keys:',x,dictA[tp].has_key(x),dictB[tp].has_key(x)
                                        print 'unique haplotypes:',unique(haps)

                                    ## Calculate chiSq value
                                    chiSq+=((contingencyTable[w][q]-expectedTable[w][q])**2)/float((expectedTable[w][q]))

                            ## Calculate chiSqdf value
                            #print poly1,poly2
                            chiSqdf=chiSq/(float((len(poly1)-1)*(len(poly2)-1)*total))

                            ## Remember chiSqdf value for this particular contingency table - saves time
                            seen[str(contingencyTable)]=chiSqdf

                            ## Print time period, segment 1 site, segment 2 site, number of polymorphisms in segment 1 at site, number of polymorphisms in segment 2 at site, chi value (not chiSq) and chiSqdf value
                            print>>f,'%s,%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],len(strains[tp]),i+1,j+1,len(poly1),len(poly2),chi,chiSqdf)
                            #print '%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],i+1,j+1,len(poly1),len(poly2),chi,chiSqdf)
                            print '%s\tseqs %s\t%s site %s\t%s site %s\t%s'%(timeline[tp],len(strains[tp]),segments[mem1],i+1,segments[mem2],j+1,chiSqdf)
        f.close()


        handle1.close()
        handle2.close()
print '\nDone!'
