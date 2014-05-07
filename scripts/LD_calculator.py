#####################################################################################
##
##  This script calculates linkage disequilibrium (LD) between sites in different alignments (nucleotide or amino acid).
##  It calculates both the Chi-Squared df and D' statistics of linkage disequilibrium.
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
##  "python LD_calculator.py" (ignore ")
##
##  To switch between nucleotide and amino acid modes comment out mode='aa' or mode='nt', respectively.
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
import itertools

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
step=40
timeline=np.arange(1980,2020,step)

## Segments to analyze                                                                       
segments=['PB1','PB2','PA','HA','NP','NA','M1','NS1']

## Taxa names to exclude from analysis
# LD analysis in manuscript did not exclude any taxa
#mixedLineage=['B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994']
mixedLineage=[]

## Select mode - nt for nucleotide alignments, aa for amino acid alignments
#mode='nt'
mode='aa'

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
            path='/Users/admin/Documents/Viral sequences/InfB reassortment/Analyses/huge LD/'
            ## Define alignment name
            seg1='%s translation alignment.nex'%(segment1)
            seg2='%s translation alignment.nex'%(segment2)
            f = open('/Users/admin/Documents/Viral sequences/InfB reassortment/Analyses/huge LD/%s_%s_Dprime.comparisonNoSelfHUGE.csv'%(segment1,segment2),'w')
            header='time point,time step size,N strains,N haplotypes usable,site 1,site 2,number of polymoprhisms at site 1,number of polymoprhisms at site 2,D,|D\'|,Chi,ChiSqDf'
            print>>f,header
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
            seen={}
            polymorphicPairs=[]

            ## assert that numbers of parsed things are the same
            assert len(alignmentA[tp])==len(alignmentB[tp])==len(dictA[tp])==len(dictB[tp])==len(strains[tp])

            ## print header
            print header.replace(',','\t')

            ## used for debugging, determines which site LD is calculated from in both alignments, 0 by default
            startFrom=0
            ## Iterate over all pairs of sites in both alignments
            for i in range(startFrom,len(alignmentA[tp][0])):
                #print timeline[tp],'site',i+1

                ## ignore invariant sites at locus 1
                if mode=='aa' and len(removeItem(unique(column(alignmentA[tp],i)),['-','?']))==1:
                    pass
                elif mode=='nt' and len(removeItem2(unique(column(alignmentA[tp],i)),['A','C','T','G']))==1:
                    pass
                else:
                    ## proceed, locus 1 is polymorphic
                    for j in range(startFrom,len(alignmentB[tp][0])):
                        ## ignore invariant sites at locus 2
                        if mode=='aa' and len(removeItem(unique(column(alignmentB[tp],j)),['-','?']))==1 or (seg1==seg2 and i==j):
                            pass
                        elif mode=='nt' and len(removeItem2(unique(column(alignmentB[tp],j)),['A','C','T','G']))==1 or (seg1==seg2 and i==j):
                            pass
                        else:
                            ## proceed, locus 2 is polymorphic, add locus 2 to list of polymorphic site pairs, should save later on time
                            polymorphicPairs.append((i,j))

            ## iterate through potentially polymorphic site pairs
            for i,j in polymorphicPairs:
                ## reconstruct haplotypes
                checkHaplotypes={strain:(dictA[tp][strain][i],dictB[tp][strain][j]) for strain in strains[tp]}

                ## filter haplotypes, keep those that do not have invalid residues
                if mode=='aa':
                    haplotypes={strain:checkHaplotypes[strain] for strain in checkHaplotypes.keys() if '-' not in checkHaplotypes[strain] and '?' not in checkHaplotypes[strain]}
                elif mode=='nt':
                    haplotypes={strain:checkHaplotypes[strain] for strain in checkHaplotypes.keys() if checkHaplotypes[strain] in ['A','C','T','G']}

                ## some polymorphic loci still have invalid residues
                haps=unique(haplotypes.values())

                ## find alleles at locus 1 and locus 2
                poly1=[haps[q][0] for q in range(len(haps))]
                poly2=[haps[q][1] for q in range(len(haps))]
                poly1=unique(poly1)
                poly2=unique(poly2)

                ## count numbers of each allele
                poly1Count=tuple(column(haplotypes.values(),0).count(x) for x in poly1)
                poly2Count=tuple(column(haplotypes.values(),1).count(x) for x in poly2)

                ## if alleles at locus 1 are associated with invalid alleles at locus 2, things won't work
                if len(poly1)>1 and len(poly2)>1:

                    ## make flat contingency table of observed haplotypes
                    a=[poly1,poly2]
                    allhaps=[x for x in itertools.product(*a)]
                    hapCount = [haplotypes.values().count(x) for x in allhaps]
                    
                    ## assert that haplotype numbers are equal to counts of each allele
                    assert sum(poly1Count)==sum(poly2Count)==sum(hapCount)
                    total=sum(hapCount)

                    ## calculate frequencies of alleles and haplotypes
                    poly1Freq=[x/float(total) for x in poly1Count]
                    poly2Freq=[x/float(total) for x in poly2Count]
                    hapFreq=[x/float(total) for x in hapCount]

##                    # debugging extras
##                    print i,j
##                    print 'alleles at site 1:',poly1,'alleles at site 2:',poly2
##                    print 'all possible haplotypes:',allhaps
##                    print 'observed haplotypes:',hapCount,'observed alleles site 1:',poly1Count,'observed alleles site 2:',poly2Count
##                    print 'haplotype frequencies:',hapFreq,'allele frequencies at site 1:',poly1Freq,'allele frequencies at site 2:',poly2Freq
##                    print 'total number of haplotypes that passed filtering:',total

                    ## check whether these haplotype frequencies were seen before
                    if seen.has_key(str(hapCount)):
                        
                        ## give answer computed earlier
                        print seen[str(hapCount)].replace(',','\t'),'(PRECOMPUTED)'
                        print>>f,seen[str(hapCount)]

                    ## calculate everything
                    else:
                        
                        ## define minor allele frequency cutoff
                        freqCutoff=0.01

                        ## proceed if minor allele frequent enough
                        if min(poly1Freq)>=freqCutoff and min(poly2Freq)>=freqCutoff:
                            
                            ## calculate D' as well if loci biallelic
                            D='NaN'
                            Dprime='NaN'
                            if len(poly1)==2 and len(poly2)==2:

                                ## calculate D
                                D=hapFreq[0]-(poly1Freq[0]*poly2Freq[0])

                                ## calculate Dmax
                                if D>=0:
                                    Dmax=min([poly1Freq[0]*(1-poly2Freq[0]),(1-poly1Freq[0])*poly2Freq[0]])
                                else:
                                    Dmax=min([(1-poly1Freq[0])*(1-poly2Freq[0]),poly1Freq[0]*poly2Freq[0]])

                                ## normalize D by Dmax
                                Dprime=np.absolute(D/float(Dmax))
                            
                            ## calculate ChiSq
                            ChiSq=0
                            for q in range(len(poly1)):
                                for w in range(len(poly2)):
                                    observed=hapCount[q*len(poly2)+w]
                                    expected=poly1Freq[q]*poly2Freq[w]*total
                                    ChiSq+=((observed-expected)**2)/float(expected)

                            ## calculate ChiSqdf
                            ChiSqdf=ChiSq/float(total*(len(poly1)-1)*(len(poly2)-1))
                            seen[str(hapCount)]='%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],step,len(strains[tp]),total,i+1,j+1,len(poly1),len(poly2),D,Dprime,ChiSq,ChiSqdf)
                            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(timeline[tp],step,len(strains[tp]),total,i+1,j+1,len(poly1),len(poly2),D,Dprime,ChiSq,ChiSqdf)
                            print>>f,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(timeline[tp],step,len(strains[tp]),total,i+1,j+1,len(poly1),len(poly2),D,Dprime,ChiSq,ChiSqdf)

        ## close output file
        f.close()

        ## close alignment files
        handle1.close()
        handle2.close()

## done
print '\nDone!'
