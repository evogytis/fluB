#####################################################################################
##
##  This script analyzes posterior sets of trees produced by BEAST.
##  It is meant to accompany the manuscript "Reassortment between influenza B lineages and the emergence of a co-adapted PB1-PB2-HA gene complex"
##  and has only been tested with trees produced specifically for it.
##  You must have python, numpy and RSPR (available at http://kiwi.cs.dal.ca/Software/RSPR) installed on your machine to use this script.
##
##  You also require at least 4 tree files to run this script.
##  It compares trees made from alignments A and B (trees A and B), then normalizes this comparison by independent analyses of alignments A and B (trees A' and B').
##  
##
##  Usage:
##  go to OPTIONS at the bottom of the script and alter as required
##  if running in batchProcessing mode point to tree files (filenames and path variables)
##  otherwise the script should prompt you to select the files manually
##
##  open command line, go to folder where the script is located and type in:
##  "python Hydra.py" (ignore ")
##  
##  The following figures in the manuscript were made using this script:
##  Figures 8, S1-S7, S11-S12
##
##
##  Gytis Dudas
##  Institute of Evolutionary Biology
##  University of Edinburgh
##  EH9 3JT
##  Edinburgh, UK
##
#####################################################################################

import re
from datetime import datetime as dt
import time
import calendar
import numpy as np
import random
import glob
import os
import sys
import subprocess

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

def makeSliceList(slices,start,end,bins):
    """
    Return a list of dates spaced out regularly within a year.
    """
    for i in np.arange(start,end,1/float(bins)):
        slices.append([i])

def mostRecent(some_list,objectContainer):
    """
    Identify the most recent tip in the tree and its date.
    """
    full=[]
    month=[]
    year=[]
    ## Summarize tip names in regular expressions
    ## Enclose the bit of the sequence name containing the date in brackets
    ## Date should be in format yyyy-mm-dd
    dateRegex='[A-Z]\/[A-Za-z\/\-\_\.0-9]+\/[0-9]+\_([0-9\-]+)'
    for i in some_list:
        C=re.search('%s'%(dateRegex),i[1])
        if C is not None:
            try:
                yrs = int(C.group(1).split('-')[0])
                mons = int(C.group(1).split('-')[1])
                dds = int(C.group(1).split('-')[2])
                full.append([C.group(),(float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))])
            except IndexError:
                try:
                    u = int(calendar.monthrange(int(C.group(1).split('-')[0]),int(C.group(1).split('-')[1]))[1])
                    l = 1
                    yrs = int(C.group(1).split('-')[0])
                    mons = int(C.group(1).split('-')[1])
                    dds=int((u+l)/float(2))
                    month.append(C.group())
                except IndexError:
                    yrs = int(C.group(1))
                    mons=1
                    dds=1
                    year.append(C.group())
    try:
        assert (len(full)+len(month)+len(year)==len(some_list)),'Warning, some dates were not captured by regex!\nReview the \"mostRecent\" function'
        for i in objectContainer:
            if isinstance(i,leaf)==True:
                if i.name in column(full,0):
                    return [i.height,column(full,1)[index(column(full,0),i.name)[0]]]
    except ValueError:
        pass 
def make_tree(data,treeObject,objectContainer,tipList):
    """
    Scans tree string looking for BEAST formatting.
    Takes a given tree object, a list to store objects in and a list of taxa names.
    NOTE - only works with trees drawn from the posterior distribution of trees.
    """
    global i
    global trait_list
    global trait_names
    global partition_name
    stored_i=None
    
    while i < len(data)-2:
        assert (stored_i != i),'\nTree string unparseable\nStopped at %s\nstring region looks like this: %s'%(data[ind],data[ind-10:ind+10])
        stored_i=i

        if data[i] == ';':
            assert (len(tipList)*2-1 == len(objectContainer)),'\nTree string has been parsed incorrectly:\nexpected number of objects in tree %s\nobjects found in the tree string: %s'%(len(tips)*2-1,len(object_list))
            if len(trait_names)==0:
                trait_names=treeObject.cur_node.traits.keys()
                print 'Trait names: ',
                for traitNames in trait_names:
                    print traitNames,
                print '\n'
            
            print 'about to break',i,len(data)
            break
        
        cerberus=re.match('\[&([A-Za-z0-9\s\.\_\,\=\"\']+)\](:|;)',data[i:i+400])
        if cerberus is not None:
            traitComment=cerberus.group(1).split(',')
            for j in traitComment:
                traitName=j.split('=')[0]
                traitValue=j.split('=')[1].strip('"')
                treeObject.cur_node.traits[traitName]=traitValue
                if traitName not in column(trait_list,0):
                    trait_list.append([traitName,[traitValue]])
                if traitValue not in trait_list[index(column(trait_list,0),traitName)[0]][1]:
                    trait_list[index(column(trait_list,0),traitName)[0]][1].append(traitValue)
                    trait_list[index(column(trait_list,0),traitName)[0]][1]=sorted(trait_list[index(column(trait_list,0),traitName)[0]][1])
            i+=len(cerberus.group())-1

        if data[i] == '(':
            treeObject.add_node(i,objectContainer)
            i+=1
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
                treeObject.add_leaf(i,cerberus.group(1),objectContainer)
                i+=len(cerberus.group(1))
            elif data[i]=='(':
                make_tree(data,treeObject,objectContainer,tipList)
                
                
        if data[i] == ',':
            i+=1
            treeObject.move_up()
            treeObject.cur_node.secondChild=True
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
                treeObject.add_leaf(i,cerberus.group(1),objectContainer)
                i+=len(cerberus.group(1))

        if data[i] == ':':
            i+=1
            comment=re.match('\[&([A-Za-z0-9\s\_\-\.\,\=]+)\]([0-9\.\-E]+)',data[i:i+400])
            if comment is not None:
                treeObject.cur_node.branch=float(comment.group(2))
                treeObject.cur_node.height=float(comment.group(2))
                com_list=[]
                com_list=comment.group(1).split(',')

                for j in com_list:
                    if 'r' in j.split('=')[0] and j.split('=')[0] not in treeObject.cur_node.traits.keys():
                        if partition_name=='':
                            partition_name=com_list[0].split('=')[0]
                            print '\nPartition name: %s'%(partition_name)
                        treeObject.cur_node.rate=float(com_list[0].split('=')[1])
                    elif 'r' in j.split('=')[0] and j.split('=')[0] in treeObject.cur_node.traits.keys():
                        treeObject.cur_node.traitRates[j.split('=')[0]]=float(j.split('=')[1])

                robustCounts=re.search('b_u_S=([0-9\.]+),S=([0-9\.]+),N=([0-9\.]+),b_u_N=([0-9\.]+)',comment.group(1))
                if robustCounts is not None:
                    treeObject.cur_node.buS=float(robustCounts.group(1))
                    treeObject.cur_node.S=float(robustCounts.group(2))
                    treeObject.cur_node.N=float(robustCounts.group(3))
                    treeObject.cur_node.buN=float(robustCounts.group(4))


                robustCounts=re.search('S=([0-9\.]+),N=([0-9\.]+)',comment.group(1))
                if robustCounts is not None:
                    treeObject.cur_node.S=float(robustCounts.group(1))
                    treeObject.cur_node.N=float(robustCounts.group(2))

                i+=len(comment.group())

        if data[i] == ')':
            treeObject.move_up()
            i+=1

################################
## end of string parsing block
################################

## Node, leaf and tree objects used to represent the tree structure.
## Nodes connect to other nodes/leaves (tips) and always have two children

class node:
    def __init__(self):
        self.rate=None
        self.branch=None
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.leftChild=None
        self.rightChild=None
        self.secondChild=False
        self.traits={}
        self.traitRates={}
        self.index=None
        self.childHeight=None
        self.numChildren=0
        self.N=None
        self.S=None
        self.buS=None
        self.buN=None
        ## contains references to all tips of this node
        self.leaves=[]

class leaf:
    def __init__(self):
        self.name=None
        self.numName=None
        self.index=None
        self.branch=0
        self.absoluteTime=None
        self.rate=None
        self.height=0
        self.parent=None
        self.traits={}
        self.traitRates={}
        self.N=None
        self.S=None
        self.buS=None
        self.buN=None


class tree:
    def __init__(self):
        self.cur_node=node()
        self.cur_node.branch=0
        self.cur_node.height=0
        self.root=self.cur_node
        self.cur_node.index='Root'

    def move_up(self):
        """
        Go to current object's parent.
        """
        node=self.cur_node
        self.cur_node=node.parent

    def add_node(self,i,objectContainer):
        """
        Attach a node to current object as a first or second child.
        i refers to the index of the bracket defining the node in the tree string and is a unique identifier.
        """
        if self.cur_node.parent==self.root:
            self.cur_node.branch=0
            self.cur_node.height=0
        new_node=node()
        new_node.index=i
        objectContainer.append(new_node)
        
        if self.cur_node.secondChild==False:
            new_node.parent=self.cur_node
            self.cur_node.leftChild=new_node
            self.cur_node=new_node
        else:
            new_node.parent=self.cur_node
            self.cur_node.rightChild=new_node
            self.cur_node=new_node
        return objectContainer

    def add_leaf(self,i,name,objectContainer):
        """
        Attach a leaf to current object as a first or second child.
        i refers to the index of the tip number in the tree string and is a unique identifier.
        """
        new_leaf=leaf()
        new_leaf.index=i
        objectContainer.append(new_leaf)
        if self.cur_node.secondChild==False:
            new_leaf.name=name
            new_leaf.parent=self.cur_node
            self.cur_node.leftChild=new_leaf
            self.cur_node=new_leaf
        else:
            new_leaf.name=name
            new_leaf.parent=self.cur_node
            self.cur_node.rightChild=new_leaf
            self.cur_node=new_leaf
        return objectContainer

    def setAbsoluteTime(self,objectContainer,height,date):
        """
        Sets the absolute time of each object in the tree given the height of the most recent tip and its date.
        """
        for i in objectContainer:
            i.absoluteTime=date-height+i.height

    def renameTips(self,tipList,objectContainer):
        """
        Changes the name attribute of leaf objects to what the sequence name is.
        """
        for i in objectContainer:
            if isinstance(i,leaf)==True:
                i.numName=i.name
                i.name=tipList[index(column(tipList,0),int(i.name))[0]][1]
        return objectContainer
    

    def allTMRCAs(self,tipList,obList):
        """
        Identifies the TMRCAs of all pairs of tips.
        """
        allTMRCA={}
        ## Create a 2D dict for all pairs of tips
        for x in range(0,len(tipList)):
            allTMRCA[tipList[x][1]]={}
            for k in range(0,len(tipList)):
                allTMRCA[tipList[x][1]][tipList[k][1]]=0

        ## Iterate over nodes
        for y in obList:
            if isinstance(y,node):
                ## For all pairs of its children...
                for k in range(0,len(y.leaves)-1):
                    for x in range(k+1,len(y.leaves)):
                        ## If a pair of tips currently have a TMRCA older than the current node's time
                        ## set it to the node's time
                        if allTMRCA[y.leaves[k]][y.leaves[x]]<=y.absoluteTime:
                            allTMRCA[y.leaves[k]][y.leaves[x]]=y.absoluteTime
                            allTMRCA[y.leaves[x]][y.leaves[k]]=y.absoluteTime
        
        return allTMRCA                       

    def TMRCA(self,objectA,objectB,obList,tipList):
        """
        Identifies the most recent common ancestor of two objects.
        Accepts either leaf or node objects.
        """
        cur_nodeA=returnTip(obList,objectA,tipList)
        cur_nodeB=returnTip(obList,objectB,tipList)
        dictA={}
        dictB={}

        while cur_nodeA.parent!=None:
            cur_nodeA=cur_nodeA.parent
            dictA[cur_nodeA.index]=cur_nodeA

        while cur_nodeB.parent!=None:
            cur_nodeB=cur_nodeB.parent
            dictB[cur_nodeB.index]=cur_nodeB
            if dictA.has_key(cur_nodeB.index):
                key=cur_nodeB.index
                break

        return dictA[key]

    def timeToTMRCA(self,tip1,tip2):
        """
        Finds the total time taken to go from one tip to another.
        Takes two leaf objects.
        """
        cur_nodeA=tip1
        cur_nodeB=tip2
        dictA={}
        dictB={}


        while cur_nodeA.parent!=None:
            cur_nodeA=cur_nodeA.parent
            dictA[cur_nodeA.index]=cur_nodeA

        while cur_nodeB.parent!=None:
            cur_nodeB=cur_nodeB.parent
            dictB[cur_nodeB.index]=cur_nodeB
            if dictA.has_key(cur_nodeB.index):
                key=cur_nodeB.index
                break

        intersect=dictA[key]

        return -2*float(intersect.absoluteTime)+float(tip1.absoluteTime)+float(tip2.absoluteTime)

    def traverse_tree(self,maxHeight):
        """
        Pre-order tree traversal.
        Required to set the height of each object in the tree, to inform each node of its descendent tips.
        """
        cur_node=None
        cur_node=self.root.leftChild

        numChildren=0
        seen=[]
        highestTip=0
        highestTip=cur_node.branch
        height=0
        height=cur_node.branch
        root=False
        while root==False:
            if isinstance(cur_node,node):
                while cur_node.leftChild.index in seen and cur_node.rightChild.index in seen:
                    if cur_node.childHeight <= highestTip:
                        cur_node.childHeight = highestTip
                    elif cur_node.childHeight > highestTip:
                        highestTip = cur_node.childHeight

                    if cur_node.parent.index=='Root':
                        root=True
                        break
                    else:
                        cur_node.parent.numChildren+=cur_node.numChildren
                        cur_node.parent.leaves+=cur_node.leaves
                        cur_node.parent.leaves=sorted(cur_node.parent.leaves)
                        cur_node.height=height
                        height-=cur_node.branch
                        cur_node=cur_node.parent

                    if cur_node.index=='Root':
                        root=True
                        break

                if cur_node.leftChild.index in seen or cur_node.index in seen:
                    height+=cur_node.rightChild.branch
                    cur_node.childHeight=highestTip
                    highestTip=0
                    cur_node=cur_node.rightChild
                    cur_node.height=height
                else:
                    height+=cur_node.leftChild.branch
                    seen.append(cur_node.index)
                    cur_node=cur_node.leftChild
                    cur_node.height=height

            if isinstance(cur_node,leaf):
                cur_node.parent.numChildren+=1
                cur_node.parent.leaves.append(cur_node.name)
                cur_node.parent.leaves=sorted(cur_node.parent.leaves)
                if height >= maxHeight:
                    maxHeight=height
                seen.append(cur_node.index)
                cur_node.height=height
                highestTip=height
                height-=cur_node.branch
                cur_node=cur_node.parent
        return maxHeight

    def returnTreeString(self):
        """
        Returns a tree string containing the topology of the tree.
        Used as input for RSPR.
        """
        outString=[]
        cur_node=self.root.leftChild
        seen=[]
        root=False
        while root==False:
            if isinstance(cur_node,node):
                while cur_node.leftChild.index in seen and cur_node.rightChild.index in seen:
                    if cur_node.parent.index=='Root':
                        root=True
                        outString.append(');')
                        break
                    else:
                        outString.append(')')
                        cur_node=cur_node.parent

                    if cur_node.index=='Root':
                        root=True
                        outString.append(');')
                        break

                if cur_node.leftChild.index in seen or cur_node.index in seen:
                    cur_node=cur_node.rightChild
                    if root!=True:
                        outString.append(',')
                else:
                    outString.append('(')
                    seen.append(cur_node.index)
                    cur_node=cur_node.leftChild

            if isinstance(cur_node,leaf):
                outString.append(str(cur_node.numName))
                seen.append(cur_node.index)
                cur_node=cur_node.parent

        return ''.join(outString)

###################################################

global object_list1
object_list1=[]
global maxHeight1
maxHeight1=0

global object_list2
object_list2=[]
global maxHeight2
maxHeight2=0

global object_list3
object_list3=[]
global maxHeight3
maxHeight3=0

global object_list4
object_list4=[]
global maxHeight4
maxHeight4=0

global l1
l1=None
global l2
l2=None
global l3
l3=None
global l4
l4=None
global i
i=0

from itertools import izip

## parses all 4 trees
def parseTreeFile(treefile1,treefile2,repfile1,repfile2,modes,additionalInfo):
    global i
    treecount=0
    print '\nUsing files:%s\t%s'%(treefile1.name.split('/')[-1],treefile2.name.split('/')[-1])
    print 'Normalizing by:%s\t%s'%(repfile1.name.split('/')[-1],repfile2.name.split('/')[-1])
    mode='_'.join(modes)
    
    output = open('/Users/admin/Documents/Viral sequences/InfB reassortment/Analyses/1000 trees/'+treefile1.name.split('_')[1]+'_'+treefile2.name.split('_')[1]+'.'+additionalInfo+'.txt','w+')
    print 'Output in:',output.name
    print 'Burnin: %s'%(burnin)

    ## Header information:
    # AB - comparison between tree A and tree B
    # AA and BB - comparisons between independent trees of the same alignment
    # Norm - either 1 or 2 is the normalized value (manuscript uses normalization 1)
    # aSPR - SPR distances
    # TMRCA - absolute TMRCA deviations
    # TT - deviations between tip-tip distances
    # Alt - TMRCA deviations (not absolute)
    labelSet1='\tAB_aSPR\tAA_aSPR\tBB_aSPR\tNormSPR_1\tNormSPR_2'
    labelSet2='\tAB_sumTMRCA\tAA_sumTMRCA\tBB_sumTMRCA\tNorm_sumTMRCA_1\tNorm_sumTMRCA_2'
    labelSet3='\tAB_meanTMRCA\tAA_meanTMRCA\tBB_meanTMRCA\tNorm_meanTMRCA_1\tNorm_meanTMRCA_2'
    labelSet4='\tAB_stdevTMRCA\tAA_stdevTMRCA\tBB_stdevTMRCA\tNorm_stdevTMRCA_1\tNorm_stdevTMRCA_2'
    labelSet5='\tsharedClades_AB\tsharedClades_AA\tsharedClades_BB\tNorm_sharedClades_1\tNorm_sharedClades_2'
    labelSet6='\tAB_sumTT\tAA_sumTT\tBB_sumTT\tNorm_sumTT_1\tNorm_sumTT_2'
    labelSet7='\tAB_meanTT\tAA_meanTT\tBB_meanTT\tNorm_meanTT_1\tNorm_meanTT_2'
    labelSet8='\tAB_stdevTT\tAA_stdevTT\tBB_stdevTT\tNorm_stdevTT_1\tNorm_stdevTT_2'
    labelSet9='\tAltAB_sumTMRCA\tAltAA_sumTMRCA\tAltBB_sumTMRCA\tNorm_sumTMRCA_1\tNorm_sumTMRCA_2'
    labelSet10='\tAltAB_meanTMRCA\tAltAA_meanTMRCA\tAltBB_meanTMRCA\tAltNorm_meanTMRCA_1\tAltNorm_meanTMRCA_2'
    labelSet11='\tAltAB_stdevTMRCA\tAltAA_stdevTMRCA\tAltBB_stdevTMRCA\tAltNorm_stdevTMRCA_1\tAltNorm_stdevTMRCA_2'
    labelSet12='\tAB_minTMRCA\tAA_minTMRCA\tBB_minTMRCA\tNorm_minTMRCA_1\tNorm_minTMRCA_2'
    labelSet13='\tAB_maxTMRCA\tAA_maxTMRCA\tBB_maxTMRCA\tNorm_maxTMRCA_1\tNorm_maxTMRCA_2'
    
    print>>output,'states%s%s%s%s%s%s%s%s%s%s%s%s%s\ttimeTaken\texact indicator'%(labelSet1,labelSet2,labelSet3,labelSet4,labelSet5,labelSet6,labelSet7,labelSet8,labelSet9,labelSet10,labelSet11,labelSet12,labelSet13)
    treeThreshold=0
    treeThresholdIncrement=1

    global trait_list
    trait_list=[]
    global trait_slices
    trait_slices=[]
    global trait_names
    trait_names=[]
    global partition_name
    partition_name=''
    global slice_list
    slice_list=[]

    processingRate=[]

    global object_list1
    object_list1=[]
    global maxHeight1
    maxHeight1=0

    global object_list2
    object_list2=[]
    global maxHeight2
    maxHeight2=0

    global object_list3
    object_list3=[]
    global maxHeight3
    maxHeight3=0

    global object_list4
    object_list4=[]
    global maxHeight4
    maxHeight4=0
    global i

    global l1
    global l2
    global r1
    global r2
    l1=None
    l2=None
    r1=None
    r2=None
    
    global tips1
    tips1=[]
    tipNum1=0

    global tips2
    tips2=[]
    tipNum2=0

    global tips3
    tips3=[]
    tipNum3=0

    global tips4
    tips4=[]
    tipNum4=0
    
    begin=dt.now()
    conditionalCounter=0

    taxonlist=False
    plate=True
    
    for line1, line2, line3, line4 in izip(treefile1, treefile2, repfile1, repfile2):
        if plate==True:
            cerberus1=re.search('Dimensions ntax\=([0-9]+)\;',line1)
            cerberus2=re.search('Dimensions ntax\=([0-9]+)\;',line2)
            cerberus3=re.search('Dimensions ntax\=([0-9]+)\;',line3)
            cerberus4=re.search('Dimensions ntax\=([0-9]+)\;',line4)
            if cerberus1 is not None and cerberus2 is not None and cerberus3 is not None and cerberus4 is not None:
                tipNum1=int(cerberus1.group(1))
                tipNum2=int(cerberus2.group(1))
                tipNum3=int(cerberus3.group(1))
                tipNum4=int(cerberus4.group(1))
                print 'N tips:',tipNum1,tipNum2,tipNum3,tipNum4
            if 'Translate' in line1 and 'Translate' in line2 and 'Translate' in line3 and 'Translate' in line4:
                taxonlist=True
        
            if taxonlist==True and ';' not in line1 and 'Translate' not in line1 and ';' not in line2 and 'Translate' not in line2 and ';' not in line3 and 'Translate' not in line3 and ';' not in line4 and 'Translate' not in line4:
                try:
                    tip_index=int(line1.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line1.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips1:
                        tips1.append([tip_index,tip])

                    tip_index=int(line2.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line2.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips2:
                        tips2.append([tip_index,tip])

                    tip_index=int(line3.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line3.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips3:
                        tips3.append([tip_index,tip])

                    tip_index=int(line4.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line4.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips4:
                        tips4.append([tip_index,tip])
                except ValueError:
                    pass

        if 'tree STATE_0' in line1 and 'tree STATE_0' in line2 and 'tree STATE_0' in line3 and 'tree STATE_0' in line4:
            plate=False
            assert (tipNum1 == len(tips1)),'Tree 1 expected number of tips: %s\nNumber of tips found: %s'%(tipNum1,len(tips1))
            assert (tipNum2 == len(tips2)),'Tree 2 expected number of tips: %s\nNumber of tips found: %s'%(tipNum2,len(tips2))
            assert (tipNum3 == len(tips3)),'Rep 1 expected number of tips: %s\nNumber of tips found: %s'%(tipNum3,len(tips3))
            assert (tipNum4 == len(tips4)),'Rep 2 expected number of tips: %s\nNumber of tips found: %s'%(tipNum4,len(tips4))
            treecount=0
            print '\nTree 1 Number of taxa: %s'%(len(tips1))
            print 'Tree 2 Number of taxa: %s'%(len(tips2))
            print 'Rep 1 number of taxa: %s'%(len(tips3))
            print 'Rep 1 number of taxa: %s'%(len(tips4))

        cerberus1=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line1)
        cerberus2=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line2)
        cerberus3=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line3)
        cerberus4=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line4)
        if cerberus1 is not None and cerberus2 is not None and cerberus3 is not None and cerberus4 is not None:
            outputcontent=''

            ## Keeps user informed about progress
            if treecount>=treeThreshold:
                timeTakenSoFar=dt.now()-begin
                timeElapsed=float(divmod(timeTakenSoFar.total_seconds(),60)[0]+(divmod(timeTakenSoFar.total_seconds(),60)[1])/float(60))
                timeRate=float(divmod(timeTakenSoFar.total_seconds(),60)[0]*60+divmod(timeTakenSoFar.total_seconds(),60)[1])/float(treecount+1)
                processingRate.append(timeRate)
                ETA=(np.mean(processingRate)*(1000-treecount))/float(60)/float(60)
                if treecount==0:
                    pass
                else:
                    sys.stdout.write('\r')
                    sys.stdout.write("[%-40s] %d%% \tTrees: %d\tTime elapsed: %.2f minutes\tETA: %.2f hours (%.4f seconds/tree)" % ('='*(treecount/25),treecount/10,treecount,timeElapsed,ETA,processingRate[-1]))
                    sys.stdout.flush()
                treeThreshold+=treeThresholdIncrement

            if int(cerberus1.group(1))==0:
                object_list1=[]
                maxHeight1=0
                l1=None
                l1=tree()

                start1=len(cerberus1.group())
                treestring1=str(line1[start1:].strip('\n'))
                i=0
                make_tree(treestring1,l1,object_list1,tips1)

                object_list2=[]
                maxHeight2=0
                l2=None
                l2=tree()

                start2=len(cerberus2.group())
                treestring2=str(line2[start2:].strip('\n'))
                i=0
                make_tree(treestring2,l2,object_list2,tips2)
                
                object_list3=[]
                maxHeight3=0
                r1=None
                r1=tree()

                start3=len(cerberus3.group())
                repstring1=str(line3[start3:].strip('\n'))
                i=0
                make_tree(repstring1,r1,object_list3,tips3)

                object_list4=[]
                maxHeight4=0
                r2=None
                r2=tree()

                start4=len(cerberus4.group())
                repstring2=str(line4[start4:].strip('\n'))
                i=0
                make_tree(repstring2,r2,object_list4,tips4)

                print '\ntree 1',id(l1)
                print 'tree 2',id(l2)
                print 'replicate 1',id(r1)
                print 'replicate 2',id(r2)
                
                #######################################################
                headerstring='states'

            if int(cerberus1.group(1)) >= burnin and int(cerberus2.group(1)) >= burnin and int(cerberus3.group(1)) >= burnin and int(cerberus4.group(1)) >= burnin:

                ## Parse all 4 (A, B, A' and B') trees
                i=0
                object_list1=[]
                maxHeight1=0
                l1=None
                l1=tree()
                start1=len(cerberus1.group())
                treestring1=str(line1[start1:])
                make_tree(treestring1,l1,object_list1,tips1)
                l1.renameTips(tips1,object_list1)
                l1.traverse_tree(maxHeight1)
                calTipHeight,calTipDate=mostRecent(tips1,object_list1)
                l1.setAbsoluteTime(object_list1,calTipHeight,calTipDate)

                i=0
                object_list2=[]
                maxHeight2=0
                l2=None
                l2=tree()
                start2=len(cerberus2.group())
                treestring2=str(line2[start2:])
                make_tree(treestring2,l2,object_list2,tips2)
                l2.renameTips(tips2,object_list2)
                l2.traverse_tree(maxHeight2)
                calTipHeight,calTipDate=mostRecent(tips2,object_list2)
                l2.setAbsoluteTime(object_list2,calTipHeight,calTipDate)

                i=0
                object_list3=[]
                maxHeight3=0
                start3=len(cerberus3.group())
                repstring1=str(line3[start3:])
                r1=None
                r1=tree()
                make_tree(repstring1,r1,object_list3,tips3)
                r1.renameTips(tips3,object_list3)
                r1.traverse_tree(maxHeight3)
                calTipHeight,calTipDate=mostRecent(tips3,object_list3)
                r1.setAbsoluteTime(object_list3,calTipHeight,calTipDate)

                i=0
                object_list4=[]
                maxHeight4=0
                start4=len(cerberus4.group())
                repstring2=str(line4[start4:])
                r2=None
                r2=tree()
                make_tree(repstring2,r2,object_list4,tips4)
                r2.renameTips(tips4,object_list4)
                r2.traverse_tree(maxHeight4)
                calTipHeight,calTipDate=mostRecent(tips4,object_list4)
                r2.setAbsoluteTime(object_list4,calTipHeight,calTipDate)

                ## Find TMRCAs for all pairs of tips
                tmrcas1=l1.allTMRCAs(tips1,object_list1)
                tmrcas2=l2.allTMRCAs(tips2,object_list2)
                tmrcas3=r1.allTMRCAs(tips3,object_list3)
                tmrcas4=r2.allTMRCAs(tips4,object_list4)

                ## will contain A:B, A:A' and B:B' TMRCA deviations
                deviationsAlt=[]
                deviations1Alt=[]
                deviations2Alt=[]

                ## will contain ABSOLUTE A:B, A:A' and B:B' TMRCA deviations
                deviations=[]
                deviations1=[]
                deviations2=[]

                tipDeviations=[]
                tipDeviations1=[]
                tipDeviations2=[]

                ## Calculates the number of clades shared between trees
                hcd=0
                for n1 in object_list1:
                    if isinstance(n1,node):
                        for n2 in object_list2:
                            if isinstance(n2,node):
                                if n1.leaves==n2.leaves:
                                    hcd+=1


                hcdAA=0
                for n1 in object_list1:
                    if isinstance(n1,node):
                        for n2 in object_list3:
                            if isinstance(n2,node):
                                if n1.leaves==n2.leaves:
                                    hcdAA+=1

                hcdBB=0
                for n1 in object_list2:
                    if isinstance(n1,node):
                        for n2 in object_list4:
                            if isinstance(n2,node):
                                if n1.leaves==n2.leaves:
                                    hcdBB+=1
                tipDates1={}
                tipDates2={}
                tipDatesr1={}
                tipDatesr2={}

                for q in object_list1:
                    if isinstance(q,leaf)==True:
                        tipDates1[q.name]=q.absoluteTime

                for q in object_list2:
                    if isinstance(q,leaf)==True:
                        tipDates2[q.name]=q.absoluteTime

                for q in object_list3:
                    if isinstance(q,leaf)==True:
                        tipDatesr1[q.name]=q.absoluteTime

                for q in object_list4:
                    if isinstance(q,leaf)==True:
                        tipDatesr2[q.name]=q.absoluteTime
                    
                
                for p in range(len(tips1)-1):
                    tipA1=tipDates1[tips1[p][1]]
                    tipA2=tipDates2[tips1[p][1]]
                    tiprA1=tipDatesr1[tips1[p][1]]
                    tiprA2=tipDatesr2[tips1[p][1]]
                    for u in range(p+1,len(tips1)):
                        
                        tipB1=tipDates1[tips1[u][1]]
                        tipB2=tipDates2[tips1[u][1]]
                        tiprB1=tipDatesr1[tips1[u][1]]
                        tiprB2=tipDatesr2[tips1[u][1]]

                        Att=tipA1+tipB1-(2*tmrcas1[tips1[p][1]][tips1[u][1]])
                        Btt=tipA2+tipB2-(2*tmrcas2[tips1[p][1]][tips1[u][1]])
                        rAtt=tiprA1+tiprB1-(2*tmrcas3[tips1[p][1]][tips1[u][1]])
                        rBtt=tiprA2+tiprB2-(2*tmrcas4[tips1[p][1]][tips1[u][1]])

                        tipDeviations.append(np.absolute(Att-Btt))
                        tipDeviations1.append(np.absolute(Att-rAtt))
                        tipDeviations2.append(np.absolute(Btt-rBtt))

                        ## TMRCA deviation
                        mdAlt=tmrcas1[tips1[p][1]][tips1[u][1]]-tmrcas2[tips1[p][1]][tips1[u][1]]
                        daaAlt=tmrcas1[tips1[p][1]][tips1[u][1]]-tmrcas3[tips1[p][1]][tips1[u][1]]
                        dbbAlt=tmrcas2[tips1[p][1]][tips1[u][1]]-tmrcas4[tips1[p][1]][tips1[u][1]]

                        deviations1Alt.append(daaAlt)
                        deviations2Alt.append(dbbAlt)
                        deviationsAlt.append(mdAlt)

                        ## ABSOLUTE TMRCA deviation
                        md=np.absolute(mdAlt)
                        daa=np.absolute(daaAlt)
                        dbb=np.absolute(dbbAlt)
                        deviations1.append(daa)
                        deviations2.append(dbb)
                        deviations.append(md)

                # normalized number of clades shared
                hcdNorm1=2*float(hcd)/(float(hcdAA)+float(hcdBB))
                hcdNorm2=float(hcd)/(float(hcdAA)*float(hcdBB))**0.5


                # sum of delta TMRCA (not absolute deviation)
                AltsumTMRCA=sum(deviationsAlt)
                AltsumTMRCA1=sum(deviations1Alt)
                AltsumTMRCA2=sum(deviations2Alt)
                AltnormSumdevA=float(AltsumTMRCA1+AltsumTMRCA2)/float(AltsumTMRCA*2)
                AltnormSumdevB=(float(AltsumTMRCA1*AltsumTMRCA2))**2/float(AltsumTMRCA)

                # mean of delta TMRCA (not absolute deviation)
                AltmeanTMRCA=np.mean(deviationsAlt)
                AltmeanTMRCA1=np.mean(deviations1Alt)
                AltmeanTMRCA2=np.mean(deviations2Alt)
                AltnormMeandevA=float(AltmeanTMRCA1+AltmeanTMRCA2)/float(AltmeanTMRCA*2)
                AltnormMeandevB=(float(AltmeanTMRCA1*AltmeanTMRCA2))**2/float(AltmeanTMRCA)

                # st dev of delta TMRCA (not absolute deviation)
                AltstTMRCA=np.std(deviationsAlt)
                AltstTMRCA1=np.std(deviations1Alt)
                AltstTMRCA2=np.std(deviations2Alt)
                AltnormStDevA=float(AltstTMRCA1+AltstTMRCA2)/float(AltstTMRCA*2)
                AltnormStDevB=(float(AltstTMRCA1*AltstTMRCA2))**2/float(AltstTMRCA)

                # min of delta TMRCA (absolute)
                minTMRCA=min(deviations)
                minTMRCA1=min(deviations1)
                minTMRCA2=min(deviations2)
                minDevA=float(minTMRCA1+minTMRCA2)/float(minTMRCA*2)
                minDevB=(float(minTMRCA1*minTMRCA2))**2/float(minTMRCA)

                # max of delta TMRCA (absolute)
                maxTMRCA=max(deviations)
                maxTMRCA1=max(deviations1)
                maxTMRCA2=max(deviations2)
                maxDevA=float(maxTMRCA1+maxTMRCA2)/float(maxTMRCA*2)
                maxDevB=(float(maxTMRCA1*maxTMRCA2))**2/float(maxTMRCA)
                
                # sum of delta TMRCA (absolute)
                sumTMRCA=sum(deviations)
                sumTMRCA1=sum(deviations1)
                sumTMRCA2=sum(deviations2)
                normSumdevA=float(sumTMRCA1+sumTMRCA2)/float(sumTMRCA*2)
                normSumdevB=(float(sumTMRCA1*sumTMRCA2))**2/float(sumTMRCA)

                ## Used in the manuscript
                # mean of delta TMRCA (absolute)
                meanTMRCA=np.mean(deviations)
                meanTMRCA1=np.mean(deviations1)
                meanTMRCA2=np.mean(deviations2)
                normMeandevA=float(meanTMRCA1+meanTMRCA2)/float(meanTMRCA*2)
                normMeandevB=(float(meanTMRCA1*meanTMRCA2))**2/float(meanTMRCA)

                # st dev of delta TMRCA
                stTMRCA=np.std(deviations)
                stTMRCA1=np.std(deviations1)
                stTMRCA2=np.std(deviations2)
                normStDevA=float(stTMRCA1+stTMRCA2)/float(stTMRCA*2)
                normStDevB=(float(stTMRCA1*stTMRCA2))**2/float(stTMRCA)

                # sum of tip-tip distances
                sumtt=sum(tipDeviations)
                sumtt1=sum(tipDeviations1)
                sumtt2=sum(tipDeviations2)
                normSumttA=float(sumtt1+sumtt2)/float(sumtt*2)
                normSumttB=(float(sumtt1*sumtt2))**2/float(sumtt)

                # mean of tip-tip distances
                meantt=np.mean(tipDeviations)
                meantt1=np.mean(tipDeviations1)
                meantt2=np.mean(tipDeviations2)
                normMeanttdevA=float(meantt1+meantt2)/float(meantt*2)
                normMeanttdevB=(float(meantt1*meantt2))**2/float(meantt)

                # st dev of tip-tip distances
                sttt=np.std(tipDeviations)
                sttt1=np.std(tipDeviations1)
                sttt2=np.std(tipDeviations2)
                normStDevttA=float(sttt1+sttt2)/float(sttt*2)
                normStDevttB=(float(sttt1*sttt2))**2/float(sttt)

                # produce tree string of each tree for analysis in RSPR
                string1=l1.returnTreeString()
                string2=l2.returnTreeString()
                stringr1=r1.returnTreeString()
                stringr2=r2.returnTreeString()

                ## For much larger trees, the -split_approx option will compute an
                ## exponential time approximation of the distance that is exact for
                ## small distances and generally within a few percent of the optimal
                ## distance otherwise.

                ## Choose approximation:
                ## 0 == exact
                ## 1 == approx_split
                ## 2 == approx
                ## 3 == don't use RSPR
                approximation=2

                ## Path to RSPR
                rsprPath='/Applications/rspr_1_2_1/'

                ## Uses RSPR to calculate SPR distance
                exactIndicator=0
                splitApproximationN=70
                if approximation==1:
                    mainCommand='echo \"%s\n%s\" | %srspr -split_approx %s'%(string1,string2,rsprPath,splitApproximationN)
                    normalize1='echo \"%s\n%s\" | %srspr -split_approx %s'%(string1,stringr1,rsprPath,splitApproximationN)
                    normalize2='echo \"%s\n%s\" | %srspr -split_approx %s'%(string2,stringr2,rsprPath,splitApproximationN)

                    oot=subprocess.check_output([mainCommand],shell=True)
                    ootCerberus=re.search('total exact drSPR\=([0-9]+)',oot.split('\n')[-2])

                    exactIndicator=1
                    outList=oot.split('\n')
                    for o in range(len(outList)):
                        if 'too large' in outList[o]:
                            exactIndicator=0

                    normOot1=subprocess.check_output([normalize1],shell=True)
                    norm1Cerberus=re.search('total exact drSPR\=([0-9]+)',normOot1.split('\n')[-2])

                    normOot2=subprocess.check_output([normalize2],shell=True)
                    norm2Cerberus=re.search('total exact drSPR\=([0-9]+)',normOot2.split('\n')[-2])

                elif approximation==2:
                    mainCommand='echo \"%s\n%s\" | %srspr -approx'%(string1,string2,rsprPath)
                    normalize1='echo \"%s\n%s\" | %srspr -approx'%(string1,stringr1,rsprPath)
                    normalize2='echo \"%s\n%s\" | %srspr -approx'%(string2,stringr2,rsprPath)

                    exactIndicator=0

                    oot=subprocess.check_output([mainCommand],shell=True)
                    ootCerberus=re.search('approx drSPR\=([0-9]+)',oot.split('\n')[-3])

                    AB_clusters2=[x for x in oot.split('\n')[-4].split(' ')[2:] if x.count(',')>0]

                    normOot1=subprocess.check_output([normalize1],shell=True)
                    norm1Cerberus=re.search('approx drSPR\=([0-9]+)',normOot1.split('\n')[-3])

                    AA_clusters2=[x for x in normOot1.split('\n')[-4].split(' ')[2:] if x.count(',')>0]

                    normOot2=subprocess.check_output([normalize2],shell=True)
                    norm2Cerberus=re.search('approx drSPR\=([0-9]+)',normOot2.split('\n')[-3])

                    BB_clusters2=[x for x in normOot2.split('\n')[-4].split(' ')[2:] if x.count(',')>0]

                elif approximation==0:
                    mainCommand='echo \"%s\n%s\" | %srspr'%(string1,string2,rsprPath)
                    normalize1='echo \"%s\n%s\" | %srspr'%(string1,stringr1,rsprPath)
                    normalize2='echo \"%s\n%s\" | %srspr'%(string2,stringr2,rsprPath)

                    exactIndicator=1
                    
                    oot=subprocess.check_output([mainCommand],shell=True)
                    ootCerberus=re.search('total exact drSPR\=([0-9]+)',oot.split('\n')[-2])

                    normOot1=subprocess.check_output([normalize1],shell=True)
                    norm1Cerberus=re.search('total exact drSPR\=([0-9]+)',normOot1.split('\n')[-2])

                    normOot2=subprocess.check_output([normalize2],shell=True)
                    norm2Cerberus=re.search('total exact drSPR\=([0-9]+)',normOot2.split('\n')[-2])
                else:
                    mm=1
                    n1=1
                    n2=1
                    normA=1
                    normB=1

                if approximation!=3:
                    if ootCerberus is None:
                        print 'Output not caught by regex:'
                        print oot.split('\n')
##                    mm=float(ootCerberus.group(1))
##                    n1=float(norm1Cerberus.group(1))
##                    n2=float(norm2Cerberus.group(1))
##                    normA=float(n1+n2)/float(mm*2)
##                    normB=float((n1*n2)**0.5)/float(mm)

                    mm=float(len(AB_clusters2))
                    n1=float(len(AA_clusters2))
                    n2=float(len(BB_clusters2))
                    normA=float(n1+n2)/float(mm*2)
                    normB=float((n1*n2)**0.5)/float(mm)

                ################################################################################
                ## Output to file
                set1='\t%s\t%s\t%s\t%s\t%s'%(int(mm),int(n1),int(n2),normA,normB)
                set2='\t%s\t%s\t%s\t%s\t%s'%(sumTMRCA,sumTMRCA1,sumTMRCA2,normSumdevA,normSumdevB)
                set3='\t%s\t%s\t%s\t%s\t%s'%(meanTMRCA,meanTMRCA1,meanTMRCA2,normMeandevA,normMeandevB)
                set4='\t%s\t%s\t%s\t%s\t%s'%(stTMRCA,stTMRCA1,stTMRCA2,normStDevA,normStDevB)
                set5='\t%s\t%s\t%s\t%s\t%s'%(hcd,hcdAA,hcdBB,hcdNorm1,hcdNorm2)
                set6='\t%s\t%s\t%s\t%s\t%s'%(sumtt,sumtt1,sumtt2,normSumttA,normSumttB)
                set7='\t%s\t%s\t%s\t%s\t%s'%(meantt,meantt1,meantt2,normMeanttdevA,normMeanttdevB)
                set8='\t%s\t%s\t%s\t%s\t%s'%(sttt,sttt1,sttt2,normStDevttA,normStDevttB)
                set9='\t%s\t%s\t%s\t%s\t%s'%(AltsumTMRCA,AltsumTMRCA1,AltsumTMRCA2,AltnormSumdevA,AltnormSumdevB)
                set10='\t%s\t%s\t%s\t%s\t%s'%(AltmeanTMRCA,AltmeanTMRCA1,AltmeanTMRCA2,AltnormMeandevA,AltnormMeandevB)
                set11='\t%s\t%s\t%s\t%s\t%s'%(AltstTMRCA,AltstTMRCA1,AltstTMRCA2,AltnormStDevA,AltnormStDevB)
                set12='\t%s\t%s\t%s\t%s\t%s'%(minTMRCA,minTMRCA1,minTMRCA2,minDevA,minDevB)
                set13='\t%s\t%s\t%s\t%s\t%s'%(maxTMRCA,maxTMRCA1,maxTMRCA2,maxDevA,maxDevB)

                print>>output,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s\t%s\t%s'%(cerberus1.group(1),set1,set2,set3,set4,set5,set6,set7,set8,set9,set10,set11,set12,set13,processingRate[-1],exactIndicator)
                
                ################################################################################
                treecount+=1
                del l1
                del l2
                del r1
                del r2

            if 'End;' in line1 and 'End;' in line2:
                print '\nEnd of files!'

    output.close()
    end=dt.now()
    timeTaken=end-begin

    print '\n\nTime taken: %s minutes %s seconds'%(divmod(timeTaken.total_seconds(),60)[0],divmod(timeTaken.total_seconds(),60)[1])
    print 'Number of trees analyzed: %s'%(treecount)
    print '\n#################################################################################################################################'

############################################################################################################
## OPTIONS
############################################################################################################

## If batch processing set to True - give path and filenames
# otherwise Tkinter will prompt you to select the files manually
batchProcessing=True

## Define path where the subsampled tree files are
path='/Users/admin/Documents/Viral sequences/InfB reassortment/Analyses/1000 trees/'

## Provide file names to iterate over
filenames=['InfB_PB1t_ALLs1.trees.subsample_1.txt','InfB_PB2t_ALLs1.trees.subsample_1.txt','InfB_PAt_ALLs1.trees.subsample_1.txt','InfB_HAt_ALLs1.trees.subsample_1.txt','InfB_NPt_ALLs1.trees.subsample_1.txt','InfB_NAt_ALLs1.trees.subsample_1.txt','InfB_M1t_ALLs1.trees.subsample_1.txt','InfB_NS1t_ALLs1.trees.subsample_1.txt']
#filenames=['InfB_PB1t_ALLs2.trees.subsample_1.txt','InfB_PB2t_ALLs2.trees.subsample_1.txt','InfB_HAt_ALLs2.trees.subsample_1.txt']
#filenames=['InfB_PB1t_ALLs3.trees.subsample_1.txt','InfB_PB2t_ALLs3.trees.subsample_1.txt','InfB_HAt_ALLs3.trees.subsample_1.txt']


## Define output file suffix
additionalInfo='Approx1NormBy2.SPRcutoff2'
modes=['']

global burnin
## Define burnin as the number of the state from which to begin analysis.
burnin=0

############################################################################################################

if batchProcessing==False:
    import Tkinter, tkFileDialog
    root=Tkinter.Tk()
    file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose the first tree file')
    if file != None:
            treefile1 = file
    root.withdraw()
    print treefile1.name
    file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose a second tree file')
    if file != None:
            treefile2 = file
    root.withdraw()
    print treefile2.name

    file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose replicate of first tree')
    if file != None:
            treefile3 = file
    root.withdraw()
    print treefile3.name
    file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose replicate of second tree')
    if file != None:
            treefile4 = file
    root.withdraw()
    print treefile4.name
    total_begin=dt.now()
    parseTreeFile(treefile1,treefile2,treefile3,treefile4,modes,additionalInfo)
    total_end=dt.now()
    total_timeTaken=total_end-total_begin
    print '\nTotal time taken: %s minutes %s seconds'%(divmod(total_timeTaken.total_seconds(),60)[0],divmod(total_timeTaken.total_seconds(),60)[1])

else:
    total_begin=dt.now()
    filecount=1

    ## Iterate through files provided
    for q in range(len(filenames)):
        for w in range(q+1,len(filenames)):
            ## Open tree A and tree B
            treefile1=open(path+filenames[q],'r')
            treefile2=open(path+filenames[w],'r')

            ## Open tree A' and tree B'
            # tree file 'InfB_PB1t_ALLs1.trees.subsample_1.txt' will be normalized by tree file 'InfB_PB1t_ALLs2.trees.subsample_1.txt'
            # and 'InfB_PB2t_ALLs1.trees.subsample_1.txt' by 'InfB_PB2t_ALLs2.trees.subsample_1.txt'
            treefile3=open(path+filenames[q].replace('ALLs1','ALLs2'),'r')
            treefile4=open(path+filenames[w].replace('ALLs1','ALLs2'),'r')
            parseTreeFile(treefile1,treefile2,treefile3,treefile4,modes,additionalInfo)

    total_end=dt.now()
    total_timeTaken=total_end-total_begin
    print '\nTotal time taken: %s minutes %s seconds'%(divmod(total_timeTaken.total_seconds(),60)[0],divmod(total_timeTaken.total_seconds(),60)[1])

######################################################

print '\n>>>>> Done! <<<<<\n'
