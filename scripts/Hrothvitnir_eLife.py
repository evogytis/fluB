#####################################################################################
##
##  This script analyzes posterior sets of trees produced by BEAST.
##  It is meant to accompany the manuscript "Reassortment between influenza B lineages and the emergence of a co-adapted PB1-PB2-HA gene complex"
##  and has only been tested with trees produced specifically for it.
##  You must have python and numpy installed on your machine to use this script.
##
##  Usage:
##  open command line, go to folder where the script is located and type in:
##  "python Hrothvitnir_eLife.py -i <<path to .trees.txt file>> -m <<mode - X, diversityOt, FstOt or stateTime>> 1> <<path to output file>>" (ignore <<, >> and ")
##  
##  The following figures in the manuscript were made using these modes:
##  Figure 3 - diversityOT
##  Figure 4 - FstOt
##  Figures 7 and S9 - stateTime
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

def mostRecent(some_list):
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
        assert (len(full)+len(month)+len(year)==len(tips)),'Warning, some dates were not captured by regex!\nReview the \"mostRecent\" function'
        for i in object_list:
            if isinstance(i,leaf)==True:
                if i.name in column(full,0):
                    return [i.height,column(full,1)[index(column(full,0),i.name)[0]]]
        
    except ValueError:
        pass
                

def make_tree(data):
    """
    Scans tree string looking for BEAST formatting.
    NOTE - only works with trees drawn from the posterior distribution of trees.
    """
    global i
    global trait_list
    global trait_names
    global partition_name
    stored_i=None
    while i < len(data):
        assert (stored_i != i),'\nTree string unparseable\nStopped at %s\nstring region looks like this: %s'%(data[i],data[i-10:i+10])
        stored_i=i
        cerberus=re.match('\[&([A-Za-z0-9\s\.\_\,\=\"\']+)\](:|;)',data[i:i+400])
        if cerberus is not None:
            traitComment=cerberus.group(1).split(',')
            for j in traitComment:
                traitName=j.split('=')[0]
                traitValue=j.split('=')[1].strip('"')
                ll.cur_node.traits[traitName]=traitValue
                if traitName not in column(trait_list,0):
                    trait_list.append([traitName,[traitValue]])
                if traitValue not in trait_list[index(column(trait_list,0),traitName)[0]][1]:
                    trait_list[index(column(trait_list,0),traitName)[0]][1].append(traitValue)
                    trait_list[index(column(trait_list,0),traitName)[0]][1]=sorted(trait_list[index(column(trait_list,0),traitName)[0]][1])
            i+=len(cerberus.group())-1

        if data[i] == '(':
            ll.add_node(i)
            i+=1
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
                ll.add_leaf(i,cerberus.group(1))
                i+=len(cerberus.group(1))
            elif data[i]=='(':
                make_tree(data)
                
        if data[i] == ',':
            i+=1
            ll.move_up()
            ll.cur_node.secondChild=True
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
                ll.add_leaf(i,cerberus.group(1))
                i+=len(cerberus.group(1))

        if data[i] == ':':
            i+=1
            comment=re.match('\[&([A-Za-z0-9\s\_\-\.\,\=]+)\]([0-9\.\-E]+)',data[i:i+400])
            if comment is not None:
                ll.cur_node.branch=float(comment.group(2))
                ll.cur_node.height=float(comment.group(2))
                com_list=[]
                com_list=comment.group(1).split(',')

                for j in com_list:
                    if 'r' in j.split('=')[0] and j.split('=')[0] not in ll.cur_node.traits.keys():
                        if partition_name=='':
                            partition_name=com_list[0].split('=')[0]
                        ll.cur_node.rate=float(com_list[0].split('=')[1])
                    elif 'r' in j.split('=')[0] and j.split('=')[0] in ll.cur_node.traits.keys():
                        ll.cur_node.traitRates[j.split('=')[0]]=float(j.split('=')[1])

                robustCounts=re.search('b_u_S=([0-9\.]+),S=([0-9\.]+),N=([0-9\.]+),b_u_N=([0-9\.]+)',comment.group(1))
                if robustCounts is not None:
                    ll.cur_node.buS=float(robustCounts.group(1))
                    ll.cur_node.S=float(robustCounts.group(2))
                    ll.cur_node.N=float(robustCounts.group(3))
                    ll.cur_node.buN=float(robustCounts.group(4))


                robustCounts=re.search('S=([0-9\.]+),N=([0-9\.]+)',comment.group(1))
                if robustCounts is not None:
                    ll.cur_node.S=float(robustCounts.group(1))
                    ll.cur_node.N=float(robustCounts.group(2))

                i+=len(comment.group())

        if data[i] == ')':
            ll.move_up()
            i+=1

        if data[i] == ';':
            assert (len(tips)*2-1 == len(object_list)),'\nTree string has been parsed incorrectly:\nexpected number of objects in tree %s\nobjects found in the tree string: %s'%(len(tips)*2-1,len(object_list))
            if len(trait_names)==0:
                trait_names=ll.cur_node.traits.keys()
            break

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
        global object_list
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

    def add_node(self,i):
        """
        Attach a node to current object as a first or second child.
        i refers to the index of the bracket defining the node in the tree string and is a unique identifier.
        """
        if self.cur_node.parent==self.root:
            self.cur_node.branch=0
            self.cur_node.height=0
        new_node=node()
        new_node.index=i
        object_list.append(new_node)
        if self.cur_node.secondChild==False:
            new_node.parent=self.cur_node
            self.cur_node.leftChild=new_node
            self.cur_node=new_node
        else:
            new_node.parent=self.cur_node
            self.cur_node.rightChild=new_node
            self.cur_node=new_node

    def add_leaf(self,i,name):
        """
        Attach a leaf to current object as a first or second child.
        i refers to the index of the tip number in the tree string and is a unique identifier.
        """
        new_leaf=leaf()
        new_leaf.index=i
        object_list.append(new_leaf)
        if self.cur_node.secondChild==False:
            new_leaf.name=name
            new_leaf.numName=name
            new_leaf.parent=self.cur_node
            self.cur_node.leftChild=new_leaf
            self.cur_node=new_leaf
        else:
            new_leaf.name=name
            new_leaf.numName=name
            new_leaf.parent=self.cur_node
            self.cur_node.rightChild=new_leaf
            self.cur_node=new_leaf

    def setAbsoluteTime(self,height,date):
        """
        Sets the absolute time of each object in the tree given the height of the most recent tip and its date.
        """
        for i in object_list:
            i.absoluteTime=date-height+i.height

    def renameTips(self):
        """
        Changes the name attribute of leaf objects to what the sequence name is.
        """
        global tips
        for i in object_list:
            if isinstance(i,leaf)==True:
                i.name=tips[index(column(tips,0),int(i.name))[0]][1]

    def TMRCA(self,objectA,objectB):
        """
        Identifies the most recent common ancestor of two objects.
        Accepts either leaf or node objects.
        """
        cur_nodeA=objectA
        cur_nodeB=objectB
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

    def traverse_tree(self):
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
                seen.append(cur_node.index)
                cur_node.height=height
                highestTip=height
                height-=cur_node.branch
                cur_node=cur_node.parent

###################################################

global object_list
object_list=[]
global maxHeight
maxHeight=0
global i
global ll
ll=None

def parseTreeFile(argv):
    """
    Parse mode argument.
    """
    try:
        opts, args = getopt.getopt(argv,"hi:m:",["ifile=","mode="])
    except getopt.GetoptError:
        print 'HrothvitnirParallel.py -i <inputfile> -m <mode>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'HrothvitnirParallel.py -i <inputfile> -m <mode>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            treefile = arg
            sys.stderr.write('tree file is: %s'%(treefile))
        elif opt in ("-m", "--mode"):
            assert arg in allModes,'\n%s mode does not exist.'%(arg)
            mode = arg

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
    
    global object_list
    object_list=[]
    global maxHeight
    maxHeight=0
    global i
    
    global ll
    ll=None
    
    makeSliceList(slice_list,StartTime,EndTime+1,Nbins)
    
    global tips
    tips=[]
    tipNum=0

    taxonlist=False
    plate=True
    for line in open(treefile,'r'):
        if plate==True:
            ## Extract useful information from the bits preceding the actual trees.
            cerberus=re.search('Dimensions ntax\=([0-9]+)\;',line)
            if cerberus is not None:
                tipNum=int(cerberus.group(1))

            if 'Translate' in line:
                taxonlist=True
        
            if taxonlist==True and ';' not in line and 'Translate' not in line:
                try:
                    ## Identifies sequence names given at the beginning of the *.trees.txt file
                    tip_index=int(line.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips:
                        tips.append([tip_index,tip])
                except ValueError:
                    pass

        if 'tree STATE_0' in line:
            plate=False
            assert (tipNum == len(tips)),'Expected number of tips: %s\nNumber of tips found: %s'%(tipNum,len(tips))
            global mostrecent

        cerberus=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line)
        if cerberus is not None:
            ## At tree state 0 insert header into output file
            if int(cerberus.group(1))==0:
                i=0
                object_list=[]
                maxHeight=0
                ll=None
                ll=tree()
                start=len(cerberus.group())
                treestring=str(line[start:])
                make_tree(treestring)
                #######################################################
                if mode=="":
                    pass
                #######################################################
                if mode=="X":
                    sys.stdout.write('states\tN objects')
                #######################################################
                if mode=='stateTime':
                    ## Define traits of interest. Here it's PB1, PB2 and HA
                    Mtrait_list=[]
                    for o in trait_list:
                        if o[0] in ['PB1sub','PB2sub','HAsub'] or o[0] in ['PB1','PB2','HA']: 
                            Mtrait_list.append(o)
                    sys.stdout.write('states\tVicCoreN\tYamCoreN\tNCN\tVicCoreS\tYamCoreS\tNCS\tVicCore_nt\tYamCore_nt\tNC_nt\tVicCoreTime\tYamCoreTime\tNCTime')
                #######################################################
                if mode=='FstOT':
                    ## Define traits of interest. Here it's PB1, PB2 and HA
                    Mtrait_list=[]
                    for o in trait_list:
                        if o[0] in ['PB1sub','PB2sub','HAsub'] or o[0] in ['PB1','PB2','HA']: 
                            Mtrait_list.append(o)

                    sys.stdout.write('states')
                    for time in slice_list:
                        for trait in Mtrait_list:
                            sys.stdout.write('\t%s_%s.tmrcaB'%(time[0],trait[0]))
                #######################################################
                if mode=='diversityOT':
                    sys.stdout.write('states')
                    for time in slice_list:
                        sys.stdout.write('\t%s_diversity'%(time[0]))
                #######################################################

            ## After burnin start processing
            if int(cerberus.group(1)) >= burnin:
                ## i is a global variable used to parse the tree string
                ## object_list contains all the objects in the tree
                ## ll is the tree object
                
                i=0
                object_list=[]
                maxHeight=0
                ll=None
                ll=tree()
                start=len(cerberus.group())
                treestring=str(line[start:])
                make_tree(treestring)

                ## Traverse the tree - sets the height of each object in the tree
                ll.traverse_tree()

                ## Rename tips so their name refers to sequence name
                ll.renameTips()

                ## Calibrate tree - find the height of the most recent tip and its date
                calTipHeight,calTipDate=mostRecent(tips)
                ll.setAbsoluteTime(calTipHeight,calTipDate)
                ################################################################################
                if mode == "":
                    pass
                ################################################################################
                if mode == "X":
                    sys.stdout.write('\n%s\t%s'%(cerberus.group(1),len(object_list)))
                ################################################################################
                if mode=='stateTime':
                    CYN=0
                    CVN=0
                    NCN=0

                    CYS=0
                    CVS=0
                    NCS=0

                    CYnt=0
                    CVnt=0
                    NCnt=0

                    CYt=0
                    CVt=0
                    NCt=0

                    ## Go through objects in the tree, collect non-synonymous, synonymous, nucleotide substitutions and total time spent under each PB1-PB2-HA constellation
                    for k in object_list:
                        if k.N is not None and k.S is not None:
                            linstate=''
                            for x in range(len(Mtrait_list)):
                                linstate+=k.traits[Mtrait_list[x][0]]
                            if linstate.count('Y')==len(Mtrait_list):
                                CYN+=k.N
                                CYS+=k.S
                                CYnt+=k.rate*k.branch
                                CYt+=k.branch
                                
                            elif linstate.count('V')==len(Mtrait_list):
                                CVN+=k.N
                                CVS+=k.S
                                CVnt+=k.rate*k.branch
                                CVt+=k.branch
                            else:
                                NCN+=k.N
                                NCS+=k.S
                                NCnt+=k.rate*k.branch
                                NCt+=k.branch

                    ## Divide substitutions by total time spent - gives the rate
                    CVNrate=float(CVN)/float(CVt)
                    CYNrate=float(CYN)/float(CYt)
                    NCNrate=float(NCN)/float(NCt)

                    CVSrate=float(CVS)/float(CVt)
                    CYSrate=float(CYS)/float(CYt)
                    NCSrate=float(NCS)/float(NCt)

                    CVntRate=float(CVnt)/float(CVt)
                    CYntRate=float(CYnt)/float(CYt)
                    NCntRate=float(NCnt)/float(NCt)

                    sys.stdout.write('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(cerberus.group(1),CVNrate,CYNrate,NCNrate,CVSrate,CYSrate,NCSrate,CVntRate,CYntRate,NCntRate,CVt,CYt,NCt))
                ################################################################################
                if mode=='FstOT':
                    ## Output matrix
                    diversityBetweenTMRCA=list(list([] for u in range(len(Mtrait_list))) for r in range(len(slice_list)))

                    ## Split objects in the tree based on when they exist
                    fractionedYears=list([] for u in range(len(slice_list)))
                    for k in object_list:
                        for x in range(len(slice_list)):
                            if k.absoluteTime>float(slice_list[x][0])>=k.parent.absoluteTime:
                                fractionedYears[x].append(k)


                    
                    done={}
                    for x in range(len(slice_list)):
                        temp=[]
                        ob=fractionedYears[x]
                        for y in range(len(Mtrait_list)):
                            ## Split objects in each time slice into Vic and Yam lineages
                            splitObjectsA,splitObjectsB=[[Obj for Obj in ob if Obj.traits[Mtrait_list[y][0]]==val] for val in Mtrait_list[y][1]]
                            if len(splitObjectsA)!=0 and len(splitObjectsB)!=0:
                                for q in splitObjectsA:
                                    for w in splitObjectsB:
                                        if done.has_key(q.index):
                                            ## For all pairwise comparisons find TMRCA of each pairing of a Vic and a Yam object
                                            if done[q.index].has_key(w.index):
                                                diversityBetweenTMRCA[x][y].append(done[q.index][w.index])
                                            else:
                                                done[q.index][w.index]=ll.TMRCA(q,w).absoluteTime
                                                diversityBetweenTMRCA[x][y].append(done[q.index][w.index])
                                        else:
                                            done[q.index]={}
                                            done[q.index][w.index]=ll.TMRCA(q,w).absoluteTime
                                            diversityBetweenTMRCA[x][y].append(done[q.index][w.index])           

                    ## Output to file
                    sys.stdout.write('\n%s'%(cerberus.group(1)))    
                    for x in range(len(slice_list)):
                        for y in range(len(Mtrait_list)):
                            sys.stdout.write('\t%s'%(np.mean(diversityBetweenTMRCA[x][y])))            
                ################################################################################
                if mode=='diversityOT':
                    ## Split objects in the tree based on when they exist
                    fractionedYears=list([] for u in range(len(slice_list)))
                    for k in object_list:
                        for x in range(len(slice_list)):
                            if k.absoluteTime>float(slice_list[x][0])>=k.parent.absoluteTime:
                                fractionedYears[x].append(k)

                    sys.stdout.write('\n%s'%(cerberus.group(1)))

                    ## For each pair of objects in each time slice find TMRCA
                    for x in range(len(slice_list)):
                        temp=[]
                        ob=fractionedYears[x]
                        for g in range(len(ob)-1):
                            for h in range(g+1,len(ob)):
                                temp.append(ll.TMRCA(ob[g],ob[h]).absoluteTime)

                        if len(temp)==0:
                            sys.stdout.write('\t0')
                        else:
                            ## Only output the oldest TMRCA from the list of TMRCAs
                            sys.stdout.write('\t%s'%(min(temp)))
                ################################################################################
            if 'End;' in line:
                print '\nEnd of file!'

######################################################

global burnin

######################################################

## Modes available for use as arguments:
# None - if no argument is given the script will simply parse the trees without any output.
# X - outputs the number of objects (nodes+leaves) in the tree at each MCMC state.
# stateTime - outputs the non-synonymous (*N), synonymous (*S) and nucleotide (*_nt) substitution rates as well as total amount of time spent under different PB1-PB2-HA constellations.
# FstOT - outputs mean pairwise diversity between V and Y labeled branches at each time point.
# diversityOT - outputs time of most recent common ancestor of all lineages existing at each time point.

allModes=['X','stateTime','FstOT','diversityOT']

## Define burnin as the number of the state from which to begin analysis.
burnin=20000000

######################################################

## Define start and end of time slices as well as how many samples each year should have.
StartTime=1984
EndTime=2013
Nbins=1

######################################################

import getopt

if __name__ == "__main__":
   parseTreeFile(sys.argv[1:])

######################################################
sys.stderr.write('\n>>>>> Done! <<<<<\n')
