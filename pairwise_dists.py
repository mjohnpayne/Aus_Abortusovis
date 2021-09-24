from Bio import SeqIO
import pandas as pd
from time import sleep as sl
import time
import itertools
from multiprocessing import Pool
from matplotlib import pyplot as plt
import seaborn as sns

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)

def snp_dist_metric(a,b):
    match = 0
    missmatch = 0
    for pos in range(len(a)):
        anuc = a[pos]
        bnuc = b[pos]
        missing = ["N","n","X","x","-"]
        if anuc in missing or bnuc in missing or anuc == bnuc:
            match +=1
        else:
            missmatch +=1
            # if missmatch >= args.max_missmatch:
            #     return missmatch
    return missmatch

def pairdist(pair,profs):
    id1 = pair[0]
    id2 = pair[1]

    ap1 = profs[id1]
    ap2 = profs[id2]

    dist = snp_dist_metric(ap1, ap2)
    return dist

def run_dist(profs,idlist,outcsv):


    start_time = time.time()
    newdf = pd.DataFrame(index=idlist,columns=idlist,dtype=str)
    newdf.index = newdf.index.astype(str)
    newdf.columns = newdf.columns.astype(str)
    # newdf.fillna(args.max_missmatch,inplace=True)

    exist = []
    for i in newdf.columns.tolist():
        if i not in exist:
            exist.append(i)
        else:
            print("{} duplicated in new df".format(i))

    pairs = itertools.combinations(idlist,r=2)
    pairs = list(pairs)
    tot = (len(idlist)*len(idlist)-len(idlist))/2
    frac = int(tot*0.01)
    if frac == 0:
        frac = 1

    c=0



    print("\nCalculating pairwise distances\nDone\t\t% done\t\ttime")
    pre = 0
    new = 0
    distlist = [pairdist(x,profs) for x in pairs]


    # distlist = run_multiprocessing(pairdist,inps,1)

    # exist = []
    # for i in newdf.columns.tolist():
    #     if i not in exist:
    #         exist.append(i)
    #     else:
    #         print("{} duplicated".format(i))


    for pos in range(len(pairs)):
        pair = pairs[pos]
        dist = distlist[pos]
        id1 = pair[0]
        id2 = pair[1]
        newdf.at[id1,id2] = dist
        newdf.at[id2,id1] = dist
        c+=1

        if c%frac == 0:
            print("{}\t\t{}%\t\t{:.3f} seconds".format(c,int((c/tot)*100),time.time() - start_time), end="\r", flush=True)

    print("{}\t\t{}%\t\t{:.3f} seconds".format(c, 100, time.time() - start_time),
          end="\r", flush=True)
    for i in idlist:
        newdf.at[i,i] = 0
    newdf.to_csv(outcsv,sep="\t")

    return newdf

snpalign = "placeholder/path/snippycore_w_westmead_aln.fasta"
inmeta = "placeholder/path/strain_metadata_for_tree_updated.xlsx"

snpaligns = SeqIO.parse(snpalign,"fasta")
snpaligns = {x.id:str(x.seq) for x in snpaligns}
ids = list(snpaligns.keys())
ids = [x for x in ids if x != "Reference"]

## First time running use below two lines
outfile = "pairwise_dist_matrix.csv"
dist_df = run_dist(snpaligns,ids,outfile)

# Subsequent runs use following line to read in precalculated pairs
# dist_df = pd.read_csv("pairwise_dist_matrix.csv",sep="\t",index_col=0)


pairs = itertools.combinations(ids,r=2)

inmeta = pd.read_excel(inmeta,sheet_name="strain_metadata_for_tree_update")
inmeta = inmeta.set_index("tree_id")

# Aus vs non Aus
inmeta["Country_Aus"] = inmeta["Country"]=="Australia"
print(inmeta["Country_Aus"])

distls = []
annotlist = []
for pair in pairs:
    dist = dist_df.loc[pair[0],pair[1]]
    type = ""
    type1 = inmeta.loc[pair[0],"Country_Aus"]
    type2 = inmeta.loc[pair[1],"Country_Aus"]
    if type1 == True and type2 == True:
        type = "Aus"
    elif type1 == False and type2 == False:
        type = "Other"
    else:
        type = "Mixed"
    distls.append(dist)
    annotlist.append(type)
histdata = [distls,annotlist]
histdata = pd.DataFrame(histdata).T
histdata.columns = ["dist","annot"]
ag = histdata.groupby(['annot']).agg({'dist':[min,max,list]})

l = histdata["dist"].to_list()


# ag contains summary of pairwise distances for three categories
print(ag)

# below plots not used
# p = sns.histplot(histdata,x="dist",hue="annot")
# p.set(xlim=(0, 50))
# plt.savefig("within_aus_pairwise_distances.pdf")
#
# p = sns.histplot(histdata,x="dist",hue="annot")
# p.set(xlim=(700, 1300))
# p.set(ylim=(0, 5))
# plt.savefig("within_non-aus_pairwise_distances.pdf")
#
# p = sns.histplot(histdata,x="dist",hue="annot")
# p.set(xlim=(22700, 22900))
# p.set(ylim=(0, 50))
# plt.savefig("between_aus_non-aus_pairwise_distances.pdf")

### clinical vs poultry

inmeta = "placeholder/path/strain_metadata_for_tree_updated.xlsx"
inmeta = inmeta[inmeta["Country"]=="Australia"]
inmeta["poultry"] = inmeta["Associated_species"]=="Chicken"

ids = list(inmeta.index)
distls = []
annotlist = []
for pair in pairs:
    if pair[1] in ids and pair[0] in ids:
        dist = dist_df.loc[pair[0],pair[1]]
        type = ""
        type1 = inmeta.loc[pair[0],"poultry"]
        type2 = inmeta.loc[pair[1],"poultry"]
        if type1 == True and type2 == True:
            type = "Poultry"
        elif type1 == False and type2 == False:
            type = "Clinical"
        else:
            type = "Mixed"
        distls.append(dist)
        annotlist.append(type)
histdata = [distls,annotlist]
histdata = pd.DataFrame(histdata).T
histdata.columns = ["dist","annot"]
ag = histdata.groupby(['annot']).agg({'dist':[min,max,list]})

for i in ag[("dist","list")]:
    print(list(i))
print(ag)


p = sns.histplot(histdata,x="dist",hue="annot",hue_order=["Clinical","Poultry","Mixed"])

plt.savefig("Aus_organism_pairwise_dists.pdf")

