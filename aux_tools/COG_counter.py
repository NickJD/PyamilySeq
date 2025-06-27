import argparse
import collections

from scipy.stats import chisquare

def anno_load(anno_in,COGs):
    count = 0
    for line in anno_in:
        if line.startswith('#') or line.startswith('\n'):
            continue
        else:
            line = line.split('\t')
            data = line[6]
            for COG in range(0, len(data)):
                COGs[data[COG]]+=1
            if "-" not in data:
                count +=1


    return COGs,count



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', action='store', dest='cog_in', required=True,
                        help='EggNog Annotations')
    options = parser.parse_args()

    COGs = collections.defaultdict(int)

    anno_in = open(options.cog_in, 'r')
    anno_COGs, count = anno_load(anno_in, COGs)
    print("Number of Seqs with COGS: " + str(count))
    print(anno_COGs)

    COG_Groups = collections.OrderedDict({"INFO":0,"CELL":0,"META":0,"POOR":0})
    COG_Groups["INFO"] = (COGs["J"]+COGs["A"]+COGs["K"]+COGs["L"]+COGs["B"])
    COG_Groups["CELL"] = (COGs["D"] + COGs["Y"] + COGs["V"] + COGs["T"] + COGs["M"]+COGs["N"]+
                          COGs["Z"]+COGs["W"]+COGs["U"]+COGs["O"])
    COG_Groups["META"] = (COGs["C"]+COGs["G"]+COGs["E"]+COGs["F"]+COGs["H"]+COGs["I"]+COGs["P"]+
                          COGs["Q"])
    COG_Groups["POOR"] = (COGs["R"]+COGs["S"])

    all = (COG_Groups["INFO"]+COG_Groups["CELL"]+COG_Groups["META"]+COG_Groups["POOR"])

    print("INFORMATION STORAGE AND PROCESSING: " + str(COG_Groups["INFO"])  + "\t" + str((COG_Groups["INFO"]/all)*100) + "\nCELLULAR PROCESSES AND SIGNALING: " +
          str(COG_Groups["CELL"])+ "\t" + str((COG_Groups["CELL"]/all)*100) + "\nMETABOLISM: " + str(COG_Groups["META"]) + "\t" + str((COG_Groups["META"]/all)*100) +
          "\nPOORLY CHARACTERIZED: " + str(COG_Groups["POOR"]) + "\t" + str((COG_Groups["POOR"]/all)*100) )

    print("\n\n#############################")
    print('J: '+str(COGs["J"])+'\nA: '+str(COGs["A"])+'\nK: '+str(COGs["K"])+'\nL: '+str(COGs["L"])+'\nB: '+str(COGs["B"]))