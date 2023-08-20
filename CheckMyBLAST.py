import pandas as pd
import argparse
import csv

parser = argparse.ArgumentParser(prog = 'CheckMyBLAST', description='Checks a BLAST results file to see if it contains problematic sequences')
parser.add_argument('BLASTRes', help='Give BLAST hits file')
parser.add_argument('ResultsName', help = "Give name of results file")
parser.add_argument('ProblemSequences', nargs = '*', help='Give both the problematic virus and non-virus sequence files (in this order)')
parser.add_argument("--Virus", action = 'store_true', help = 'Flags problematic non-virus subjects')
parser.add_argument("--NonVirus", action = 'store_true', help = 'Flags problematic non-virus subjects')
parser.add_argument("--Strict", action = 'store_true', help= "Removes any hits including problematic subjects")
args = parser.parse_args()


problemVirus = args.ProblemSequences[0]
problemNonVirus = args.ProblemSequences[1]

# read in problematic sequence file and get list of ids
def getiddict(file, type):
    data = pd.read_csv(file, sep='\t')
    # get unique accessions
    uniqueproblem = data.drop_duplicates(subset=[f"{type}_accession"])
    # sets accession to index
    outdict = uniqueproblem.set_index(f"{type}_accession")
    # converts to dictionary where flag and action can be access as keys
    outdict = outdict[["Flag", "Action"]].to_dict('index')
    return outdict


# get problematic virus ids
if args.Virus == True:
    problemseqs = getiddict(problemVirus, "Query")

# get problematic nonvirus ids
if args.NonVirus == True:
    problemseqs = getiddict(problemNonVirus, "Subject")

if args.NonVirus == True and args.Virus == True:
    problemseqsQ = getiddict(args.ProblemSequences, "Query")
    problemseqS = getiddict(args.ProblemSequences, "Subject")
    problemseqs = problemseqS | problemseqsQ

matches = []
toremove = []
toflag = []
# check if ids are in Blast results file

# for the problematic sequences
for id in problemseqs.keys():
    # open the blast file and the summary file (to write down the warning info)
    with open(args.BLASTRes, "r") as hitfile, open("WarningSummary.txt", "a") as summary:
        # read through file and check if it contains a match with problematic id
            for line in hitfile:
                # checks if line contains match with problematic ids
                if id in line:
                    # if action is to remove, add line to remove list
                    if problemseqs[id]["Action"] == 'Remove':
                        print(f"Warning, remove due to {problemseqs[id]['Flag']}")
                        summary.write(f"Warning, remove due to {problemseqs[id]['Flag']} : {line} contains {id}")
                        toremove.append(line)
                    # if action is to flag, add line to flag list
                    if problemseqs[id]["Action"] == 'Flag':
                        print(f"Warning, {problemseqs[id]['Flag']}")
                        summary.write(f"Warning, {problemseqs[id]['Flag']: {line} contains {id}}")
                        toflag.append(line)

# makes a results file with all the lines containing problematic sequences removed
if args.Strict == True:
    with open(args.BLASTRes, "r") as hitfile, open(args.ResultsName, "w") as output:
        for line in hitfile:
            # does not write any lines which are in toremove
            if line not in toremove:
                output.write(line)






