import pandas as pd
import numpy as np

refseq_hits = pd.read_csv('RefseqHits.csv', sep ='\t')
# get unique
#get count of number unique subject accessions total
refuni = refseq_hits['Subject_accession'].unique()
refunilist = refuni.tolist()
refuni = pd.DataFrame(refuni)
refuni.to_csv('allrefseq_sub_id.csv', sep=',',index=False, header=False)


# get all unique ids from all hits > 95
all_hits = pd.read_csv('AllHits.csv', sep='\t')
all_hits95 = all_hits.query('Percentage >95')
all_hits95 = all_hits95['Subject_accession'].unique()

# remove refseq hits to avoid BLAST redundant searches
# converet all hits back to dataframe, now only has unique subject ids
allsubuni = pd.DataFrame(all_hits95)
allsubuni.rename(columns={0: 'Subject_accession'}, inplace=True)
# checks if any subejct ids match ids in refseq hits >95
allsubuni['MatchRefseq'] = allsubuni['Subject_accession'].isin(refunilist)
# gets all sub hits which are not in refseq hits
allsubuni_uni = allsubuni.query('MatchRefseq == False')
# get list of all hits
all_sub_acclist = allsubuni_uni['Subject_accession'].unique().tolist()

# split list into 3 parts
splitlist = np.array_split(all_sub_acclist, 3)

part1 = splitlist[0]
part2 = splitlist[1]
part3 = splitlist[2]

# save parts to dataframe and csv
part_one = pd.DataFrame(part1)
part_one.to_csv('95allsubs_part1.csv', sep=',', index=False, header=False)
part_two = pd.DataFrame(part2)
part_two.to_csv('95allsubs_part2.csv', sep=',', index=False, header=False)
part_three = pd.DataFrame(part3)
part_three.to_csv('95allsubs_part3.csv', sep=',', index=False, header=False)



