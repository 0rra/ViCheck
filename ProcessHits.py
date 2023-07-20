#import logger
import argparse
import sqlite3
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from ete3 import NCBITaxa
import seaborn as sb
import time

def get_tax(data, subquery, desrank):
    ncbi = NCBITaxa()
    if subquery == 'subject':
        taxids = data['Subject_taxid'].tolist()
    if subquery == 'query':
        taxids = data['Query_taxid'].tolist()

    def get_rank_code(taxid, desired_ranks):
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        lineage2ranks = ncbi.get_rank(names)
        #print(lineage2ranks)
        # switches round taxid and rank order and repackages in dictionary
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        # formats with rank:
        return [ranks2lineage.get(rank, '<not present>') for rank in desired_ranks]

    def get_rank_info(rankname):
        taxinfo = {}
        # desired_ranks = ['kingdom','family','species']
        desired_ranks = [rankname]
        for taxid in taxids:
            out = []
            ranks = get_rank_code(taxid, desired_ranks)
            for i in ranks:
                if i != '<not present>':
                    ranklist = list(ncbi.get_taxid_translator([i]).values())[0]
                    taxinfo.update({taxid: ranklist})
                    #print(ranklist)
                else:
                    ranklist = 'Not present'
                    taxinfo.update({taxid: ranklist})
        return taxinfo


    rankinfo = get_rank_info(desrank)
    if subquery == 'subject':
        data[f'Subject_{desrank}'] = data.Subject_taxid.map(rankinfo)
    if subquery == 'query':
        data[f'Query_{desrank}'] = data.Query_taxid.map(rankinfo)
    return data

# add seq contamination info
def contamcat(data, subquery):
    conditions = [
        (data[f'{subquery}_start'] <100), (data[f'{subquery}_length'] - data[f'{subquery}_end'] < 100), (data['Alignment_length'] > (data['Subject_length'] * 0.5))
    ]
    choices= ['StartContamination', 'EndContamination', 'BigOverlap']
    data['Contamination'] = np.select(conditions, choices, default='Unclassified')
    return data

# given results file name: read in chunks
def results_annotation(resultsfile, filter):
    # read in file, add column names
    resultsread = pd.read_csv(resultsfile, sep="\t", chunksize=1000, header=None, names=range(17))
    resultsdf = pd.concat(resultsread)

    resultsdf.columns = ['Query', 'Subject', 'Percentage', 'Alignment_length', 'Mismatch', 'Gap_Open', 'Query_start',
                         'Query_end', 'Subject_start', 'Subject_end', 'Evalue', 'Bitscore', 'Query_length',
                         'Subject_length', 'QCOVS', 'QCOVUS', 'QCOVHSP']



    # split query and subject ids
    resultsdf[['Query_accession', 'Query_taxid', 'Query_division', 'Query_name']] = resultsdf['Query'].str.split('|',
                                                                                                                 expand=True)

    if filter == 'virus':
        # remove phages
        resultsdf = resultsdf[resultsdf['Query_division'] == 'VIRUS']

    if filter == 'all':
        pass

    resultsdf[['Subject_accession', 'Subject_taxid', 'Subject_division', 'Subject_name']] = resultsdf[
        'Subject'].str.split('|', expand=True)

    # add taxonomy information
    resultsdf = get_tax(resultsdf, 'subject', 'family')
    resultsdf = get_tax(resultsdf, 'subject', 'species')
    resultsdf = get_tax(resultsdf, 'query', 'family')
    resultsdf = get_tax(resultsdf, 'query', 'species')

    # basic filter
    resultsdf = resultsdf.query(f'Percentage > 85')
    resultsdf = resultsdf.query(f'Alignment_length > 100')

    # flag retrovirus
    retrolist = []
    with open('retroviruses.txt') as retrofile:
        for row in retrofile:
            row = row.strip('\n')
            retrolist.append(row)
    resultsdf['Retro'] = resultsdf['Query_taxid'].isin(retrolist)
    if filter == 'virus':
        resultsdf = resultsdf[~resultsdf['Query_taxid'].isin(retrolist)]

    master = pd.read_csv('ICTV_Master.csv', sep=',')
    dsDNAent = master[master['Genome Composition'].str.contains("dsDNA")]
    dsDNAfam = dsDNAent['Family'].unique()
    resultsdf['dsDNA'] = resultsdf['Query_family'].isin(dsDNAfam)
    # add contamination info
    resultsdf = contamcat(resultsdf, 'Query')

    return resultsdf

# refseq splitting section
VIRUS = results_annotation('refseq_hits.txt', 'virus')
#NORETRO = VIRUS[VIRUS['Retro'] == False]
allref = results_annotation('refseq_hits.txt', 'all')
PHAGE = allref[allref['Query_division'] == 'PHAGE']
RETRO = allref[allref['Retro'] == True]

# for other 2 datasets
comp = results_annotation('complete_hits.txt', 'virus')
part = results_annotation('partial_hits.txt', 'virus')

def add_db(hitsdf,tablename):

    def make_table(table):
        connection = sqlite3.connect('hits2.db')
        cur = connection.cursor()
        createerror = False
        try:
            cur.execute(f'CREATE TABLE {tablename} (Query, Subject, Percentage, Alignment_length, Mismatch, Gap_Open, Query_start, Query_end, Subject_start, Subject_end, Evalue, Bitscore, Query_length, Subject_length, QCOVS, QCOVUS, QCOVHSP, Query_accession, Query_taxid, Query_division, Query_name, Subject_accession, Subject_taxid, Subject_division, Subject_name, Subject_family, Subject_species, Query_family, Query_species, Retro, dsDNA, Contamination);')
        except sqlite3.OperationalError:
            connection.rollback()
            createerror = True
        if not createerror:
            connection.commit()
        cur.close()
        connection.close()

    def adddata(data, table):
        connection = sqlite3.connect('hits2.db')
        n = 1000
        data = [data[i:i + n] for i in range(0, len(data), n)]
        for part in data:
            # add chunk to sql
            part.to_sql(table, connection, if_exists='append', index=False)
            connection.commit()
        connection.close()

    make_table(tablename)
    adddata(hitsdf, tablename)

    start = hitsdf[(hitsdf['Contamination'] == 'StartContamination')]
    end = hitsdf[(hitsdf['Contamination'] == 'EndContamination')]
    make_table(f'Start{tablename}')
    make_table(f'End{tablename}')
    adddata(start, f'Start{tablename}')
    adddata(end, f'End{tablename}')
    # save all of these to csv too
    hitsdf.to_csv(f'{tablename}.csv', sep='\t', header=True)
    start.to_csv(f'Start{tablename}.csv', sep='\t', header=True)
    end.to_csv(f'End{tablename}.csv', sep='\t', header=True)

add_db(VIRUS, 'NonRetroResults')
add_db(allref, 'RefseqResults')
add_db(PHAGE, 'PhageResults2')
add_db(RETRO, 'RetroResults2')
add_db(part, 'PartialResults')
add_db(comp, 'CompleteResults')

