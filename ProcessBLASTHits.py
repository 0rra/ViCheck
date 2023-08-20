import argparse
import pandas as pd
import os
from GetTax import get_tax
import sqlite3
import numpy as np

# script to process BLAST hit result result files, adding information, saving to csv and adding to database
# command line interface for process blast hits
parser = argparse.ArgumentParser(prog = 'ProcessBLASTHits', description='Processes BLAST result files, adds taxonomy and position of potential contamination')
parser.add_argument('BlastFile', help='Give blastn results file')
parser.add_argument('ResultsName', help = "Give name of results csv file and name of results table (do not include file extension)")
parser.add_argument('--RetrovirusFile', default = "retroviruses.txt", help='List of retroviruses needed to filter retro and non retro')
# optional arguments
parser.add_argument('--OriginDatabase', help="Give name of BLAST database used to add column containing database info e.g. RefSeq. Ideal if merging datasets")
parser.add_argument("--AddDb", nargs='?', const="Hits.db", help = "Give database name to add results to")
args = parser.parse_args()


class AnnotateBLAST:
    def __init__(self, resfile):
        # reads in file as dataframe in chunks
        resultsread = pd.read_csv(resfile, sep="\t", chunksize=1000, header=None, names=range(17))
        resultsdf = pd.concat(resultsread)
        self.df = resultsdf

    # calculates position of alignment from query start and end
    def distancefrom(self):
        data = self.df
        data['DistanceStart'] = data['Query_start'] - 1
        data['DistanceEnd'] = data['Query_length'] - data['Query_end']
        return data

    # categorises if potential contamination is at sequence start or end
    def contamcat(self, subquery):
        data = self.df
        conditions = [
            (data[f'{subquery}_start'] <100), (data[f'{subquery}_length'] - data[f'{subquery}_end'] < 100)
        ]
        choices= ['Start', 'End']
        data[f'{subquery}_contamination'] = np.select(conditions, choices, default='Unclassified')
        return data

    # dataframe processing and annotation
    def results_annotation(self):
        # adds column names
        self.df.columns = ['Query', 'Subject', 'Percentage', 'Alignment_length', 'Mismatch', 'Gap_Open', 'Query_start',
                             'Query_end', 'Subject_start', 'Subject_end', 'Evalue', 'Bitscore', 'Query_length',
                             'Subject_length', 'QCOVS', 'QCOVUS', 'QCOVHSP']

        # split query and subject information columns into separate columns
        self.df[['Query_accession', 'Query_taxid', 'Query_division', 'Query_name']] = self.df[
            'Query'].str.split('|',expand=True)

        self.df[['Subject_accession', 'Subject_taxid', 'Subject_division', 'Subject_name']] = self.df[
            'Subject'].str.split('|', expand=True)

        # add family and species taxonomy information for subject and queries
        self.df = get_tax(self.df, 'subject', 'family')
        self.df = get_tax(self.df, 'subject', 'species')
        self.df = get_tax(self.df, 'query', 'family')
        self.df = get_tax(self.df, 'query', 'species')

        #flag retroviruses
        retrolist = []
        with open(args.RetrovirusFile) as retrofile:
            print(args.RetrovirusFile)
            for row in retrofile:
                row = row.strip('\n')
                retrolist.append(row)
        self.df['Retro'] = self.df['Query_taxid'].isin(retrolist)

        # add contamination info
        self.df = self.contamcat('Query')
        self.df = self.contamcat('Subject')
        self.df = self.distancefrom()

        return self.df

# checking if blast results file exists
if not os.path.isfile(args.BlastFile):
    raise SystemExit(1)
    print('Please give a blast results file')
else:
    # processes the results file with the annotate BLAST class
    sample = AnnotateBLAST(args.BlastFile).results_annotation()


# check if origin database argument given
if args.OriginDatabase == None:
    sample["Origin"] = 'NA'
# if given, annotate origin database
if not args.OriginDatabase == None:
    sample["Origin"] = args.OriginDatabase

# add results to sql database
class AddToDb:
    def __init__(self, dbpath):
        self.db = dbpath

    # makes results table in database
    def make_table(self, table):
        connection = sqlite3.connect(self.db)
        cur = connection.cursor()
        try:
            cur.execute(f'CREATE TABLE {table} (Query, Subject, Percentage, Alignment_length, Mismatch, Gap_Open, Query_start, Query_end, Subject_start, Subject_end, Evalue, Bitscore, Query_length, Subject_length, QCOVS, QCOVUS, QCOVHSP, Query_accession, Query_taxid, Query_division, Query_name, Subject_accession, Subject_taxid, Subject_division, Subject_name, Subject_family, Subject_species, Query_family, Query_species, Retro, dsDNA, Query_contamination, Subject_contamination, DistanceStart, DistanceEnd, Origin, SubjectCheckCode, SubjectCheckResult);')
            connection.commit()
        except sqlite3.OperationalError:
            connection.rollback()
            print('Table already exists')
        cur.close()
        connection.close()

    def adddata(self, table, data):
        connection = sqlite3.connect(self.db)
        # splits up the database 1000 rows at a time
        data = [data[i:i + 1000] for i in range(0, len(data), 1000)]
        # adds data to database in chunks to maintain performance when adding large results files
        for chunk in data:
            chunk.to_sql(table, connection, if_exists='append', index=False)
            connection.commit()
        connection.close()
        print('added to database!')

# if option to add to db used
if not args.AddDb == None:
    # makes table for results and adds results to table
    usedb = AddToDb(args.AddDb)
    usedb.make_table(args.ResultsName)
    usedb.adddata(args.ResultsName, sample)

# writes processed blast results to csv
sample.to_csv(f'{args.ResultsName}.csv', sep='\t', header=True)
print('writing to csv')

