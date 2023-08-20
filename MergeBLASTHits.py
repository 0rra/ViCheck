import argparse
import pandas as pd
import os
import sqlite3
#from ProcessBLASTHits_c import AddToDb


parser = argparse.ArgumentParser(prog = 'MergeBLASTHits', description='Merge blast hit results files')
parser.add_argument('Processedfile', help='Give processed BLAST hits file')
parser.add_argument('ResultsName', help = "Give name of merged results file / name of merged results table")
parser.add_argument('FilesToMerge', nargs='*', help='Give files to merge with processed file, must give at least one file')
parser.add_argument("--AddDb", nargs='?', const="Hits.db", help = "Give database name to add results to")
args = parser.parse_args()


if not os.path.isfile(args.Processedfile):
    print('Please give a blast results file')
    raise SystemExit(1)
else:
    print('Reading files')
    sample = pd.read_csv(args.Processedfile, sep="\t")

    dataframelist = [sample]
# check if datasets to merge have similar structure

for i in args.FilesToMerge:
    df = pd.read_csv(i, sep="\t")
    dataframelist.append(df)

#merge dfs and reset indexes
mergedfile = pd.concat(dataframelist, ignore_index=True)
print(f'Merging {args.FilesToMerge} and {args.Processedfile}')
# drop unnamed column
mergedfile = mergedfile.drop(labels=['Unnamed: 0'], axis=1)
#mergedfile = mergedfile.drop("Unnamed: 0", axis=1)
mergedfile.to_csv(f'{args.ResultsName}.csv', sep="\t")

class AddToDb:
    def __init__(self, dbpath):
        self.db = dbpath

    def make_table(self, table):
        connection = sqlite3.connect(self.db)
        cur = connection.cursor()
        try:
            cur.execute(
                f'CREATE TABLE {table} (Query, Subject, Percentage, Alignment_length, Mismatch, Gap_Open, Query_start, Query_end, Subject_start, Subject_end, Evalue, Bitscore, Query_length, Subject_length, QCOVS, QCOVUS, QCOVHSP, Query_accession, Query_taxid, Query_division, Query_name, Subject_accession, Subject_taxid, Subject_division, Subject_name, Subject_family, Subject_species, Query_family, Query_species, Retro, dsDNA, Query_contamination, Subject_contamination, DistanceStart, DistanceEnd, Origin, SubjectCheckCode, SubjectCheckResult);')
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

# if argument to add to database given
if not args.AddDb == None:
    # makes table for results and adds results to table
    usedb = AddToDb(args.AddDb)
    usedb.make_table(args.ResultsName)
    usedb.adddata(args.ResultsName, mergedfile)
if args.AddDb == None:
    print('Not adding to db')
