import os.path
from GetTax import get_tax
import sqlite3
import pandas as pd
pd.options.mode.chained_assignment = None


# function which processes the subject hit results, adding column info, taxonomy info
def process_subhits(data):
    # add in column names
    subhitdata = data.copy()
    subhitdata.columns = ['Query', 'Subject', 'Percentage', 'Alignment_length', 'Mismatch', 'Gap_open', 'Query_start',
                          'Query_end', 'Subject_start', 'Subject_end', 'Evalue', 'QCOVS', 'Subject_group',
                          'Subject_title', 'Subject_taxid', 'Subject_kingdom', 'Subject_species', 'sallseqid']
    # removing | from end of query/subject
    subhitdata['Query'] = subhitdata['Query'].str.strip('|')
    subhitdata['Subject'] = subhitdata['Subject'].str.strip('|')
    # split up query and subject into constituent columns
    subhitdata[['gi1', 'Query_gi', 'Query_origindb', 'Query_accession']] = subhitdata['Query'].str.split('|',
                                                                                                         expand=True)
    subhitdata[['gi2', 'Subject_gi', 'Subject_origindb', 'Subject_accession', 'spare']] = subhitdata[
        'Subject'].str.split('|', expand=True)
    # calculate query length and subject length
    subhitdata['Query_length'] = subhitdata['Query_end'] - subhitdata['Query_start']
    subhitdata['Subject_length'] = subhitdata['Subject_end'] - subhitdata['Subject_start']
    # remove unwanted columns
    subhits = subhitdata.drop(labels=['spare', 'gi1', 'gi2', 'Query_gi', 'Subject_gi', 'sallseqid'], axis=1)
    # re order columns
    subhits.rename(columns={'staxids': 'Subject_taxid'}, inplace=True)
    # annotating subject with tax info
    subhits = get_tax(subhits, 'subject', 'family')

    # re-ordering and renaming columns, adding prefix 'Sub' to distinguish from primary results data
    colorder = ['Query', 'Subject', 'Percentage', 'Alignment_length', 'Mismatch', 'Gap_open', 'Query_start',
                'Query_end',
                'Subject_start', 'Subject_end', 'Evalue', 'Query_length', 'Subject_length', 'QCOVS', 'Query_accession',
                'Subject_accession', 'Subject_taxid', 'Subject_kingdom', 'Subject_group', 'Subject_title',
                'Subject_family', 'Subject_species', 'Query_origindb', 'Subject_origindb']
    # add sub prefix to column names
    colnames = ['Sub_' + name for name in colorder]
    subhits = subhits[colorder]
    subhits.columns = colnames
    return subhits


# add to database function, makes different tables depending on type of input specified
def add_db(hitsdf, tablename, hittype):
    if hittype == 'sub':
        createquery = f'CREATE TABLE {tablename} (Sub_Query, Sub_Subject, Sub_Percentage, Sub_Alignment_length, Sub_Mismatch, Sub_Gap_open, Sub_Query_start, Sub_Query_end, Sub_Subject_start, Sub_Subject_end, Sub_Evalue, Sub_Query_length, Sub_Subject_length, Sub_QCOVS, Sub_Query_accession, Sub_Subject_accession, Sub_Subject_taxid, Sub_Subject_kingdom, Sub_Subject_group, Sub_Subject_title, Sub_Subject_family, Sub_Subject_species, Sub_Query_origindb, Sub_Subject_origindb), FOREIGN KEY (Sub_Query_accession) REFERENCES AllHits(Subject_accession);'

    if hittype != 'sub':
        createquery = f'CREATE TABLE {tablename} (Query, Subject, Percentage, Alignment_length, Mismatch, Gap_Open, Query_start, Query_end, Subject_start, Subject_end, Evalue, Bitscore, Query_length, Subject_length, QCOVS, QCOVUS, QCOVHSP, Query_accession, Query_taxid, Query_division, Query_name, Subject_accession, Subject_taxid, Subject_division, Subject_name, Subject_family, Subject_species, Query_family, Query_species, Retro, dsDNA, Query_contamination, Subject_contamination, DistanceStart, DistanceEnd, Origin, SubjectCheckCode, SubjectCheckResult, PercentageOverlap);'

    # makes table for results
    def make_table(table, create):
        connection = sqlite3.connect('Hits.db')
        cur = connection.cursor()
        try:
            cur.execute(create)
            connection.commit()
        except sqlite3.OperationalError:
            connection.rollback()
        cur.close()
        connection.close()

    # adds data to table
    def adddata(data, table):
        connection = sqlite3.connect('Hits.db')
        # adds data to database in chunks
        data = [data[i:i + 1000] for i in range(0, len(data), 1000)]
        for part in data:
            part.to_sql(table, connection, if_exists='append', index=False)
            connection.commit()
        connection.close()
    # makes table with create query depending on selected results type (subject or not subject)
    make_table(tablename, createquery)
    adddata(hitsdf, tablename)
    # saves to csv
    hitsdf.to_csv(f'{tablename}.csv', sep='\t', header=True)
    print('added to database!')


# for all refseq sub hits
if not os.path.isfile('AllRefseqSubjectHits.csv'):
    allrefhits = pd.read_csv('allrefsub_hits.txt', sep='\t', header=None)
    allrefsubhits = process_subhits(allrefhits)
    add_db(allrefsubhits, 'AllRefseqSubjectHits', 'sub')

# for 95 all sub hits
if not os.path.isfile('95AllSubjectHits.csv'):
    allhit_1 = pd.read_csv('95allsubs_part1_hits.txt', sep='\t', header=None)
    allhit_2 = pd.read_csv('95allsubs_part2_hits.txt', sep='\t', header=None)
    allhit_3 = pd.read_csv('95allsubs_part3_hits.txt', sep='\t', header=None)
    # join results together
    allsubhits = pd.concat([allhit_1, allhit_2, allhit_3], ignore_index=True)
    # process all sub hits
    allsubhits = process_subhits(allsubhits)
    # add to database
    add_db(allsubhits, '95AllSubjectHits', 'sub')


# adding sub info to results files


# list of words to use for filter function
EVEwords = ['endogenous', 'leukemia', 'sarcoma', 'viral', 'proto-oncogene', 'oncogene', 'leukosis', 'tumour',
            'tumor', 'cytoma', 'carcinoma', 'virus', 'HERV']
unannowords = ['uncharacterized LOC', 'uncharacterized protein', 'predicted gene']

# filter by keyword count function
def filterkey_count(keywords, data):
    # counts those which include key word
    keycount = len(data[data['Sub_Subject_title'].str.contains('|'.join(keywords))])
    return keycount

# function to analyse subject hits and add results to original virus blast hits dataframe
def add_sub_info(subfile, hitfile):
    # using list of present subject ids to filter refseq db
    subhits_no_dupe = subfile.drop_duplicates(subset=['Sub_Query_accession'])
    sub_acclist = subhits_no_dupe['Sub_Query_accession'].tolist()
    # checking if all queries have a match with the list of subject results ids
    hitdf = hitfile[hitfile['Subject_accession'].isin(sub_acclist)]

    # function which is used to evaluate each subject in the blast hits dataframe
    def subjecttrust(index):
        # gets row from the blast hits dataframe and gets sub accession and query accession
        hitdfrow = hitdf.loc[hitdf['Unnamed: 0'] == index]
        sub_accession = hitdfrow['Subject_accession'].values[0]
        query_accession = hitdfrow['Query_accession'].values[0]
        # subsets subject hits results by the accession
        subsetsub = subfile[subfile['Sub_Query_accession'].str.contains(sub_accession)]
        # removes potential false positive by excluding original virus query accession
        subsetsub = subsetsub[~subsetsub['Sub_Subject_accession'].str.contains(query_accession)]

        # evaluates the subjects in the subset of subject hit results
        def subject_eval(subsetdata):
            # get number of subject hits
            # subset hits by percentage identity
            subsetdata = subsetdata.query('Sub_Percentage > 95')
            # get total number of high similarity subjects
            subhitnumber = len(subsetdata)

            # find if sequences overlap in order to evaluate subjects
            def subjectoverlap(subindex):
                overlapout = ''
                # find where the subject start and ends for the original subject-virus hit
                subject_start = int(hitdfrow['Subject_start'].values[0])
                subject_end = int(hitdfrow['Subject_end'].values[0])
                # find where the subject hit start and ends are
                subrow = subsetdata.loc[subsetdata['Unnamed: 0'] == subindex]
                subQStart = subrow['Sub_Query_start'].values[0]
                subQEnd = subrow['Sub_Query_end'].values[0]
                if subject_end < subject_start and subQEnd < subQStart:
                    # print('both reverse here', hitdfrow['Subject_accession'].values[0])
                    pass
                elif subject_end < subject_start:
                    # print('subject is reverse align with query', hitdfrow['Subject_accession'].values[0])
                    subject_end = int(hitdfrow['Subject_start'].values[0])
                    subject_start = int(hitdfrow['Subject_end'].values[0])
                # get the size of subject sequence, and subject hit sequence
                sub = range(subject_start, subject_end)
                subq = range(subQStart, subQEnd)
                # find if the subject sequence and subject hit sequence overlap
                overlaprange = len(list(range(max(sub[0], subq[0]), min(sub[-1], subq[-1]) + 1)))
                if overlaprange == 0:
                    overlapout = 'no overlap'
                if overlaprange > 0:
                    overlapout = 'overlap'
                # return overlap information
                return overlapout

            subsetdata['Overlap'] = subsetdata['Unnamed: 0'].map(subjectoverlap)

            # find out how many of the subject hits overlap with the original virus-subject hit
            subsetdata = subsetdata[subsetdata['Overlap'] == 'overlap']
            overlapnumber = len(subsetdata)

            # analyse the subset of subject hits which overlapped with the original virus-subject hits
            # get number of viruses
            subvirusnumber = len(subsetdata[subsetdata['Sub_Subject_group'].str.contains('viruses')])
            bacteriacount = len(subsetdata[subsetdata['Sub_Subject_group'].str.contains('bacteria')])
            subartnumber = len(subsetdata[subsetdata['Sub_Subject_group'].str.contains('other sequences')])

            # subset to exclude virus
            novirusdata = subsetdata[~subsetdata['Sub_Subject_group'].str.contains('viruses')]
            # get number of EVEs
            subEVEnumber = filterkey_count(EVEwords, novirusdata)
            # using unanno because lack of annotation more important for non-virus sequences
            unannonumber = filterkey_count(unannowords, novirusdata)

            # add all counts to a dictionary to store all the count info
            results = {}
            results.update({'Total': subhitnumber})
            results.update({'Overlap': overlapnumber})
            results.update({'Virus': subvirusnumber})
            results.update({'EVE': subEVEnumber})
            results.update({'Unannotated': unannonumber})
            results.update({'Artificial': subartnumber})
            results.update({'Bacteria': bacteriacount})
            return results

        out = subject_eval(subsetsub)
        return out
    # adds subject check code, containing count info of each subject hit category back to virus-subject df
    hitdf['SubjectCheckCode'] = hitdf['Unnamed: 0'].map(subjecttrust)

    # function which evaluates the counts in the subject code and assigns a category to the subject
    def subjectscore(index):
        code = hitdf.loc[hitdf['Unnamed: 0'] == index, 'SubjectCheckCode'].values[0]
        hitdfrow = hitdf.loc[hitdf['Unnamed: 0'] == index]
        # default assignment is no result
        out = 'No_result'
        # virus contamination assignment
        if code['Virus'] > 0 and code['EVE'] == 0:
            out = 'Virus_contam'
        # EVE contamination assingment
        if code['EVE'] > 0:
            if code['EVE'] + code['Virus'] + code['Unannotated'] > code['Total'] * 0.5:
                # if eve < 2/ low proportion of hits, check if unnanoated + virus
                out = 'EVE'
        # Unannotated assignment
        if code['Unannotated'] > code['Virus'] and code['Unannotated'] > code['EVE']:
            unannowords = ['uncharacterized_LOC', 'uncharacterized_protein', 'predicted_gene']
            # assigns only when the subject name in the blast results is also unannotated
            if len(hitdfrow[hitdfrow['Subject_name'].str.contains('|'.join(unannowords))]) == 1:
                if code['Unannotated'] >= code['Total'] * 0.5:
                    out = 'Unannotated'
        # artificial contmaination assignment
        if code['Artificial'] > 1 and all(
                x < code['Artificial'] for x in (code['Virus'], code['EVE'], code['Unannotated'])):
            out = 'Artificial_contam'
        # bacterial contamiantion assignment
        # checks that subject is not bacteria before seeing if it has problematic hits with bacteria
        if len(hitdfrow[hitdfrow['Subject_division'].str.contains('BACTERIA')]) == 0:
            if code['Bacteria'] > 0 and all(
                    x < code['Bacteria'] for x in (code['Virus'], code['EVE'], code['Artificial'])):
                out = 'Bacteria_contam'
        return out
    # adds subject check result ot virus-subject dataframe
    hitdf['SubjectCheckResult'] = hitdf['Unnamed: 0'].map(subjectscore)

    # calculate percentage of many of the subject hits overlap with the original virus-subject hit
    def get_percent_overlap(index):
        code = hitdf.loc[hitdf['Unnamed: 0'] == index, 'SubjectCheckCode'].values[0]
        hitdfrow = hitdf.loc[hitdf['Unnamed: 0'] == index]
        percentageoverlap = code['Overlap'] / code['Total'] * 100
        return percentageoverlap
    # adds percentage overlap info to virus-subject dataframe
    hitdf["PercentageOverlap"] = hitdf["Unnamed: 0"].map(get_percent_overlap)

    return hitdf

# checks if processed subject hit files exist and coressponding checked blast results do not exist
if os.path.isfile('AllRefseqSubjectHits.csv') and not os.path.isfile('AllRefseqSubChecked.csv'):
    print('found all refseq sub file!')
    allrefsubhits = pd.read_csv('AllRefseqSubjectHits.csv', sep='\t')
    # read in results file
    allrefseq = pd.read_csv('RefseqHits.csv', sep='\t')
    # adds subject annotation info to blast results file
    annorefseq = add_sub_info(allrefsubhits, allrefseq)
    # remove the unnamed:0 column, remove unwanted columns
    annorefseq = annorefseq.drop(labels=['Unnamed: 0'], axis=1)
    # convert dict column to string
    annorefseq['SubjectCheckCode'] = annorefseq['SubjectCheckCode'].astype(str)
    # adds to database
    add_db(annorefseq, 'AllRefseqSubChecked', 'nosub')

# same process as above, but with 95% identity blast results
if os.path.isfile('95AllSubjectHits.csv') and not os.path.isfile('95AllHitsSubChecked.csv'):
    print('found 95 all subject hits file!')
    all_hits = pd.read_csv('AllHits.csv', sep='\t')
    all_hits = all_hits.query('Percentage >95')
    allsubhits = pd.read_csv('95AllSubjectHits.csv', sep='\t')
    annoallseq_95 = add_sub_info(allsubhits, all_hits)
    annoallseq_95 = annoallseq_95.drop(labels=['Unnamed: 0'], axis=1)
    # convert dict column to string
    annoallseq_95['SubjectCheckCode'] = annoallseq_95['SubjectCheckCode'].astype(str)
    add_db(annoallseq_95, '95AllHitsSubChecked', 'nosub')

# combining all subject hit results and adding to dataframe
if not os.path.isfile('AllSubjectHits.csv'):
    print('found all subject hits')
    allrefsubhits = pd.read_csv('AllRefseqSubjectHits.csv', sep='\t')
    allsubhits95 = pd.read_csv('95AllSubjectHits.csv', sep='\t')
    # merge together
    merged = pd.concat([allrefsubhits, allsubhits95], ignore_index=True)
    # remove "index" column to accurately drop duplicate rows
    merged = merged.drop('Unnamed: 0', axis=1)
    merged = merged.drop_duplicates().reset_index(drop=True)
    add_db(merged, 'AllSubjectHits', 'sub')
    merged.to_csv('AllSubjectHits.csv', sep='\t', header=True)

#  using all subject hits to annotate the all hits blast results (can annotate more hits with less than 95% identity)
if os.path.isfile('AllSubjectHits.csv') and not os.path.isfile('AllHitsSubChecked.csv'):
    print('found all subject hits, now sub checking')
    allsubjecthits = pd.read_csv('AllSubjectHits.csv', sep='\t')
    # run this with all hits
    all_hits = pd.read_csv('AllHits.csv', sep='\t')
    annoallseq = add_sub_info(allsubjecthits, all_hits)
    print('preparing to add to db!')
    annoallseq = annoallseq.drop(labels=['Unnamed: 0'], axis=1)
    # convert dict column to string
    annoallseq['SubjectCheckCode'] = annoallseq['SubjectCheckCode'].astype(str)
    add_db(annoallseq, 'AllHitsSubChecked', 'nosub')
