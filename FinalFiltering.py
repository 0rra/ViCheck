import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb
import plotly.io as pio
pio.renderers.default = 'png'

# read in all results and all checked results
allhits = pd.read_csv("AllHits.csv", sep="\t")
checkedhits = pd.read_csv("AllHitsSubChecked.csv", sep="\t")

# get unique subject accessions in allhits
subid = allhits['Subject_accession'].drop_duplicates().tolist()
hitsub = checkedhits['Subject_accession'].drop_duplicates().tolist()
uniquesubs = list(set(subid) - set(hitsub))
# subset all hits by unique subs
leftoverhits = allhits[allhits["Subject_accession"].isin(uniquesubs)]

# add leftover hits to all hits
everything = pd.concat([checkedhits, leftoverhits], ignore_index=True)

# function used to flag contamination/problematic category to each hit
def flagcontamtype(index):
    hitdfrow = everything.loc[everything['Unnamed: 0'] == index]
    # flagging monkeypox and suncus bax motif homolog hit
    if hitdfrow["Query_species"].values[0] == "Monkeypox virus" and hitdfrow["Subject_accession"].values[0] == "XM_049782773.1":
        return "Mpox homolog"
    # flagging subject annotation results
    if hitdfrow["SubjectCheckResult"].values[0] == "Bacteria_contam":
        return "Bacteria_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "Virus_contam":
        return "Virus_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "Artificial_contam":
        return "Artificial_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "EVE":
        return "EVE"
    # flagging ribosomal contam hits
    ribowords = ['ribosomal_RNA', 'rRNA']
    if hitdfrow["Subject_name"].str.contains('|'.join(ribowords)).any():
        return "Ribosome_contam"
    #if hitdfrow["SubjectCheckResult"].values[0] == "Unannotated":
    #    return "Unannotated"
    # flagging mag hits
    if hitdfrow["Query_name"].str.contains('|'.join(["_MAG_","_MAG:"])).any():
        return "MAG"
    # flagging retroviruses
    if hitdfrow["Retro"].values[0]:
        return "Retrovirus"
    else:
        return 'Unknown'

# flagging based on contamination categories identified so far
everything['Flag'] = everything['Unnamed: 0'].map(flagcontamtype)

# plotting the counts of each flagged category
def plottypes(dataa, order):
    plt.figure(figsize=(6, 3))
    sb.set_palette(sb.color_palette("husl"))
    # option to order bars (choice depended on flag function used)
    if order == 'yes':
        plot = sb.countplot(data=dataa, y='Flag', order =['Unknown', 'Mpox homolog', 'EVE','Virus_contam', 'Bacteria_contam','Artificial_contam','Ribosome_contam', 'MAG'])
    if order != 'yes':
        plot = sb.countplot(data=dataa, y='Flag', order =['Unknown', 'Homolog pair', 'Unannotated','EVE','Virus_contam', 'Bacteria_contam','Artificial_contam','Ribosome_contam','HighPercentageHit','Polydnaviriformidae','Geminiviridae','Mimiviridae','MAG'])
    plot.bar_label(plot.containers[0])
    plt.tight_layout()
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    fig2 = plot.get_figure()
    fig2.show()
plottypes(everything, 'yes')


# updated flag everything function, including all identified homology hits, and other identified groups
def flagmorecontam(index):
    hitdfrow = everything.loc[everything['Unnamed: 0'] == index]
    # flagging homolog pairs
    # flags all bax hits in poxviridae family
    if hitdfrow["Query_family"].values[0] == "Poxviridae" and hitdfrow["Subject_accession"].values[0] == "XM_049782773.1":
        return "Homolog pair"
    if hitdfrow["Query_family"].values[0] == "Poxviridae" and hitdfrow["Subject_name"].str.contains('|'.join(["ribonucleoside_reductase", "ribonucleotide_reductase","ribonucleoside-diphosphate"])).any():
        return "Homolog pair"
    if hitdfrow["Query_family"].values[0] == "Herpesviridae" and hitdfrow["Subject_name"].str.contains("thymidylate").any():
        return "Homolog pair"
    if hitdfrow["Query_family"].values[0] == "Iridoviridae" and hitdfrow["Subject_name"].str.contains("thymidylate").any():
        return "Homolog pair"
    if hitdfrow["Query_family"].values[0] == "Flaviviridae" and hitdfrow["Subject_name"].str.contains("DnaJ_heat_shock_protein").any():
        return "Homolog pair"
    if hitdfrow["Query_family"].values[0] == "Flaviviridae" and hitdfrow["Subject_name"].str.contains('|'.join(["NEDD8", "ubiquitin","ribosomal_protein", "GABA", "GATE", "light_chain_3"])).any() and hitdfrow["Subject_family"].values[0] == "Bovidae":
        return "Homolog pair"
    # flagging mimiviridae hits
    if hitdfrow["Query_family"].values[0] == "Mimiviridae" and hitdfrow["Subject_species"].values[0] == "Daktulosphaira vitifoliae":
        return "Mimiviridae"
    # flagging subject annotation results
    if hitdfrow["SubjectCheckResult"].values[0] == "Bacteria_contam":
        return "Bacteria_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "Virus_contam":
        return "Virus_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "Artificial_contam":
        return "Artificial_contam"
    if hitdfrow["SubjectCheckResult"].values[0] == "EVE":
        return "EVE"
    # flagging subject annotation result (unannotated/ ambiguous hit)
    if hitdfrow["SubjectCheckResult"].values[0] == "Unannotated":
        return "Unannotated"
    # flagging ribosomes
    ribowords = ['ribosomal_RNA', 'rRNA']
    if hitdfrow["Subject_name"].str.contains('|'.join(ribowords)).any():
        return "Ribosome_contam"
    # flagging mags
    if hitdfrow["Query_name"].str.contains('|'.join(["_MAG_","_MAG:"])).any():
        return "MAG"
    # flag retrovirus
    if hitdfrow["Retro"].values[0]:
        return "Retrovirus"
    # flagging families with EVE-like relationships
    if hitdfrow["Query_family"].values[0] == "Polydnaviriformidae" and hitdfrow["Subject_family"].values[0] == "Braconidae":
        return "Polydnaviriformidae"
    if hitdfrow["Query_family"].values[0] == "Geminiviridae" and hitdfrow["Subject_family"].values[0] == "Solanaceae":
        return "Geminiviridae"
    # flagging remaining high percentage, possible contaminated virus hits
    if hitdfrow["Percentage"].values[0] > 95 and hitdfrow["Query_contamination"].isin(['Start', 'End']).any():
        return "HighPercentageHit"
    else:
        return 'Unknown'

everything['Flag'] = everything['Unnamed: 0'].map(flagmorecontam)

plottypes(everything, 'no')

# split into problematic virus and problematic subject
problemvirus = everything[everything["Flag"].isin(["Homolog pair", "MAG", "Ribosome_contam", "HighPercentageHit"])]
problemsubject = everything[everything["Flag"].isin(["EVE","Bacteria_contam","Virus_contam","Artificial_contam"])]
def addaction(data):
    def suggestaction(index):
        hitdfrow = data.loc[data['Unnamed: 0'] == index]
        if hitdfrow["Flag"].values[0] == "Homolog pair":
            return "Flagged"
        if hitdfrow["Flag"].values[0] == "MAG" and hitdfrow["Percentage"].values[0] > 95:
            return "Remove"
        if hitdfrow["Flag"].values[0] == "MAG" and hitdfrow["Percentage"].values[0] <= 95:
            return "Flagged"
        if hitdfrow["Flag"].values[0] == "Ribosome_contam":
            return "Remove"
        if hitdfrow["Flag"].values[0] == "HighPercentageHit":
            return "Remove"
        if hitdfrow["Flag"].values[0] == "EVE":
            return "Flagged"
        if hitdfrow["Flag"].values[0] == "Virus_contam":
            return "Flagged"
        if hitdfrow["Flag"].values[0] == "Bacteria_contam":
            return "Remove"
        if hitdfrow["Flag"].values[0] == "Artificial_contam":
            return "Remove"
    data["Action"] = data["Unnamed: 0"].map(suggestaction)
    return data

problemvirus = addaction(problemvirus)
# save problematic virus to csv
problemvirus = problemvirus.drop(labels=['Unnamed: 0'], axis=1)
# drop columns except query accession, subject accession, query taxid, subject_taxid, origin, subject check code, subjectcheck result, flag, action
selectcolumns = ["Query_accession","Subject_accession","Query_taxid","Subject_taxid","Origin","SubjectCheckCode","SubjectCheckResult","Flag","Action"]
problemvirus = problemvirus[selectcolumns]
problemvirus.to_csv('ProblematicVirus.csv', sep = "\t", header=True)

problemsubject = addaction(problemsubject)
# save problematic subjects to csv
problemsubject = problemsubject.drop(labels=['Unnamed: 0'], axis=1)
problemsubject = problemsubject[selectcolumns]
problemsubject.to_csv('ProblematicNonVirus.csv', sep="\t", header=True)



# investigating the leftover hits unknown:

uncat = everything[everything["Flag"] == "Unknown"]
# plot percentage
def alignlength(dataa, bins, title):
    #plt.figure(1)
    alignplot = sb.histplot(data = dataa, x = 'Alignment_length', binwidth=bins, element='step')
    #plt.yscale('log')
    plt.title(title)
    plt.tight_layout()
    alignplot.spines['right'].set_visible(False)
    alignplot.spines['top'].set_visible(False)
    fig3 = alignplot.get_figure()
    plt.show()

alignlength(uncat, 200, 'Uncatergorised alignment lengths')

# plotting bar to see query and subejct combinations in leftover hits

def divisioncombobar(dataa):
    combo = dataa.groupby(['Query_division', 'Subject_division']).size().reset_index()
    combo.rename(columns={0: 'Count'}, inplace=True)
    # to make labels combine query and subject columns
    combo['Combination'] = combo[['Query_division', 'Subject_division']].agg('-'.join, axis=1)
    combo = combo.sort_values("Count", ascending=False)
    comboplot = sb.barplot(data=combo, y='Combination', x='Count', palette="crest")
    comboplot.bar_label(comboplot.containers[0])
    comboplot.spines['right'].set_visible(False)
    comboplot.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.show()

divisioncombobar(uncat)

# check percentage identity of leftover hits
def percentident(dataa):
    #plt.figure(1)
    percentplot = sb.histplot(data=dataa, x='Percentage', element='step')
    #plt.yscale('log')
    plt.tight_layout()
    percentplot.spines['right'].set_visible(False)
    percentplot.spines['top'].set_visible(False)
    fig3 = percentplot.get_figure()
    plt.show()

percentident(uncat)
# filter high percentage
highpercent = uncat[uncat["Percentage"] >= 95]
lowpercent = uncat[uncat["Percentage"] < 95]
# see division combos for remaining high percent hits
divisioncombobar(highpercent)

# get largest families in leftover hits (high percentage and low percentage hits)
def famcount(dataa):
    famcount = dataa.groupby(['Query_species']).size().reset_index()
    famcount.rename(columns={0: 'Count'}, inplace=True)
    print(f"Hits with family not present: {dataa['Subject_species'].str.contains('Not present').sum()}")
    famcount = famcount.sort_values('Count', ascending= False)
    return famcount

highfam = famcount(highpercent)
lowfam = famcount(lowpercent)


