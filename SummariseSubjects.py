import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb
import plotly.io as pio
pio.renderers.default = 'png'


checkedhits = pd.read_csv('AllHitsSubChecked.csv', sep='\t')

# make query type column for retrovirus, non-retrovirus, phage
def retromap(index):
    hitdfrow = checkedhits.loc[checkedhits['Unnamed: 0'] == index]
    if hitdfrow["Query_division"].values[0] == "VIRUS" and hitdfrow["Retro"].values[0]:
        return "Retrovirus"
    if hitdfrow["Query_division"].values[0] == "VIRUS" and not hitdfrow["Retro"].values[0]:
        return "Non-Retrovirus"
    if hitdfrow["Query_division"].values[0] == "PHAGE":
        return "Phage"
# uses mapping function to go through each row in dataframe and identify virus type
checkedhits["Query_type"] = checkedhits["Unnamed: 0"].map(retromap)

# grouped bar plot of the assigned subject catergories counts for retroviruses, phages, and non-retroviruses
def plot_groupedSubCat(countdata):
    colors = ["grey", "silver", "violet", "red", "dodgerblue", "gold", "#30738d"]
    sb.set_palette(sb.color_palette(colors))
    plt.rcParams['font.size'] = 15.5
    plt.figure(figsize=(12, 7))
    plot = sb.countplot(data=countdata, x='Query_type', hue="SubjectCheckResult", hue_order=['No_result','Unannotated','EVE','Virus_contam','Bacteria_contam','Artificial_contam'])
    for i in plot.containers:
        plot.bar_label(i)
    plot.set_yscale("log")
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.show()
plot_groupedSubCat(checkedhits)

#re set colour palette
sb.set_palette(sb.color_palette("tab10"))

# subset results by subject category and do break down analysis of each type
# analysing bacteria contam hits, how many are with phages
bacteriahits = checkedhits[checkedhits['SubjectCheckResult'] == 'Bacteria_contam']
bactphage = bacteriahits[bacteriahits['Query_division'] == 'PHAGE']

# division plot
def familycombo(dataa):
    combo = dataa.groupby(['Query_family', 'Subject_family']).size().reset_index()
    combo.rename(columns={0: 'Count'}, inplace=True)
    # .reset_index().rename(columns={0:'count'})
    print(combo.head())
    # to make labels combine query and subject columns
    combo['Combination'] = combo[['Query_family', 'Subject_family']].agg('-'.join, axis=1)
    # plot combination
    comboplot = sb.barplot(data=combo, y='Combination', x='Count')
    comboplot.bar_label(comboplot.containers[0])
    plt.tight_layout()
    plt.show()

# family plot of top subjecs
def famcount(dataa):
    famcount = dataa.groupby(['Subject_species']).size().reset_index()
    famcount.rename(columns={0: 'Count'}, inplace=True)
    famcount = famcount.sort_values('Count', ascending= False)
    return famcount

def get_famplot(dataa, titlename):
    plot = sb.barplot(data=dataa, y='Subject_species', x='Count')
    plot.bar_label(plot.containers[0])
    plt.title(titlename)
    plt.tight_layout()
    fig2 = plot.get_figure()
    fig2.show()

# top families of subjects which are contaminated with bacteria and hitting phages
topfamiliesbactphage = famcount(bactphage)
# plot top families
get_famplot(topfamiliesbactphage[0:6], 'Top bacteria contam families with phages')

# analysis of virus contaminted subjects
virusits = checkedhits[checkedhits['SubjectCheckResult'] == 'Virus_contam']
topfamiliesvi= famcount(virusits)
get_famplot(topfamiliesvi[0:6], 'Top virus contaminated subject families')

# penaeus monodon and parvovirus example, investigation into results
prawnhits = checkedhits[checkedhits['Subject_family'] == 'Penaeidae']

def alignlength(dataa, bins, title):
    #plt.figure(1)
    contamplot = sb.histplot(data = dataa, x = 'Alignment_length', binwidth=bins, kde=True, element = "step", multiple = "stack", hue = 'SubjectCheckResult')
    #plt.xscale('log')
    plt.title(title)
    plt.tight_layout()
    fig3 = contamplot.get_figure()
    plt.show()

#alignlength(prawnhits, 50, 'Comparison of alignment lenghts of prawn hits')
#longhits = prawnhits[prawnhits['Alignment_length'] > 200]
#alignlength(longhits, 500, 'Long alignment lengths')
#shorterhits = prawnhits[prawnhits['Alignment_length'] < 200]
#alignlength(shorterhits, 10, 'Short alignment lengths')
#print(len(shorterhits)/ len(prawnhits) *100)
#print(len(longhits))
#print(len(prawnhits[prawnhits['Alignment_length'] < 500])/ len(prawnhits) *100)

# investigation into artificial contamination
artefacts = checkedhits[checkedhits['SubjectCheckResult'] == 'Artificial_contam']
topfamiliesart = famcount(artefacts)
get_famplot(topfamiliesart[0:6], 'Top artificial contamination families')

# top species combo
def spcombocount(dataa):
    spcount = dataa.groupby(['Subject_species', 'Query_species']).size().reset_index()
    spcount.rename(columns={0: 'Count'}, inplace=True)
    spcount = spcount.sort_values('Count', ascending= False)
    return spcount
# find top virus and subject combinations with artefact contamination
topcombo = spcombocount(artefacts)

# biggest EVE species
eve = checkedhits[checkedhits['SubjectCheckResult'] == 'EVE']
topfamilieseve = famcount(eve)
get_famplot(topfamilieseve[0:6], 'Top EVE families')

# biggest unannotated species
unanno = checkedhits[checkedhits['SubjectCheckResult'] == 'Unannotated']
topfamiliesunanno = famcount(unanno)
get_famplot(topfamiliesunanno[0:6], 'Top Unannotated families')


# biggest no result species
nores = checkedhits[checkedhits['SubjectCheckResult'] == 'No_result']
topfamiliesnores = famcount(nores)
get_famplot(topfamiliesnores[0:6], 'Top no result families')

# compare all hits subject checked with all hits to see how many of the hits lack subject annotation
allhits = pd.read_csv("AllHits.csv", sep = '\t')
subid = allhits['Subject_accession'].drop_duplicates().tolist()
hitsub = checkedhits['Subject_accession'].drop_duplicates().tolist()
uniquesubs = list(set(subid) - set(hitsub))
# subset all hits by unique subs
leftoverhits = allhits[allhits["Subject_accession"].isin(uniquesubs)]

# plot percentage identity and alignment
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

def percentident(dataa, multi):
    #plt.figure(1)
    if multi == 'yes':
        percentplot = sb.histplot(data = dataa, x = 'Percentage', element='step', hue = 'Origin',
                                  palette=['limegreen','blue'], binwidth = 1)
    if multi == 'no':
        percentplot = sb.histplot(data=dataa, x='Percentage', element='step')
    #plt.yscale('log')
    plt.tight_layout()
    percentplot.spines['right'].set_visible(False)
    percentplot.spines['top'].set_visible(False)
    fig3 = percentplot.get_figure()
    plt.show()

def contamplot(data, fignum):
    plt.figure(fignum)
    conttlot = sb.countplot(data=data, x = 'Subject_contamination', order =['Unclassified','Start','End'])
    conttlot.bar_label(conttlot.containers[0])
    plt.tight_layout()
    plt.show()

# plots to look at leftover hits

percentident(leftoverhits, 'no')
alignlength(leftoverhits,200,"leftover hits alignment length")
contamplot(leftoverhits, 4)

