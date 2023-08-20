import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
pio.renderers.default = 'png'


# set file name here
file = 'AllHits.csv'
dataname = 'All'
hitdf = pd.read_csv(file, sep="\t")

hitdffil = hitdf.query('Percentage >95')


# gets proportions of different virus counts in dataset
def virus_prop(data, name):
    Name = name
    total = len(data)
    # count number of phage
    phageno = sum(data['Query_division'] == 'PHAGE')
    virusno = sum(data['Query_division'] == 'VIRUS')
    Retrono = sum(data['Retro'] == True)
    Nonretro = virusno - sum(data['Retro'] == True)
    dict = {'Name': Name, 'Total': total, 'Phage': phageno, 'Virus': virusno, 'Retrovirus': Retrono, 'Non-retrovirus': Nonretro}
    return dict


# gets virus counts for genbank and refseq sequences
gendict = virus_prop(hitdf[hitdf["Origin"] == "GenBank"], 'Genbank')
refseqdict = virus_prop(hitdf[hitdf["Origin"] == "RefSeq"], 'RefSeq')
totaldict = virus_prop(hitdf, 'All')

dictlist = [refseqdict, gendict, totaldict]
countdf = pd.DataFrame(dictlist)


# makes stacked barchart of virus types
def plot_virusprop(data):
    data["Phage %"] = round(data["Phage"] / data["Total"] * 100)
    data["Retrovirus %"] = round(data["Retrovirus"] / data["Total"] * 100)
    data["Non-Retrovirus %"] = round(data["Non-retrovirus"] / data["Total"] * 100)

    fig = go.Figure()
    fig.add_trace(go.Bar(y=data["Phage %"], x= data["Name"], name = "Phage %", text = data["Phage %"],
                         marker=dict(color='#FE4C44',line=dict(color='#FE4C44', width=0.05))))
    fig.add_trace(go.Bar(y=data["Retrovirus %"], x= data["Name"], name = "Retrovirus %", text= data["Retrovirus %"],
                         marker=dict(color='#fFD500',line=dict(color='#fFD500', width=0.05))))
    fig.add_trace(go.Bar(y=data["Non-Retrovirus %"], x= data["Name"], name = "Non-Retrovirus %", text = data["Non-Retrovirus %"],
                         marker=dict(color='#52a4ee',line=dict(color='#52a4ee', width=0.05))))
    fig.update_layout(yaxis=dict(title_text="%"), width=650, height=640,plot_bgcolor='white', barmode='stack')
    fig.update_layout(uniformtext_minsize=20, uniformtext_mode='hide', font=dict(size=24, color = 'black'))
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', tickfont=dict(size=21))
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', tickfont=dict(size=21))
    fig.update_layout(title={'text': "BLAST Hits",'y': 0.03 ,'x': 0.35,'xanchor': 'center', 'font' :dict(size=20)})
    pio.write_image(fig, file='bar.png', format='png', width=650, height=640, scale=5)
    #fig.show(renderer = "browser")


plot_virusprop(countdf)


# histogram  of alignment lengths
def alignlength(dataa, bins, title):
    #plt.figure(1)
    plt.figure(figsize=(8, 6))
    plt.rcParams['font.size'] = 20
    alignplot = sb.histplot(data = dataa, x = 'Alignment_length', binwidth=bins, element='step')
    #plt.yscale('log')
    plt.title(title)
    plt.tight_layout()

    fig3 = alignplot.get_figure()
    plt.show()

# plot all alignments
alignlength(hitdf, 200, 'All alignment lengths')
# long alignment
longhits = hitdf[hitdf['Alignment_length'] > 2500]
alignlength(longhits, 500, 'Long alignment lengths')
shorterhits = hitdf[hitdf['Alignment_length'] < 2500]
print(len(shorterhits)/ len(hitdf) *100)
print(len(longhits))
print(len(hitdf[hitdf['Alignment_length'] < 500])/ len(hitdf) *100)


# counts of query family making hits
def famcount(dataa):
    famcount = dataa.groupby(['Query_family']).size().reset_index()
    famcount.rename(columns={0: 'Count'}, inplace=True)
    print(f"Hits with family not present: {dataa['Query_family'].str.contains('Not present').sum()}")
    famcount = famcount.sort_values('Count', ascending= False)
    return famcount


# plots top families
def get_famplot(dataa, titlename):
    colors = ["#30738d", "#d07ba7", "#7bb8d0"]
    sb.set_palette(sb.color_palette(colors))
    plot = sb.barplot(data=dataa, y='Query_family', x='Count')
    plot.bar_label(plot.containers[0])
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.title(titlename)
    plt.tight_layout()
    fig2 = plot.get_figure()
    fig2.show()

topfamilies = famcount(hitdf)
# plot top 7 families in all hit databse
get_famplot(topfamilies[0:6], 'Top families - All Hits')

# plot high percentage hit families
topfamilies95 = famcount(hitdf[hitdf["Percentage"] >95])

# get top families in refseq and genbank hits
reftopfam = famcount(hitdf[hitdf["Origin"] == "RefSeq"])
gentopfam = famcount(hitdf[hitdf["Origin"] == "GenBank"])

# bar plot of query vs subject division combinations
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

# pie chart to show the subject divisions of the virus hits and the phage hits
def divisioncombo(dataa):
    # grouping data by subject and query division combinations to get counts of most common combinations
    combo = dataa.groupby(['Query_division', 'Subject_division']).size().reset_index()
    combo.rename(columns={0: 'Count'}, inplace=True)
    # setting colours for subject divisoin type
    colourmap = {"MAMMAL": "orange", "INVERTEBRATE": "#b8e9ff", "BACTERIA": "#262626", "RODENT": "#ff7863", "PRIMATES": "#6c00a3",
                 "VERTEBRATE": "#d390ff", "PLANTS_&_FUNGI": "#3cab89"}
    combo["Colour"] = combo["Subject_division"].apply(lambda x: colourmap.get(x))

    # plot two pie charts together
    fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'domain'}, {'type': 'domain'}]])
    # plot phage and subject pie chart
    fig.add_trace(go.Pie(labels=combo[combo["Query_division"] =="PHAGE"]['Subject_division'],
                         values=combo[combo["Query_division"] =="PHAGE"]["Count"],name="PHAGE", sort =True), 1, 1)
    # plot virus and subject pie chart
    fig.add_trace(go.Pie(labels=combo[combo["Query_division"] =="VIRUS"]['Subject_division'],
                         values=combo[combo["Query_division"] == "VIRUS"]["Count"], name="VIRUS", sort = True), 1, 2)
    # update the colours of the piechart
    fig.update_traces(marker=dict(colors=combo["Colour"]))
    # change the size of font and colour to black
    fig.update_layout(uniformtext_minsize=25, uniformtext_mode='hide', legend=dict(font=dict(size=25, color='black')))
    # show percentage of subject divisions in total virus/phage hits
    fig.update_traces(textinfo ='percent')
    # change plot size and add titles to each pie chart
    fig.update_layout(width=1100, height=900, annotations=[dict(text="Phage", x = 0.19, y=0.95, font_size = 25, showarrow=False), dict(text='Virus', x=0.81, y=0.95, font_size = 25, showarrow=False)])
    #fig.show(renderer = "browser")
    #save image
    pio.write_image(fig, file='pie.png', format='png', width=1100, height=900, scale=5)

divisioncombo(hitdf)


# histogram for percentage identity
def percentident(dataa, titlename):
    #plt.figure(1)
    plt.rcParams['font.size'] = 17
    plt.figure(figsize=(6, 6))
    percentplot = sb.histplot(data=dataa, x='Percentage', element='step')
    #plt.yscale('log')
    plt.title(titlename)
    plt.tight_layout()
    percentplot.spines['right'].set_visible(False)
    percentplot.spines['top'].set_visible(False)

    fig3 = percentplot.get_figure()
    plt.show()

# percentage identities of the refseq and genbank hits
percentident(hitdf[hitdf['Origin']=='RefSeq'], 'RefSeq hits')
percentident(hitdf[hitdf['Origin']=='GenBank'], 'GenBank NT hits')
print(hitdf['Percentage'].mean())
refseq = hitdf[hitdf['Origin']=='RefSeq']
genbank =hitdf[hitdf['Origin']=='GenBank']
print(refseq['Percentage'].mean())
print(genbank['Percentage'].mean())
print(len(refseq[refseq['Percentage'] >90])/ len(refseq) * 100)
print(len(genbank[genbank['Percentage'] >90])/ len(genbank) * 100)

## 90 peak investigation
# filtering dataset
genbankpeak = genbank[(genbank['Percentage'] > 89.5) & (genbank['Percentage'] < 90.5)]

percentident(genbankpeak, 'GenBank 90% peak')

# genbank peak species vs species, see what the hits are between
def spcount(dataa):
    speciescount = dataa.groupby(['Query_species']).size().reset_index()
    speciescount.rename(columns={0: 'Count'}, inplace=True)
    print(f"Hits with species not present: {dataa['Query_species'].str.contains('Not present').sum()}")
    speciescount = speciescount.sort_values('Count', ascending= False)
    return speciescount

def get_spplot(dataa, titlename):
    plt.figure(figsize=(8, 6))
    plot = sb.barplot(data=dataa, y='Query_species', x='Count')
    plot.bar_label(plot.containers[0])
    plt.title(titlename)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.tight_layout()
    fig2 = plot.get_figure()
    fig2.show()
# find main virus species in peak
toppeak = spcount(genbankpeak)
get_spplot(toppeak[0:6], 'GenBank Peak virus species')

# plot virus and subject species combinations
def speciescombo(dataa):
    plt.figure(figsize=(8, 6))
    combo = dataa.groupby(['Query_species', 'Subject_species']).size().reset_index()
    combo.rename(columns={0: 'Count'}, inplace=True)
    # .reset_index().rename(columns={0:'count'})
    print(combo.head())
    # to make labels combine query and subject columns
    combo['Combination'] = combo[['Query_species', 'Subject_species']].agg('-'.join, axis=1)
    # plot combination
    comboplot = sb.barplot(data=combo, y='Combination', x='Count')
    comboplot.bar_label(comboplot.containers[0])
    comboplot.spines['right'].set_visible(False)
    comboplot.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.show()

# find what subject species is hitting with monkeypox virus
speciescombo(genbankpeak[genbankpeak['Query_species'] == 'Monkeypox virus'])

# filtered genbank nt results, removing monkey pox virus
filtergen = genbank[genbank['Query_species'] != 'Monkeypox virus']
# plotting percentage identity plot after removing monkeypox
percentident(filtergen, 'GenBank NT hits filtered')

# bar plot of contamination position type
def contamplot(data, fignum):
    plt.figure(fignum)
    plt.figure(figsize=(6, 6))
    conttlot = sb.countplot(data=data, x = 'Query_contamination', order =['Unclassified','Start','End'])
    conttlot.bar_label(conttlot.containers[0])
    plt.tight_layout()
    conttlot.spines['right'].set_visible(False)
    conttlot.spines['top'].set_visible(False)
    plt.show()

contamplot(hitdf, 5)
contamplot(hitdffil,6)

# ribosomal hits
allhits = pd.read_csv('AllHits.csv', sep='\t')
ribowords = ['ribosomal_RNA', 'rRNA']

def filterkey_data(keywords, data, column):
    #includes for those which include key word
    keydata = data[data[column].str.contains('|'.join(keywords))]
    return keydata

ribohits = filterkey_data(ribowords, allhits, 'Subject_name')
# histogram
alignlength(ribohits, 100, 'Distribution of alignment lengths for rRNA subject hits')
contamplot(ribohits, 5)
percentident(ribohits, '')
allhits95 = pd.read_csv('95AllHitsSubChecked.csv', sep='\t')
ribohits95 = filterkey_data(ribowords, allhits95, 'Subject_name')
topribohits = ribohits95[ribohits95['SubjectCheckResult'].isin(['No_result','Unannotated'])]
contamplot(topribohits, 5)

# identifying EVES
# list of words to filter
EVEwords = ['endogenous', 'leukemia', 'sarcoma', 'viral', 'proto-oncogene', 'oncogene', 'leukosis', 'tumour', 'tumor',
            'cytoma', 'carcinoma', 'virus', 'HERV']
unannowords = ['uncharacterized LOC', 'uncharacterized protein', 'predicted gene']

evehits = filterkey_data(EVEwords, allhits, 'Subject_name')
unannohits = filterkey_data(unannowords, allhits, 'Subject_name')

# human pathogenic virus hits
priorityvirus = pd.read_csv('WHOpriorityvirus.txt', sep='\t')
priorityfamily = priorityvirus['Family'].tolist()
# subset hits if they have priority family
priorityhits = allhits[allhits['Query_family'].isin(priorityfamily)]
# identify ribosome hits in pathogenic viruses
priortiyribohits = filterkey_data(ribowords, priorityhits, 'Subject_name')
# how many unique species names
uniquespnames = len(priortiyribohits["Query_species"].unique())
# found mag sequences
MAGs = filterkey_data(["_MAG_", "_MAG:"], allhits, "Query_name")


