import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sb
import os
from GetTax import get_tax
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = 'png'

# function to flag retroviruses
def flag_retro(data, retrovirusefile):
    retrolist = []
    with open(retrovirusefile) as retrofile:
        for row in retrofile:
            row = row.strip('\n')
            retrolist.append(row)
    data['Retro'] = data['Query_taxid'].isin(retrolist)
    return data

# processes fasta headers of refseq and genbank viral databases
def process_head(file):
    headdfchunks = pd.read_csv(file, sep='\t', header=None, chunksize = 8000)
    result = pd.DataFrame()
    # processes in chunks due to the large size of some of the files
    for chunk in headdfchunks:
        chunk.rename(columns={0: 'Query'}, inplace=True)
        chunk['Query_accession'] = chunk.Query.str.split('|', expand=True)[0]
        print(chunk['Query_accession'])
        chunk['Query_taxid'] = chunk.Query.str.split('|', expand=True)[1]
        chunk['Query_division'] = chunk.Query.str.split('|', expand=True)[2]
        chunk['Query_name'] = chunk.Query.str.split('|', expand=True)[3]
        chunk['Query_accession'] = chunk['Query_accession'].str.replace('>', '')
        # flags retroviruses
        chunk = flag_retro(chunk, 'retroviruses.txt')
        # adds family taxonomy
        chunk = get_tax(chunk, 'query', 'family')
        result = pd.concat([chunk, result], ignore_index=True)
    return result

# checks if processed files exist, if not processes
if not os.path.isfile('Process_refseq_head.csv'):
    refhead = process_head('refseq_head.txt')
    refhead.to_csv(f'Process_refseq_head.csv', sep='\t', header=True)

if not os.path.isfile('Process_partial_head.csv'):
    parthead = process_head('partial_head.txt')
    parthead.to_csv(f'Process_partial_head.csv', sep='\t', header=True)

if not os.path.isfile('Process_complete_head.csv'):
    comphead = process_head('complete_head.txt')
    comphead.to_csv('Process_complete_head.csv', sep='\t', header=True)


# function to count virus types, returns counts in dictioanry
def virus_prop(data, name):
    Name = name
    total = len(data)
    phageno = sum(data['Query_division'] == 'PHAGE')
    virusno = sum(data['Query_division'] == 'VIRUS')
    Retrono = sum(data['Retro'] == True)
    Nonretro = virusno - sum(data['Retro'] == True)
    dict = {'Name': Name, 'Total': total, 'Phage': phageno, 'Virus': virusno, 'Retrovirus': Retrono, 'Non-retrovirus': Nonretro}
    return dict

# reads processed file and gets virus type count dict
if os.path.isfile('Process_refseq_head.csv'):
    refseq = pd.read_csv('Process_refseq_head.csv', sep='\t')
    # total number of sequences
    refseqdict = virus_prop(refseq, 'RefSeq')

# merges the genbank nt files and gets virus type count dict
if os.path.isfile('Process_complete_head.csv') and os.path.isfile('Process_partial_head.csv'):
    complete = pd.read_csv('Process_complete_head.csv', sep='\t')
    partial = pd.read_csv('Process_partial_head.csv', sep='\t')
    GenBank = pd.concat([complete, partial], ignore_index=True)
    gendict = virus_prop(GenBank, 'GenBank')

# combines all datasets together
if os.path.isfile('Process_complete_head.csv') and os.path.isfile('Process_partial_head.csv') and os.path.isfile('Process_refseq_head.csv'):

    allseqs = pd.concat([refseq, partial, complete], ignore_index=True)

    # flags data type
    refseqdict["Type"] = "Query FASTA"
    gendict["Type"] = "Query FASTA"
    # calculated virus prop for all seqs in the query fasta
    totaldict = virus_prop(allseqs, 'All')
    totaldict["Type"] = "Query FASTA"
    dictlist = [refseqdict, gendict, totaldict]
    countdf = pd.DataFrame(dictlist)


    # function to plot virus proportion in refseq and genbank database virus sequences
    def plot_virusprop(data):
        data["Phage %"] = round(data["Phage"] / data["Total"] * 100)
        data["Retrovirus %"] = round(data["Retrovirus"] / data["Total"] * 100)
        data["Non-Retrovirus %"] = round(data["Non-retrovirus"] / data["Total"] * 100)

        # plot stacked proportion bar charts for refseq, genbank and total
        fig = go.Figure()
        fig.add_trace(go.Bar(y=data["Phage %"], x= data["Name"], name = "Phage %", text = data["Phage %"],
                             marker=dict(color='#FE4C44',line=dict(color='#FE4C44', width=0.05))))
        fig.add_trace(go.Bar(y=data["Retrovirus %"], x= data["Name"], name = "Retrovirus %", text= data["Retrovirus %"],
                             marker=dict(color='#fFD500',line=dict(color='#fFD500', width=0.05))))
        fig.add_trace(go.Bar(y=data["Non-Retrovirus %"], x= data["Name"], name = "Non-Retrovirus %", text = data["Non-Retrovirus %"],
                             marker=dict(color='#52a4ee',line=dict(color='#52a4ee', width=0.05))))
        # set font size, axis lines, title
        fig.update_layout(yaxis=dict(title_text="%"), width=650, height=640,plot_bgcolor='white', barmode='stack')
        fig.update_layout(uniformtext_minsize=20, uniformtext_mode='hide', font=dict(size=24, color = 'black'))
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black', tickfont=dict(size=21))
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black', tickfont=dict(size=21))
        fig.update_layout(title={'text': "Query FASTA",'y': 0.03 ,'x': 0.35,'xanchor': 'center', 'font' :dict(size=20)})
        # save plot image
        pio.write_image(fig, file='bar.png', format='png', width=650, height=640, scale=5)
        #fig.show(renderer="browser")

    plot_virusprop(countdf)

    # get retrieve total sequences for each catergoy from countdf (dropping virus/phage/retrovirus info)
    dbcounts = countdf[["Name", "Total", "Type"]].copy()
    # drop all
    dbcounts.drop(dbcounts[dbcounts["Name"] == "All"].index, inplace = True)

    # reading in all hits to make comparison with blast results hits
    file = 'AllHits.csv'
    hitdf = pd.read_csv(file, sep="\t")
    # getting virus proportions from blast results hits
    gendict = virus_prop(hitdf[hitdf["Origin"] == "GenBank"], 'GenBank')
    refseqdict = virus_prop(hitdf[hitdf["Origin"] == "RefSeq"], 'RefSeq')
    gendict["Type"] = "Hit"
    refseqdict["Type"] = "Hit"
    dictlist = [refseqdict, gendict]
    hitcountdf = pd.DataFrame(dictlist)
    hitcounts = hitcountdf[["Name", "Total", "Type"]].copy()

    # combining blast hit counts and the refseq+genbank database sequence counts
    seqcounts = pd.concat([dbcounts, hitcounts], ignore_index=True)

    # plot counts for the blast hits and refseq+genbank database sequences
    def bar_plot(countdata):
        plt.figure(figsize=(8,3.5))
        plt.rcParams['font.size'] = 15
        plot = sb.barplot(data= countdata, x='Total', y='Name', hue="Type", palette=["#f2559f", "#0fb045"])
        plot.bar_label(plot.containers[0], fmt='{:,.0f}')
        plot.bar_label(plot.containers[1], fmt='{:,.0f}')
        plot.set_xscale("log")
        plot.spines['right'].set_visible(False)
        plot.spines['top'].set_visible(False)
        plt.tight_layout()
        plt.show()
    bar_plot(seqcounts)

    # finding largest virus families in the datasets
    def famcount(dataa):
        famcount = dataa.groupby(['Query_family']).size().reset_index()
        famcount.rename(columns={0: 'Count'}, inplace=True)
        #print(f"Hits with family not present: {dataa['Query_family'].str.contains('Not present').sum()}")
        famcount = famcount.sort_values('Count', ascending= False)
        return famcount

    def get_famplot(dataa, titlename):
        plt.figure(figsize=(10, 6))
        plot = sb.barplot(data=dataa, y='Query_family', x='Count', palette="crest")
        plot.bar_label(plot.containers[0])
        plt.title(titlename)
        plt.tight_layout()
        plot.spines['right'].set_visible(False)
        plot.spines['top'].set_visible(False)
        fig2 = plot.get_figure()
        fig2.show()

    # find largest families for all sequences in genbank and refseq viral databases
    allfam = famcount(allseqs)
    topallfam = allfam[0:9]
    get_famplot(topallfam, 'Top Families - GenBank and RefSeq Query FASTAs')

    # find largest families in all seqs in refseq viral database
    reffam = famcount(refseq)
    topreffam = reffam[0:6]
    get_famplot(topreffam, 'Top Families - RefSeq FASTA')

    # find largest families in all seqs in genbank viral database
    genfam = famcount(GenBank)
    topgenfam = genfam[0:6]
    get_famplot(topgenfam, 'Top Families - GenBank FASTA')