import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap  # Import ListedColormap
from matplotlib.colors import ListedColormap
import InputOutput_class as io

def bar_plot(taxa:pd.DataFrame, DF_bars:pd.DataFrame, Graph_folder_path:str, set_threshold:int):
    colors = np.vstack((plt.cm.tab20(np.linspace(0, 1, 20)),plt.cm.tab20b(np.linspace(0, 1, 20))))
    combined_cmap = ListedColormap(colors)

    bars_transpose = DF_bars.transpose()
    bars_transpose.rename(columns=taxa.loc[:,"ASVs_genus"], inplace=True)
    
    fig = plt.figure()
    fig.set_figheight(15)
    fig.set_figwidth(10)
    # Replace values below the threshold with NaN and sum them up for the "Other" category
    bars_transpose_filtered = bars_transpose.where(bars_transpose >= set_threshold, 0)
    bars_transpose_filtered['Other'] = bars_transpose.where(bars_transpose < set_threshold, 0).sum(axis=1)
    bars_transpose_filtered = bars_transpose_filtered.loc[:, (bars_transpose_filtered > 0).any(axis=0)]
        
    bars_transpose_filtered.plot(kind='bar', stacked=True, colormap=combined_cmap)

    plt.xticks(rotation=90)
    plt.legend(
    loc='center left', 
        bbox_to_anchor=(1.05, 0.5), 
        fancybox=True, 
        ncol=4, 
        labelspacing=1.5, 
        frameon=True,  # Add frame around legend
        framealpha=1,  # Make the frame opaque
        edgecolor='black',  # Frame color
        fontsize='large'  # Increase font       
    )
    plt.xlabel("Sample")
    plt.ylabel("ASV Abundance %")
    plt.title("Releative Abundance >1% ASVs")
    graph_path = Graph_folder_path + 'stacked_bar_plot.png'
    plt.savefig(graph_path, dpi=250, bbox_inches='tight')
    print("Stacked barplot saved in:\t\t", graph_path)                
    plt.clf()
    
    return True
      
def heatmap_plot(taxa:pd.DataFrame, DF_heatmap:pd.DataFrame, Graph_folder_path:str, set_threshold:int):
    fold_change_gradient = LinearSegmentedColormap.from_list('my_gradient', [
            (0.000, (0.259, 1.000, 0.357)),
            (0.450, (0.125, 0.510, 0.149)),
            (0.500, (0.161, 0.161, 0.161)),
            (0.550, (0.569, 0.137, 0.137)),
            (1.000, (0.980, 0.322, 0.322))
    ])# colors from green to black to red
    sample_list = []
    sample_list = DF_heatmap.columns.to_list()
    
    different_samples = {}
    """
    Con lo scopo di rappresentare i valori graficamente come visto nella relazione del prof
    il programma è stato modificato per rappresentare i dati in modo simile
    il for successivo si occupa di identificare i diversi campioni presenti nel dataset
    tra forsu e biob e di contare quanti ce ne sono, servirà per il fold change
    """
    for sample in sample_list:
        if "_" in sample:
            name = sample.split("_", 1)[0]
        else:
            name = sample
        
        if name in different_samples:
            different_samples[name] += 1
        else:
            different_samples[name] = 1
            
    fig = plt.figure(figsize=(10, 14))
    
    ch_graph_path = Graph_folder_path + 'comparative_heatmap_plot.png' #fold change directory
    h_graph_path = Graph_folder_path + 'heatmap_plot.png' #heatmap (absolute value) directory
    
    DF_heatmap_bool_threshold = (DF_heatmap > set_threshold).any(axis=1)
    DF_heatmap_threshold = DF_heatmap.loc[DF_heatmap_bool_threshold]
    taxa_threshold = taxa.loc[DF_heatmap_bool_threshold]
    
    DF_heatmap_threshold = DF_heatmap_threshold.copy()
    DF_heatmap_threshold.rename(index= taxa_threshold.loc[:,"ASVs_genus"], inplace=True)
    
    heatmap = pd.DataFrame(DF_heatmap_threshold).to_numpy()
    DF_heatmap_log2fc = np.zeros((len(DF_heatmap_threshold.index), len(DF_heatmap_threshold.columns)))
    DF_heatmap_log2fc= np.log2((heatmap+ 1) / (heatmap[:,0][:, np.newaxis]+ 1))

    #per BIOB e FORSU
    """
    Questa porzione di codice fa fede alla rappresentazione dei dati nella heatmap della relazione del prof
    di conseguenza è stata modificata per dare la stessa rappresentazione e quindi in condizioni diverse
    potrebbe rompersi
    """ 
    """min = 0
    value = 0
    num = 0
    if len(different_samples) > 1:
        value = 1
    for key in different_samples:
        value = value + different_samples[key] + value
        if value > len(DF_heatmap_threshold.columns): # serve per analizzare prima biob poi forsu nel vs
            value = len(DF_heatmap_threshold.columns)
            num = -1
        DF_heatmap_log2fc[:,min:value] = np.log2((heatmap[:, min:value] + 1) / (heatmap[:, min+num][:, np.newaxis]+ 1))
        min = min + different_samples[key] + value"""

    fold_change= pd.DataFrame(DF_heatmap_log2fc, index=DF_heatmap_threshold.index, columns=DF_heatmap_threshold.columns)
    fold_change.drop(fold_change.columns[0], axis=1, inplace=True) #rimuove la prima colonna
    #fold_change.drop(fold_change.columns[0], axis=1, inplace=True)
    
    sns.heatmap(fold_change,  yticklabels=True, cmap=fold_change_gradient, vmin=-2, vmax=2) # heatmap per il fold change
    plt.xticks(rotation=90)
    plt.xlabel("Sample", fontsize=20)
    plt.ylabel("ASVs", fontsize=20)
    plt.title("Fold Change", fontsize=22)
    plt.savefig(ch_graph_path, dpi=250, bbox_inches='tight')
    print("Comparative heatmap plot saved in:\t", ch_graph_path)          
    plt.clf()
    
    
    sns.heatmap(DF_heatmap_threshold, cmap='gnuplot', vmin=0, vmax=20) # heatmap per l'abbondanza assoluta
    plt.xticks(rotation=90)
    plt.xlabel("Sample", fontsize=20)
    plt.ylabel("ASVs", fontsize=20)
    plt.title("Relative Abundance >0.5%", fontsize=22)
    plt.savefig(h_graph_path, dpi=250, bbox_inches='tight')
    print("Relative Abundance heatmap plot saved in: ", ch_graph_path)          
    plt.clf()

    return True

def ASVs_table(Graph_folder_path:str, DF_all_sample_abundance:pd.DataFrame, DF_legend:pd.DataFrame, set_threshold:int):
    
    
    ASVs = DF_all_sample_abundance.iloc[:,0:7]
    RA_aboundance = DF_all_sample_abundance.loc[:,DF_legend.loc[:,"Samples"]]
    
    sample_list = []
    sample_list = DF_legend.loc[:,"Samples"].to_list()
    
    different_samples = {}

    for sample in sample_list:
        if "_" in sample:
            name = sample.split("_", 1)[0]
        else:
            name = sample
        
        if name in different_samples:
            different_samples[name] += 1
        else:
            different_samples[name] = 1
            
    RA_bool_threshold = (RA_aboundance > set_threshold).any(axis=1)
    RA_aboundance_05 = DF_all_sample_abundance.loc[RA_bool_threshold]
    ASVs = ASVs.loc[RA_bool_threshold]
    RA_samples = RA_aboundance_05.loc[:,DF_legend.loc[:,"Samples"]]
    min = 0
    value = 0
    for key in different_samples:
        value = value + different_samples[key]
        ASVs.loc[:, key] = RA_samples.iloc[:,min:value].mean(axis=1).round(2)
        min  = min + different_samples[key]
    
    
    fig, ax = plt.subplots(figsize=(18, 20))  # Adjust figsize as needed

    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False)
    ax.set_frame_on(False)

    table = ax.table(cellText=ASVs.values, colLabels=ASVs.columns, cellLoc='center', loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(28)
    table.scale(2, 4)
    colors = ['#f0f0f0', '#d9d9d9']
      
    for i in range(len(ASVs)):
        color = colors[i % len(colors)]
        for j in range(len(ASVs.columns)):
            table[(i + 1, j)].set_facecolor(color)
    
    plt.savefig(Graph_folder_path + 'ASVs.png', bbox_inches='tight', dpi=330)
    plt.close()
    print("ASVs table saved in:\t\t\t", Graph_folder_path + 'ASVs.png')
    ASVs = ASVs.fillna("NA")
    Phylum_grouped = ASVs.groupby("Phylum").sum()
    start_column = len(Phylum_grouped.columns) - len(different_samples)
    
    all_average_phy = Phylum_grouped.iloc[:, start_column:].mean(axis=1).round(2)
    others = 100 - all_average_phy.sum().round(2)
    all_average_phy = pd.concat([all_average_phy, pd.Series(others, index=['Others'])], ignore_index=False)    

    colors = sns.color_palette('pastel')[0:20]
   
    patches, texts = plt.pie(all_average_phy, colors=colors, startangle=90, radius=1.2)
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(all_average_phy.index, all_average_phy)]
    sort_legend = True
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, all_average_phy),
                                            key=lambda x: x[2],
                                            reverse=True))

    plt.legend(patches, labels, loc='lower left',  bbox_to_anchor=(0, -0.7), ncol =1, fontsize=16)
    plt.title("Phylum Abundance", fontsize=22)  
    plt.savefig(Graph_folder_path + 'Phylum_pie_chart.png', bbox_inches='tight', dpi=300)
    
    ############################################################
    
    Genera_grouped = ASVs.groupby("Genus").sum()
    start_column = len(Genera_grouped.columns) - len(different_samples)
    
    all_average_ge = Genera_grouped.iloc[:, start_column:].mean(axis=1).round(2)
    others = 100 - all_average_ge.sum().round(2)
    all_average_ge = pd.concat([all_average_ge, pd.Series(others, index=['Others'])], ignore_index=False)    

    
    cmap = sns.color_palette("Spectral", as_cmap=True)
    colors = [cmap(i) for i in np.linspace(0, 1, len(all_average_ge))]

   
    patches, texts = plt.pie(all_average_ge, colors=colors, startangle=90, radius=1.2)
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(all_average_ge.index, all_average_ge)]
    sort_legend = True
    if sort_legend:
        patches, labels, dummy=  zip(*sorted(zip(patches, labels, all_average_ge),
                                            key=lambda x: x[2],
                                            reverse=True))

    plt.legend(patches, labels, loc='center',  bbox_to_anchor=(0.5, -0.6), ncol = 2, fontsize=16)
    plt.title("Genus Abundance", fontsize=22)
    plt.savefig(Graph_folder_path + 'Genus_pie_chart.png', bbox_inches='tight', dpi=300)

    return True

def main(Graph_folder_path:str, Legend_csv_path:str, all_sample_aboundance_path:str):
    
    set_threshold = 1 #SET VALUE THRESHOLD

    taxa_id = ["ASVs_genus","Domain","Phylum","Class","Order","Family","Genus"]
    DF_all_sample_abundance = io.ReadFiles.readToDataframe(all_sample_aboundance_path, None)
    DF_legend = io.ReadFiles.readToDataframe(Legend_csv_path, None)

    sample_index = []
    taxa = pd.DataFrame(DF_all_sample_abundance.loc[:,taxa_id])
    for pos in range(0, len(DF_legend.loc[:, "Samples"])):
        if DF_legend.loc[pos, "Order"] == 0:
            sample_index.append(DF_legend.loc[pos, "Samples"])

    DF_ASV = pd.DataFrame(DF_all_sample_abundance.loc[:,sample_index])
    
    Graph_out = bar_plot(taxa,DF_ASV,Graph_folder_path, set_threshold)
    Heatmap_out = heatmap_plot(taxa,DF_ASV,Graph_folder_path, set_threshold)
    ASVs_out = ASVs_table(Graph_folder_path, DF_all_sample_abundance, DF_legend, set_threshold)

    if Heatmap_out and Graph_out and ASVs_out == True:
        return True
    else:
        return False
