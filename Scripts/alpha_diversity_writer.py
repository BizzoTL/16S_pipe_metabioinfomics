import pandas as pd
import matplotlib.pyplot as plt
import InputOutput_class as io

"""
# info
utilizzando i file 'rarefaction_observe.csv' 'rarefaction_shannon.csv' 'OUT_16S_Amplicon_Analysis.Xlsx/Legend' e 'group_sample_size.csv'
genera la seconda pagina del file excel, Stat page contenente:
1.una tabella con le colonne di 'quality_info' e 'alpha_diversity'
2.una tabella con la media di ASVs e Shannon per campione
2.un grafico con la media di ASVs per campione

# paths
QualityInfo_path  = '/Users/lucas/16S_amplicon_analysis/phyloseq/quality_info.csv'
AlphaDiversity_path  = '/Users/lucas/16S_amplicon_analysis/phyloseq/alpha_diversity.csv'
Xlsx_output = "/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx"
Sample_size_path  = "/Users/lucas/16S_amplicon_analysis/group_sample_size.csv"
"""

def rarefaction_observe_table(Header:list, DF_output:pd.DataFrame, DF_from_file:pd.DataFrame, legend_xlxs:pd.DataFrame, field:str):
    DF_FirstRow = pd.DataFrame(columns = Header)
    DF_output = pd.concat([DF_output, DF_FirstRow], axis=1)
    pos = 0
    
    
    legend_xlxs_sort = legend_xlxs.sort_values(by=['NGS_ID'], axis=0)
    DF_output_sort = DF_output.reindex(legend_xlxs_sort.index)
    # perché serve? praticamente su R i file vengono ordinati in base all'NGS_ID in ordine *alfabetico* e non in base all'ordine di campionamento
    # cambio l'ordine della leggenda in alfabetico per così poter leggere rarefaction_observe.csv e rarefaction_shannon.csv senza dover riordinarli
    
    for y in range(0, len(DF_output_sort['Samples'])):
        for x in range(1, len(DF_output_sort.columns)):
            if pos == len(DF_from_file['sample']):
                break
            DF_output_sort.loc[y, DF_output_sort.columns[x]] = DF_from_file.loc[pos, field]
            pos += 1
            
    DF_output_sorted = DF_output_sort.sort_index()

    return DF_output_sorted

def average_table(DF_output:pd.DataFrame, DF_table:pd.DataFrame, legend_xlxs:pd.DataFrame, field:str):
        average_row_value = []
        for y in range(0, len(legend_xlxs['Samples'])):
            avg = DF_table.loc[y,:]
            avg = avg.drop('Samples')
            average_row_value.append(round(avg.mean(),2))
        DF_output[field] = average_row_value
        
        return DF_output

def rarefaction_graph(DF_AlDiv_NumbASV:pd.DataFrame, DF_legend:pd.DataFrame, Graph_folder_path:str):
    
    # inizializzazione grafico
    plt.figure(figsize=(10, 16))
    plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.xlabel('Number of Reads')
    plt.ylabel('ASVs')
    plt.title('dada2')
    
    # plot con colori diversi per ogni campione
    num = 0
    n = 0
    marker = ['o','s','^','D','*','P','X','H']
 
    for x in range(0,len(DF_legend['Samples'])):
        if x % 20 == 0 and x != 0:
            n += 1
        col2 =DF_AlDiv_NumbASV.loc[x,:].drop('Samples') # Y AXIS
        
        if DF_legend['Order'][x] == 0:
            dot_color = 'C'+str(num) # colore diverso per ogni campione
            label_name = DF_legend['Samples'][x]
            num += 1
            
        else:
            label_name = None

        #print("DF_AlDiv_NumbASV", DF_AlDiv_NumbASV.loc[x,"Samples"], "label_name", label_name, "marker", marker[n])
        plt.plot(col2, color=dot_color, label=label_name, linewidth=1.5, alpha=0.8, marker=marker[n])
        
    plt.legend(loc='right', bbox_to_anchor=(1.35, 0.5),fancybox=True)
    graph_path  = Graph_folder_path  + 'rarefaction_curve.png'
    plt.savefig(graph_path , dpi=550, bbox_inches='tight')
    plt.clf()
    
    return graph_path
    
    
def main(Xlsx_output:str, Legend_csv_path:str, Graph_folder_path:str, Rarefaction_obsv_path:str, Shannon_entropy_path:str):
    rarefaction_observe = io.ReadFiles.readToDataframe(Rarefaction_obsv_path ,None) # legge il file rarefaction_observe.csv
    shannon_entropy = io.ReadFiles.readToDataframe(Shannon_entropy_path ,None) # legge il file rarefacrion_shannon.csv
    DF_legend = io.ReadFiles.readToDataframe(Legend_csv_path, None) # legge il file "Legend.csv"

    DF_AlDiv_NumbASV = pd.DataFrame()
    DF_AlDiv_NumbASV['Samples'] = DF_legend['Samples']
    DF_AlDiv_Sahannon_Entropy = pd.DataFrame()
    DF_AlDiv_Sahannon_Entropy['Samples'] = DF_legend['Samples']

    Header = []
    for pos in range (0,len(rarefaction_observe['sample'])):
        if rarefaction_observe['sample'][pos] == DF_legend['NGS_ID'][0]:
            Header.append(rarefaction_observe.loc[pos,'readsNums'])

    DF_AlDiv_NumbASV = rarefaction_observe_table(Header, DF_AlDiv_NumbASV, rarefaction_observe, DF_legend, 'value')
    DF_Al_divShannon_Entropy = rarefaction_observe_table(Header, DF_AlDiv_Sahannon_Entropy, shannon_entropy, DF_legend, 'value')
    
    #print(DF_AlDiv_NumbASV.iloc[63:71,:])
    
    AVG_DF_AlDiv_NumbASV = pd.DataFrame()
    AVG_DF_Al_divShannon_Entropy= pd.DataFrame()

    AVG_DF_AlDiv_NumbASV = average_table(AVG_DF_AlDiv_NumbASV, DF_AlDiv_NumbASV, DF_legend, 'Alpha Diversity - Number of ASVs')
    AVG_DF_Al_divShannon_Entropy = average_table(AVG_DF_Al_divShannon_Entropy, DF_Al_divShannon_Entropy, DF_legend, 'Alpha Diversity - Shannon index')
    graph_path = rarefaction_graph(DF_AlDiv_NumbASV, DF_legend, Graph_folder_path)# crea il grafico
    
    # scrittura su file
    DF_AlDiv_NumbASV_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Alpha', DF_AlDiv_NumbASV, 1, 0)
    y_coordinate = len(DF_AlDiv_NumbASV) + 4
    DF_Al_divShannon_Entropy_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Alpha', DF_Al_divShannon_Entropy, y_coordinate, 0)
    x_coordinate = len(DF_AlDiv_NumbASV.columns)+3
    AVG_DF_AlDiv_NumbASV_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Alpha', AVG_DF_AlDiv_NumbASV, 1, x_coordinate)
    AVG_DF_Al_divShannon_Entropy_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Alpha', AVG_DF_Al_divShannon_Entropy, y_coordinate, x_coordinate)
    
    Graph_out = io.WriteFiles.writeGraphs(graph_path , Xlsx_output, 'Alpha', 'BG4')
    
    if DF_AlDiv_NumbASV_out and DF_Al_divShannon_Entropy_out and AVG_DF_AlDiv_NumbASV_out and AVG_DF_Al_divShannon_Entropy_out and Graph_out == True:
        return True
    else:
        return False
