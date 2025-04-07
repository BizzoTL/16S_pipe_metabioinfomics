import pandas as pd
import matplotlib.pyplot as plt
import InputOutput_class as io

"""
# info
utilizzando i file 'beta_diversity_values.csv', 'beta_diversity_vectors.csv', 'Legend.csv' e 'group_sample_size.csv'
genera la terza pagina del file excel, Beta page contenente:
1.  una tabella (VECTOR_table) composta alla prima riga dei valori di beta_diversity_values.csv (in particolare la colonna 'relative_eig') 
    convertiti in percentuale e rinominati in 'Bray-Curtis', alla seconda riga i nomi delle colonne di beta_diversity_vectors.csv 
2.  una tabella (PCoA_2D_table) con le prime 2 colonne di VECTOR_table per essere utilizzata in un grafico
3.  un grafico composto dai valori numeri di PCoA_2D_table, colorati in base al Sample specificato tramite 'Legend.csv'

#input directories
    Group_size_path = '/Users/lucas/16S_amplicon_analysis/group_sample_size.csv'
    Beta_VECTOR_path= '/Users/lucas/16S_amplicon_analysis/phyloseq/beta_diversity_vectors.csv'
    BetaDiversity_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/beta_diversity_values.csv'
#output directories
    Xlsx_path = '/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx'
    Graph_folder_path = '/Users/lucas/16S_amplicon_analysis/graphs'
"""
def graphScatter(Graph_folder_path:str ,PCoA_2D_table:pd.DataFrame, DF_Legend:pd.DataFrame):
    # inizializzo il grafico
    plt.figure(figsize=(15, 15))
    plt.axhline(0,color='black', linewidth=0.8)
    plt.axvline(0,color='black', linewidth=0.8)
    plt.grid(color='gray', axis='y', linewidth=0.4)
    plt.xlabel('PCo1', weight='bold', fontsize=25)
    plt.ylabel('PCo2', weight='bold', fontsize=25)
    plt.title('PCoA 2D',weight='bold', fontsize=30)
    offset = 0.003
    # creo lo scatter plot
    num = 0
    m = 0
    sep = '_'
    cmap = plt.get_cmap('tab20c')
    marker = ['^','D','o','x','s']
    
    PCoA_2D_table.reset_index(drop=True, inplace=True)
    for pos in range(0,len(PCoA_2D_table)):

        if DF_Legend.loc[pos,'Order'] == 0: # -1 per allineare l'indice di DF_Legend con quello di PCoA_2D_table
            color = cmap(num%20)
            num += 1
            if num%20==0:
                m += 1  # cambio il marker ogni 20 campioni
            
            for val in range(pos+1, len(DF_Legend['Samples'])):
                if DF_Legend.loc[val,'Order'] ==  0:  
                    break  
            
            PCoA_2D_table_X = PCoA_2D_table.iloc[pos:val, 1].mean()# calcola la media dei campioni
            PCoA_2D_table_Y = PCoA_2D_table.iloc[pos:val, 2].mean()# calcola la media dei campioni
            
            legend_name = DF_Legend['Samples'][pos] # nomi dei campioni
            
            if set(DF_Legend['Samples'][pos]) == {"_"}: # se il nome del campione contiene un '_' lo splitta e prende solo la seconda parte es BIOB_22/06/2020 -> 22/06/2020
                plot_name = DF_Legend['Samples'][pos].split(sep, 1)[1]
            else:
                plot_name = DF_Legend['Samples'][pos]
                
            plt.scatter(PCoA_2D_table_X, PCoA_2D_table_Y, s=300, alpha=0.8, color=color, label=legend_name, marker=marker[m])
            plt.text(PCoA_2D_table_X + offset, PCoA_2D_table_Y + offset, plot_name, fontsize=9,  bbox = dict(facecolor = 'gray', alpha = 0.1, edgecolor='black', boxstyle='round,pad=0.1'))
            
    
    plt.legend(
    loc='center left', 
    bbox_to_anchor=(1.05, 0.5), 
    fancybox=True, 
    ncol=1, 
    labelspacing=1.5, 
    prop={'weight': 'bold'},  # Make legend text bold
    frameon=True,  # Add frame around legend
    framealpha=1,  # Make the frame opaque
    edgecolor='black',  # Frame color
    fontsize='large'  # Increase font size
    )
    graph_path = Graph_folder_path + 'PCoA_2D.png'
    plt.savefig(graph_path, dpi=250, bbox_inches='tight')
    print("PCoA 2D plot saved in:\t\t\t", graph_path)
    plt.clf()
    return graph_path

def main(Xlsx_output:str, DF_legend_path:str, Graph_folder_path:str, BetaDiversity_VALUES_path:str, BetaDiversity_VECTOR_path:str):
    DF_VALUES_beta_diversity = io.ReadFiles.readToDataframe(BetaDiversity_VALUES_path, None)
    DF_VECTOR_beta_diversity = io.ReadFiles.readToDataframe(BetaDiversity_VECTOR_path, None)
    DF_Legend = io.ReadFiles.readToDataframe(DF_legend_path, None)

    """
    Genera la prima colonna della tabella contenente i valori di beta_diversity_vectors.csv
    questa prima colonna è composta dai valori di beta_diversity_values.csv
    convertiti in percentuale e sottoposti a trimming (rimozione dei valori negativi)
    successivamente rinominata in 'Bray-Curtis'
    """

    relative_eig = pd.DataFrame(DF_VALUES_beta_diversity['Relative_eig']) # seleziona la colonna 'Relative_eig' di beta_diversity_values.csv
    relative_eig = round(relative_eig*100,2)
    relative_eig = relative_eig[relative_eig['Relative_eig'] > 0] # rimuove i valori negativi 'trimming' dei valori
    relative_eig.rename(columns={'Relative_eig':'0' }, inplace=True)
    relative_eig = relative_eig.transpose()
    relative_eig.insert(0, '0', 'Bray-Curtis')
    relative_eig.columns = range(0,len(relative_eig.columns))
    # relative_eig rappresenta la prima riga futura VECTOR_table, contenenti le % di beta_diversity_values.csv

    DF_VECTOR_beta_diversity.insert(1, 'Samples', DF_Legend['Samples'])  # aggiunta colonna 'Samples' a DF_VECTOR_beta_diversity
    for x in range(2, len(DF_VECTOR_beta_diversity.columns)):
        DF_VECTOR_beta_diversity.rename(columns={DF_VECTOR_beta_diversity.columns[x]: 'PCo' + str(x-1)}, inplace=True)

    DF_VECTOR_beta_diversity =DF_VECTOR_beta_diversity.drop(DF_VECTOR_beta_diversity.columns[0], axis=1) 
    #rimuovo la colonna degli NGS ID siccome ho aggiunto la colonna 'Samples'

    DF_header = pd.DataFrame(data=DF_VECTOR_beta_diversity.columns)
    DF_header = DF_header.transpose()
    VECTOR_table = pd.concat([relative_eig, DF_header ], axis=0) # creo la tabella VECTOR_table
    #prendo gli indici di colonna di DF_VECTOR_beta_diversity e li aggiungo come seconda riga della tabella VECTOR_table

    beta_diversity = pd.DataFrame(data=DF_VECTOR_beta_diversity)
    beta_diversity.columns = range(len(beta_diversity.columns)) # rinomino le colonne di beta_diversity, che derivano da DF_VECTOR_beta_diversity in numeri
    VECTOR_table = pd.concat([VECTOR_table, beta_diversity], axis=0)
    #VECTOR_table.index = range(0,len(VECTOR_table)) # rinomino gli indici(righe) di VECTOR_table
    VECTOR_table = VECTOR_table.reset_index(drop=True)

    PCoA_2D_table = pd.DataFrame(VECTOR_table.iloc[:,0:3]) # seleziono le prime 3 colonne di VECTOR_table e le salvo in DF_PCoA_2D
    graph_path = graphScatter(Graph_folder_path, PCoA_2D_table.iloc[2:len(PCoA_2D_table.index),:], DF_Legend) # creo il grafico PCoA_2D
    VECTOR_table_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Beta', VECTOR_table, 1, 1)
    y_coordinate = len(VECTOR_table)+3
    PCoA_2D_table_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Beta', PCoA_2D_table, y_coordinate, 1)
    Graph_out = io.WriteFiles.writeGraphs(graph_path, Xlsx_output, 'Beta', 'V2')
    
    if VECTOR_table_out and PCoA_2D_table_out and Graph_out == True:
        return True
    else:
        return False
