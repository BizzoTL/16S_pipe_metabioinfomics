import pandas as pd
import InputOutput_class as io

"""
#info
il programma modifica 'Stat' page riempendo le colonne mancanti come 'assigned to ASVs', 'perc of assigned to ASVs','ASVs >0.01%', 'ASVs >0.5%'
sulla base dei dati contenuti in 'Raw_Data' sotto le colonne di 'Abundance of NGS_ID'

#paths
    Xlsx_output = '/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx'
    TotalCSV_path = '/Users/lucas/16S_amplicon_analysis/total_page.csv'
    Legend_path = '/Users/lucas/16S_amplicon_analysis/legend.csv'
"""

def main(Xlsx_output:str, Legend_path:str, Filtered_rows_path:str, Total_rows_path:str, input_folder_path:str):

    raw_data_path = input_folder_path + 'Raw_Data.csv'
    stat_path = input_folder_path + 'Stat.csv'
    raw_data_page = io.ReadFiles.readToDataframe(raw_data_path, None)
    stat_page = io.ReadFiles.readToDataframe(stat_path, None)
    legend = io.ReadFiles.readToDataframe(Legend_path, None)

    # somma totale dei Raw Reads (escludendo l'ultimo valore che è la media)
    NGSID_ASVs = []

    for pos in range(0, len(legend)): #se hai un'errore qua, controlla raw_data_page_writer! linea 14 del codice e cambia tra 7-8
        column = 'Abundance of ' + legend.loc[pos, 'NGS_ID']
        NGSID_ASVs.append(raw_data_page[column].sum())
        
    total_ASVs = sum(NGSID_ASVs)

    for pos in range(0, len(stat_page)):
        if pd.isna(stat_page.loc[pos, 'Samples']) == False:
            stat_page.loc[pos, 'assigned to ASVs'] = NGSID_ASVs[pos]
            stat_page.loc[pos, 'perc of assigned to ASVs'] = round(NGSID_ASVs[pos]/total_ASVs, 3)
            
            if stat_page.loc[pos, 'perc of assigned to ASVs'] > 0.5:
                stat_page.loc[pos, 'ASVs >0.5%'] = 1
                stat_page.loc[pos, 'ASVs >0.01%'] = 1
            else:
                stat_page.loc[pos, 'ASVs >0.5%'] = 0
                if stat_page.loc[pos, 'perc of assigned to ASVs'] > 0.01:
                    stat_page.loc[pos, 'ASVs >0.01%'] = 1
                else:
                    stat_page.loc[pos, 'ASVs >0.01%'] = 0


    Header = ('ASVs_genus', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'ASVs')
    total_page = pd.DataFrame(columns=Header)
    raw_data_page = raw_data_page.sort_values(by = 'Combined abundance', ascending = False)
    total_page['ASVs_genus'] = raw_data_page['ID']
    total_page['Domain'] = raw_data_page['Kingdom']
    total_page['Phylum'] = raw_data_page['Phylum']
    total_page['Class'] = raw_data_page['Class']
    total_page['Order'] = raw_data_page['Order']
    total_page['Family'] = raw_data_page['Family']
    total_page['Genus'] = raw_data_page['Genus']
    total_page['ASVs'] = raw_data_page['ID']

    for pos in range(0, len(legend)):
        column = 'Abundance of ' + legend.loc[pos, 'NGS_ID']
        total_page[column] = raw_data_page[column]
        total_page.rename(columns={column: legend.loc[pos, 'NGS_ID']}, inplace=True)


    DF_ASV_aboundance = round(total_page[legend.loc[:, 'NGS_ID']]/NGSID_ASVs[:]*100, 2)
    for pos in range(0, len(DF_ASV_aboundance.columns)):
        DF_ASV_aboundance.rename(columns={DF_ASV_aboundance.columns[pos]: DF_ASV_aboundance.columns[pos] + '%'}, inplace=True)
    
    total_page = pd.concat([total_page, DF_ASV_aboundance], axis=1)
    all_samples =pd.concat([total_page.iloc[:,0:7], DF_ASV_aboundance], axis=1)
    """
    for pos in range(0, len(legend)): # calcolo abbondanza relativa delle specie (in ASVs) in ogni campione (NGS_ID) e la salvo in una nuova colonna
        total_page[legend.loc[pos, 'NGS_ID'] + '%'] = round(total_page[legend.loc[pos, 'NGS_ID']]/NGSID_ASVs[pos]*100, 2) 
    """
    
    DF_percentages001 = pd.DataFrame(total_page.iloc[:,len(Header)+len(legend):])
    DF_percentages005 = pd.DataFrame(total_page.iloc[:,len(Header)+len(legend):])
    total_page['AVG'] = round(total_page.iloc[:,len(Header)+len(legend):].mean(axis=1), 2)  # calcolo la media delle abbondanze relative (header + il primo blocco di di NGS_ID)

    DF_percentages001 = DF_percentages001> 0.01 # convert bool table into 0/1 values
    DF_percentages001 = DF_percentages001.astype(int)
    DF_percentages001['>0.01%'] = DF_percentages001.sum(axis=1)

    DF_percentages005 = DF_percentages005> 0.5 # convert bool table into 0/1 values
    DF_percentages005 = DF_percentages005.astype(int)
    DF_percentages005['>0.5%'] = DF_percentages005.sum(axis=1)


    total_page = pd.concat([total_page, DF_percentages001], axis=1)
    total_page = pd.concat([total_page, DF_percentages005], axis=1)

    filtered_rows = total_page[total_page['>0.5%'] > 0]
    filtered_rows = filtered_rows.iloc[:,0:len(Header)+len(legend)]
    
    io.WriteFiles.writeToCsv(Total_rows_path, all_samples)
    io.WriteFiles.writeToCsv(Filtered_rows_path, filtered_rows) # scrive la tabella 'total_page' nel file csv, utile per il la scrittura del file successivo
    total_out = io.WriteFiles.writeToXlsx(Xlsx_output,'Total', total_page, 0, 0) # scrive la tabella 'Total' nel file excel
    
    if total_out == True:
        return True
    else:
        return False
