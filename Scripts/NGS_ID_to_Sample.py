import pandas as pd
import  InputOutput_class as io

def main(Xlsx_output, Legend_csv_path, Xlsx_input):

    df = io.ReadFiles.readToDataframe(Xlsx_input, 'Foglio1')
    Legend_Header = ["Position", "NGS_ID", "Samples", "Replicates", "Order"] #Ordine tabella output
    DF_Legend = pd.DataFrame(columns=Legend_Header)


    pos = 0
    order = 0
    for sample in df["Nome"]:
        sample_id = sample[:-2]
        DF_Legend.loc[pos, "Position"] = pos+1
    
        NGS_ID = df.loc[pos, "Nome file fastq (forward reads)"]
        sep = '_'
        DF_Legend.loc[pos, "NGS_ID"] = NGS_ID.split(sep, 1)[0] #Splitto il nome del file per ottenere solo l'ID (rimuovo fastq)
        
        DF_Legend.loc[pos, "Samples"] = df.loc[pos, "Nome"][:-2]
        DF_Legend.loc[pos, "Replicates"] = df.loc[pos, "Replica"]
        if pos != 0 and DF_Legend.loc[pos-1, "Samples"] == sample_id: #Pos != 0 per evitare errore di out of range
            order += 1
        else:
            order = 0
        DF_Legend.loc[pos, "Order"] = order
        pos += 1
            
    legend_xlsx_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Legend', DF_Legend, 0, 0) # scrive la leggenda in Legend.xlsx
    legend_csv_out = io.WriteFiles.writeToCsv(Legend_csv_path, DF_Legend) # scrive la leggenda in Legend.csv

    if legend_xlsx_out and legend_csv_out == True: # verifico se i file sono stati scritti correttamente
        return True
    else:
        return False
