import pandas as pd
import os.path as os
import openpyxl as op


class ReadFiles:
    
    def readToDataframe(path:str, name:str): # apre un file e restituisce un dataframe
        if path.endswith('.csv'): # 1 caso: file csv
            file = pd.read_csv(path, sep=',', header=0)
            DF_file = pd.DataFrame(file)
            if 'Unnamed: 0' in DF_file.columns: # rinomina la colonna 'Unnamed: 0' in 'NGS_ID' per i file 'quality_info.csv' e 'alpha_diversity.csv'
                DF_file.rename(columns={'Unnamed: 0': 'NGS_ID'}, inplace=True)
        if path.endswith('.xlsx'): # 2 caso: file excel
            file = pd.read_excel(path, sheet_name=name, header=0)
            DF_file = pd.DataFrame(file)
        return DF_file
    

class WriteFiles:
    
    def writeToCsv(path:str, tb:pd.DataFrame):
        tb.to_csv(path, index=False)
        print("CSV file saved in:\t\t\t", path)
        return True
    
    def writeToXlsx(path:str, sheet:str, tb:pd.DataFrame, y_coordinate:int, x_coordinate:int):
        if os.isfile(path):
            # se il file esiste, aggiungi la pagina 'append mode' e 'overlay' per aggiungere elementi
            with pd.ExcelWriter(path, engine='openpyxl', mode='a', if_sheet_exists='overlay' ) as writer:
                tb.to_excel(writer, sheet_name=sheet, index=False, startrow=y_coordinate, startcol=x_coordinate)
                print(sheet, "element saved in:\t\t\t", path)
                return True
        else:
            # se il file non esiste, crea il file 'write mode'
            with pd.ExcelWriter(path, engine='openpyxl', mode='w') as writer:
                tb.to_excel(writer, sheet_name=sheet, index=False, startrow=y_coordinate, startcol=x_coordinate)
                # coordinate x y calcolate dal programma che chiama la funzione
                print(sheet, "page created in:\t\t\t", path)
                return True            
    
    def writeGraphs(graph_path:str, path:str, name:str, cell_pos:str):
        wb = op.load_workbook(path)
        ws = wb[name]
        img = op.drawing.image.Image(graph_path)
        img.anchor = cell_pos # valore arbitrario dato dal programma che chiama la funzione
        img.width = 1000
        img.height = 600
        ws.add_image(img)
        wb.save(path)
        print("Graph saved in",name,"page at:\t\t", path)
        return True
