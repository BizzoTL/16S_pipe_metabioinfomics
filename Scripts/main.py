import time
import NGS_ID_to_Sample
import stat_page_writer
import alpha_diversity_writer
import beta_diversity_writer
import raw_data_page_writer
import total_page_writer
import total_order_refined_page_writer
import data_visualization

"""
*librerie usate*
pandas
numpy
matplotlib.pyplot
openpyxl
InputOutput_class (scritta da me)
os.path
"""
# prima di avviarlo una seconda volta, cancellare il file di output per un runtime ridotto
start_time = time.time()

input_folder_path = '/home/baccus/Desktop/16S_Zebra/dada2_results/'  # folder dove sono posizionati i file di input
output_folder_path = '/home/baccus/Desktop/16S_Zebra/result_graphs/'  # folder dove verrà posizionato il file di output
Graph_folder_path = output_folder_path + 'graphs/' #creare cartella nested in output_folder_path

#  output directories
OUT_16SAmpliconAnalysis_path = output_folder_path + '16S_Amplicon_Analysis.xlsx'
Legend_csv_path = output_folder_path + 'Legend.csv'
Filtered_rows_path = output_folder_path + 'filtered_rows_total_page.csv'
Total_rows_path = output_folder_path + 'total_rows.csv'
Sample_aboundance_path = output_folder_path + 'sample_aboundance_table.csv'
all_sample_aboundance_path = output_folder_path + 'all_sample_aboundance.csv'
# input directories
Xlsx_input = input_folder_path + 'Campioni.xlsx' 

# preso dal folder input
QualityInfo_path = input_folder_path + 'quality_info.csv'
AlphaDiversity_path = input_folder_path + 'alpha_diversity.csv'
BetaDiversity_VALUES_path = input_folder_path + 'beta_diversity_values.csv'
BetaDiversity_VECTOR_path = input_folder_path + 'beta_diversity_vectors.csv'
Raw_data_path = input_folder_path + 'final_table.csv'
Rarefaction_obsv_path = input_folder_path + 'rarefaction_observe.csv'
Shannon_entropy_path = input_folder_path + 'rarefaction_shannon.csv'


print("main.py started\n")

print("NGS_ID_to_Sample.py started")
NGS_IDtoSample_message = NGS_ID_to_Sample.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Xlsx_input)
print(NGS_IDtoSample_message,'\n')

print("stat_page_writer.py started")
Stat_page_message = stat_page_writer.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Graph_folder_path, QualityInfo_path, AlphaDiversity_path, output_folder_path)
print(Stat_page_message,'\n')

print("alpha_diversity_writer.py started")
Alpha_page_message = alpha_diversity_writer.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Graph_folder_path, Rarefaction_obsv_path, Shannon_entropy_path)
print(Alpha_page_message,'\n')

print("beta_diversity_writer.py started")
Beta_page_message = beta_diversity_writer.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Graph_folder_path, BetaDiversity_VALUES_path, BetaDiversity_VECTOR_path)
print(Beta_page_message,'\n')

print("raw_data_page_writer.py started")
Raw_data_page_message = raw_data_page_writer.main(OUT_16SAmpliconAnalysis_path, Raw_data_path, output_folder_path)
print(Raw_data_page_message,'\n')

print("total_page_writer.py started")
Total_page_message = total_page_writer.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Filtered_rows_path, Total_rows_path, output_folder_path)
print(Total_page_message,'\n')

print("total_order_page_writer.py started")
Total_order_page_message = total_order_refined_page_writer.main(OUT_16SAmpliconAnalysis_path, Legend_csv_path, Filtered_rows_path, Total_rows_path, Sample_aboundance_path, all_sample_aboundance_path)
print(Total_order_page_message,'\n')

print("data_visualization.py started")
Data_visualization_message = data_visualization.main(Graph_folder_path, Legend_csv_path, all_sample_aboundance_path)
print(Data_visualization_message,'\n')

print("main.py ended\nruntime:", (time.time() - start_time), "seconds")

"""
Questo è il main del programma, che si occupa di chiamare le funzioni principali 

# l'ordine di esecuzione delle funzioni è importante
#NGS_ID_to_Sample.main()  
    genera la prima pagina del file excel, 'Legend', con la leggenda NGS_ID to Sample, 
    che serve essere letto dai successivi script, minore ridondanza di chiamate
    dello script avendo il file scritto a cui i file successivi possono fare riferimento
#stat_page_writer.main()
    genera la seconda pagina del file excel, 'Stat',
    contenente i dati di quality_info.csv e alpha_diversity.csv
#alpha_diversity_writer.main()
    genera la terza pagina del file excel, 'Alpha', 
    contenente i dati di rarefaction_observe.csv e rarefaction_shannon.csv
    e il grafico 'dada2.png'
#beta_diversity_writer.main()
    genera la quarta pagina del file excel, 'Beta', 
    contenente i dati di beta_diversity_values.csv e beta_diversity_vectors.csv
    e il grafico 'PCoA_2D.png'
#raw_data_page_writer.main()
    genera la quinta pagina del file excel, 'Raw_Data', 
    contenente i dati di final_table.csv
#total_page_writer.main()
    genera la sesta pagina del file excel, 'Total' 
    genera il file 'filtered_rows_total_page.csv' necessario per total_order_refined_page_writer.main()
#total_order_refined_page_writer.main()
    genera la settima pagina del file excel, 'Total_Order', e l'ottava pagina 'Refined'

#total_order_page_writer.main()

#/Users/lucas/16S_amplicon_analysis/pipeline_code/NGS_ID_to_Sample.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/stat_page_writer.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/alpha_diversity_writer.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/beta_diversity_writer.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/raw_data_page_writer.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/total_page_writer.py
#/Users/lucas/16S_amplicon_analysis/pipeline_code/total_order_refined_page_writer.py
"""

#Cappelletti Lucas, lucas.cappelletti@studenti.unipd.it lucas.cappelletti02@gmail.com
