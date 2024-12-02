##############################################################################################################
#Import Dependency 
##############################################################################################################
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
import io
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde
from shiny import App, ui, render, reactive
from shinywidgets import output_widget
from shinywidgets import output_widget, render_widget
from google.cloud import storage
from scipy.stats import pearsonr  # Import pearsonr for calculating Pearson correlation

#Taigapy can not be used for Shiny Deployment, test only 
from taigapy import default_tc3 as tc
from taigapy.client_v3 import UploadedFile, LocalFormat, TaigaReference




##############################################################################################################
#Datasets Used 
##############################################################################################################

##############################################################################################################

ParalogScreenQCReport = tc.get(name='paralogscreen08232024-3db1',file='ParalogScreenQCReport')
ParalogCommonEssentialControlsFromAVANA = tc.get(name='paralogv2-library-files-d8b3', version=20, file='ParalogCommonEssentialControlsFromAVANA')
ParalogNonessentialControlsFromAVANA = tc.get(name='paralogv2-library-files-d8b3', version=20, file='ParalogNonessentialControlsFromAVANA')
ParalogFullLfcGeneSeq = tc.get(name='paralogscreen08232024-3db1', version=11, file='ParalogFullLfcGeneSeq')

common_ess = ParalogCommonEssentialControlsFromAVANA["Symbol"]
none_ess = ParalogNonessentialControlsFromAVANA["Symbol"]
paralog_list = ['OUMS23', 'SNU308']

AvanaLogfoldChange = tc.get(name='paralogscreen08232024-3db1', version=11, file='AvanaLogfoldChange')
print("AvanaLogfoldChange: ", AvanaLogfoldChange)
AvanaLogfoldChange.set_index("Gene", inplace=True)
print("AvanaLogfoldChange2: ", AvanaLogfoldChange )
AvanaLogfoldChange_GPP = AvanaLogfoldChange
print("AvanaLogfoldChange_GPP: ", AvanaLogfoldChange_GPP)
ParalogFullLfcGeneScreen = tc.get(name='paralogscreen08232024-3db1', version=11, file='ParalogFullLfcGeneScreen')

ParalogFullLfcGeneScreen.set_index("GuideTargetSymbol", inplace=True)

ParalogFullLfcGeneScreen_GPP = ParalogFullLfcGeneScreen[paralog_list]
print("ParalogFullLfcGeneScreen_GPP: ", ParalogFullLfcGeneScreen_GPP)

common_indices = ParalogFullLfcGeneScreen_GPP.index.intersection(AvanaLogfoldChange_GPP.index)
print("length of common indices", len(common_indices))
ParalogFullLfcGeneScreen_GPP = ParalogFullLfcGeneScreen_GPP.loc[common_indices]
AvanaLogfoldChange_GPP = AvanaLogfoldChange_GPP.loc[common_indices]
print("AvanaLogfoldChange_GPP : ", AvanaLogfoldChange_GPP)


##############################################################################################################
# # Load static data
# ParalogFullGuideMap = pd.read_csv('Data/paralogscreen06272024_v9-paralogfullguidemap.csv')
# ParalogFullLfcGeneSeq = pd.read_csv('Data/paralogfulllfcgeneseq.csv')
# ParalogSequenceMap = pd.read_csv('Data/paralogscreen06272024_v9-paralogsequencemap.csv')
# ParalogSequenceQCReport = pd.read_csv('Data/paralogscreen06272024_v9-paralogsequenceqcreport.csv')
# ParalogCommonEssentialControlsFromAVANA = pd.read_csv('Data/paralogv2-library-files_v20-paralogcommonessentialcontrolsfromavana.csv')
# ParalogNonessentialControlsFromAVANA = pd.read_csv('Data/paralogv2-library-files_v20-paralognonessentialcontrolsfromavana.csv')
ParalogGeneEffect = pd.read_csv('Data/ParalogGeneEffect.csv')
# #AvanaLFC = pd.read_csv('Data/AvanaLogfoldChange.csv')
# # Load local counted plates data
counted_plates_sequence_df = pd.read_csv('Data/counted plates - sequenceid-samplename (1).csv')
counted_plates_sheet1_df = pd.read_csv('Data/counted plates - Sheet1 (counted plates).csv')

# LFC and norm counts
ParalogFullGuideMap = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullGuideMap')
ParalogFullLfcGeneScreen = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullLfcGeneScreen')
ParalogFullLfcGeneSeq = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullLfcGeneSeq')
#ParalogFullLfcSgrnaScreen = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullLfcSgrnaScreen')  # download_to_cache for raw
#ParalogFullLfcSgrnaSeq = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullLfcSgrnaSeq')  # download_to_cache for raw
ParalogFullNormCountsSgrnaSeq = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogFullNormCountsSgrnaSeq')  # download_to_cache for raw
#ParalogGuideAgreementScore = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogGuideAgreementScore')
#ParalogRawCounts = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogRawCounts')
#ParalogRawCountsAllpDNAWells = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogRawCountsAllpDNAWells')
#ParalogScreenMap = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogScreenMap')
#ParalogScreenQCReport = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogScreenQCReport')
ParalogSequenceMap = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogSequenceMap')
ParalogSequenceQCReport = tc.get(name='paralogscreen06062024-938c', version=5, file='ParalogSequenceQCReport')
ParalogCommonEssentialControlsFromAVANA = tc.get(name='paralogv2-library-files-d8b3', version=17, file='ParalogCommonEssentialControlsFromAVANA')
ParalogNonessentialControlsFromAVANA = tc.get(name='paralogv2-library-files-d8b3', version=17, file='ParalogNonessentialControlsFromAVANA')




common_ess = ParalogCommonEssentialControlsFromAVANA["Symbol"]
none_ess = ParalogNonessentialControlsFromAVANA["Symbol"]
counted_plates = pd.read_csv('/Users/wangyura/Desktop/shiny2/Data/counted plates - sequenceid-samplename (1).csv')
celllines = set([line for line, seq_id in zip(counted_plates["CellLine"], counted_plates["SequenceID"]) if "SEL000" in seq_id])
# avana_list = [line for line in celllines if line not in ['BxPC3', 'CCLFPEDS0001T(PEDS005TSUSP)', 'CCLFPEDS0001T (PEDS005TSUSP)', 'COLO205', 'Hs766T', 'PANC0327', 'SNUC5', 'HCT15', 'LS123', 'SSP25', 'SNUC2A', 'CL34', 'SW1710', 'RBE', 'ICC12', 'LS1034', 'PANC0504']]
# paralog_list = [line for line in avana_list if line not in ['CORL23', 'PANC1', 'LU65']]
#avana_list = list(celllines)  # Convert set to list for consistency
avana_list = [line for line in celllines if line not in ['BxPC3', 'CCLFPEDS0001T(PEDS005TSUSP)', 'CCLFPEDS0001T (PEDS005TSUSP)', 'COLO205', 'Hs766T', 'PANC0327', 'SNUC5','LU65', 'PANC1', 'CORL23','HCT15', 'LS123', 'SSP25', 'SNUC2A', 'CL34', 'SW1710', 'RBE', 'ICC12', 'LS1034', 'PANC0504']]
paralog_list = avana_list  # Use the same list for paralog_list
#AvanaLogfoldChange = tc.get(name='internal-24q2-3719', version=87, file='ScreenNaiveGeneScore')  # download_to_cache for raw
AvanaLogfoldChange = tc.get(name='paralogscreen08232024-3db1', version=11, file='AvanaLogfoldChange')
print("AvanaLogfoldChange: ", AvanaLogfoldChange)
AvanaLogfoldChange.set_index("Gene", inplace=True)
print("AvanaLogfoldChange2: ", AvanaLogfoldChange )
AvanaLogfoldChange_GPP = AvanaLogfoldChange
print("AvanaLogfoldChange_GPP: ", AvanaLogfoldChange_GPP)
ParalogFullLfcGeneScreen = tc.get(name='paralogscreen08232024-3db1', version=11, file='ParalogFullLfcGeneScreen')

ParalogFullLfcGeneScreen.set_index("GuideTargetSymbol", inplace=True)

ParalogFullLfcGeneScreen_GPP = ParalogFullLfcGeneScreen[paralog_list]

common_indices = ParalogFullLfcGeneScreen_GPP.index.intersection(AvanaLogfoldChange_GPP.index)
print("length of common indices", len(common_indices))
ParalogFullLfcGeneScreen_GPP = ParalogFullLfcGeneScreen_GPP.loc[common_indices]
AvanaLogfoldChange_GPP = AvanaLogfoldChange_GPP.loc[common_indices]

print("AvanaLogfoldChange_GPP : ", AvanaLogfoldChange_GPP)

##############################################################################################################
#Google Cloud Set Up 
##############################################################################################################

# Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the correct path
google_key_path = 'Data/paralog-sequencing-447e4f3ab1a7.json'
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = google_key_path

# Print statement to check if the key file exists
if os.path.exists(google_key_path):
    print(f"Google Cloud key found at: {google_key_path}")
else:
    print(f"Google Cloud key not found at: {google_key_path}")



def load_csv_from_gcs(bucket_name, file_path, delimiter=','):
    try:
        # Test if the environment variable is correctly set
        google_cred_path = os.getenv("GOOGLE_APPLICATION_CREDENTIALS")
        if google_cred_path:
            print(f"GOOGLE_APPLICATION_CREDENTIALS is set to: {google_cred_path}")
        else:
            print("GOOGLE_APPLICATION_CREDENTIALS is not set.")
        
        # Initialize the Google Cloud storage client
        storage_client = storage.Client()
        print("Google Cloud Storage client initialized successfully.")

        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(file_path)
        print(f"Attempting to download data from: gs://{bucket_name}/{file_path}")
        
        # Download the file as text
        data = blob.download_as_text()
        print(f"Successfully loaded data from {bucket_name}/{file_path}")

        # Determine the delimiter and load the data into a DataFrame
        if '\t' in data:
            data = pd.read_csv(io.StringIO(data), delimiter='\t')
        else:
            data = pd.read_csv(io.StringIO(data), delimiter=delimiter)
        
        return data
    except Exception as e:
        print(f"Error loading data from {bucket_name}/{file_path}: {e}")
        return pd.DataFrame()

def get_platemap_path(plate_id):
    path_row = counted_plates_sheet1_df[counted_plates_sheet1_df['PlateID'] == plate_id]
    if not path_row.empty:
        return path_row['PlateMap'].values[0]
    return None

def parse_gcs_path(gcs_path):
    # Ensure the path starts with "gs://"
    if gcs_path.startswith("gs://"):
        # Remove the "gs://" prefix
        gcs_path = gcs_path[5:]
        # Split into bucket name and file path
        parts = gcs_path.split("/", 1)
        
        if len(parts) == 2:
            bucket_name = parts[0]
            file_path = parts[1]
            return bucket_name, file_path
        
    # If the path doesn't start with "gs://" or the split fails, return None
    return None, None

def get_metrics_directory_path(plate_id):
    return f'gs://paralog-sequencing-files/Processed_Data/{plate_id}_counts/metrics_{plate_id}.csv'

def get_merge_directory_path(plate_id):
    # Retrieve the path from the 'PlateMap' column in the DataFrame based on the plate_id
    path_row = counted_plates_sheet1_df[counted_plates_sheet1_df['PlateID'] == plate_id]
    
    if not path_row.empty:
        # Get the path from the 'PlateMap' column
        plate_map_path = path_row['PlateMap'].values[0]
        
        # Check if the path ends with '.csv'
        if plate_map_path.endswith('.csv'):
            return plate_map_path
        else:
            print(f"Path for PlateID {plate_id} does not end with .csv: {plate_map_path}")
            return None  # Or handle it as needed (e.g., log a warning)
    
    return None  # Return None if no path is found

def list_files_in_gcs_directory(bucket_name, directory_path):
    try:
        storage_client = storage.Client()
        print("Google Cloud Storage client initialized successfully for listing files.")
        
        bucket = storage_client.bucket(bucket_name)
        print(f"Attempting to list files in directory: gs://{bucket_name}/{directory_path}")

        blobs = bucket.list_blobs(prefix=directory_path)
        files = [blob.name for blob in blobs if blob.name != directory_path]  # Avoid listing the directory itself
        print(f"Files in directory {directory_path}: {files}")
        return files
    except Exception as e:
        print(f"Error listing files in {bucket_name}/{directory_path}: {e}")
        return []
    
def process_metrics_data(metrics_data):
    if 'both_matched' in metrics_data.columns and 'matched_pairs' in metrics_data.columns:
        metrics_data['mismatched'] = metrics_data['both_matched'] - metrics_data['matched_pairs']
    else:
        metrics_data['mismatched'] = 0

    if 'total_reads' in metrics_data.columns and 'left_only' in metrics_data.columns and 'right_only' in metrics_data.columns and 'both_matched' in metrics_data.columns:
        metrics_data['none'] = metrics_data['total_reads'] - metrics_data['left_only'] - metrics_data['right_only'] - metrics_data['both_matched']
    else:
        metrics_data['none'] = 0

    return metrics_data

counted_plates_sequence_df.rename(columns={
    'SequenceID': 'PlateID',
    'SampleName': 'SequenceID',
    'CellLine': 'CellLineID'
}, inplace=True)

ParalogFullLfcGeneSeq.set_index('GuideTargetSymbol', inplace=True)


##############################################################################################################
#Shiny UI
##############################################################################################################

app_ui = ui.page_fluid(
    ui.input_action_button("show", "Show details of data"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h2("PlateID and CellLineID Selector"),
            ui.input_radio_buttons(
                "select_type",
                "Select by:",
                choices={"PlateID": "PlateID", "CellLineID": "CellLineID"},
                selected="PlateID"
            ),
            ui.input_text("input_value", "Enter value:"),
            ui.output_text_verbatim("output"),
            #ui.output_text_verbatim("output_avana")
        ),
        ui.layout_column_wrap(
            ui.accordion(
                ui.accordion_panel(
                    "Bar Plot",
                    ui.input_slider("percent_threshold", "Percentage Threshold", 0.0, 1.0, 0.1, step=0.01),
                    output_widget("generate_bar_plots"),
                ),
                ui.accordion_panel(
                    "Scatter Plot",
                    ui.card(
                        ui.popover(
                            "Options",
                            ui.input_radio_buttons(
                                "scatter_plot_data",  # Unique ID for Scatter Plot
                                "Choose data:",
                                ["ParalogFullLfcGeneSeq", "ParalogGeneEffect"],
                                selected="ParalogFullLfcGeneSeq",
                                inline=True,
                            ),
                            title="Choose data",
                        ),
                        class_="d-flex justify-content-between align-items-center",
                    ),
                    output_widget("scatter_matrix_plot"),
                    full_screen=True,
                ),
                ui.accordion_panel(
                    "Avana Scatter Plot",
                    ui.card(
                        ui.popover(
                            "Options",
                            ui.input_radio_buttons(
                                "avana_scatter_plot_data",  # Unique ID for Avana Scatter Plot
                                "Choose data source:",
                                choices=["Avana", "ParalogFullLfcGeneSeq"],
                                selected="Avana",
                                inline=True,
                            ),
                            title="Data Source Selection",
                        ),
                        class_="d-flex justify-content-between align-items-center",
                    ),
                    output_widget("Avana_scatter_plot"),
                    full_screen=True,
                ),
                ui.accordion_panel(
                    "Histogram Plot",
                    ui.card(
                        ui.card_header(
                            "Histogram Plot",
                            ui.popover(
                                "Options",
                                ui.input_radio_buttons(
                                    "histogram_color",
                                    "Color by:",
                                    ["CommonEssentials", "Nonessentials"],
                                    selected="CommonEssentials",
                                    inline=True,
                                ),
                                title="Add a color variable",
                            ),
                            class_="d-flex justify-content-between align-items-center",
                        ),
                        output_widget("histogram_widget"),
                        full_screen=True,
                    ),
                ),
                id="acc",
                open="Bar Plot",
          
                )
    
        ),
         title="Paralog Visualization"
    )

)

##############################################################################################################
#Shiny Server
##############################################################################################################

def server(input, output, session):
    @reactive.effect
    @reactive.event(input.show)
    def _():
        m = ui.modal(  
            "This is where we put the description message",  
            title="Detail messages ",  
            easy_close=True,  
        )  
        ui.modal_show(m) 
##############################################################################################################
# filtered_data
##############################################################################################################

    @reactive.Calc
    def filtered_data():
        selected_type = input.select_type()
        input_value = input.input_value().strip() if input.input_value().strip() else 'SEL00073'

        if not input_value:
            return pd.DataFrame()

        if selected_type == "PlateID":
            filtered = counted_plates_sequence_df[counted_plates_sequence_df['PlateID'].str.contains(input_value, na=False)]
        elif selected_type == "CellLineID":
            filtered = counted_plates_sequence_df[counted_plates_sequence_df['CellLineID'].str.contains(input_value, na=False)]
        else:
            return pd.DataFrame()
        return filtered
    
    
##############################################################################################################
# filtered_data_avana
##############################################################################################################
    
    # @reactive.Calc
    # def filtered_data_avana():
    #     selected_type = input.select_type()
    #     input_value = input.input_value().strip() if input.input_value().strip() else 'SEL00073'

    #     if not input_value:
    #         return pd.DataFrame()

    #     if selected_type == "PlateID":
    #         # If searching by PlateID, first find the corresponding CellLineIDs
    #         filtered_by_plate_id = counted_plates_sequence_df[counted_plates_sequence_df['PlateID'].str.contains(input_value, na=False)]
    #         cell_line_ids = filtered_by_plate_id['CellLineID'].unique()
    #         filtered = counted_plates_sequence_df[counted_plates_sequence_df['CellLineID'].isin(cell_line_ids)]
    #     elif selected_type == "CellLineID":
    #         filtered = counted_plates_sequence_df[counted_plates_sequence_df['CellLineID'].str.contains(input_value, na=False)]
    #     else:
    #         return pd.DataFrame()

    #     # Additional filtering for Avana plot (same as before)
    #     celllines = set(filtered['CellLineID'].unique())
    #     avana_list = [line for line in celllines if line not in ['BxPC3', 'CCLFPEDS0001T(PEDS005TSUSP)', 'CCLFPEDS0001T (PEDS005TSUSP)', 'COLO205', 'Hs766T', 'PANC0327', 'SNUC5','LU65', 'PANC1', 'CORL23','HCT15', 'LS123', 'SSP25', 'SNUC2A', 'CL34', 'SW1710', 'RBE', 'ICC12', 'LS1034', 'PANC0504']]
    #     return filtered[filtered['CellLineID'].isin(avana_list)]
        
##############################################################################################################
# sequence_used 
##############################################################################################################

    @reactive.Calc
    def sequence_used():
        data = filtered_data()
        if data.empty:
            return []
        return data['SequenceID'].unique().tolist()
    print("sequence_used: ", sequence_used)
    
##############################################################################################################
# sequence_used_avana
##############################################################################################################

    # @reactive.Calc
    # def sequence_used_avana():
    #     data = filtered_data_avana()  # Use the Avana-specific filtered data
    #     if data.empty:
    #         return []
    #     return data['SequenceID'].unique().tolist()
    
##############################################################################################################
# text output for scatter 
##############################################################################################################

    @output
    @render.text
    def output():
        data = filtered_data()
        
        if data.empty:
            return "No matches found."
        else:
            result = data[['PlateID', 'SequenceID', 'CellLineID']].to_string(index=False)
            return result
        
##############################################################################################################
# Bar Plot Logic 
##############################################################################################################

    @render_widget
    def generate_bar_plots():
        print("Starting the generate_bar_plots function...")
        
        data = filtered_data()

        if data.empty:
            print("No filtered data available")
            return None
        
        print(f"Filtered data contains {len(data)} rows.")

        # Get the unique Plate IDs from the filtered data
        plate_ids = data['PlateID'].unique()
        print(f"Unique Plate IDs: {plate_ids}")
        if len(plate_ids) == 0:
            
            return None

        # Assuming we are using the first Plate ID for simplicity
        plate_id = plate_ids[0]

        # Get the full path to the metrics file
        metrics_file_path = get_metrics_directory_path(plate_id)
        print(f"Metrics file path: {metrics_file_path}")

        if metrics_file_path:
            # Parse the GCS path into bucket name and file path
            bucket_name, file_path = parse_gcs_path(metrics_file_path)
            print(f"Parsed GCS path - Bucket: {bucket_name}, File path: {file_path}")


            # Load the CSV from GCS using the parsed bucket name and file path
            filtered_sheet1_data = load_csv_from_gcs(bucket_name, file_path, delimiter=',')
            print("Loaded metrics data from GCS.")

        if filtered_sheet1_data.empty:
            print("No data returned from GCS for the metrics file.")
            return None

        print(f"Metrics data preview (first 5 rows):\n{filtered_sheet1_data.head()}")



        # Load the CSV from GCS using the parsed bucket name and file path
        filtered_sheet1_data = load_csv_from_gcs(bucket_name, file_path, delimiter=',')
        # Your original data cleaning and processing logic
        if isinstance(filtered_sheet1_data, str):
            filtered_sheet1_data = pd.read_csv(io.StringIO(filtered_sheet1_data), index_col=None)

        filtered_sheet1_data.rename(columns={filtered_sheet1_data.columns[0]: 'SampleIndex'}, inplace=True)
        filtered_sheet1_data.columns = filtered_sheet1_data.columns.str.strip().str.replace(r'\\', '', regex=True)

        # Clean object columns by removing backslashes
        for column in filtered_sheet1_data.select_dtypes(include=['object']).columns:
            filtered_sheet1_data[column] = filtered_sheet1_data[column].str.replace(r'\\', '', regex=True)

        metrics_data = filtered_sheet1_data

        # Convert numeric columns to floats
        numeric_columns = ['total_reads', 'left_only', 'right_only', 'both_matched', 'matched_pairs', 'avg_mismatch_lib_corr']
        metrics_data[numeric_columns] = metrics_data[numeric_columns].astype(float)

        # Calculate 'mismatched' and 'none' columns
        metrics_data['mismatched'] = metrics_data['both_matched'] - metrics_data['matched_pairs']
        metrics_data['none'] = metrics_data['total_reads'] - metrics_data['left_only'] - metrics_data['right_only'] - metrics_data['both_matched']

        # Reshape the data for plotting
        df_plot = (
            metrics_data
            .drop(['avg_mismatch_lib_corr', 'both_matched'], axis=1, errors='ignore')
            .melt(id_vars=['SampleIndex', 'total_reads'], var_name='type', value_name='counts')
        )

        df_plot['percentage'] = df_plot['counts'] / df_plot['total_reads']
        print(f"Data for plotting after melt:\n{df_plot.head()}")

        # Ensure all required types are present in the DataFrame
        required_types = ['matched_pairs', 'left_only', 'right_only', 'mismatched', 'none']
        existing_types = df_plot['type'].unique()

        for required_type in required_types:
            if required_type not in existing_types:
                missing_type_df = df_plot[['SampleIndex', 'total_reads']].drop_duplicates()
                missing_type_df['type'] = required_type
                missing_type_df['counts'] = 0
                missing_type_df['percentage'] = 0
                df_plot = pd.concat([df_plot, missing_type_df], ignore_index=True)

        print(f"Data for plotting after ensuring all types are present:\n{df_plot.head()}")
        # Assuming get_merge_directory_path(plate_id) returns the full path to the merge CSV file
        merge_file_path = get_merge_directory_path(plate_id)

        # Parse the GCS path into bucket name and file path
        bucket_name, file_path = parse_gcs_path(merge_file_path)

        # Load the CSV directly using the full path
        if merge_file_path:
            merge_data = load_csv_from_gcs(bucket_name, file_path, delimiter=',')
            
            # Handle the case where the data is returned as a string
            if isinstance(merge_data, str):
                merge_data = pd.read_csv(io.StringIO(merge_data), index_col=None)
                print("Loaded merge data from GCS.")
                

            # Clean the column names in the merge_data
            merge_data.columns = merge_data.columns.str.strip().str.replace(r'\\', '', regex=True)
            print(f"Merge data preview (first 5 rows):\n{merge_data.head()}")

            # Merge the data with df_plot
            df_plot = merge_data.merge(df_plot, how='left', left_on='IndexBarcode1', right_on='SampleIndex')
            print(f"Data after merging:\n{df_plot.head()}")

            # Filter the data based on the percentage threshold from the UI input
            filtered_df = df_plot[df_plot['percentage'] >= input.percent_threshold()]
            print(f"Filtered data for bar plot (first 5 rows):\n{filtered_df.head()}")

            # Create the bar plot using Plotly Express
            bar_plot = px.bar(
                filtered_df, x='SampleIndex', y='percentage', color='type', barmode='stack',
                labels={'percentage': 'Percentage (%)'},
                hover_data={'SampleName': True},  # Assuming 'SampleName' exists in merge_data
                color_discrete_sequence=px.colors.sequential.Plasma_r,
                category_orders={'type': ['matched_pairs', 'left_only', 'right_only', 'mismatched', 'none']}
            )

            return bar_plot
        else:
            print("No merge file found.")
            return None
    



##############################################################################################################
# Scatter Plot
##############################################################################################################
    @render_widget
    def scatter_matrix_plot():
    
        sequences = sequence_used()
        print("Sequences used for scatter plot: ", sequences)
        if not sequences:
            return None

        # Use the unique input ID for Scatter Plot
        data_source = input.scatter_plot_data()
        if data_source == "ParalogFullLfcGeneSeq":
            lfc_used = ParalogFullLfcGeneSeq.loc[:, sequences].copy()
        elif data_source == "ParalogGeneEffect":
            lfc_used = ParalogGeneEffect.loc[:, sequences].copy()
        else:
            lfc_used = pd.DataFrame()

        if lfc_used.empty:
            return None

        lfc_used['TargetType'] = 'Others'
        lfc_used.loc[lfc_used.index.isin(ParalogCommonEssentialControlsFromAVANA['DepmapSymbol']), 'TargetType'] = 'CommonEssentials'
        lfc_used.loc[lfc_used.index.isin(ParalogNonessentialControlsFromAVANA['DepmapSymbol']), 'TargetType'] = 'Nonessentials'
        lfc_used.loc[lfc_used.index.str.contains('_'), 'TargetType'] = 'pair'

        color_map = {
            'CommonEssentials': 'red',
            'Nonessentials': 'blue',
            'pair': 'green',
            'Others': 'grey'
        }

        fig = px.scatter_matrix(
            lfc_used,
            dimensions=sequences,
            width=750,
            height=650,
            color='TargetType',
            color_discrete_map=color_map,
            title='Scatter Matrix',
            labels={col: col for col in sequences},
            hover_data={'index': lfc_used.index}
        )

        fig.update_traces(marker=dict(size=2))
        return fig

##############################################################################################################
# Scatter Plot with Avana
##############################################################################################################
    @render_widget
    def Avana_scatter_plot():
        # Ensure we have sequences to use
        sequences = ['OUMS23', 'SNU308']
        if not sequences:
            print("No sequences available for plotting.")
            return None

        # Filter and align data
        common_index = AvanaLogfoldChange_GPP.index.intersection(ParalogFullLfcGeneScreen_GPP.index)
        if common_index.empty:
            print("No common indices available for plotting.")
            return None

        lfc_avana = AvanaLogfoldChange_GPP.loc[common_index, sequences]
        lfc_paralog = ParalogFullLfcGeneScreen_GPP.loc[common_index, sequences]

        # Prepare combined DataFrame for plotting
        combined_data = pd.DataFrame({
            'AvanaLFC': lfc_avana.stack(),
            'ParalogLFC': lfc_paralog.stack()
        }).reset_index()
        combined_data.columns = ['Gene', 'CellLine', 'AvanaLFC', 'ParalogLFC']

        # Assign TargetType
        combined_data['TargetType'] = 'Others'
        combined_data.loc[combined_data['Gene'].isin(common_ess), 'TargetType'] = 'CommonEssentials'
        combined_data.loc[combined_data['Gene'].isin(none_ess), 'TargetType'] = 'Nonessentials'

        # Define color mapping
        color_map = {
            'CommonEssentials': 'red',
            'Nonessentials': 'blue',
            'Others': 'grey'
        }

        # Create scatter plot
        fig = px.scatter(
            combined_data,
            x='AvanaLFC',
            y='ParalogLFC',
            color='TargetType',
            color_discrete_map=color_map,
            facet_col='CellLine',  # Facet by cell line
            facet_col_wrap=2,      # Wrap after 2 columns
            title="Avana vs Paralog Scatter Plot",
            labels={
                'AvanaLFC': "Avana mean LFC",
                'ParalogLFC': "Mean LFC in our screen"
            },
            width=900,
            height=600
        )

        # Simplify facet titles
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

        # Add reference lines
        for shape in [
            dict(type="line", x0=-6, y0=-6, x1=2, y1=2, line=dict(color="black", width=1)),
            dict(type="line", x0=-6, y0=0, x1=2, y1=0, line=dict(color="gray", width=1, dash="dash")),
            dict(type="line", x0=0, y0=-6, x1=0, y1=2, line=dict(color="gray", width=1, dash="dash"))
        ]:
            fig.add_shape(**shape)

        fig.update_traces(marker=dict(size=6, opacity=0.7))

        # Return the single Plotly figure
        return fig

    # @render_widget
    # def Avana_scatter_plot():
    #     data = filtered_data()

    #     if data.empty:
    #         print("No filtered data available")
    #         return None
        
    #     sequences = data['CellLineID'].unique()
    
    #     if len(sequences)==0:
    #         print("No sequences found in the data")
    #         return None
        
    
    #     print("Sequence used: ", sequences)
    #     # Use the unique input ID for Avana Scatter Plot
    #     data_source = input.scatter_plot_data()
    #     if data_source == "Avana":
    #         lfc_used = AvanaLogfoldChange_GPP.loc[:, sequences].copy()  # Use Avana data
    #     elif data_source == "ParalogFullLfcGeneSeq":
    #         lfc_used = ParalogFullLfcGeneSeq.loc[:, sequences].copy()
    #     else:
    #         lfc_used = pd.DataFrame()

    #     if lfc_used.empty:
    #         return None
    #             # Create DataFrame for plotting, merging with Paralog data
    #     df = lfc_used.merge(ParalogFullLfcGeneScreen_GPP, left_index=True, right_index=True, how='inner', suffixes=('_avana', '_paralog'))
        
    #             # Create DataFrame for plotting

    #     df.loc[df.index.isin(common_ess), 'TargetType'] = 'CommonEssentials'
    #     df.loc[df.index.isin(none_ess), 'TargetType'] = 'Nonessentials'
    #     df.loc[df.index.str.contains('_'), 'TargetType'] = 'pair'  # If you have pairs

    #     color_map = {
    #         'CommonEssentials': 'red',
    #         'Nonessentials': 'blue',
    #         'pair': 'green',
    #         'Others': 'grey'
    #     }
    # # Create separate figures for each screen
    #     figs = []
    #     for screen in sequences:
    #         fig = px.scatter(
    #             df,
    #             x=screen,#screen
    #             y=screen,#paralog
    #             color='TargetType',
    #             color_discrete_map=color_map,
    #             width=750,
    #             height=650,
    #             title=f'Cell Line: {screen}',
    #             labels={'x': "Avana mean LFC", 'y': "Mean LFC in our screen"},
    #             hover_data={'index': df.index}
    #         )

    #     # fig = px.scatter(
    #     #     df,
    #     #     x='Avana_LFC',
    #     #     y='Paralog_LFC',
    #     #     color='TargetType',
    #     #     color_discrete_map=color_map,
    #     #     width=750,
    #     #     height=650,
    #     #     title=f'Cell Line: {input.cellline()}',  # Dynamic title
    #     #     labels={'Avana_LFC': "Avana mean LFC", 'Paralog_LFC': "Mean LFC in our screen"},
    #     #     hover_data={'index': df.index}
    #     # )

    #     # Add the rectangle and square annotations
    #     fig.add_shape(
    #         type="rect",
    #         x0=-6, y0=-6, x1=-1.5, y1=-1.5,
    #         line=dict(color="gray", width=1.5),
    #         fillcolor="#fdecec",
    #         opacity=0.5,
    #         layer="below"
    #     )

    #     # Add diagonal line
    #     fig.add_shape(
    #         type="line",
    #         x0=-6, y0=-6, x1=2, y1=2,
    #         line=dict(color="black", width=1),
    #         layer="below"
    #     )

    #     # Add dashed lines at x=0 and y=0
    #     fig.add_shape(
    #         type="line",
    #         x0=-6, y0=0, x1=2, y1=0,
    #         line=dict(color="gray", width=1, dash="dash"),
    #         layer="below"
    #     )
    #     fig.add_shape(
    #         type="line",
    #         x0=0, y0=-6, x1=0, y1=2,
    #         line=dict(color="gray", width=1, dash="dash"),
    #         layer="below"
    #     )

    #     # Calculate Pearson correlation for each screen
    #     # for screen in sequences:
    #     #     common_index = df.index
    #     #     if len(common_index) > 0:
    #     #         correlation, _ = pearsonr(df[f'{screen}_avana'], df[f'{screen}_paralog'])
    #     #         fig.add_annotation(
    #     #             x=-3, y=1 - 0.2 * sequences.index(screen),  # Adjust y-position for each annotation
    #     #             text=f"r ({screen})= {correlation:.3f}",
    #     #             showarrow=False,
    #     #             font=dict(size=14)
    #     #         )
    #     # Calculate Pearson correlation for this screen
    #     common_index = df.index
    #     if len(common_index) > 0:
    #         correlation, _ = pearsonr(df[screen], df[screen])
    #         fig.add_annotation(
    #             x=-3, y=1,
    #             text=f"r = {correlation:.3f}",
    #             showarrow=False,
    #             font=dict(size=14)
    #         )

    #     fig.update_traces(marker=dict(size=2))
    #     figs.append(fig)
    #         # Return a subplot with all the figures
    #     subplot_fig = make_subplots(rows=len(figs), cols=1)
    #     for i, fig in enumerate(figs):
    #         for trace in fig['data']:
    #             subplot_fig.add_trace(trace, row=i+1, col=1)
        

    #     return subplot_fig

        # fig.update_traces(marker=dict(size=2))
        # return fig

    #     lfc_compare_plot = lfc_used.merge(avana_lfc_gene_used, left_index=True, right_index=True)
    #     lfc_compare_plot['TargetType'] = 'Others'
    #     lfc_compare_plot.loc[lfc_used.index.isin(ParalogCommonEssentialControlsFromAVANA['DepmapSymbol']), 'TargetType'] = 'CommonEssentials'
    #     lfc_compare_plot.loc[lfc_used.index.isin(ParalogNonessentialControlsFromAVANA['DepmapSymbol']), 'TargetType'] = 'Nonessentials'
    #     lfc_compare_plot.loc[lfc_used.index.str.contains('_'), 'TargetType'] = 'pair'
   
    #     color_map = {
    #         'CommonEssentials': 'red',
    #         'Nonessentials': 'blue',
    #         'pair': 'green',
    #         'Others': 'grey'
    #     }

    #     fig = px.scatter_matrix(
    #         lfc_used,
    #         dimensions=sequences,
    #         width=750,
    #         height=650,
    #         color='TargetType',
    #         color_discrete_map=color_map,
    #         title='Scatter Matrix',
    #         labels={col: col for col in sequences},
    #         hover_data={'index': lfc_compare_plot.index}
    #     )

    #     fig.update_traces(marker=dict(size=2))
    #     return fig
##############################################################################################################
# Histograph Plot Logic 
##############################################################################################################

    @render_widget
    def histogram_widget():
        sequences = sequence_used()
        if not sequences:
            return None

        lfc_used = ParalogFullLfcGeneSeq.loc[:, sequences].copy()
        lfc1 = lfc_used.loc[ParalogCommonEssentialControlsFromAVANA['DepmapSymbol'], :]
        lfc1['TargetType'] = 'CommonEssentials'

        lfc2 = lfc_used.loc[ParalogNonessentialControlsFromAVANA['DepmapSymbol'], :]
        lfc2['TargetType'] = 'Nonessentials'
        lfc_plot = pd.concat([lfc1, lfc2])

        if sequences:
            print("LFC Plot Data:", lfc_plot.head())  # Debug print
            fig = plot_3d_ridge_plot(lfc_plot, sequences)  # Use all sequences
            return fig

##############################################################################################################
# Histograph Plot 
##############################################################################################################

def plot_3d_ridge_plot(df, samples):
    fig = go.Figure()

    for idx, sample in enumerate(samples):
        for target in df['TargetType'].unique():
            subset = df[df['TargetType'] == target]
            color = colors[target]

            # Add KDE for the target type
            kde = gaussian_kde(subset[sample])
            x = np.linspace(subset[sample].min(), subset[sample].max(), 1000)
            y = kde(x)
            z = np.ones_like(x) * idx

            print(f"Sample: {sample}, Target: {target}, X: {x[:5]}, Y: {y[:5]}, Z: [0, 0, 0, 0, 0]")  # Debug print

            fig.add_trace(go.Scatter3d(
                x=x, y=z, z=y,
                mode='lines',
                line=dict(color=color),
                name=f'{target} - {sample}'
            ))

            # Calculate and print the correlation coefficient
            for other_sample in samples:
                if sample != other_sample:
                    correlation = np.corrcoef(df[sample], df[other_sample])[0][1]
                    print(f"Correlation between {sample} and {other_sample}: {correlation}")

    fig.update_layout(
        title='3D Ridge Plot by TargetType',
        scene=dict(
            xaxis_title='Value',
            yaxis_title='Samples',
            zaxis_title='Density'
        ),
        width=800,
        height=600
    )

    return fig

#ParalogFullLfcGeneSeq, used for the 3D 

# Define colors for each TargetType
colors = {
    'CommonEssentials': 'red',
    'Nonessentials': 'blue'
}

##############################################################################################################
# Run App
##############################################################################################################
app = App(app_ui, server)

if __name__ == "__main__":
    app.run()
