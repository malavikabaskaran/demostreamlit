import streamlit as st
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from io import StringIO
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Streamlit app
st.title('Protein Data Analysis App')

# Retrieve protein data using Uniprot ID
def get_protein_data(uniprot_id):
    # Make API request to retrieve protein data
    protein_data = requests.get(f'https://www.uniprot.org/uniprot/{uniprot_id}.txt')
    return protein_data.text

# Function to analyze protein sequence
def analyze_protein_sequence(sequence):
    protein = ProteinAnalysis(sequence)
    return protein

# Fetch protein-protein interaction network
def fetch_interaction_network(uniprot_id):
    url = f"https://string-db.org/api/tsv/network?identifiers={uniprot_id}"
    response = requests.get(url)
    if response.ok:
        network_data = response.text
        network_df = pd.read_csv(StringIO(network_data), sep='\t')
        return network_df
    else:
        return None

# Main section of the app
def main():
    # Sidebar for user input
    st.sidebar.title('Protein Data Explorer:')
    option = st.sidebar.radio('Select Option', ['Uniprot ID', 'Protein Sequence'])
    user_input = None
    if option == 'Uniprot ID':
        user_input = st.sidebar.text_input('Enter Uniprot ID')
    elif option == 'Protein Sequence':
        user_input = st.sidebar.text_area('Enter Protein Sequence')

    st.sidebar.markdown('---')
    st.sidebar.text("Created by:")
    st.sidebar.text("Malavika Baskaran - B22EC0069")
    st.sidebar.text("Noor Hannani Syamimi Binti Mohd Suffian - A21EC0104")
    st.sidebar.text("Maathuree A/P Veerabalan - A21EC0051")

    if user_input:
        # Check if input is Uniprot ID or Protein Sequence
        if option == 'Protein Sequence':
            # Analyze Protein Sequence
            analyzed_data = analyze_protein_sequence(user_input)
            st.subheader('Protein Analysis Results:')
            st.write("Molecular Weight:", analyzed_data.molecular_weight())
            st.write("Amino Acid Composition:", analyzed_data.count_amino_acids())
            st.write("Aromaticity:", analyzed_data.aromaticity())
            st.write("Instability Index:", analyzed_data.instability_index())
            st.write("Isoelectric Point:", analyzed_data.isoelectric_point())
            st.write("Secondary Structure Fraction:", analyzed_data.secondary_structure_fraction())

        elif option == 'Uniprot ID':
            # Retrieve protein data using Uniprot ID
            protein_data = get_protein_data(user_input)
            st.subheader('Protein Characteristics:')
            st.code(protein_data)

            # Fetch and display protein-protein interaction network
            interaction_network = fetch_interaction_network(user_input)
            if interaction_network is not None:
                st.write("Protein-Protein Interaction Network:")
                st.write(interaction_network)
                ppi_graph = nx.from_pandas_edgelist(interaction_network, "preferredName_A", "preferredName_B")
                pos = nx.spring_layout(ppi_graph)  # Define node positions
                fig, ax = plt.subplots()  # Create a new figure
                nx.draw(ppi_graph, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=8, ax=ax)
                st.pyplot(fig)  # Display the graph in Streamlit with the specified figure
            else:
                st.error("Failed to fetch interaction network. Please try again later.")

if __name__ == '__main__':
    main()
