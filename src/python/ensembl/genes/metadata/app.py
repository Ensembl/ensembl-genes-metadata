import streamlit as st
import GB_metadata_reporting
import Fetch_annotations


assemblies = st.Page("GB_metadata_reporting.py", title="Assembly registry", icon=":material/genetics:")
annotations = st.Page("Fetch_annotations.py", title="Annotations", icon=":material/search_check_2:")

pg = st.navigation([assemblies,annotations])


# Local page config
st.set_page_config(
    page_title="Genebuild Metadata",
    page_icon=":material/genetics:",
    layout="wide",
    menu_items={
        'Get Help': 'https://github.com/Ensembl/ensembl-genes-metadata/tree/dev/gb_metadata_handling',
        'Report a bug': "mailto:lazar@ebi.ac.uk",
        'About': "# This is the Genebuild metadata reporting system. This is an *extremely* cool app!"
    }
)

st.logo("app_utils/embl_ebi_logo.svg", size="large")


# Read the CSS file
with open('app_utils/styles.css') as f:
    css = f.read()
st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)

pg.run()

# Function to validate BioProject ID
def is_valid_bioproject_id(bioproject_id):
    return True if bioproject_id.startswith("PRJNA") else False


# Mapping for more readable metrics
metric_mapping = {
    "gc_percent": "GC Percent",
    "total_sequence_length": "Total Sequence Length (bp)",
    "contig_n50": "Contig N50",
    "number_of_contigs": "Number of Contigs",
    "number_of_scaffolds": "Number of Scaffolds",
    "scaffold_n50": "Scaffold N50",
    "genome_coverage": "Genome Coverage (X)",
    "bioproject_id": "Bioproject ID",
    "asm_level": "Assembly Level",
    "gca_chain": "GCA Chain",
    "gca_version": "GCA Version",
    "asm_type": "Assembly Type",
    "release_date": "Release Date",
    "taxon_id": "Taxon ID"
}

def rename_metrics(df):
    """Renames the columns in the dataframe to more readable versions."""
    return df.rename(columns=metric_mapping)


def app():
    """Main function."""

    st.title("Genebuild Metadata Reporting")
    with st.container():
        st.caption("The assembly registry starts at _2019-01-01_. Apply filters below to query registry.")

        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            with st.popover("Bioproject ID", use_container_width=True, icon=":material/tactic:"):
                bioproject_ids = st.text_area(" Enter BioProject ID", height=68,
                                              placeholder="e.g., PRJNA391427, PRJNA607328")
        with col2:
            with st.popover("Release Date", use_container_width=True, icon=":material/calendar_month:"):
                release_date = st.text_input("Enter Release Date. Format: YYYY-MM-DD", placeholder="2019-01-01", max_chars=10)
        with col3:
            with st.popover("Taxon ID", use_container_width=True, icon=":material/badge:"):
                taxon_id = st.text_input("Enter Taxon ID", placeholder="40674",
                                         value="")

            # If taxon_id is provided, try to convert it to a number
            if taxon_id:
                try:
                    taxon_id = float(taxon_id)
                except ValueError:
                    st.error("Please enter a valid number for Taxon ID.")

        with col4:
            with st.popover("Assembly Metrics", use_container_width=True, icon=":material/query_stats:"):
                # Display filter options immediately
                excluded_metrics = ["bioproject_id", "gca_chain", "gca_version", "release_date", "taxon_id"]

                selected_friendly_metrics = st.pills(
                    "Select metrics to use as filters.",
                    [value for key, value in metric_mapping.items() if key not in excluded_metrics],
                    selection_mode="multi", label_visibility= "collapsed")

                selected_metrics = [key for key, value in metric_mapping.items() if value in selected_friendly_metrics]

                metric_thresholds = {}
                asm_level = None
                asm_type = None

                if selected_metrics:
                    metric_thresholds = {
                        metric: st.number_input(f"Set threshold for {metric_mapping[metric]}", min_value=0, value=None, step=10)
                        for metric in selected_metrics if metric not in ["asm_level", "asm_type"]
                    }

                    asm_level = st.pills("Select Assembly Level", ["Contig", "Scaffold", "Chromosome", "Complete genome"],
                                                 selection_mode="multi") if "asm_level" in selected_metrics else None

                    asm_type = st.pills("Select Assembly Type",
                                                ["haploid", "alternate-pseudohaplotype", "unresolved-diploid",
                                                 "haploid-with-alt-loci", "diploid"],
                                                selection_mode="multi") if "asm_type" in selected_metrics else None

        # Add a checkbox for 'Reference genome'
        with col5:
            with st.popover("Reference", use_container_width=True, icon=":material/check_circle:" ):
                add_reference_genome = st.toggle(
                        "Check if GCA is reference genome. This function uses the NCBI API and will break if used for large queries.",
                        value=False)

    if "assemblies_data" not in st.session_state:
        st.session_state.assemblies_data = None
        st.session_state.summary_data = None
        st.session_state.info_data = None
        st.session_state.gca_list = None

    bioproject_id_list = [id.strip() for id in bioproject_ids.split(",") if id.strip()]
    invalid_bioproject_ids = [id for id in bioproject_id_list if not is_valid_bioproject_id(id)]

    if invalid_bioproject_ids:
        st.error(
            f"The following BioProject IDs are invalid: {', '.join(invalid_bioproject_ids)}. Please enter valid BioProject IDs.",
            icon=":material/error:")
        return  # Stop further processing

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Keep the button enabled but prevent execution if no BioProject ID is provided
    if st.button("Get Assemblies", use_container_width=True, type="primary"):
        st.session_state.assemblies_data = None
        st.session_state.summary_data = None
        st.session_state.info_data = None
        st.session_state.gca_list = None
        status_placeholder = st.empty()

        with status_placeholder.status("Fetching assemblies... Please wait."):
            df, summary_df, info_result, gca_list = GB_metadata_reporting.get_filtered_assemblies(
                bioproject_id_list, metric_thresholds, all_metrics, asm_level, asm_type, release_date, taxon_id
            )

        status_placeholder.empty()

        if not isinstance(df, str):
            st.session_state.assemblies_data = rename_metrics(df)
            st.session_state.summary_data = rename_metrics(summary_df)
            st.session_state.info_data = info_result
            st.session_state.gca_list = gca_list

            # Apply the 'is_reference_genome' function if checkbox is checked
            if add_reference_genome:
                st.session_state.info_data["Reference genome"] = st.session_state.info_data["GCA"].apply(
                    GB_metadata_reporting.is_reference_genome)
        else:
            st.session_state.assemblies_data = df

    if st.session_state.assemblies_data is None:
        st.image("app_utils/background.svg", use_container_width=True)
    else:
        if isinstance(st.session_state.assemblies_data, str):
            st.error(st.session_state.assemblies_data, icon=":material/error:")
        else:
            st.divider()
            st.header("Filtered Assemblies")
            st.dataframe(st.session_state.assemblies_data, hide_index=True)
            st.download_button("Download Filtered Assemblies",
                               data=st.session_state.assemblies_data.to_csv(index=False),
                               file_name="filtered_assemblies.csv", mime="text/csv")

            st.divider()

            st.header("Summary Statistics")
            st.dataframe(st.session_state.summary_data)
            st.download_button("Download Summary Statistics", data=st.session_state.summary_data.to_csv(index=False),
                               file_name="summary_statistics.csv", mime="text/csv")

            st.divider()

            st.header("Assembly Info")
            st.dataframe(st.session_state.info_data, hide_index=True, column_config={
                "Release date": st.column_config.DateColumn(), "Lowest taxon ID": st.column_config.NumberColumn(format="%d"),
                "Species taxon ID":st.column_config.NumberColumn(format="%d"),})
            st.download_button("Download Assembly Info", data=st.session_state.info_data.to_csv(index=False),
                               file_name="assembly_info.csv", mime="text/csv")

            st.divider()

            st.header("Download GCA List")
            st.download_button("Download GCA List", data=st.session_state.gca_list.to_csv(index=False, header=False),
                               file_name="gca_list.csv", mime="text/csv")


if __name__ == "__main__":
    app()
