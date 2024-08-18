# WGS-DSD_VEP_viewer
# Variant Effect Prediction Shiny App

## Overview

The **VEP Viewer** Shiny app provides an interactive platform for visualizing and analyzing the effects of genetic variants on gene expression, based on various filtering options. The app allows users to explore the variants table containing information about variants, their associated regions (e.g., TADs, INTERVAL_IDs), and other genomic features. Users can apply multiple filters, view detailed summaries of the data, and visualize the relationship between original and variant scores.

## Key Features

- **Interactive Filters**: Users can filter the dataset by `Index`, `TAD`, and `Interval_ID` to narrow down the variants of interest.
- **Dynamic Data Visualization**: The app generates scatter plots to visualize the relationship between original and variant scores for the selected variants.
- **Detailed Data Table**: A summary of the filtered data is displayed in a table, allowing users to explore the dataset in more detail.

## Dataset

The app uses a dataset (`qualityDSD_variants_effect_prediction.tsv` can be found in 'WGS_on_DSD/data/pipeline_outputs/variants_with_layers/2024-06-16/.') containing variant information. The dataset includes columns such as:

- `variant_id`: Unique identifier for each variant.
- `CHROM`, `POS`, `REF`, `ALT`: Genomic location and sequence information.
- `AF_popmax`: Maximum allele frequency in the population.
- `geneHancer`: GeneHancer information.
- `TAD`: Topologically Associating Domain.
- `distance_from_nearest_DSD_gene`: Distance from the nearest DSD gene.
- `INTERVAL_ID`: Interval ID.
- `TFBS_delta`: Change in Transcription Factor Binding Site.
- `total_probands`, `probands_names`: Information about probands.
- `healthy_members`, `healthy_names`: Information about healthy members.
- `conservation`, `conservation_4way`: Conservation scores.
- `in_exon`: Whether the variant is in an exon.
- `gene_name_of_exon`: Gene name of the exon.
- `gonad_exon`: Information about gonad-specific exons.
- `contains_human_cells`, `contains_mouse_cells`: According to ATAC-SEQ data.

## Installation and Usage

### Prerequisites

Ensure that you have the following installed:

- Python 3.x
- `pandas`, `plotly`, `seaborn`, `shiny`, and `shinywidgets` packages.

### Setup

1. Place the dataset (`qualityDSD_variants_effect_prediction.tsv`) in the same directory as the Shiny app script.
2. Install the required Python packages using pip:

   ```bash
   pip install pandas plotly seaborn shiny shinywidgets
   ```
3. Run the Shiny app:

   ```bash
   shiny run --reload --launch-browser app.py
   ```

### Usage
1. Open the app in your browser.
2. Use the checkbox group to select filters (e.g., Index, TAD, Interval_ID).
3. Select specific values from the dropdown menus for each filter.
4. The summary data table will update to reflect the filtered dataset.
5. Select specific row from data table.
6. View the scatter plot for selected variant to compare original and variant scores.
