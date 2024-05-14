from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_widget

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from handle_df import get_scatter_data, handle_original_df

DF_FILE = "qualityDSD_variants_effect_prediction.tsv"
COLUMNS = ["variant_id","CHROM", "POS", "REF", "ALT", 'AF_popmax', 'geneHancer', 'TAD', "distance_from_nearest_DSD_gene",
            "INTERVAL_ID", "total_probands", "probands_names", "healthy_members", "healthy_names",
            "conservation", "conservation_4way", "in_exon", "gene_name_of_exon", "gonad_exon",
            "contains_human_cells", "contains_mouse_cells"]

def display_df(df):
    df_display = df[COLUMNS].copy()
    df_display["REF"] = df_display["REF"].apply(lambda x: x[:15] + "..." if len(x) > 15 else x)
    df_display["ALT"] = df_display["ALT"].apply(lambda x: x[:15] + "..." if len(x) > 15 else x)
    df_display["gonad_exon"] = df_display["gonad_exon"].astype(str).apply(lambda x: x[:15] + "..." if len(x) > 15 else x)
    return df_display[COLUMNS]

df = pd.read_csv(DF_FILE, sep='\t', encoding='latin-1')
#reset index and change columns types
df = handle_original_df(df)
palette = sns.color_palette("pastel6").as_hex()

app_ui = ui.page_fillable(
    {"title": "Variant Effect Prediction"},
    ui.layout_columns(
        ui.card(
            ui.input_checkbox_group("filter_options", label="Choose filter(s)", choices=["Index", "TAD", "Interval_ID"], inline=True),
            ui.input_select("index_filter", label="Filter by index", choices=["All"] + list(df['variant_id'].astype(str).unique())),
            ui.input_select("tad_filter", label="Filter by TAD", choices=["All"] + sorted(df['TAD'].unique())),
            ui.input_select("interval_filter", label="Filter by interval ID", choices=["All"] + list(df['INTERVAL_ID'].unique())),
            height="800px"
        ),
        ui.card(output_widget("plot_scatter"), height="400px"),
        ui.card(ui.output_data_frame("summary_data"), height="400px"),
        
        col_widths=[3, 9, 12]
        #row_heights=[4, 2],
    ),
)
def server(input, output, session):
    def apply_filters():
        filters = input.filter_options()
        filtered_df = df.copy()
        if "Index" in filters:
            index_filter = input.index_filter()
            if index_filter != "All":
                filtered_df = filtered_df[filtered_df['variant_id'] == int(index_filter)]
        if "TAD" in filters:
            tad_filter = input.tad_filter()
            if tad_filter != "All":
                filtered_df = filtered_df[filtered_df['TAD'] == tad_filter]
        if "Interval_ID" in filters:
            interval_filter = input.interval_filter()
            if interval_filter != "All":
                filtered_df = filtered_df[filtered_df['INTERVAL_ID'] == interval_filter]
        return filtered_df

    @render.data_frame
    def summary_data():
        filtered_df = apply_filters()
        return render.DataGrid(display_df(filtered_df), selection_mode="rows")

    @reactive.calc
    def filtered_df():
        row = summary_data.data_view(selected=True)
        return row.index

    @render_widget
    def plot_scatter():
        selected_rows = filtered_df()
        if len(selected_rows) == 0:
            fig = go.Figure()
            fig.add_annotation(text="Please select a variant", showarrow=False, font=dict(size=20), x=0.5, y=0.5)
            fig.update_xaxes(type="linear", range=[0, 1], showticklabels=False)  # Set x-axis type to linear
            fig.update_yaxes(type="linear", range=[0, 1], showticklabels=False)  # Set y-axis type to linear
            fig.update_layout(
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
            )
            return fig
        else:
            data = get_scatter_data(df.iloc[selected_rows])
            p = px.scatter(data, x="original_score", y="variant_score", color="TF",
                       text="TF", title="Original vs Variant Scores", color_discrete_sequence=palette)
            p.add_shape(type="line",
                    x0=0, y0=0,
                    x1=1, y1=1,
                    line=dict(color="gray", width=2, dash="dash"))
            p.update_xaxes(type="linear", range=[-0.01, 1.01])  # Set x-axis type to linear
            p.update_yaxes(type="linear", range=[-0.01, 1.01])  # Set y-axis type to linear
            p.layout.update(showlegend=True)
            return p

app = App(app_ui, server)
