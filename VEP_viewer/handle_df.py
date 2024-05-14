import pandas as pd
import re
import numpy as np

def extract_cols_data(row, pattern, delim):
    #find columns that match pattern
    cols = row.index.to_list()
    names = [s for s in cols if pattern.match(s)]
    #keep only columns with scores
    row = row[names]
    # #split rows with multiple values
    row = row.apply(lambda x: x.split('_') if isinstance(x, str) else x)
    row = row.explode()
    row = row.to_frame().transpose()
    data = pd.DataFrame(row)
    data = data.transpose()
    data.reset_index(inplace=True)
    # rename columns
    score_col = "score" + delim
    data.columns = ["tf", score_col]
    data = data[data[score_col].notna()]
    data["tf"] = data["tf"].str.split("_")
    data['tf'] = data['tf'].apply(lambda x: x[0])
    data = data.replace("NA", 0)
    return data

def get_scatter_data(row):
    row = row.squeeze(axis=0)
    # Regular expression pattern to match strings ending with '_original'
    orig_pattern = re.compile(r'.*_original$')
    #get original score data
    orig_data = extract_cols_data(row, orig_pattern, "_o")
    # Regular expression pattern to match strings ending with '_mutated'
    mut_pattern = re.compile(r'.*_mutated$')
    # get mutation score data
    mut_data = extract_cols_data(row, mut_pattern, "_m")
    #merge original and mutation data
    data = pd.concat([orig_data, mut_data['score_m']], axis=1)
    data.columns = ["TF", "original_score", "variant_score"]
    data = data.replace(np.nan, 0)
    return data

def handle_original_df(df):
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'variant_id'}, inplace=True)
    df['variant_id'] = df['variant_id'].apply(lambda x: x+1)
    df['contains_human_cells'] = df['contains_human_cells'].astype(str)
    df['contains_mouse_cells'] = df['contains_mouse_cells'].astype(str)
    df['distance_from_nearest_DSD_gene'] = df['distance_from_nearest_DSD_gene'].apply(lambda x: "{:,.0f}".format(float(x)))
    return df
