# %% ---------------------------- Imports ---------------------------- %%
from pathlib import Path
import pandas as pd

# %% ------------------------- Main Function ------------------------- %%
pathway_dir = Path("./brite_table")
all_tables = list(pathway_dir.glob("*.tsv"))

pathways = pd.read_excel("./pombe_kegg_brite.xlsx", index_col=[0]).reset_index()
pathways["Code"] = "spo" + pathways["Code"].astype(str).str.zfill(5)
    

dfs = []
for table in all_tables:
    df = pd.read_csv(table, sep="\t", header=None)
    dfs.append(df)

concat_df = pd.concat(dfs, ignore_index=True)
concat_df.columns = ["Name"] + ["Level_" + str(i) for i in range(1, len(concat_df.columns))]
concat_df = concat_df[concat_df["Name"].notna()]
concat_df = concat_df.ffill(axis=1)
# %%
concat_df = concat_df.merge(
    pathways,
    left_on="Level_1",
    right_on="Code",
    how="left"
).drop(columns=["Code"])

# %%
concat_df = concat_df[
    ["Name", "Category", "Description"] + [col for col in concat_df.columns if col.startswith("Level_")]
]

# concat_df = concat_df.query("Level_1 != 'spo00001'")
# %%
concat_df.to_csv("combined_brite_table.tsv", sep="\t", index=False)
# %%
