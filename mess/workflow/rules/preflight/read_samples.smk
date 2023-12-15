if os.path.isfile(INPUT):
    files = INPUT
else:
    files = list_files(INPUT, "tsv")
# concat tables into one and add sample columns
dfs = []
for file in files:
    df = pd.read_csv(file, sep="\t")
    dfs.append(df)
    try:
        samples = list(set(df["sample"]))
    except KeyError:
        sample = os.path.basename(file).split(".")[0]
        df["sample"] = [sample] * len(df)

DFS = pd.concat(dfs).reset_index(drop=True).sort_values("sample")
