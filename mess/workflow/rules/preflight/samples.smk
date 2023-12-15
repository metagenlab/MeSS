replicates = list(range(1, config["replicates"] + 1))

samples = []
if os.path.isfile(config["input_path"]):
    files = config["input_path"]
else:
    files = list_files(config["input_path"], "tsv")

for file in files:
    df = pd.read_csv(file, sep="\t")
    try:
        samples.append(list(set(df["sample"])))
    except KeyError:
        samples.append([os.path.basename(file).split(".")[0]])

samples = list(chain.from_iterable(samples))
samples_rep = list(product(samples, replicates))
renamed_samples = [
    "-".join([os.path.basename(sample[0]).split(".")[0], str(sample[1])])
    for sample in samples_rep
]
