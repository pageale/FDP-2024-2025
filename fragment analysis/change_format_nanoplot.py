import pandas as pd
import glob
from pathlib import Path

nanoplot_info_path = "/home/alessandra/Projects/FDP/fragment_distribution/nanoplot_nanopore/nanopore_downsampled_mapped/nanoplot/*/"
file_list = glob.glob(nanoplot_info_path + "/*.tsv.gz")
output_path = "/home/alessandra/Projects/FDP/fragment_distribution/samples_similar_cov/nanopore/"


for f in file_list:
    df = pd.read_csv(f, compression='gzip', sep='\t')
    name = Path(f).stem
    name = name.replace(".sorted.aligned.bam_NanoPlot-data.tsv", "") 
    print(name)
    histogram = df['aligned_lengths'].value_counts().sort_index()
    histogram_df = histogram.reset_index()
    histogram_df.columns = ['insert_size', 'All_Reads.fr_count']

    out_file = f"{output_path}/{name}_aligned_length_histogram.metrics"

    with open(out_file, 'w') as out:
        histogram_df.to_csv(out, sep='\t', index=False)

    print(f"Saved histogram for {name} to {out_file}")

