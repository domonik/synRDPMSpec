from RAPDOR.datastructures import RAPDORData
from multiprocessing import set_start_method

set_start_method("spawn", force=True)
import pandas as pd

def main(input, output, threads):

    df = pd.read_csv(input.intensities, sep="\t")
    design = pd.read_csv(input.design, sep="\t")
    rbpmdata = RAPDORData(df=df, design=design, logbase=2, min_replicates=3)
    rbpmdata.run_preprocessing("Jensen-Shannon-Distance", kernel=3, impute=True, impute_quantile=0.95)
    rbpmdata.calc_all_scores()
    rbpmdata.rank_table(["ANOSIM R", "Mean Distance"], ascending=(False, False))
    rbpmdata.calc_anosim_p_value(999, threads=threads, mode="global")
    rbpmdata.to_json(output.json)
    rbpmdata.export_csv(output.tsv, sep="\t")
    del rbpmdata


if __name__ == '__main__':
    main(snakemake.input, snakemake.output, snakemake.threads)