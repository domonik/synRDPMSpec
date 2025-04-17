from RAPDOR.datastructures import RAPDORData
from multiprocessing import set_start_method

set_start_method("spawn", force=True)
import pandas as pd

def main(input, output, pval, impute):

    df = pd.read_csv(input.intensities, sep="\t")
    design = pd.read_csv(input.design, sep="\t")
    rbpmdata = RAPDORData(df=df, design=design, logbase=2, min_replicates=3)
    rbpmdata.run_preprocessing("Jensen-Shannon-Distance", kernel=3, impute=impute, impute_quantile=0.95)
    rbpmdata.calc_all_scores()
    rbpmdata.rank_table(["ANOSIM R", "Mean Distance"], ascending=(False, False))
    if pval:
        rbpmdata.calc_anosim_p_value(999, threads=1, mode="global")
    rbpmdata.pca()
    rbpmdata.to_json(output.json)
    rbpmdata.export_csv(output.tsv, sep="\t")
    del rbpmdata


if __name__ == '__main__':
    main(snakemake.input, snakemake.output, snakemake.params.pval, snakemake.params.impute)