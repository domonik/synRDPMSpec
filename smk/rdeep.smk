
include: "rapdor.smk"

rule prepareForRDeep:
    input:
        file = config["intensities"],
        design= config["design"],
        rdpdata = rules.run_RAPDOR.output.json
    output:
        file = "Pipeline/RDeep/preparedIntensities.csv"
    run:
        import pandas as pd
        from RAPDOR.datastructures import RAPDORData
        import re
        import numpy as np
        df = pd.read_csv(input.file, sep="\t")
        design = pd.read_csv(input.design, sep="\t")
        with open(input.rdpdata) as handle:
            rap_data = RAPDORData.from_json(handle.read())


        df = df.loc[:, df.columns.isin(list(design.Name) + ["id"])]


        # Define a custom sorting key
        def custom_sort_key(column_name):
            parts = column_name.split('-')
            frac = int(parts[1])
            rep = int(parts[0][-1])
            treat = parts[0][:-1]
            return (frac, rep, treat)


        # Sort columns based on the custom sorting key
        sorted_columns = sorted(df.columns[1:] ,key=custom_sort_key)

        # Reorder DataFrame columns
        df = df[["id"] + sorted_columns]
        df.columns = ["name"] + list(df.columns[1:])
        df.iloc[:, 1:] = np.power(2, df.iloc[:, 1:])
        df = df.fillna(0)

        df.to_csv(output.file, sep=";", index=False)
        print(df)



rule runRDeep:
    input:
        raw_masspec_rdeep = rules.prepareForRDeep.output
    output:
        outfile = "Pipeline/RDeep/MS_Analysis_Shifts.csv",
        normalized_counts = "Pipeline/RDeep/normalized_rdeep_counts.tsv",
    benchmark: repeat("Data/benchmarks/rdeepbenchmark.txt", config["benchmark_repeats"]),
    conda: "../envs/rdeep.yml"
    script:
        "../Rscripts/MS_Statistical_Analysis.R"



rule runRAPDORonRDeePNorm:
    input:
        file = rules.runRDeep.output.normalized_counts,
        design = config["design"],
        orig_df = config["intensities"]
    output:
        json = "Pipeline/RDeep/RDeepNormalized.json"
    run:
        import pandas as pd
        from RAPDOR.datastructures import RAPDORData
        intensities = pd.read_csv(input.file, sep="\t", index_col=None)
        orig_intensities = pd.read_csv(input.orig_df, sep="\t")


        intensities.columns = intensities.columns.str.replace(".", "-")
        intensities["id"] = intensities.index
        intensities = intensities[["id"] + [col for col in intensities.columns if col != "id"]]
        orig_intensities = orig_intensities[["id", "Gene", "old_locus_tag"]]
        intensities = intensities.merge(orig_intensities, on="id")

        design = pd.read_csv(input.design ,sep="\t")
        rapdor_data = RAPDORData(df=intensities, design=design)
        rapdor_data.normalize_and_get_distances(method="Jensen-Shannon-Distance",kernel=3)
        rapdor_data.calc_all_scores()

        rapdor_data.rank_table(["ANOSIM R", "Mean Distance"],ascending=(False, False))

        rapdor_data.to_json(output.json)


rule plotSingleProteinRDeeP:
    input:
        file = rules.runRDeep.output.outfile,
    output:
        file = "Pipeline/RDeep/singleplot.png"
    conda: "../envs/rdeep.yml"
    script:
        "../Rscripts/plotSingleRDeeP.R"