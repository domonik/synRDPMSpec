
rule prepareForRDeep:
    input:
        file = "Data/synIntensities.tsv",
        design= "Data/synDesign.tsv"
    output:
        file = "Pipeline/RDeep/preparedIntensities.csv"
    run:
        import pandas as pd
        import re
        import numpy as np
        df = pd.read_csv(input.file, sep="\t")
        design = pd.read_csv(input.design, sep="\t")


        def rename_columns(column_name):
            match_free = re.match(r'^LFQ intensity RNAse_free_(\d+)_BR(\d+)$',column_name)
            match_rnase = re.match(r'^LFQ intensity RNAse_(\d+)_BR(\d+)$',column_name)

            if match_free:
                return f'ctrl{match_free.group(2)}-{match_free.group(1)}'
            elif match_rnase:
                return f'rnase{match_rnase.group(2)}-{match_rnase.group(1)}'
            else:
                return column_name
        df = df.loc[:, df.columns.isin(list(design.Name) + ["id"])]
        df = df.rename(columns=rename_columns)


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
        df.iloc[:, 1:] = np.square(df.iloc[:, 1:])
        df = df.fillna(0)

        df.to_csv(output.file, sep=";", index=False)
        print(df)



rule runRDeep:
    input:
        raw_masspec_rdeep = rules.prepareForRDeep.output
    output:
        outfile = "Pipeline/RDeep/MS_Analysis_Shifts.csv",
        normalized_counts = "Pipeline/RDeep/normalized_rdeep_counts.tsv"
    benchmark: repeat("Data/benchmarks/rdeepbenchmark.txt", 10),
    conda: "../envs/rdeep.yml"
    script:
        "../Rscripts/MS_Statistical_Analysis.R"
