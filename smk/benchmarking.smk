

rule meanRuntimeAndRAM:
    input:
        file = rules.run_RAPDORForBenchmark.benchmark,
        synthetic = rules.run_on_synthetic_data.benchmark,
        anosim = rules.ANOSIM_on_synthetic_data.benchmark,
        rdeep = rules.runRDeep.benchmark,
        rdeeporiginal = rules.originalRDeeP.benchmark,
        rapdorOnRdeep = rules.runRAPDORonRDeeP.benchmark,
    output:
        file ="Pipeline/benchmarks/meanStats.tsv",
        source_data ="Pipeline/Paper/SourceData/runtime.tsv"
    run:
        import pandas as pd
        import numpy as np
        d = {
            "RealData": (input.file, 3, False),
            "SyntheticData": (input.synthetic, 9, False),
            "SyntheticAnsoim": (input.anosim, 9, True),
            "RdeepData": (input.rdeep, 3, False),
            "RDeePonRDeeP": (input.rdeeporiginal, 3, False),
            "RAPDORonRDeeP": (input.rapdorOnRdeep, 3, False),
        }
        dfs = []
        source_data = []
        for key, (file, replicates, anosim) in d.items():
            df = pd.read_csv(file, sep="\t")

            stats = df[["s", "max_uss"]].mean(axis=0)
            for q in [0.25, 0.5, 0.75]:
                quantile = np.quantile(df[["s", "max_uss"]], q=q, axis=0)
                stats[f"q{q}_s"] = quantile[0]
                stats[f"q{q}_max_uss"] = quantile[1]
            stats["Tool"] = "RAPDOR" if key not in ("RdeepData", "RDeePonRDeeP") else "R-DeeP"
            stats["dataset"] = "Synechocystis GradR" if key not in ("RAPDORonRDeeP", "RDeePonRDeeP") else "HeLa R-DeeP"
            stats["stdev_s"] = np.std(df["s"])
            stats["stdev_max_uss"] = np.std(df["max_uss"])
            stats["ANOSIM"] =  anosim
            stats["Replicates per group"] = replicates
            stats["Name"] = key
            df["Tool"] = "RAPDOR" if key not in ("RdeepData", "RDeePonRDeeP") else "R-DeeP"
            df["dataset"] = "Synechocystis GradR" if key not in ("RAPDORonRDeeP", "RDeePonRDeeP") else "HeLa R-DeeP"
            df["Replicates per group"] = replicates


            source_data.append(df)
            dfs.append(stats)
        df = pd.concat(dfs, axis=1 , ignore_index=True).transpose()
        source_data = pd.concat(source_data, axis=0 , ignore_index=True)
        df.to_csv(output.file, sep="\t", index=False)
        source_data.to_csv(output.source_data, sep="\t", index=False)