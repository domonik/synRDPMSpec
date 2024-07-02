include: "data_prep.smk"

rule run_on_synthetic_data:
    input:
        intensities=rules.multiplyData.output.multipleIntensities,
        design=rules.multiplyData.output.multipleDesign
    output:
        tsv = "Pipeline/RAPDORAnalysis/analyzedFake.tsv",
        json="Pipeline/RAPDORAnalysis/analyzedFake.json",
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmarkSynthetic.csv", config["benchmark_repeats"])
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RAPDORData(df=df,design=design,logbase=2)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance",kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.rank_table(["ANOSIM R", "Mean Distance"],ascending=(False, False))
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv,sep="\t")

rule ANOSIM_on_synthetic_data:
    input:
        intensities=rules.multiplyData.output.multipleIntensities,
        design=rules.multiplyData.output.multipleDesign
    output:
        tsv = "Pipeline/RAPDORAnalysis/analyzedFakeANOSIM.tsv",
        json="Pipeline/RAPDORAnalysis/analyzedFakeANOSIM.json",
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmarkANOSIMSynthetic.csv",10)
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RAPDORData(df=df,design=design,logbase=2)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance",kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.calc_anosim_p_value(permutations=999, mode="local", threads=1)
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv,sep="\t")


rule run_RAPDOR:
    input:
        intensities = rules.prepareinitialData.output.sanitized_df,
        design = rules.prepareinitialData.output.design
    output:
        tsv = "Pipeline/RAPDORAnalysis/GradRData.tsv",
        json = "Pipeline/RAPDORAnalysis/GradRData.json"
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmark.csv", config["benchmark_repeats"])
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities, sep="\t")
        design = pd.read_csv(input.design, sep="\t")
        rbpmdata = RAPDORData(df=df, design=design, logbase=2)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance", kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.rank_table(["ANOSIM R", "Mean Distance"], ascending=(False, False))
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv, sep="\t")
        del rbpmdata


rule runIdentifierOnSalmonella:
    input:
        intensities = rules.prepareSalmonellaData.output.file,
        design = rules.prepareSalmonellaData.output.file2,
    output:
        json = "Pipeline/Salmonella/SalmonellaAnalyzed.json"
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RAPDORData(df=df,design=design)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance", kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.rank_table(["Mean Distance"], ascending=(False,))
        rbpmdata.to_json(output.json)



rule runOnSynData:
    input:
        intensities=rules.prepareSynConditionedMS.output.intensities,
        design=rules.prepareSynConditionedMS.output.design,
    output:
        json = "Pipeline/ConditionedSynechocystis/analyzed/{condition}_Analyzed.json",
        tsv = "Pipeline/ConditionedSynechocystis/analyzed/{condition}_Analyzed.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RAPDORData(df=df,design=design)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance")
        rbpmdata.calc_all_scores()
        rbpmdata.rank_table(["ANOSIM R", "Mean Distance"],ascending=(False, False))
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv,sep="\t")

