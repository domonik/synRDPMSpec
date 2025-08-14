include: "data_prep.smk"

rule run_on_synthetic_data:
    input:
        intensities=rules.multiplyData.output.multipleIntensities,
        design=rules.multiplyData.output.multipleDesign
    output:
        tsv = "Pipeline/RAPDORAnalysis/analyzedFake.tsv",
        json="Pipeline/RAPDORAnalysis/analyzedFake.json",
    threads: 999
    params:
        impute = False,
        pval = False
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmarkSynthetic.csv", config["benchmark_repeats"])
    script: "../pyfunctions/runRAPDOR.py"


rule ANOSIM_on_synthetic_data:
    input:
        intensities=rules.multiplyData.output.multipleIntensities,
        design=rules.multiplyData.output.multipleDesign
    output:
        tsv = "Pipeline/RAPDORAnalysis/analyzedFakeANOSIM.tsv",
        json="Pipeline/RAPDORAnalysis/analyzedFakeANOSIM.json",
    threads: 999
    params:
        impute = False,
        pval = True
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmarkANOSIMSynthetic.csv", config["benchmark_repeats"])
    script: "../pyfunctions/runRAPDOR.py"


rule run_RAPDOR:
    input:
        intensities = rules.prepareinitialData.output.sanitized_df,
        design = rules.prepareinitialData.output.design
    output:
        tsv = "Pipeline/RAPDORAnalysis/GradRData.tsv",
        json = "Pipeline/RAPDORAnalysis/GradRData.json"
    params:
        impute = True,
        pval = False
    threads: 999
    script: "../pyfunctions/runRAPDOR.py"

rule run_RAPDORForBenchmark:
    input:
        intensities = rules.prepareinitialData.output.sanitized_df,
        design = rules.prepareinitialData.output.design
    output:
        tsv = "Pipeline/RAPDORAnalysis/BenchmarkGradRData.tsv",
        json = "Pipeline/RAPDORAnalysis/BenchmarkGradRData.json"
    params:
        impute = False,
        pval = False
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/rapdorBenchmark.csv", config["benchmark_repeats"])
    script: "../pyfunctions/runRAPDOR.py"


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

rule exportImputedValues:
    input:
        json = rules.run_RAPDOR.output.json,
    output:
        csv = "Pipeline/RDeep/preparedIntensities.csv"
    script: "../pyfunctions/exportRAPDOR.py"

