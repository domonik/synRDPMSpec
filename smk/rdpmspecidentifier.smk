include: "data_prep.smk"

rule run_on_synthetic_data:
    input:
        intensities="Pipeline/synIntensitiesFake.tsv",
        design="Pipeline/synDesignFake.tsv"
    output:
        tsv = "Pipeline/RDPMSpecAnalysis/analyzedFake.tsv",
        json="Pipeline/RDPMSpecAnalysis/analyzedFake.json",
    benchmark:
        repeat("Pipeline/benchmarks/rdpmspecBenchmarkSynthetic.csv",10)
    run:
        from RDPMSpecIdentifier.datastructures import RDPMSpecData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RDPMSpecData(df=df,design=design,logbase=2)
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
        tsv = "Pipeline/RDPMSpecAnalysis/analyzedFakeANOSIM.tsv",
        json="Pipeline/RDPMSpecAnalysis/analyzedFakeANOSIM.json",
    benchmark:
        repeat("Pipeline/benchmarks/rdpmspecBenchmarkANOSIMSynthetic.csv",10)
    run:
        from RDPMSpecIdentifier.datastructures import RDPMSpecData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RDPMSpecData(df=df,design=design,logbase=2)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance",kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.calc_anosim_p_value(permutations=999, mode="local", threads=1)
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv,sep="\t")


rule run_rdpmspecidentifier:
    input:
        intensities = "Data/synIntensities.tsv",
        design = "Data/synDesign.tsv"
    output:
        tsv = "Pipeline/RDPMSpecAnalysis/analyzed.tsv",
        json = "Pipeline/RDPMSpecAnalysis/analyzed.json"
    benchmark:
        repeat("Pipeline/benchmarks/rdpmspecBenchmark.csv", 10)
    run:
        from RDPMSpecIdentifier.datastructures import RDPMSpecData
        import pandas as pd

        df = pd.read_csv(input.intensities, sep="\t")
        design = pd.read_csv(input.design, sep="\t")
        rbpmdata = RDPMSpecData(df=df, design=design, logbase=2)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance", kernel=3)
        rbpmdata.calc_all_scores()

        rbpmdata.rank_table(["ANOSIM R", "Mean Distance"], ascending=(False, False))
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv, sep="\t")
