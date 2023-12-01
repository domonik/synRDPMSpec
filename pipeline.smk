

rule all:
    input:
        file = "Pipeline/images/VennDiagram.svg",
        file2="Pipeline/benchmarks/meanStats.tsv",
        multi= "Pipeline/analyzedFake.tsv",
        single= "Pipeline/analyzed.tsv",
        rdeep =  "Pipeline/RDeep/normalized_rdeep_counts.tsv"

rule multiplyData:
    input:
        intensities = "Data/synIntensities.tsv",
        design = "Data/synDesign.tsv"
    output:
        multipleIntensities = "Pipeline/synIntensitiesFake.tsv",
        multipleDesign = "Pipeline/synDesignFake.tsv",
    run:
        import pandas as pd
        intensities = pd.read_csv(input.intensities, sep="\t")
        design = pd.read_csv(input.design, sep="\t")
        print(intensities.columns)
        names = []
        rnases = []
        fractions = []
        replicates = []
        newcols = {"index": intensities.iloc[:, 0]}
        for x in intensities.columns:
            if x.startswith("LFQ intensity RNAse_") and "peptides" not in x:
                print(x)
                d = x.split("LFQ intensity RNAse_")[1].split("_")
                print(d)
                if len(d) == 3:
                    RNase = True
                    d = d[1:]
                else:
                    RNase = False
                print(d)
                frac = int(d[0])
                for nr in range(3):
                    br = d[1] + f".{nr}"
                    fakedata = intensities[x].sample(frac=1).values
                    name = x + f".{nr}"

                    names.append(name)
                    rnases.append(RNase)
                    fractions.append(frac)
                    replicates.append(br)
                    newcols[name] = fakedata

        design = {"Name": names, "RNase": rnases, "Fraction": fractions, "Replicate": replicates}
        df = pd.DataFrame(design)
        df.to_csv(output.multipleDesign, sep="\t", index=False)
        intensities = pd.DataFrame(newcols)
        intensities.to_csv(output.multipleIntensities, sep="\t", index=False)

rule run_on_synthetic_data:
    input:
        intensities="Pipeline/synIntensitiesFake.tsv",
        design="Pipeline/synDesignFake.tsv"
    output:
        tsv = "Pipeline/analyzedFake.tsv",
        json="Pipeline/analyzedFake.json",
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
        intensities="Pipeline/synIntensitiesFake.tsv",
        design="Pipeline/synDesignFake.tsv"
    output:
        tsv = "Pipeline/analyzedFakeANOSIM.tsv",
        json="Pipeline/analyzedFakeANOSIM.json",
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
        tsv = "Pipeline/analyzed.tsv",
        json = "Pipeline/analyzed.json"
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
        rbpmdata.calc_anosim_p_value(999, 1, distance_cutoff=0.1)
        rbpmdata.to_json(output.json)
        rbpmdata.export_csv(output.tsv, sep="\t")


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
        print(df)
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
    conda: "envs/rdeep.yml"
    script:
        "Rscripts/MS_Statistical_Analysis.R"

rule meanRuntime:
    input:
        file = rules.run_rdpmspecidentifier.benchmark,
        synthetic = rules.run_on_synthetic_data.benchmark,
        anosim = rules.ANOSIM_on_synthetic_data.benchmark,
        rdeep = rules.runRDeep.benchmark
    output:
        file ="Pipeline/benchmarks/meanStats.tsv"
    run:
        import pandas as pd
        import numpy as np
        d = {
            "RealData": (input.file, 3, False),
            "SyntheticData": (input.synthetic, 9, False),
            "SyntheticAnsoim": (input.anosim, 9, True),
            "RdeepData": (input.rdeep, 3, False)
        }
        dfs = []
        for key, (file, replicates, anosim) in d.items():
            df = pd.read_csv(file, sep="\t")
            stats = df[["s", "max_uss"]].mean(axis=0)
            for q in [0.25, 0.5, 0.75]:
                quantile = np.quantile(df[["s", "max_uss"]], q=q, axis=0)
                stats[f"q{q}_s"] = quantile[0]
                stats[f"q{q}_max_uss"] = quantile[1]
            stats["ANSOIM"] =  anosim
            stats["Replicates per group"] = replicates
            stats["Name"] = key
            dfs.append(stats)
        df = pd.concat(dfs, axis=1 , ignore_index=True).transpose()
        df.to_csv(output.file, sep="\t", index=False)


rule joinSVMandGradR:
    input:
        svm = "Data/svm.tsv",
        gradR = "Pipeline/analyzed.tsv"
    output:
        svg = "Pipeline/images/VennDiagram.svg",
        html = "Pipeline/images/VennDiagram.html",
        svg_ribo= "Pipeline/images/VennDiagramRiboIncluded.svg",
        html_ribo="Pipeline/images/VennDiagramRiboIncluded.html",
        joined = "Pipeline/joined.tsv",
        intersection = "Pipeline/intersect.tsv",
    run:
        import pandas as pd
        import plotly.express as px
        import plotly.graph_objects as go
        from scipy.stats import spearmanr
        import numpy as np
        from RDPMSpecIdentifier.plots import COLOR_SCHEMES
        import plotly.io as pio
        import copy
        from pyfunctions.plotlyVenn import venn_to_plotly

        NUMBER_PROTEINS = 200
        pio.templates["FFDefault"] = copy.deepcopy(pio.templates["plotly_white"])

        pio.templates["FFDefault"].update(
            {
                "layout": {
                    # e.g. you want to change the background to transparent
                    "paper_bgcolor": "rgba(255,255,255,1)",
                    "plot_bgcolor": " rgba(0,0,0,0)",
                    "font": dict(color="black"),
                    "xaxis": dict(linecolor="black",showline=False),
                    "yaxis": dict(linecolor="black",showline=False),
                    "coloraxis": dict(colorbar=dict(outlinewidth=1,outlinecolor="black"))
                }
            }
        )

        template = pio.templates["FFDefault"]
        svm = pd.read_csv(input.svm, sep="\t", decimal=",")
        rdpm = pd.read_csv(input.gradR, sep="\t")
        svm = svm.rename({"Gene": "old_locus_tag", "Score": "SVMScore"}, axis=1)
        svm = svm.drop_duplicates(subset=["old_locus_tag"])


        result_df = pd.merge(svm, rdpm, on="old_locus_tag")
        result_df = result_df.sort_values(by="SVMScore", ascending=False)
        result_df["SVMRank"] = np.arange(1, len(result_df)+1)
        result_df.to_csv(output.joined, sep="\t", index=False)
        result_df["Gene"] = result_df["Gene"].fillna("None")
        no_ribosomes = result_df[~result_df['Gene'].str.contains('rpl|rps|Rpl23|Rps')]
        no_ribosomes = result_df[~result_df['ProteinFunction'].str.contains('ribosomal protein')]
        #top200SVM = set(no_ribosomes.iloc[0:200]["old_locus_tag"])
        top200SVM = set(no_ribosomes[(no_ribosomes["SVMRank"] <= NUMBER_PROTEINS) & (no_ribosomes["Prediction"] == "RNA-binding protein")]["old_locus_tag"])
        top200SVMribo = set(result_df[(result_df["SVMRank"] <= NUMBER_PROTEINS) & (result_df["Prediction"] == "RNA-binding protein")]["old_locus_tag"])
        no_ribosomes = no_ribosomes.sort_values(by="Rank", ascending=True)

        #top200RDPM = set(no_ribosomes.iloc[0:200]["old_locus_tag"])
        top200RDPM = set(no_ribosomes[no_ribosomes["Rank"] <= NUMBER_PROTEINS]["old_locus_tag"])
        top200RDPMribo = set(result_df[result_df["Rank"] <= NUMBER_PROTEINS]["old_locus_tag"])
        duplicates = no_ribosomes[no_ribosomes.duplicated('old_locus_tag',keep=False)]
        print(duplicates['old_locus_tag'].unique())

        print(len(top200RDPM))
        intersect = top200SVM.intersection(top200RDPM)
        fig = venn_to_plotly(
            (top200SVM, top200RDPM),
            L_labels=("TripepSVM", "RDPMSpecidentifier"),
            L_color=COLOR_SCHEMES["Dolphin"])
        fig.update_layout(font=dict(size=18), template=template, xaxis=dict(showgrid=False, zeroline=False), yaxis=dict(showgrid=False, zeroline=False))
        fig.write_image(output.svg)
        fig.write_html(output.html)
        intersect = no_ribosomes[no_ribosomes["old_locus_tag"].isin(intersect)]
        intersect.to_csv(output.intersection, sep="\t", index=False)

        fig = venn_to_plotly(
            (top200SVMribo, top200RDPMribo),
            L_labels=("TripepSVM", "RDPMSpecidentifier"),
            L_color=COLOR_SCHEMES["Dolphin"])
        fig.update_layout(font=dict(size=18),template=template,xaxis=dict(showgrid=False,zeroline=False),yaxis=dict(showgrid=False,zeroline=False))
        fig.write_image(output.svg_ribo)
        fig.write_html(output.html_ribo)



