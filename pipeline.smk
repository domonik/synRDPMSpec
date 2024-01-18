
include: "smk/rdeep.smk"
include: "smk/rapdor.smk"
include: "smk/data_prep.smk"
include: "smk/postProcessingAndPlots.smk"

rule all:
    input:
        #rdeep =  rules.runRDeep.output,
        rdpmspec = rules.plotMeanDistribution.output,
        p = rules.plotBarcodePlot.output,

        #joined = expand(rules.plotVennDiagramm.output, ribo=[False, True, "only"]),
        allvenns = rules.plotAllVenns.output,
        #distribution = expand(rules.plotDistribution.output, distributionids=config["distributions"].keys()),
        bubble_plot = expand(rules.createBubblePlot.output, highlight=["overlapping"]),
        #rna_binding = rules.extractGORNABinding.output,
        #enrich = rules.GOEnrichment.output,
        benchmark = rules.plotRuntime.output,
        #salmonella = rules.runIdentifierOnSalmonella.output,
        syn_cond = expand(rules.runOnSynData.output, condition=["COLD", "HEAT", "DARK", "N", "Fe"]),
        overlapping_data = rules.extract_overlapping_proteins.output,
        fig4 = rules.createFigure4.output,
        table1 = rules.createTable1.output,
        figs2 = rules.plotConditionedSynechochoColdRibo.output







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
        from RAPDOR.plots import COLOR_SCHEMES
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



