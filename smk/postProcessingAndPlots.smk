import os.path
import pandas as pd
import numpy as np

include: "rapdor.smk"
include: "benchmarking.smk"
include: "runRAPDORonRDeePData.smk"

from RAPDOR.plots import DEFAULT_TEMPLATE, COLOR_SCHEMES


DEFAULT_COLORS = COLOR_SCHEMES["Dolphin"]

DEFAULT_TEMPLATE = DEFAULT_TEMPLATE.update(
    layout=dict(
        font=config["fonts"]["default"],
        xaxis=dict(
            title=dict(font=config["fonts"]["axis"]),
            tickfont=config["fonts"]["axis_ticks"],

        ),
        yaxis=dict(
            title=dict(font=config["fonts"]["axis"]),
            tickfont=config["fonts"]["axis_ticks"],
        ),
        legend=dict(font=config["fonts"]["legend"]),
        legend2=dict(font=config["fonts"]["legend"]),
        annotationdefaults=dict(font=config["fonts"]["annotations"]),
        margin=config["margin"],
        coloraxis=dict(colorbar=dict(tickfont=config["fonts"]["legend"]))

    )

)

rule plotIntensitiesDistribution:
    input:
        json = rules.run_RAPDOR.output.json
    output:
        html = "Pipeline/plots/IntensitiesFractionDistribution.html"
    run:
        from RAPDOR.datastructures import RAPDORData
        from plotly.subplots import make_subplots
        import plotly.graph_objs as go
        rapdor = RAPDORData.from_file(input.json)
        array = np.log2(rapdor.array)
        x = np.tile(np.arange(array.shape[-1]) + 1, array.shape[0])

        fig = make_subplots(rows=array.shape[-2], cols=1)
        for i in range(array.shape[-2]):
            y = array[:, i, :]
            _sum = np.log2(rapdor.array[:, i, :].sum(axis=1))
            _sum = _sum[_sum > 0]
            y = y.flatten()
            indices = np.argwhere(y > 0)
            fig.add_trace(
                go.Box(
                    x=x[indices].flatten(),
                    y=y[indices].flatten(),
                    width=0.9
                ),
                row=i+1,
                col=1
            )
            fig.add_trace(
                go.Box(
                    x=x[indices].flatten(),
                    y=y[indices].flatten(),
                    width=0.9
                ),
                row=i + 1,
                col=1
            )
            fig.add_trace(
                go.Box(
                    x=[array.shape[-1] + 1 for _ in _sum],
                    y=_sum
                ),
                row=i + 1,
                col=1
            )
        fig.show()



rule fitMultiGaussian:
    input:
        json = rules.run_RAPDOR.output.json,
    output:
        file = "Pipeline/gaussianFit/table.tsv"
    threads: 8
    run:
        from pyfunctions.helpers import determineFits
        df = determineFits(input.json, nr_cores=threads)
        df.to_csv(output.file)


rule pca:
    input:
        file=rules.run_RAPDOR.output.json,
    output:
        svg = "Pipeline/plots/QC_PCA.svg",
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        import math
        from RAPDOR.plots import plot_sample_pca, COLOR_SCHEMES
        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        dolphin = COLOR_SCHEMES["Dolphin"]
        fig = plot_sample_pca(data, plot_dims=(1, 2), ntop=0.2, summarize_fractions=True, colors=dolphin)
        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            margin=dict(b=50, l=30, r=30, t=20),
            height=config["QCPlot"]["B"]["height"],
            width=config["width"] // 2
        )
        fig.write_image(output.svg)

rule sampleCorrelation:
    input:
        file=rules.run_RAPDOR.output.json,
    output:
        svg = "Pipeline/plots/QC_Correlation.svg",
        png = "Pipeline/plots/QC_Correlation.png",
        corr= "Pipeline/plots/QC_Correlation_onlyReplicate.svg",
        corr_png= "Pipeline/plots/QC_Correlation_onlyReplicate.png",
        corr_html= "Pipeline/plots/QC_Correlation_onlyReplicate.html",
        histo = "Pipeline/plots/QC_CorrHistogram.svg",
    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import plot_sample_correlation, COLOR_SCHEMES, plot_sample_histogram, plot_sample_pca
        data = RAPDORData.from_file(input.file)
        dolphin = list(COLOR_SCHEMES["Dolphin"])
        dolphin.insert(1, "white")
        fig = plot_sample_correlation(
            data,
            method="pearson",
            summarize_fractions=False,
            use_raw=False,
            highlight_replicates=False,
            ntop=None,
            colors=dolphin
        )
        fig.update_layout(template=DEFAULT_TEMPLATE, margin=dict(b=50, l=50, r=50, t=50), height=config["QCPlot"]["A"]["height"])
        fig.update_xaxes(showticklabels=False)
        fig.write_image(output.svg)
        fig.write_image(output.png)
        fig = plot_sample_correlation(
            data,
            method="spearman",
            summarize_fractions=True,
            use_raw=False,
            highlight_replicates=False,
            show_values=True,
            ntop=None,
            colors=dolphin
        )
        fig.update_layout(template=DEFAULT_TEMPLATE,margin=dict(b=70,l=50,r=50,t=25),
            height=config["QCPlot"]["A"]["height"],
            width=config["width"] // 2,
        )
        x_vals = [trace.x for trace in fig.data if hasattr(trace,'x')][0]
        new_ticktext = [x.split(' - ')[-1] for x in x_vals]
        fig.update_yaxes(showticklabels=True, range=(-.5, 5.5), ticktext=new_ticktext, tickvals=x_vals, tickmode='array',)
        fig.update_xaxes(showticklabels=True, range=(-.5, 5.5), ticktext=new_ticktext, tickvals=x_vals, tickmode='array',)

        fig.update_traces(colorbar=dict(thickness=10))
        y = -0.1
        fig.add_shape(type="line", x0=0, x1=2, y0=y, y1=y, yref="paper", line_color="black")
        fig.add_shape(type="line", x0=3, x1=5, y0=y, y1=y, yref="paper", line_color="black")
        fig.add_annotation(
            yref="paper",
            text="Control",
            x=1,
            y=y,
            yanchor="top",
            showarrow=False
        )
        fig.add_annotation(
            yref="paper",
            text="RNase",
            x=4,
            y=y,
            yanchor="top",
            showarrow=False
        )
        fig.add_shape(type="line", y0=0, y1=2, x0=-0.15, x1=-0.15, xref="paper",line_color="black")
        fig.add_shape(type="line", y0=3, y1=5, x0=-0.15, x1=-0.15, xref="paper",line_color="black")
        fig.add_annotation(
            xref="paper",
            text="Control",
            y=1,
            x=-0.15,
            xanchor="right",
            showarrow=False,
            textangle = 270

        )
        fig.add_annotation(
            xref="paper",
            text="RNase",
            y=4,
            x=-0.15,
            xanchor="right",
            showarrow=False,
            textangle=270
        )
        fig.write_image(output.corr)
        fig.write_image(output.corr_png)
        fig.write_html(output.corr_html)
        fig = plot_sample_histogram(rapdordata=data,method="spearman",colors=COLOR_SCHEMES["Dolphin"])
        fig.update_layout(template=DEFAULT_TEMPLATE,margin=dict(b=50,l=70,r=100,t=60),
            height=config["QCPlot"]["C"]["height"],
            width=config["width"],
        )
        fig.update_annotations(font=config["fonts"]["annotations"])
        fig.write_image(output.histo)


rule plotSumofIntensities:
    input:
        json = rules.run_RAPDOR.output.json,
    output:
        html = "Pipeline/plots/SumOfIntensities.html"
    script: "../pyfunctions/sumOfIntensities.py"



rule QCRDeePData:
    input:
        rdeep = rules.runRAPDORonRDeeP.output.json,
    output:
        svg = "Pipeline/plots/QC_Correlation_RDeeP.svg",
        png = "Pipeline/plots/QC_Correlation_RDeeP.png",
        corr = "Pipeline/plots/QC_Correlation_onlyReplicate_RDeeP.svg",
        corr_png = "Pipeline/plots/QC_Correlation_onlyReplicate_RDeeP.png",
        histo = "Pipeline/plots/QC_CorrHistogram_RDeeP.svg",
        pca = "Pipeline/plots/QC_PCA_RDeeP.svg",
    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import plot_sample_correlation, COLOR_SCHEMES, plot_sample_histogram, plot_sample_pca

        data = RAPDORData.from_file(input.rdeep)
        dolphin = list(COLOR_SCHEMES["Dolphin"])
        dolphin.insert(1,"white")
        fig = plot_sample_correlation(
            data,
            method="pearson",
            summarize_fractions=False,
            use_raw=False,
            highlight_replicates=False,
            ntop=None,
            colors=dolphin
        )
        fig.update_layout(template=DEFAULT_TEMPLATE,margin=dict(b=50,l=50,r=50,t=50),height=config["QCPlot"]["A"][
            "height"])
        fig.update_xaxes(showticklabels=False)
        fig.write_image(output.svg)
        fig.write_image(output.png)
        fig = plot_sample_histogram(rapdordata=data,method="spearman",colors=COLOR_SCHEMES["Dolphin"])
        fig.update_layout(template=DEFAULT_TEMPLATE,margin=dict(b=70,l=70,r=100,t=100),height=config["QCPlot"]["A"][
            "height"])
        fig.write_image(output.histo)
        fig = plot_sample_pca(data, plot_dims=(1, 2),ntop=0.2, summarize_fractions=True, colors=COLOR_SCHEMES["Dolphin"])
        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            margin=dict(b=50,l=70,r=200,t=20),
            height=config["QCPlot"]["B"]["height"],
            width=config["width"]
        )
        fig.write_image(output.pca)
        fig = plot_sample_correlation(
            data,
            method="spearman",
            summarize_fractions=True,
            use_raw=False,
            highlight_replicates=False,
            show_values=True,
            ntop=None,
            colors=dolphin
        )
        fig.update_layout(template=DEFAULT_TEMPLATE,margin=dict(b=50,l=50,r=50,t=50),
            height=config["QCPlot"]["A"]["height"],
            width=config["width"] // 2
        )
        fig.update_xaxes(showticklabels=False)
        fig.write_image(output.corr)
        fig.write_image(output.corr_png)




rule joinAnalysisMethods:
    input:
        svm = "Data/svm.tsv",
        rapdor = rules.run_RAPDOR.output.tsv,
        rdeep = rules.runRDeep.output.outfile,
    output:
        file = "Pipeline/postProcessing/joinedTable.tsv"
    run:
        import pandas as pd
        import numpy as np
        svm_table = pd.read_csv(input.svm, sep="\t",  decimal=",")
        svm_table["SVM Rank"] = np.arange(1,len(svm_table) + 1)

        rapdor_table = pd.read_csv(input.rapdor, sep="\t")
        rdeep_table = pd.read_csv(input.rdeep, sep=";")

        indices = list(rdeep_table.columns[0:1]) + list(rdeep_table.columns[-14:])
        rdeep_table = rdeep_table[indices]
        rdeep_table = rdeep_table[rdeep_table["nb_shift"] > 0]
        rdeep_table = rdeep_table.sort_values(by=["rnase_peak_p_value", "ctrl_peak_p_value"])
        rdeep_table = rdeep_table.drop_duplicates("protein_name")
        rdeep_table = rdeep_table.rename({"protein_name": "RAPDORid"}, axis=1)

        rapdor_table = rapdor_table[["RAPDORid"] + [col for col in rapdor_table.columns if col != "RAPDORid"]]

        rapdor_table = rapdor_table.merge(rdeep_table, on="RAPDORid", how="left")

        svm_table = svm_table.rename({"Gene": "old_locus_tag", "Score": "SVM Score"}, axis=1)
        svm_table["SVM RNA-binding"] = svm_table["Prediction"] == "RNA-binding protein"

        svm_table = svm_table[["old_locus_tag", "SVM Score", "SVM RNA-binding", "SVM Rank"]]
        svm_table = svm_table.sort_values(by="SVM Score", ascending=False)
        svm_table = svm_table.drop_duplicates("old_locus_tag")

        rapdor_table = rapdor_table.merge(svm_table, on="old_locus_tag", how="left")


        rapdor_table = rapdor_table.sort_values(by="SVM Score",ascending=False)
        rapdor_table["Gene"] = rapdor_table["Gene"].fillna("None")
        rapdor_table["Gene"] = rapdor_table["Gene"].replace("None", pd.NA)

        rapdor_table["RDeeP significant"] = (rapdor_table["ctrl_peak_p_value"] <= 0.05) | (rapdor_table["rnase_peak_p_value"] <= 0.05 )

        rapdor_table = rapdor_table.sort_values(by="Rank", ascending=False)

        rapdor_table.to_csv(output.file ,sep="\t", index=False)

rule extract_tophits:
    input:
        tsv = rules.joinAnalysisMethods.output.file
    output:
        tsv = "Pipeline/Paper/unused/Tables/overlappingProteins.tsv",
        tsv2 = "Pipeline/Paper/unused/Tables/overlappingProteins_no_ribosomes.tsv",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")

        svm_candidates = [300, 643, 421]
        top200SVM = set(
            df[df["RAPDORid"].isin(svm_candidates)]["old_locus_tag"]
        )
        rdp_candidates = config["bubble_plot"]["locus_tag"]
        topRDPM = set(
            df[df["old_locus_tag"].isin(rdp_candidates)]["old_locus_tag"]
        )
        df = df[
            [
                "RAPDORid",
                "old_locus_tag",
                "Gene",
                "ProteinFunction",
                "ANOSIM R",
                "Mean Distance",
                "Rank",
                "SVM Rank",
                "SVM Score",
                "RDeeP significant",
                "Control expected shift",
                "RNase expected shift",
                "relative fraction shift",
                "ribosomal protein",
                "contains empty replicate",
            ]
        ]
        all = topRDPM
        all = df[df["old_locus_tag"].isin(all)]
        all.to_csv(output.tsv,sep="\t",index=False)
        all = all[all["ribosomal protein"] == False]
        all.to_csv(output.tsv2,sep="\t",index=False)


rule createTable1andTableS1:
    input:
        tsv = rules.extract_tophits.output.tsv2,
        tsv2 = rules.joinAnalysisMethods.output.file
    output:
        tsv = "Pipeline/Paper/Tables/Table1.tsv",
        xlsx = "Pipeline/Paper/Supplementary/Tables/TableS5.xlsx",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")
        df = df[["old_locus_tag", "Gene", "Rank", "Mean Distance"]]
        df = df.rename({"old_locus_tag": "locus tag"}, axis=1)
        df = df.sort_values(by="Rank")
        df.loc[df.Gene.str.contains("sll|slr", na=False), "Gene"] = ""
        df.to_csv(output.tsv, sep="\t", index=False, float_format='%.2f')
        df = pd.read_csv(input.tsv2, sep="\t")

        df = df[[
                "RAPDORid",
                "old_locus_tag",
                "Gene",
                "ProteinFunction",
                "ANOSIM R",
                "Mean Distance",
                "Rank",
                "SVM Rank",
                "SVM Score",
                'SVM RNA-binding',
                "RDeeP significant",
                "Control expected shift",
                "RNase expected shift",
                "relative fraction shift",
                "ribosomal protein",
                "contains empty replicate",
            ]]
        df = df.sort_values(by="Rank")
        df.to_excel(output.xlsx, index=False)



rule plotVennDiagramm:
        input:
            tsv = rules.joinAnalysisMethods.output.file
        output:
            html = "Pipeline/plots/VennDiagramm_set{set}.html",
            svg = "Pipeline/plots/VennDiagramm_set{set}.svg",
            tsv = "Pipeline/postProcessing/topCandidates_set{set}.tsv",
            json = "Pipeline/plots/VennDiagramm_set{set}.json",
        run:
            from RAPDOR.plots import COLOR_SCHEMES, _color_to_calpha
            import pandas as pd
            import numpy as np
            from pyfunctions.plotlyVenn import venn_to_plotly
            result_df = pd.read_csv(input.tsv, sep="\t")
            if wildcards.set == "others":
                result_df = result_df[~result_df["ribosomal protein"]]
            elif wildcards.set == "ribosomes":
                result_df = result_df[result_df["ribosomal protein"]]
            elif wildcards.set == "rna_binding":
                result_df = result_df[~result_df["ribosomal protein"]]
                result_df = result_df[result_df["GO RNA-binding"] == True]

            top200SVM = set(
                result_df[result_df["SVM RNA-binding"] == True]["old_locus_tag"]
            )
            topRDPM = set(result_df[result_df["ANOSIM R"] >= config["anosimCutoff"]]["old_locus_tag"])

            topRDeep = set(
                result_df[(result_df["rnase_peak_p_value"] <= 0.05) | (result_df["ctrl_peak_p_value"] <= 0.05)]["old_locus_tag"]
            )

            all = topRDeep.union(topRDPM).union(top200SVM)
            all = result_df[result_df["old_locus_tag"].isin(all)]
            all = all[["RAPDORid", "old_locus_tag", "Gene",  "Mean Distance", "relative fraction shift",  "ANOSIM R", "RDeeP significant", "SVM RNA-binding", "Rank", "ribosomal protein"]]
            all.to_csv(output.tsv, sep="\t", index=False)
            colors = list(COLOR_SCHEMES["Dolphin"]) + [COLOR_SCHEMES["Viking"][0]]
            colors = [_color_to_calpha(col, alpha=0.6) for col in colors]
            fig = venn_to_plotly(L_sets=(topRDPM, top200SVM, topRDeep), L_labels=("RAPDOR", "TripepSVM", "RDeeP"), L_color=colors)
            fig.update_layout(template=DEFAULT_TEMPLATE)
            fig.update_layout(
                xaxis=dict(showgrid=False,zeroline=False, showline=False),
                yaxis=dict(showgrid=False,zeroline=False, showline=False)
            )
            fig.add_annotation(
                text=f"n: {len(result_df)}",
                yanchor="top",
                xanchor="right",
                xref="x domain",
                yref="y domain",
                x=1,
                y=1,
                showarrow=False,
            )


            fig.write_html(output.html)
            fig.write_image(output.svg)
            fig.write_json(output.json)

rule prepareforEnrichment:
    input:
        tsv = rules.joinAnalysisMethods.output.file
    output:
        file = "Pipeline/GO/rankedTable.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")
        df = df[~df["old_locus_tag"].isna()]
        df = df[["old_locus_tag", "Rank"]]
        df = df.drop_duplicates("old_locus_tag")
        df.to_csv(output.file, sep="\t", index=False)

rule GOEnrichment:
    input:
        annotation_db = rules.generateOrgDB.output.annotation_db,
        ranked_file = rules.prepareforEnrichment.output.file
    output:
        enriched = "Pipeline/GO/GOEnrichment.tsv"
    script:
        "../Rscripts/goEnrichment.R"



rule rplaDistribution:
    input:
        file = rules.run_RAPDOR.output.json,
        ctrl_img = "Data/Westernblot_ctrl.png",
        rnase_img = "Data/Westernblot_rnase.png",
    output:
        svg="Pipeline/Paper/Subfigures/rplAFigure.svg",
        tsv="Pipeline/Paper/SourceData/F3A.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        from PIL import Image
        from RAPDOR.plots import plot_distribution, plot_barcode_plot, COLOR_SCHEMES
        from plotly.subplots import make_subplots
        import plotly.graph_objs as go

        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        rpla = data.df[data.df["Gene"] == "rplA"]["RAPDORid"]
        subdata = data.norm_array[data[rpla][0]]
        row_heights = [0.8, 0.1, 0.1]
        fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.01, row_heights=row_heights)
        idesign = data.internal_design_matrix[data.internal_design_matrix["Replicate"] == "BR1"]
        fig1 = plot_distribution(subdata, data.internal_design_matrix, offset=data.state.kernel_size // 2, colors=COLOR_SCHEMES["Dolphin"], show_outliers=False)
        source_data = pd.DataFrame({"Treatment": ["Control"] * 3 + ["RNase"] * 3, "Replicate": list(range(1, 4)) * 2 })
        source_data.loc[:, [f"rel. Intensity Fraction {x+2}" for x in range(subdata.shape[1])]] = subdata
        source_data.to_csv(output.tsv, sep="\t", index=False)
        #fig2 = plot_barcode_plot(subdata2, idesign, colors=COLOR_SCHEMES["Dolphin"], fractions=data.fractions, scale_max=False)
        fig.update_yaxes(fig1.layout.yaxis, row=1)
        for trace in fig1["data"]:
            t = "Control" if trace.legend == "legend" else "RNase"
            trace.update(
                line=dict(width=2),
                marker=dict(size=3),
                showlegend=True,
                legend="legend",
                legendgroup=t,
                legendgrouptitle=dict(text=t)
            )
            fig.add_trace(trace, row=1, col=1)
        # for i_idx, trace in enumerate(fig2["data"]):
        #     v = row_heights[1] + row_heights[2]
        #
        #     trace.colorbar.update(
        #         len=v,
        #         yref="paper",
        #         y=row_heights[0] - v * (i_idx // 2),
        #         yanchor="top",
        #         nticks=5,
        #         x=1. + 0.05 if i_idx % 2 else 1.,
        #         showticklabels=True,
        #         thickness=0.05,
        #         thicknessmode="fraction",
        #         ticklabelposition="inside",
        #         title=dict(text="iBAQ BR1" if i_idx == 0 else "-", font=dict(size=12, color="rgba(0,0,0,0)" if i_idx else None)),
        #         tickfont=config["fonts"]["legend"],
        #         outlinewidth=1,
        #         outlinecolor="black"
        #
        #     )
        #     fig.add_trace(trace,row=2 + i_idx, col=1)

        for row in [2, 3]:
            fig.update_yaxes(range=[0, 1], row=row)
            cat = ["Control" if row == 2 else "RNase"]
            fig.update_yaxes(row=row,showgrid=False, autorange="reversed", categoryorder='total ascending', categoryarray=cat)

            fig.add_trace(
                go.Bar(
                    x=[20],
                    y = cat,
                    marker=dict(color= "rgba(0,0,0,0)"),
                    orientation="h",
                    showlegend=False
                ),
                row=row, col=1
            )
        for row in range(2, 4):
            fig.update_yaxes(row=row,showgrid=False)
            fig.update_xaxes(row=row,showgrid=False)


        ctrl_western = Image.open(input.ctrl_img)
        rnase_western = Image.open(input.rnase_img)
        fig.add_layout_image(dict(
            source=ctrl_western,
            x=0,
            y=0.5,
            xref="x2 domain",
            yref="y2 domain",
            xanchor="left",
            yanchor="middle",

        )
        )
        fig.add_layout_image(dict(
            source=rnase_western,
            x=0,
            y=0.5,
            xref="x3 domain",
            yref="y3 domain",
            xanchor="left",
            yanchor="middle",
        )
        )
        fig.update_layout_images(dict(
            sizex=1,
            sizey=1,
            sizing="fill"

        ))
        fig.update_xaxes(dtick=1)
        fig.update_layout(template=DEFAULT_TEMPLATE, height=config["F3"]["A"]["height"], width=config["width"] // 2)
        fig.update_yaxes(row=4, showgrid=False)

        fig.update_yaxes( row=5, showgrid=False)
        fig.update_xaxes(title = "Fraction", row=3, range=[0.5, 20.5])
        fig.update_layout(legend=dict(tracegroupgap=0, orientation="v", x=1.005, y=1, xanchor="left", yanchor="top"))
        fig.add_annotation(
            text="Anti-RplA",
            xref="paper",
            yref="paper",
            x=1.005,
            y=row_heights[-1],
            yanchor="middle",
            xanchor="left",
            showarrow=False
        )
        fig.write_image(output.svg)


rule plotTopHitDistributions:
    input:
        tophits = rules.extract_tophits.output.tsv2,
        json = rules.run_RAPDOR.output.json
    output:
        svg = "Pipeline/Paper/Subfigures/Distribution{distribution}.svg",
        source_data = "Pipeline/Paper/SourceData/Distributions{distribution}.tsv"
    run:
        import pandas as pd
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        import math
        from RAPDOR.plots import plot_protein_distributions
        dist_config = config["distributions"][wildcards.distribution]
        with open(input.json) as handle:
            rapdor_data = RAPDORData.from_json(handle.read())
        df = rapdor_data.df
        topids = df[df["old_locus_tag"].isin(dist_config["locus_tags"] )]
        topids = topids.sort_values(["Rank"])
        topids = topids["RAPDORid"]

        source_data = rapdor_data.df.loc[rapdor_data.df["RAPDORid"].isin(topids), ["Gene"]]
        indices = source_data.index

        source_data = source_data.loc[source_data.index.repeat(6)].reset_index(drop=True)
        source_data["type"] = (["Control"] * 3 + ["RNase"] * 3) * (source_data.shape[0] // 6)

        subdata = rapdor_data.norm_array[indices]
        rsubdata = subdata.reshape(source_data.shape[0],-1)

        source_data.loc[:, [f"rel. Intensity Fraction {x + 2}" for x in range(rsubdata.shape[1])]] = rsubdata
        source_data.to_csv(output.source_data,sep="\t",index=False)

        fig = plot_protein_distributions(
            topids, rapdordata=rapdor_data, colors=COLOR_SCHEMES["Dolphin"],
            plot_type="zoomed", column_widths=[0.7, 0.3], horizontal_spacing=0.075, title_col="Protein name", vertical_spacing=dist_config["vertical_spacing"]
        )
        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=dist_config["height"])
        fig.update_layout(
            legend=dict(font=config["fonts"]["legend"], y=1.015),
            legend2=dict(font=config["fonts"]["legend"], y=dist_config["legend2_y"]),
        )
        fig.update_layout(margin=dict(l=40, r=30, b=37))
        fig.update_annotations(font=config["fonts"]["annotations"])
        for annotation in fig.layout.annotations:
            if annotation.text not in ("rel. Protein Intensities", "Zoom to strongest shift", "Fraction"):
                annotation.update(x=1.1, xanchor="right")
            elif annotation.text == "Fraction":
                annotation.update(yshift=-15)
            elif annotation.text == "rel. Protein Intensities":
                annotation.update(x=.03)
            if annotation.text == "RaiA/LrtA":
                annotation.text = "LrtA"
            if annotation.text == "Sll1371":
                annotation.text = "SyCrp1"
        fig.update_layout(
            legend=dict(title=dict(font=dict(size=config["fonts"]["default"]["size"]))),
            legend2=dict(title=dict(font=dict(size=config["fonts"]["default"]["size"]))),
        )
        fig.update_traces(line=dict(width=2),
                marker=dict(size=3),)
        fig.update_yaxes(nticks=2, col=2)
        fig.update_yaxes(nticks=3, col=1)
        fig.write_image(output.svg)

rule plotMeanDistribution:
    input:
        file= rules.run_RAPDOR.output.json,
    output:
        svg="Pipeline/Paper/Subfigures/MeanDistributions.svg",
        tsv="Pipeline/Paper/SourceData/F3B.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        import math
        from RAPDOR.plots import multi_means_and_histo, COLOR_SCHEMES, DEFAULT_TEMPLATE, plot_distance_histo

        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        ids = data.df[data.df["small_ribosomal"] == True]["RAPDORid"]
        ids2 = data.df[data.df["large_ribosomal"] == True]["RAPDORid"]
        ids3 = data.df[data.df["photosystem"] == True]["RAPDORid"]
        ids4 = data.df[data.df["RNA ploymerase"] == True]["RAPDORid"]
        joined_ids = set(pd.concat([ids, ids2, ids3, ids4]).tolist())
        source_data = data.df.loc[data.df["RAPDORid"].isin(joined_ids), ["Gene", "Mean Distance", "ANOSIM R", "small_ribosomal", "large_ribosomal", "photosystem"]]
        indices = source_data.index
        source_data = source_data.loc[source_data.index.repeat(2)].reset_index(drop=True)
        source_data["type"] = ["Control", "RNase"] * (len(source_data) // 2)
        means = data._treatment_means[:, indices]
        means2 = means.swapaxes(0, 1).reshape(150,-1)
        source_data.loc[:, [f"rel. Intensity Fraction {x+2}" for x in range(means2.shape[1])]] = means2
        source_data.to_csv(output.tsv, sep="\t", index=False)
        d = {"small subunit": ids, "large subunit": ids2, "photosystem": ids3}
        fig = multi_means_and_histo(d, data, colors=COLOR_SCHEMES["Dolphin"] + COLOR_SCHEMES["Viking"], row_heights=[0.2, 0.2, 0.6], vertical_spacing=0.08)
        fig.update_yaxes(row=3, range=[0, 0.6])
        fig.update_yaxes(row=1, dtick=0.4)
        fig.update_traces(
            marker_color="lightgrey",
            marker_line=dict(width=0.5,color='black'),
            row=1
        )
        fig.update_traces(
            marker_color="lightgrey",
            marker_line=dict(width=0.5,color='black'),
            row=2
        )

        fig.update_layout(
            height=config["F3"]["B"]["height"],
            width=config["width"] // 2,
            template=DEFAULT_TEMPLATE,
            margin=dict(
                b=30,
                t=20,
                r=30,
                l=30
            ),
        )
        fig.update_annotations(
            dict(font=config["fonts"]["axis"])
        )
        fig.write_image(output.svg)


rule plotBarcodePlot:
    input:
        file = rules.joinAnalysisMethods.output.file,
        rapdor = rules.run_RAPDOR.output.json
    output:
        svg = "Pipeline/Paper/Subfigures/ribosomalIdentifier.svg",
        html = "Pipeline/Paper/Subfigures/ribosomalIdentifier.html",
        tsv="Pipeline/Paper/SourceData/F3C.tsv"
    run:
        import plotly.graph_objs as go
        import pandas as pd
        import numpy as np
        from RAPDOR.plots import rank_plot
        from RAPDOR.datastructures import RAPDORData

        with open(input.rapdor) as handle:
            data = RAPDORData.from_json(handle.read())


        df = pd.read_csv(input.file, sep="\t")
        #df = df[df["ribosomal protein"] == True]
        df = df.sort_values(by="Rank", ascending=True)
        df = df[~df["ANOSIM R"].isna()]
        x = list(range(int(df["Rank"].min()), int(df.Rank.max())+1))
        small_subunit = df[df["small_ribosomal"] == True]["RAPDORid"]
        large_subunit = df[df["large_ribosomal"] == True]["RAPDORid"]
        photosystem = df[df["photosystem"] == True]["RAPDORid"]

        joined_ids = set(pd.concat([small_subunit, large_subunit, photosystem]).tolist())
        source_data = data.df.loc[
            data.df["RAPDORid"].isin(joined_ids), ["Gene", "Rank", "small_ribosomal",
                                                   "large_ribosomal", "photosystem"]]
        source_data.to_csv(output.tsv,sep="\t",index=False)

        colors = (COLOR_SCHEMES["Dolphin"][0], COLOR_SCHEMES["Dolphin"][1], COLOR_SCHEMES["Viking"][0])

        d = {"small subunit": small_subunit, "large subunit": large_subunit, "photosystem": photosystem}
        fig = rank_plot(d, data, colors, orientation="h", triangles="inside", tri_x=17, tri_y=0.2)
        fig.update_shapes(showlegend=False)
        fig.layout.annotations = None

        fig.update_layout(
            width=config["width"] // 2, height=config["F3"]["C"]["height"],
            margin=dict(
                b=30,
                t=20,
                r=30,
                l=30
            ),
            legend=dict(bgcolor = 'rgba(0,0,0,0)', )

        )





        fig.write_image(output.svg)
        fig.write_html(output.html)

rule plotAllVenns:
    input:
        files = expand(rules.plotVennDiagramm.output.json, set=["others", "ribosomes", "rna_binding"])
    output:
        svg = "Pipeline/Paper/Figure4.svg"
    run:
        import plotly.graph_objs as go
        import plotly.subplots as sp
        import plotly.io as pio
        import string
        mapping = {input.files[1]: "ribosomal proteins", input.files[0]: "non ribosomal proteins", input.files[2]: "GO RNA-binding"}
        multi_fig = sp.make_subplots(rows=1,cols=3)
        annos = string.ascii_uppercase
        for idx, file in enumerate(input.files, 1):
            fig = pio.read_json(file)
            for iidx, shape in enumerate(fig['layout']['shapes']):
                shape["xref"] = f"x{idx}" if idx > 1 else "x"
                shape["yref"] = f"y{idx}" if idx > 1 else "y"
                multi_fig.add_shape(shape)
                color = shape["fillcolor"]
                if idx == 1 and color != "rgba(255,0,0,0)":
                    multi_fig.add_trace(
                        go.Scatter(
                            x=[None],
                            y=[None],
                            mode="markers",
                            name=fig['layout']['annotations'][iidx]["text"],
                            marker=dict(
                                color=color,
                                size=900,
                                # line=dict(
                                #     color="black",
                                #     width=3
                                # ),
                            ),
                            showlegend=True

                        )
                    )
            for annotation in fig['layout']['annotations'][3:]:
                if annotation["xref"] == "x":
                    annotation["xref"] = f"x{idx}" if idx > 1 else "x"
                    annotation["yref"] = f"y{idx}" if idx > 1 else "y"
                else:
                    annotation["xref"] = f"x{idx} domain" if idx > 1 else "x domain"
                    annotation["yref"] = f"y{idx} domain" if idx > 1 else "y domain"
                multi_fig.add_annotation(annotation)

            multi_fig.add_annotation(
                text=f"{mapping[file]}",
                xref=f"x{idx} domain" if idx > 1 else "x domain",
                xanchor="center",
                x=0.5,
                yref=f"y{idx} domain" if idx > 1 else "y domain",
                yanchor="bottom",
                y=0,
                font=dict(size=12),
                showarrow=False

            )

            multi_fig.update_xaxes(fig["layout"]["xaxis"], row=1, col=idx)
            multi_fig.update_yaxes(fig["layout"]["yaxis"], row=1, col=idx)
            multi_fig.update_yaxes(scaleanchor="x" if idx == 1 else f"x{idx}", scaleratio=1, col=idx)

        multi_fig.update_layout(template=DEFAULT_TEMPLATE)
        multi_fig.update_layout(
            legend={'itemsizing': 'trace', "orientation": "h", "yanchor": "bottom", "y": 1.01},
            width=config["width"], height=250, font=config["fonts"]["legend"]
        )
        for annotation in multi_fig.layout.annotations:
            if annotation.text.startswith("<b>"):
                pass
            elif annotation.text.startswith("n:"):
                annotation.update(font=config["fonts"]["annotations"])
            else:
                annotation.update(font=config["fonts"]["annotations"], xanchor="center", yanchor="middle")
                if annotation.text == "0":
                    if annotation.xref == "x":
                        annotation.update(showarrow=True, ay=-0.5, ax=0.5, axref=annotation.xref, ayref=annotation.yref, arrowcolor='black')
                if annotation.text == "1":
                    if annotation.xref == "x":
                        annotation.update(showarrow=True, ay=-0.3, ax=0.7, axref=annotation.xref, ayref=annotation.yref, arrowcolor='black')
                if annotation.text == "4":
                    annotation.update(x=annotation.x - 0.1, y=annotation.y - 0.1)
                if annotation.text == "8":
                    annotation.update(x=annotation.x + 0.05)
                if annotation.text == "12":
                    annotation.update(y=annotation.y + 0.05)


        multi_fig.write_image(output.svg)


rule plotRuntime:
    input:
        all_benchmark = rules.meanRuntimeAndRAM.output.file
    output:
        html = "Pipeline/Paper/Supplementary/Figures/html/benchmark.html",
        svg = "Pipeline/Paper/Supplementary/Figures/FigureS7.svg"
    run:
        import plotly.graph_objs as go
        import plotly.express as px
        import pandas as pd
        import plotly.subplots as sp
        import string

        annos = string.ascii_uppercase
        df = pd.read_csv(input.all_benchmark, sep="\t")
        datasets = ["Synechocystis GradR", "HeLa R-DeeP"]
        fig = sp.make_subplots(rows=2,cols=2,
            x_title="Replicates per group",
            horizontal_spacing=0.025,
            vertical_spacing=0.05,
            column_titles=datasets,
            shared_yaxes=True,
        )
        names = ["s", "max_uss"]
        df["Replicates per group"] = df["Replicates per group"].astype(str)
        for col in range(1, 3):
            ds = datasets[col-1]
            sdf = df[df["dataset"] == ds]
            for row in range(1, 3):
                name = names[row-1]

                fig.add_trace(
                    go.Bar(
                        x=sdf[(sdf["Tool"] == "R-DeeP")]["Replicates per group"],
                        y=sdf[(sdf["Tool"] == "R-DeeP")][name],
                        error_y=dict(
                            type='data',# value of error bar given in data coordinates
                            array=sdf[(sdf["Tool"] == "R-DeeP")][f"stdev_{name}"],
                            visible=True,
                            color="black"
                        ),
                        legendgroup="R-DeeP",
                        legendgrouptitle=dict(text="R-DeeP"),
                        name="R-DeeP",
                        marker_color=COLOR_SCHEMES["Viking"][0],
                        marker_line=dict(width=2,color='black'),
                        showlegend=True if row == 1 and col == 1 else False,

                    ),
                    col=col, row=row
                )


                fig.add_trace(
                    go.Bar(
                        x=sdf[(sdf["ANOSIM"] == False) & (sdf["Tool"] == "RAPDOR")]["Replicates per group"],
                        y=sdf[(sdf["ANOSIM"] == False) & (sdf["Tool"] == "RAPDOR")][name],
                        error_y=dict(
                            type='data',# value of error bar given in data coordinates
                            array=sdf[(sdf["ANOSIM"] == False) & (df["Tool"] == "RAPDOR")][f"stdev_{name}"],
                            visible=True,
                            color="black"
                        ),

                        legendgroup="RAPDOR",
                        legendgrouptitle=dict(text="RAPDOR"),
                        marker_color=COLOR_SCHEMES["Dolphin"][0],
                        marker_line=dict(width=2,color='black'),
                        showlegend=True if row == 1 and col == 1 else False,


                        name="Ranking"
                    ),
                    col=col, row=row
                )
                fig.add_trace(
                    go.Bar(
                        x=sdf[(sdf["ANOSIM"] == True) & (sdf["Tool"] == "RAPDOR")]["Replicates per group"],
                        y=sdf[(sdf["ANOSIM"] == True) & (sdf["Tool"] == "RAPDOR")][name],
                        error_y=dict(
                            type='data',# value of error bar given in data coordinates
                            array=sdf[(sdf["ANOSIM"] == True) & (sdf["Tool"] == "RAPDOR")][f"stdev_{name}"],
                            visible=True,
                            color="black"
                        ),
                        legendgroup="RAPDOR",
                        legendgrouptitle=dict(text="RAPDOR"),
                        marker_color=COLOR_SCHEMES["Dolphin"][1],
                        marker_line=dict(width=2,color='black'),
                        showlegend=True if row == 1 and col == 1 else False,


                        name="ANOSIM"
                    ),
                    col=col,row=row

                )
                # fig.add_annotation(
                #     text=f"<b>{annos[col - 1]}</b>",
                #     xref=f"x{col} domain" if col > 1 else "x domain",
                #     xanchor="left",
                #     x=-0.4,
                #     yref=f"y{col} domain" if col > 1 else "y domain",
                #     yanchor="top",
                #     y=1.1,
                #     font=dict(size=config["multipanel_font_size"]),
                #     showarrow=False
                #
                # )

        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=300,
            margin=config["margin"]

        )
        fig.update_yaxes(type="log", row=1)
        fig.update_yaxes(title=dict(text="Average Runtime [s]"), row=1, col=1)
        fig.update_yaxes(title=dict(text="Memory [Mb]"), row=2, col=1)
        fig.update_annotations(
            dict(font=config["fonts"]["axis"])
        )

        #fig.update_xaxes(type='linear')
        #fig.update_xaxes(col=2, range=(3-4.5, 3+4.5))
        fig.update_layout(bargap=0.2,bargroupgap=0.1)
        fig.update_xaxes(col=2,range=(-1.5, 1.5))
        fig.update_xaxes(showticklabels=False, row=1)
        fig.write_image(output.svg)
        fig.write_html(output.html)


rule createBubblePlot:
    input:
        file=rules.run_RAPDOR.output.json,
        joined=rules.extract_tophits.output.tsv2,
    output:
        svg = "Pipeline/Paper/Subfigures/DimensionReduction_ids{highlight}.svg",
        html = "Pipeline/plots/DimensionReduction_ids{highlight}.html",
        json = "Pipeline/plots/DimensionReduction_ids{highlight}.json",
        source_data = "Pipeline/Paper/SourceData/BubblePlot_{highlight}.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        import math
        from RAPDOR.plots import plot_dimension_reduction
        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        data.calc_distribution_features()
        if wildcards.highlight == "overlapping":
            overlapping = pd.read_csv(input.joined, sep="\t")
            ids = data.df[data.df["old_locus_tag"].isin(overlapping["old_locus_tag"])]["RAPDORid"]
        elif wildcards.highlight == "topHits":
            df = data.df
            ids = df[df["old_locus_tag"].isin(config["bubble_plot"]["locus_tag"])]["RAPDORid"]
        else:
            raise ValueError("Not supported")
        embedding = data.current_embedding
        rel_dist_change = embedding[:, 1]
        data.df['Gene'] = data.df.apply(lambda row: row['old_locus_tag'] if pd.isnull(row['Gene']) or row['Gene'] == '' else row['Gene'], axis=1)
        source_data = data.df.loc[:, ["Gene", "Mean Distance", "relative fraction shift",]]
        source_data["relative distribution change"] = rel_dist_change
        source_data.to_csv(output.source_data, sep="\t", index=False)
        fig = plot_dimension_reduction(
            rapdordata=data,
            colors=COLOR_SCHEMES["Dolphin"],
            highlight=ids,
            title_col="Protein name",
            legend_spread=0.125,
            legend_start=0.175
        )

        fig.layout.update(
            template=DEFAULT_TEMPLATE,
            width=config["width"],
            height=config["height"],
            legend=dict(visible=False),
            margin=config["margin"]
        )
        fig.update_xaxes(dtick=5)
        fig.update_annotations(
            font=config["fonts"]["legend"],
        )
        for annotation in fig.layout.annotations:
            if annotation.text == "0.9":
                annotation.update(x=annotation.x + 0.01)
            elif annotation.text == "PheS":
                annotation.update(showarrow=True,ay=-.1,ax=-2,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "GroEL1":
                annotation.update(showarrow=True,ay=-.3,ax=-10,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "CphA":
                annotation.update(showarrow=True,ax=-17, ay=0.7, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Tsf":
                annotation.update(showarrow=True,ay=.1,ax=annotation.x,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "FtsH2":
                annotation.update(showarrow=True,ay=-.1,ax=-8,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "CysC":
                annotation.update(showarrow=True,ay=.075,ax=1,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RpoB":
                annotation.update(showarrow=True,ay=annotation.y,ax=0, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RbpA":
                annotation.update(showarrow=True,ay=1,ax=-7.5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0147":
                annotation.update(showarrow=True,ay=1.25,ax=-6, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr1143":
                annotation.update(showarrow=True,ay=1.5,ax=-4, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RaiA/LrtA":
                annotation.update(showarrow=True,ay=-1.25,ax=-19, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1388":
                annotation.update(showarrow=True,ax=-11,ay=.9, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1961":
                annotation.update(showarrow=True,ax=-18.5 ,ay=.75, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll0921":
                annotation.update(showarrow=True,ax=-18,ay=-1.6, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Pgm":
                annotation.update(showarrow=True,ax=0, ay=1.1, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1967":
                annotation.update(showarrow=True,ax=-5,ay=-2, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0782":
                annotation.update(showarrow=True,ax=0.5,ay=-.7, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "QueF":
                annotation.update(showarrow=True,ax=-8.5,ay=-1, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "SecE":
                annotation.update(showarrow=True,ax=-15,ay=1,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Ffh":
                annotation.update(showarrow=True,ax=14,ay=1.5,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll0284":
                annotation.update(showarrow=True,ax=5,ay=-1.25,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "AroB":
                annotation.update(showarrow=True,ax=2,ay=.75,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0678":
                annotation.update(showarrow=True,ax=0,ay=-1.43,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll7087":
                annotation.update(showarrow=True,ax=-15,ay=-2,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RpoC1":
                annotation.update(showarrow=True,ax=-7,ay=-1.25,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "ChlI":
                annotation.update(showarrow=True,ax=-19,ay=-.4,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Ssr3189":
                annotation.update(showarrow=True,ax=-11,ay=-1.5,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Ssl2245":
                annotation.update(showarrow=True,ax=-1.5,ay=1.27,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1371":
                annotation.update(text="SyCrp1", showarrow=True,ax=-3,ay=-1.75,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1315":
                annotation.update(showarrow=True, ax=-13,  ay=1.25,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0670":
                annotation.update(showarrow=True,ax=-8.5 ,ay=.75, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "BioU":
                    annotation.update(showarrow=True,ax=7 ,ay=-1.5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            #fig.update_xaxes(range=[-6.5, -1], row=2)
        #fig.update_yaxes(range=[-.75, .75], row=2)
        fig.write_image(output.svg)
        fig.write_html(output.html)
        fig.write_json(output.json)


rule createFigure5:
    input:
        bubble = expand(rules.createBubblePlot.output.svg, highlight="topHits")
    output:
        svg = "Pipeline/Paper/Figure6.svg"
    shell:
        """
        cp {input.bubble[0]} {output.svg}
        """

rule cpRedistSynecho:
    input:
        json = rules.runOnSynData.output.json
    output:
        json = "Pipeline/Paper/Supplementary/JSON/RAPDORforSynechocystisabioticStress{condition}.json"
    shell:
        "cp {input.json} {output.json}"



rule plotConditionedSynechochoColdRibo:
    input:
        json = expand(rules.runOnSynData.output.json, condition="COLD")
    output:
        svg = "Pipeline/Paper/unused/FigureS2.svg",
        svg2 = "Pipeline/Paper/unused/Figures/FigureS6.svg",
        html2 = "Pipeline/Paper/unused/Figures/html/FigureS6.html",
        tsv = "Pipeline/Paper/unused/Tables/TableS6.xlsx",
        #json = "Pipeline/Paper/Supplementary/JSON/RAPDORforSynechocystisColdShock.json",
    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import multi_means_and_histo, rank_plot
        file = input.json[0]
        with open(file) as handle:
            data = RAPDORData.from_json(handle.read())
        r_to_membrane = data.df[(data.df["Gene"].str.contains("rps|rpl")) & (data.df["position strongest shift"] == "Membrane")]["RAPDORid"]
        names = [f"rpl{id}" for id in [1, 13, 18, 22, 23, 28, 29, 33,]]
        names += [f"rps{id}" for id in [10, 20, 7, 8]]
        identified_inpub = data.df[data.df["Gene"].isin(names)]["RAPDORid"]
        d = {"ribosomal & membrane": r_to_membrane}
        fig = multi_means_and_histo(d,data,colors=COLOR_SCHEMES["Dolphin"] + COLOR_SCHEMES["Viking"])
        fig.update_traces(
            marker_color="lightgrey",
            marker_line=dict(width=0.5,color='black'),
            row=1
        )
        fig.update_traces(
            marker_color="lightgrey",
            marker_line=dict(width=0.5,color='black'),
            row=2
        )
        fig.write_image(output.svg)
        fig = rank_plot(d, rapdordata=data, colors = COLOR_SCHEMES["Dolphin"] + COLOR_SCHEMES["Viking"])
        fig.layout.shapes = None
        fig.layout.annotations = None
        fig.update_layout(width=config["width"])
        fig.update_layout(height=config["height"] / 2)

        fig.write_image(output.svg2)
        fig.write_html(output.html2)
        data.df["mentioned in wang"] = data.df["Gene"].isin(names)
        data.export_csv(output.tsv, sep="\t")
        data: RAPDORData
        df = data.extra_df.drop(["id"],axis=1)
        df.to_excel(output.tsv, index=False)


rule analyzeProteinAbundance:
    input:
        json=rules.run_RAPDOR.output.json,
    output:
        svg = "Pipeline/Paper/Subfigures/AbundanceDistribution.svg",

    run:
        from RAPDOR.datastructures import RAPDORData
        from scipy.stats import pearsonr
        import numpy as np
        import plotly.graph_objs as go
        from sklearn.decomposition import PCA

        file = input.json
        with open(file) as handle:
            data = RAPDORData.from_json(handle.read())
        array = data.norm_array[:-1]
        #array = np.swapaxes(array, 1, 2)
        cov = np.cov(array[0], rowvar=False)
        covariance_matrices = np.empty((array.shape[0], array.shape[-1], array.shape[-1]))
        for i in range(array.shape[0]):
            # Get the data for the current experiment
            experiment_data = array[i]

            # Calculate the covariance matrix for the current experiment
            covariance_matrix = np.cov(experiment_data,rowvar=False)

            # Store the covariance matrix in the array
            covariance_matrices[i] = covariance_matrix
        print(cov.shape)
        print(covariance_matrices.shape)
        print(covariance_matrices[0] == cov)
        dets = np.linalg.det(covariance_matrices)
        print(np.nanargmax(dets))
        print(dets[np.nanargmax(dets)])
        sorted_indices = np.argsort(dets)[::-1]

        # Filter out NaN values
        non_nan_indices = sorted_indices[~np.isnan(dets[sorted_indices])]

        # Get the top 20 indices
        ind = non_nan_indices[:100]
        print(ind)
        print(dets[ind])
        print(data.df.loc[ind][["RAPDORid", "Rank", "Gene"]])
        fig = go.Figure()
        print(array.shape)
        array = array[ind]
        shape = array.shape
        array = array.reshape(-1, array.shape[0])

        pca = PCA(n_components=2)
        new_x = pca.fit_transform(array)
        new_x = new_x.reshape(shape[1], shape[2], -1)
        print(new_x.shape)
        print(pca.explained_variance_ratio_)
        for fraction in range(new_x.shape[1]):
            #name = f"{row['Treatment']}-{row['Replicate']}"
            for p in range(2):
                i = "control" if p == 0 else "RNase"
                name = f"fraction{fraction}- {i}"

                embedding = new_x[p*3:(1+p)*3, fraction, :]
                fig.add_trace(
                    go.Scatter(
                        x=embedding[:, 0],
                        y=embedding[:, 1],
                        #z=[new_x[idx, 2]],
                        name=name,
                        marker=dict(size=6),
                        mode="markers",

                    )
                )
        fig.show()
        exit()

rule plotSpearmans:
    input:
        json=rules.run_RAPDOR.output.json,
    output:
        svg = "Pipeline/Paper/Subfigures/SampleSpearman.svg",
        svg2 = "Pipeline/Paper/Subfigures/SampleJSD.svg",

    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import plot_sample_histogram
        import numpy as np
        import plotly.graph_objs as go
        from scipy.stats import spearmanr
        from plotly.subplots import make_subplots

        file = input.json
        with open(file) as handle:
            data = RAPDORData.from_json(handle.read())
        fig = plot_sample_histogram(data,method="spearman",  colors=DEFAULT_COLORS)
        fig.update_layout(template=DEFAULT_TEMPLATE, margin=dict(b=70, l=70, r=70, t=70))
        fig.update_annotations(font=config["fonts"]["axis"])

        fig.write_image(output.svg)
        fig.update_traces(
            marker_line=dict(color="black", width=0.1)
        )
        fig = plot_sample_histogram(data, method="jsd", colors=DEFAULT_COLORS)
        fig.update_traces(
            marker_line=dict(color="black", width=0.1)
        )
        fig.update_layout(template=DEFAULT_TEMPLATE, margin=dict(b=70, l=70, r=70, t=70))
        fig.update_annotations(font=config["fonts"]["axis"])

        fig.write_image(output.svg2)


                # Store the coefficient in the matrix

        # names = [f"{row['Treatment']} - {row['Replicate']}" for idx, row in data.internal_design_matrix.iterrows()]
        # fig = go.Figure(data=go.Heatmap(
        #     z=pearson_matrix,
        #     x = names,
        #     y=names,
        # ))
        # fig.show()

rule combineFigure3:
    input:
        rpla = rules.rplaDistribution.output.svg,
        means = rules.plotMeanDistribution.output.svg,
        barcode = rules.plotBarcodePlot.output.svg
    output:
        svg = "Pipeline/Paper/Figure3.svg",

    run:
        from svgutils.compose import Figure, Panel, SVG, Text
        b_y = config["F3"]["A"]["height"]
        c_y = config["F3"]["A"]["height"] + config["F3"]["B"]["height"]
        ges_y = config["F3"]["A"]["height"] + config["F3"]["B"]["height"] + config["F3"]["C"]["height"]
        f = Figure("624px",f"{ges_y}px",
            Panel(
                SVG(input.rpla),
                Text("A",2,15,size=config["multipanel_font_size"],weight='bold', font="Arial")

            ),
            Panel(
                SVG(input.means),
                Text("B",2,2,size=config["multipanel_font_size"],weight='bold', font="Arial")
            ).move(0,b_y),
            Panel(
                SVG(input.barcode),
                Text("C",2,2,size=config["multipanel_font_size"],weight='bold', font="Arial")
            ).move(0,c_y)
        )
        svg_string = f.tostr()
        svg_string = svg_string.decode().replace("encoding='ASCII'", "encoding='utf-8'")
        with open(output.svg, "w") as handle:
            handle.write(svg_string)

rule detectedProteinTable:
    input:
        json = rules.run_RAPDOR.output.json
    output:
        file2 = "Pipeline/Paper/Supplementary/Tables/TableS4.xlsx"
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        from scipy.stats import pearsonr
        import pandas as pd
        import plotly.graph_objs as go
        from matplotlib.colors import LinearSegmentedColormap

        file = input.json
        with open(file) as handle:
            data = RAPDORData.from_json(handle.read())
        p = data.array
        mask = np.all(np.isnan(p), axis=-1)
        z = p
        z[mask] = 1
        idx = np.nanargmax(p, axis=-1) + 1
        idx[mask] = -1

        i = np.nansum(p > 0, axis=0).reshape(p.shape[-1], -1)#
        counts = np.empty(i.shape)
        ks = data.state.kernel_size // 2
        ks = 0

        for x in data.fractions:
            count_along_axis = np.count_nonzero(idx == x, axis=0)
            counts[x-1] = count_along_axis


        row_names = data.fractions[ks:-ks ] if ks != 0 else data.fractions
        column_names = [f"{row['Treatment']}-{row['Replicate']}" for _, row in data.internal_design_matrix.iterrows()]
        df = pd.DataFrame(i, columns=column_names, index=row_names)
        df2 = pd.DataFrame(counts, columns=column_names, index=row_names)
        corr = np.empty((df.shape[1], df.shape[1]))
        for idx1, col1 in enumerate(df.columns):
            for idx2, col2 in enumerate(df.columns):
                a1 = df[col1]
                a2 = df[col2]
                pearson_r, _ = pearsonr(a1, a2)
                corr[idx1, idx2] = pearson_r if idx1 != idx2 else np.nan
        fig = go.Figure(data=go.Heatmap(
            z=corr,
            x=df.columns,
            y=df.columns
        ))
        df.index.name = "Fraction"
        #df.to_csv(output.file, sep="\t")
        colors = ["#008e97", "#fc4c02"]

        # Create the custom color map
        cmap = LinearSegmentedColormap.from_list('custom_cmap',colors)


        def apply_gradient(column):
            return column.style.background_gradient(cmap=cmap)


        df2 = df2.style.background_gradient(cmap=cmap, axis=0)
        df.loc['Sum'] = df.sum(axis=0)
        df2.to_excel(output.file2)
        #df2.loc['Sum'] = df2.sum(axis=0)

rule calcANOSIMDistribution:
    input:
        json=rules.run_RAPDOR.output.json,

    output:
        gradr_np="Pipeline/analyzed/ANOSIMDistributionGradR.np",
    threads: 10
    run:
        from RAPDOR.datastructures import RAPDORData
        import multiprocessing

        multiprocessing.set_start_method(method="spawn", force=True)
        gradr_data = input.json
        iter_data = [
            (gradr_data, 3, output.gradr_np),

        ]
        for plot_id, data in enumerate(iter_data):
            data, samples, outfile = data
            data = RAPDORData.from_file(data)
            data: RAPDORData
            print("starting")
            _, distribution = data.calc_anosim_p_value(-1,threads=threads, mode="global")
            print(data.df["global ANOSIM adj p-Value"])
            indices = data.df[data.df["min replicates per group"] == samples]["id"].to_numpy()
            del data
            distribution = distribution[:, indices].flatten()
            with open(outfile, "wb") as handle:
                np.save(handle, distribution)
            del distribution
            print("saved")

rule anosimEGF:
    input:
        json = rules.AnalyzeNatureWithRAPDOR.output.json
    output:
        np = "Pipeline/analyzed/ANOSIMDistribution{experiment}.np",
        json = "Pipeline/analyzed/ANOSIMDistribution{experiment}.json",
    threads: 10
    run:
        from RAPDOR.datastructures import RAPDORData
        import multiprocessing

        multiprocessing.set_start_method(method="spawn", force=True)
        data, samples, outfile = (input.json, 4, output.np)
        data = RAPDORData.from_file(data)
        data: RAPDORData
        print("starting")
        _, distribution = data.calc_anosim_p_value(-1,threads=threads, mode="global")
        print(data.df["global ANOSIM adj p-Value"])
        indices = data.df[data.df["min replicates per group"] == samples]["id"].to_numpy()
        distribution = distribution[:, indices].flatten()
        with open(outfile,"wb") as handle:
            np.save(handle,distribution)
        data._anosim_distribution = None
        del distribution
        data.to_json(output.json)
        print("saved")



rule plotANOSIMRDistribution:
    input:
        np=rules.calcANOSIMDistribution.output,
        json=rules.calcANOSIMDistribution.input,
        npegf = expand(rules.anosimEGF.output.np, experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
        jsonegf = expand(rules.AnalyzeNatureWithRAPDOR.output.json, experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
        rdeppdistribution = rules.runRAPDORonRDeeP.output.RDistribution,
        rdeppjson = rules.runRAPDORonRDeeP.output.json
    output:
        svg = "Pipeline/Paper/Supplementary/Figures/FigureS5.svg",
        source_gradr = "Pipeline/Paper/SourceData/gradr.tsv",
        source_rdeep = "Pipeline/Paper/SourceData/rdeep.tsv",
        source_egf2 = "Pipeline/Paper/SourceData/egf2.tsv",
        source_egf8 = "Pipeline/Paper/SourceData/egf8.tsv",
        source_egf20 = "Pipeline/Paper/SourceData/egf20.tsv",
        source_egf90 = "Pipeline/Paper/SourceData/egf90.tsv",
    threads: 1
    run:
        from RAPDOR.datastructures import RAPDORData
        import plotly.graph_objs as go
        from plotly.subplots import make_subplots
        ip_files = list(input.np) + [input.rdeppdistribution] + list(input.npegf)
        json_files = list(input.json) + [input.rdeppjson] + list(input.jsonegf)
        ip_files = zip(ip_files, json_files)
        row_titles = ["GradR", "RDeeP", "2 min", "8 min", "20 min", "90 min"]
        fig = make_subplots(
            rows=len(json_files),
            shared_xaxes=True,
            x_title="ANOSIM R",
            y_title="Probability",
            row_titles=row_titles,
            vertical_spacing=0.02
        )
        syn_ys = [anno["y"] for anno in fig.layout.annotations[0:2]]
        syn_y = np.mean(syn_ys)
        fig.update_layout(font=config["fonts"]["default"])
        egf_ys = [anno["y"] for anno in fig.layout.annotations[2:6]]
        egf_y = np.mean(egf_ys)
        # fig.add_annotation(
        #     text="GradR",
        #     xref="x2 domain",
        #     yref="paper",
        #     yanchor="middle",
        #     xanchor="left",
        #     textangle=90,
        #     x=1.15,
        #     font=dict(size=16,
        #     ),
        #     y=syn_y
        # )
        x_0 = 1.1
        # fig.add_shape(type="line",
        #     x0=x_0,y0=syn_ys[0] + 0.025,x1=x_0,y1=syn_ys[-1] - 0.025,xref="x2 domain",yref="paper",line_color="black"
        # )
        fig.add_annotation(
            text="HeLa EGF",
            xref="x2 domain",
            yref="paper",
            yanchor="middle",
            xanchor="left",
            textangle=90,
            x=1.15,
            font=dict(size=16),
            y=egf_y
        )
        fig.add_shape(type="line",
            x0=x_0,y0=egf_ys[0],x1=x_0,y1=egf_ys[-1],xref="x2 domain",yref="paper",line_color="black"
        )
        source_data_out = [
            output.source_gradr,
            output.source_rdeep,
            output.source_egf2,
            output.source_egf8,
            output.source_egf20,
            output.source_egf90,
        ]
        for plot_id, (data, json_file) in enumerate(ip_files):

            with open(data, "rb") as handle:
                distribution = np.load(handle)
            rapdor_data = RAPDORData.from_file(json_file)

            source_df = rapdor_data.df
            title = row_titles[plot_id]
            if title == "GradR":
                name = "Gene"
                rep = 3
                treat = "RNase"
                add = 2
            elif title == "RDeeP":
                name = "Protein"
                rep = 3
                treat = "RNase"
                add = 2

            else:
                name = "Gene.names"
                rep = 4
                treat = "EGF"
                add = 1

            source_data = rapdor_data.df.loc[:, [name]]


            source_data = source_data.loc[source_data.index.repeat(2*rep)].reset_index(drop=True)
            source_data["type"] = (["Control"] * rep + [treat] * rep) * (source_data.shape[0] // (2*rep))

            subdata = rapdor_data.norm_array
            rsubdata = subdata.reshape(source_data.shape[0],-1)

            source_data.loc[:, [f"rel. Intensity Fraction {x + add}" for x in range(rsubdata.shape[1])]] = rsubdata
            source_data.to_csv(source_data_out[plot_id],sep="\t",index=False)


            original_dist = rapdor_data.df["ANOSIM R"]

            y, x = np.histogram(
                distribution,
                bins=np.linspace(-1, 1, int(np.floor(2/0.01)))
            )
            y = y / np.sum(y)
            percentile = np.nanquantile(distribution, 0.95)

            fig.add_trace(go.Bar(
                x=x,
                y=y,
                showlegend=False,
                marker=dict(
                    line=dict(color="black", width=1),
                    color=DEFAULT_COLORS[0],
                ),
            ), row=plot_id+1, col=1)
            fig.add_vline(
                x=percentile,
                annotation_text=f"{percentile:.3f}",
                annotation_position="top left",
                row=plot_id + 1,col=1
            )
            y, x = np.histogram(
                original_dist,
                bins=np.linspace(-1,1,int(np.floor(2 / 0.01)))
            )
            y = y / np.sum(y)

            # fig.add_trace(go.Bar(
            #     x=x,
            #     y=y,
            #     showlegend=False,
            #     marker=dict(
            #         line=dict(color="black",width=1),
            #         color=DEFAULT_COLORS[1],
            #         opacity=0.25
            #     ),
            # ),row=plot_id + 1,col=1)

        fig.update_xaxes(range=[-1, 1])
        fig.update_layout(template=DEFAULT_TEMPLATE)
        fig.update_layout(width=config["width"])
        fig.update_layout(height=config["height"]*1.5)
        fig.update_layout(margin=dict(
            b=50,
            t=20,
            r=70,
            l=70,
        ))
        fig.update_layout(barmode="overlay", bargap=0)
        fig.write_image(output.svg)

rule calcMobilityScore:
    input:
        json = rules.anosimEGF.output.json,
    output:
        json = "Pipeline/Paper/Supplementary/JSON/HeLaEGFTreatment_{experiment}.json",
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        from scipy.stats import combine_pvalues
        from scipy.stats import  ttest_ind
        from statsmodels.stats.multitest import multipletests
        import plotly.graph_objs as go

        data = RAPDORData.from_file(input.json)
        print(data.norm_array.shape)
        positionwise_mobility = data.norm_array
        indices = data.internal_design_matrix.groupby("Treatment",group_keys=True,observed=False).apply(lambda x: list(x.index))
        print(data._treatment_means.shape)
        mobility = np.mean(np.abs(data._treatment_means[0] - data._treatment_means[1]), axis=-1)
        max_ctrl = np.argmax((data._treatment_means[0] - data._treatment_means[1]), axis=-1)
        max_treat = np.argmax((data._treatment_means[1] - data._treatment_means[0]), axis=-1)
        data.df["max ctrl"] = max_ctrl
        data.df["max treat"] = max_treat
        data.df["mobility score"] = mobility
        data.to_json(output.json)


rule plotJSDMScoreCorrelation:
    input:
        jsons = expand(rules.calcMobilityScore.output.json, experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)])
    output:
        svg = "Pipeline/Paper/Supplementary/Figures/FigureS9.svg",
        source_data = "Pipeline/Paper/SourceData/FS9.tsv"
    run:
        import plotly.graph_objs as go
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import _color_to_calpha
        from plotly.subplots import make_subplots
        from scipy.stats import pearsonr
        import numpy as np
        rows = cols = 2
        titles = [f"{x} min" for x in (2, 8, 20, 90)]
        color = _color_to_calpha(DEFAULT_COLORS[0],0.3)
        fig = make_subplots(rows =rows, cols = cols, y_title="Mobility Score", x_title="Jensen-Shannon Distance", subplot_titles=titles, vertical_spacing=0.15)
        dfs = []
        for idx, data in enumerate(input.jsons):
            data = RAPDORData.from_file(data)
            x = data.df["Mean Distance"]
            y = data.df["mobility score"]
            source_data = pd.DataFrame()
            source_data["Gene.names"] = data.df["Gene.names"]

            source_data[f"Jensen-Shannon Distance - {titles[idx]}"] = x
            source_data[f"mobility score - {titles[idx]}"] = y
            dfs.append(source_data)
            corr, pval = pearsonr(x, y)
            row = idx // cols
            col = idx % cols
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    showlegend=False,
                    mode="markers",
                    marker=dict(color=color),

                ),
                row=row+1, col=col+1
            )
            domain = str(idx+1) if idx > 0 else ""
            fig.add_annotation(
                text=f"Pearson R: {corr:.2f}",
                x=1,
                y=0,
                xref=f"x{domain} domain",
                yref=f"y{domain} domain",
                xanchor="right",
                yanchor="bottom",
                showarrow=False,
                bgcolor="#dbdbdb",
                borderpad=4,

            )
        fig.update_layout(template=DEFAULT_TEMPLATE)
        fig.update_layout(
            width=config["width"],
            height=config["height"]
        )
        fig.write_image(output.svg)
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.merge(df, on="Gene.names")
        merged_df.to_csv(output.source_data, sep="\t", index=False)

rule prepareForLimma:
    input:
        json = rules.calcMobilityScore.output.json
    output:
        tsv = temporary("Pipeline/NatureSpatial/limma/prep{experiment}FR{fraction}.tsv")
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        data = RAPDORData.from_file(input.json)
        fraction = int(wildcards.fraction) - 1
        proteins = data.df[(data.df["max ctrl"] == fraction) | (data.df["max treat"] == fraction)]["RAPDORid"]

        d_fraction = data.norm_array[:, :, fraction]
        d_ctrl = d_fraction[:, data.indices[0]]
        d_treat = d_fraction[:, data.indices[1]]
        ctrl_columns = ['CTRL1', 'CTRL2', 'CTRL3', 'CTRL4']
        treat_columns = ['EGF1', 'EGF2', 'EGF3', 'EGF4']
        df_ctrl = pd.DataFrame(d_ctrl, columns=ctrl_columns)
        df_treat = pd.DataFrame(d_treat, columns=treat_columns)
        df = pd.concat([data.df["RAPDORid"], df_ctrl, df_treat], axis=1)
        ctrl_filter = np.sum((~(df.loc[:, ctrl_columns] == 0)), axis=1)
        treat_filter = np.sum((~(df.loc[:, treat_columns] == 0)), axis=1)
        filter = ((ctrl_filter >= 3) | (treat_filter >= 3))
        df = df[filter]
        df = df[df["RAPDORid"].isin(proteins)]

        df.to_csv(output.tsv, sep="\t", index=False)
        #df.to_csv(output.tsv)


rule joinRAPDORNature:
    input:
        jsons = expand(rules.anosimEGF.output.json,experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
    output:
        file ="FOOOO.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        dfs = []
        d = (2, 8, 20, 90)
        for idx, data in enumerate(input.jsons):
            data = RAPDORData.from_file(data)
            df = data.df
            if idx == 0:
                df = df[["Gene.names", "global ANOSIM adj p-Value", "Mean Distance"]]
            else:
                df = df[["global ANOSIM adj p-Value", "Mean Distance"]]
            df = df.rename({"global ANOSIM adj p-Value": f"pval {d[idx]}"}, axis=1)
            df[f"{d[idx]} sig"] = df[f"pval {d[idx]}"] <= 0.1

            dfs.append(df)

        mdf = pd.concat(dfs, axis=1)
        mdf["s"] = mdf[[f"{d[idx]} sig" for idx, _ in enumerate(input.jsons)]].sum(axis=1)
        mdf = mdf[mdf["s"] >= 3]
        print(mdf)



rule runLimma:
    input:
        tsv = rules.prepareForLimma.output.tsv,
    output:
        tsv = "Pipeline/NatureSpatial/limma/pvals{experiment}FR{fraction}.tsv"
    conda: "../envs/limma.yml"
    script:
        "../Rscripts/limmaR.R"


rule collectLimmaResults:
    input:
        json=rules.calcMobilityScore.output.json,
        fraction_data= expand(rules.runLimma.output.tsv, fraction=list(range(1, 7)), allow_missing=True)
    output:
        json = "Pipeline/NatureSpatial/limma/final{experiment}.json"
    run:
        from RAPDOR.datastructures import RAPDORData
        from scipy.stats import combine_pvalues
        import numpy as np
        from statsmodels.stats.multitest import multipletests
        data = RAPDORData.from_file(input.json)
        frac_data = []
        for idx, file in enumerate(input.fraction_data):
            df = pd.read_csv(file, sep="\t" )
            df = df[["adj.P.Val"]]
            df = df.rename({"adj.P.Val": f"p.adj Fraction{idx+1}"}, axis=1)
            frac_data.append(df)
            t = "WDR75"
        frac_data = pd.concat(frac_data, axis=1)
        max_ctrl = data.df["max ctrl"]
        max_treat = data.df["max treat"]
        pvals = []
        pv1s = []
        pv2s = []
        for idx, _ in enumerate(max_ctrl):
            protein = data.df.loc[idx, "RAPDORid"]
            gene = data.df
            m_ctrl = max_ctrl[idx]
            m_treat = max_treat[idx]
            if protein in frac_data.index:
                pval1 = frac_data.loc[protein].iloc[m_ctrl]
                pval2 = frac_data.loc[protein].iloc[m_treat]
            else:
                pval1 = np.nan
                pval2 = np.nan
            if np.isnan(pval1):
                pval1 = 1
            if np.isnan(pval2):
                pval2 = 1
            comb = combine_pvalues([pval1, pval2])
            pv1s.append(pval1)
            pv2s.append(pval2)
            pvals.append(comb.pvalue)
        pvals = np.asarray(pvals)
        mask = np.isnan(pvals)
        _, pvals[~mask], _, _ = multipletests(pvals[~mask],method="fdr_bh")
        data.df["limma adj p-Value"] = pvals
        data.df["ctrl p-value"] = pv1s
        data.df["treat p-value"] = pv2s
        data.to_json(output.json)




rule plotEGFHeLa:
    input:
        jsons = expand(rules.anosimEGF.output.json,experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
    output:
        svg = "Pipeline/Paper/Subfigures/HeLaPlot.svg",
        html = "Pipeline/Paper/Subfigures/HeLaPlot.html",
        dist_svg = "Pipeline/Paper/Subfigures/HeLaPlotDistribution.svg",
        source_data = "Pipeline/Paper/SourceData/F9A.tsv",
        source_datab = "Pipeline/Paper/SourceData/F9B.tsv",
    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import plot_protein_distributions, _color_to_calpha
        from pyfunctions.plotlyVenn import venn_to_plotly
        import plotly.graph_objs as go
        from plotly.subplots import make_subplots
        novel = ["PRELID1", "CHUK", "RNF220", "ELL2", "GRB2", "SHC1", "CBL"]
        identified = {
            0: ["TXNL4B", "DNAJB2", "BDP1", "EDEM1", "PDE8A", "PRELID1", "CHUK", "STBD1", "GRB2", "ZNF131", "UNC80", "SLC25A33", "SHC1", "CIAO1", "CBL", "JMJD4",],
            1: ["RNF220", "CHUK", "SFT2D2", "PRELID1", "MYH10", "TBC1D10B", "NBEAL1", "FANCG", "KCTD18", "CES3", "GTF2IRD1", "GRB2", "ELL2", "CEMIP2", "CBL", "CST3", "SHC1"],
            2: ["MYH10", "CHUK", "BDP1", "PWWP3A", "WRAP73", "NBEAL1", "MAEA", "PYROXD2", "TEDC1", "CSNK1G2", "SIRT7", "RNF220", "CES3", "FIG4", "GORASP2", "ACP6", "FAM208B", "UNC80", "NDUFA11"],
            3: ["EGR1", "P3H4", "GATAD2B", "NSMCE3", "IRF2BP2", "MBNL1", "RBM7", "CD276"]
        }
        rows = 4
        cols = 2
        titles = [f"{x} min" for x in (2, 8, 20, 90)]
        mapping = (2, 8, 20, 90)
        fig = make_subplots(
            rows =rows,
            cols = cols,
            row_titles=titles,
            y_title="log<sub>10</sub>(p-value)",
            vertical_spacing=0.05,
            horizontal_spacing=0.05,
            column_widths=[0.7, 0.3]
        )
        dreg = []
        dfs = []
        for idx, data in enumerate(input.jsons):
            row = idx // cols
            col = idx % cols
            plot_ref = (idx + 1) * 2 - 1
            mins = mapping[idx]

            data = RAPDORData.from_file(data)
            found = identified[idx] if idx in identified else novel
            f  = data.df["Gene.names"].isin(found)
            check = list(data.df[f]["Gene.names"])

            df = data.df
            df["-log10(p-value)"] = -1 * np.log10(df["global ANOSIM adj p-Value"])
            difreg = df[(df["global ANOSIM adj p-Value"] <= 0.05) & (df["Mean Distance"] >= 0.2)]
            difreg = set(difreg["Gene.names"])
            dreg.append(difreg)

            if len(check) != len(found):

                for p in found:
                    if p not in check:
                        print(p)
                raise ValueError
            fig.add_trace(
                go.Scatter(
                    y=data.df[~f & (data.df["-log10(p-value)"] > 1)]["-log10(p-value)"],
                    x=data.df[~f & (data.df["-log10(p-value)"] > 1)]["Mean Distance"],
                    mode="markers",
                    text=data.df[~f]["Gene.names"], marker=dict(color=DEFAULT_COLORS[0]),
                    showlegend=True if idx == 0 else False,
                    name="RAPDOR p-value \u2264 0.1",

                ), row = idx+1, col=1
            )
            ccolor = _color_to_calpha(DEFAULT_COLORS[0], 0.3)
            fig.add_trace(
                go.Scatter(
                    y=data.df[~f & (data.df["-log10(p-value)"] <= 1)]["-log10(p-value)"],
                    x=data.df[~f & (data.df["-log10(p-value)"] <= 1)]["Mean Distance"],
                    mode="markers",
                    text=data.df[~f]["Gene.names"],marker=dict(color=ccolor),
                    showlegend=True if idx == 0 else False,
                    name="RAPDOR p-value > 0.1",

                ),row=idx + 1,col=1
            )

            fig.add_trace(
                go.Scatter(
                    y=data.df[f]["-log10(p-value)"],
                    x=data.df[f]["Mean Distance"],
                    mode="markers",
                    name="shift in Martinez-Val et al.",
                    showlegend=True if idx == 0 else False,
                    text=data.df[f]["Gene.names"],marker=dict(color=DEFAULT_COLORS[1])
                ),row=idx + 1,col=1
            )
            fig.add_hline(
                y=-1 * np.log10(0.1),
                row=idx+1, col=1, line=dict(dash='dot')
            )
            rapdor_identifed = set(data.df[data.df["-log10(p-value)"] > 1]["RAPDORid"].tolist())
            original_identifed = set(data.df[f]["RAPDORid"].tolist())
            venn = venn_to_plotly(L_sets=(rapdor_identifed, original_identifed), L_labels=("RAPDOR", "Martinez-Val et. al"), L_color=DEFAULT_COLORS)
            for iidx, shape in enumerate(venn['layout']['shapes']):
                shape["xref"] = f"x{plot_ref+1}"
                shape["yref"] = f"y{plot_ref+1}"
                fig.add_shape(shape)
            for annotation in venn['layout']['annotations'][2:]:
                if annotation["xref"] == "x":
                    annotation["xref"] = f"x{plot_ref+1}"
                    annotation["yref"] = f"y{plot_ref+1}"
                else:
                    annotation["xref"] = f"x{plot_ref+1} domain"
                    annotation["yref"] = f"y{plot_ref+1} domain"
                fig.add_annotation(annotation)


            fig.update_xaxes(venn["layout"]["xaxis"],row=idx+1,col=2)
            fig.update_yaxes(venn["layout"]["yaxis"],row=idx+1,col=2)
            fig.update_yaxes(scaleanchor=f"x{plot_ref+1}",scaleratio=1,col=2, row=idx+1)
            #fig.update_xaxes(scaleanchor=f"y{plot_ref+1}", scaleratio=1,col=2, row=idx+1)

            highlight = {"ITSN1": (0.05, 0.5), "MITF": (0., 1.45), "FOXJ3": (0.5, 0.6)} if idx > 1 else {"ITSN1":(0, 1.1) , "MITF": (0., 1.4), "FOXJ3": (0.2, 1.55),  "GRB2": (.5, .6), "CBL": (.2, .55), "SHC1":(0.12, 0.8) }
            if idx == 3:
                #highlight["ITSN1"] = (0.11, 0.3)
                highlight["FOXJ3"] = (0.5, 1.35)
            sdf = df[df["Gene.names"].isin(highlight)]
            for (_, row) in sdf.iterrows():
                y = row["-log10(p-value)"]
                x = row["Mean Distance"]
                ax, ay = highlight[row["Gene.names"]]
                anno = dict(
                    text=row["Gene.names"],
                    x=x,
                    y=y,
                    xanchor="center",
                    yanchor="middle",
                    showarrow=True,
                    xref="x" if idx == 0 else f"x{plot_ref}",
                    yref="y" if idx == 0 else f"y{plot_ref}",
                    ax=ax,
                    ay=ay,
                    axref="x" if idx == 0 else f"x{plot_ref}",
                    ayref="y" if idx == 0 else f"y{plot_ref}",

                )
                fig.add_annotation(
                    anno
                )
            df = df.loc[:, ["Gene.names", "-log10(p-value)", "Mean Distance"]]
            df = df.rename({"-log10(p-value)": f"-log10(p-value) {mins} min", "Mean Distance": f"Jensen-Shannon distance {mins} min"}, axis=1)
            df[f"shift in Martines-Val et al. at {mins} min"] = df["Gene.names"].isin(found)
            dfs.append(df)
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.merge(df,on="Gene.names")
        merged_df.to_csv(output.source_data, sep="\t", index=False)

        data = RAPDORData.from_file(input.jsons[0])
        sdf = data.df
        ids = sdf[sdf["Gene.names"].isin(["FOXJ3", "MITF", "ITSN1"])]["RAPDORid"]

        source_data = data.df.loc[data.df["RAPDORid"].isin(ids), ["Gene.names"]]
        indices = source_data.index

        source_data = source_data.loc[source_data.index.repeat(8)].reset_index(drop=True)
        source_data["type"] = (["Control"] * 4 + ["EGF 2 min"] * 4) * (source_data.shape[0] // 8)

        subdata = data.norm_array[indices]
        rsubdata = subdata.reshape(source_data.shape[0],-1)

        source_data.loc[:, [f"rel. Intensity (FR {x+1})" for x in range(rsubdata.shape[1])]] = rsubdata
        source_data.to_csv(output.source_datab,sep="\t",index=False)


        dist_fig = plot_protein_distributions(ids,data,mode="bar",colors=DEFAULT_COLORS,barmode="overlay", vertical_spacing=0.03)
        dist_fig.update_traces(legend="legend1")
        dist_fig.data[1].update(showlegend=False)
        dist_fig.data[3].update(showlegend=False)
        dist_fig.data[0].update(name="Control")
        dist_fig.data[2].update(name="EGF")
        dist_fig.update_traces(marker=dict(size=5),selector=dict(mode='markers'))
        # dist_fig.update_layout(
        #     barmode="overlay",
        #     bargroupgap=0,
        #     bargap=0,
        # )
        # dist_fig.update_xaxes(
        #     tickvals=list(range(len(data.fractions))),
        #     ticktext=[val.replace(" ", "<br>").replace("<br>&<br>", " &<br>") for val in data.fractions],
        #     tickmode="array"
        # )
        fig.update_layout(
            legend=dict(
                orientation="h",
                x=0.1,
                y=0,
                xref="container",
                yref="container",
                xanchor="left",
                yanchor="bottom",
            )
        )
        dist_fig.update_layout(
            template=DEFAULT_TEMPLATE,
            legend=dict(
                orientation="h",
                x=0.1,
                y=0,
                xref="container",
                yref="container",
                xanchor="left",
                yanchor="bottom",
            )
        )
        dist_fig.update_layout(width=config["width"])
        dist_fig.update_layout(height=config["HeLaPlot"]["B"]["height"])
        dist_fig.update_layout(margin=config["margin"])

        fig.update_layout(template=DEFAULT_TEMPLATE)
        fig.update_layout(width=config["width"])
        fig.update_layout(height=config["HeLaPlot"]["A"]["height"])
        fig.update_layout(margin=config["margin"])

        fig.update_annotations(
            dict(font=config["fonts"]["axis"])
        )
        fig.update_xaxes(
            showgrid=False, zeroline=False, showline=False,
            col=2
        )
        fig.update_yaxes(
            showgrid=False,zeroline=False, showline=False,
            col=2
        )
        fig.update_xaxes(
            title="Jensen-Shannon distance",
            row=4, col=1
        )
        dist_fig.update_annotations(
            dict(font=config["fonts"]["axis"])
        )
        dist_fig.update_layout(font=config["fonts"]["default"], legend=dict(font=config["fonts"]["legend"], title=None))
        dist_fig.update_layout(margin=config["margin"])
        dist_fig.update_layout(margin=dict(b=10, t=10))
        for annotation in dist_fig.layout.annotations:
            if annotation.text == "Fraction":
                annotation.update(y=annotation.y - 0.04)
        fig.write_image(output.svg)
        fig.write_html(output.html)
        dist_fig.write_image(output.dist_svg)



rule joinHeLaPlot:
    input:
        svg1 = rules.plotEGFHeLa.output.svg,
        svg2 = rules.plotEGFHeLa.output.dist_svg,
    output:
        svg = "Pipeline/Paper/Figure9.svg",
    run:
        from svgutils.compose import Figure, Panel, SVG, Text

        b_y = config["HeLaPlot"]["A"]["height"]
        ges_y = b_y + config["HeLaPlot"]["B"]["height"]
        f = Figure("624px",f"{ges_y}px",
            Panel(
                SVG(input.svg1),
                Text("A",2,15,size=config["multipanel_font_size"],weight='bold',font="Arial")

            ),
            Panel(
                SVG(input.svg2),
                Text("B",2,2,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ).move(0,b_y),
        )
        svg_string = f.tostr()
        svg_string = svg_string.decode().replace("encoding='ASCII'","encoding='utf-8'")
        with open(output.svg,"w") as handle:
            handle.write(svg_string)

rule theoreticalDistinctRValues:
    input:
        json=rules.run_RAPDOR.output.json,
    output:
        tsv = "Pipeline/Paper/unused/NrDistinctRValues.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        import pandas as pd
        import math
        data = RAPDORData.from_file(input.json)

        def distinct_perms(n):
            return math.factorial(2 * n) / (2 * np.square(math.factorial(n)))

        nr_protein = len(data.df[~data.df["contains empty replicate"].to_numpy()])
        df = {
            "Replicates per group": [],
            "local distinct R values": [],
            "global distinct R values": []
        }
        for n in range(2, 7):
            df["Replicates per group"].append(n)
            perms = distinct_perms(n)
            df["local distinct R values"].append(perms)
            df["global distinct R values"].append(perms * nr_protein)
        df = pd.DataFrame(df)
        df.to_csv(output.tsv, sep="\t", index=False)


rule produceArtificialSoftArgMaxExample:
    output:
        svg = "Pipeline/Paper/unused/artificialSoftArgMaxExample.svg"
    run:
        import numpy as np
        from scipy.special import rel_entr
        import plotly.graph_objs as go
        from plotly.subplots import make_subplots

        ctrl_mean = np.asarray([0.02, 0.02, 0.02, 0.02, 0, 0,0.25, 0.42, 0.25, 0., 0])

        treatment_mean = np.asarray([0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.16, 0.14, 0.1, 0])
        assert np.allclose(ctrl_mean.sum(), 1), ctrl_mean.sum()
        assert np.allclose(treatment_mean.sum(), 1), treatment_mean.sum()
        mid = .5 * (ctrl_mean + treatment_mean)
        rel1 = rel_entr(ctrl_mean,mid)
        rel2 = rel_entr(treatment_mean,mid)
        fig = go.Figure()
        x_range = np.arange(0, len(ctrl_mean))
        fig.add_trace(
            go.Scatter(
                y=ctrl_mean,
                x=x_range,
                name="Control mean",
                line=dict(color=DEFAULT_COLORS[0])
            )
        )
        fig.add_trace(
            go.Scatter(
                y=treatment_mean,
                x=x_range,
                name="Treatment mean",
                line=dict(color=DEFAULT_COLORS[1])

            )
        )

        for idx, beta in enumerate((10, 1000)):
            softmax1 = ((np.exp(beta * rel1)) / np.nansum(np.exp(beta * rel1),axis=-1,keepdims=True))
            softmax2 = ((np.exp(beta * rel2)) / np.nansum(np.exp(beta * rel2),axis=-1,keepdims=True))
            r1 = np.nansum(softmax1 * x_range, axis=-1)
            r2 = np.nansum(softmax2 * x_range, axis=-1)
            fig.add_vline(
                x=r1,
                #line=dict(color=DEFAULT_COLORS[idx]),
                name="r1"
            )
            fig.add_vline(
                x=r2,
                #line=dict(color=DEFAULT_COLORS[idx]),
                name=""

            )
            fig.add_annotation(
                x=r2,
                y=0.25,
                ax=r1,
                ay=0.25,
                axref="x",
                ayref="y",
                showarrow=True,
                arrowwidth=3,
                arrowhead=3,
                text=""
            )
            fig.add_annotation(
                x=(r2 + r1) / 2,
                y=0.26,
                axref="x",
                ayref="y",
                yanchor="bottom",
                text=f"<b>RPS</b><br>(\u03b2= {beta})",
                showarrow=False
            )
        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            width=config["width"],
            height=300,
            font=config["fonts"]["default"]
        )
        fig.write_image(output.svg)


rule GOTermEnrichmentMouseNature:
    input:
        file = rules.AnalyzeNatureWithRAPDOR.output.tsv
    conda: "../envs/ClusterProfiler.yml"
    threads: 10
    output:
        enriched = directory("Pipeline/NatureSpatial/GO/{experiment}GOEnrichment/",)
    script: "../Rscripts/AnalyzeGoTermsMouse.R"


rule KEGGEnrichmentMouseNature:
    input:
        file = rules.AnalyzeNatureWithRAPDOR.output.tsv
    conda: "../envs/ClusterProfiler.yml"
    threads: 10
    output:
        enriched = directory("Pipeline/NatureSpatial/KEGG/{experiment}GOEnrichment/",)
    script: "../Rscripts/AnalyzeKEGG.R"

rule plotKEGGEnrichment:
    input:
        file = expand(rules.KEGGEnrichmentMouseNature.output.enriched, experiment="egf_liver")
    output:
        svg = "Pipeline/Paper/Subfigures/KEGGEnrichment.svg"
    run:
        from pyfunctions.helpers import enrichment_plot_from_cp_table
        df_file = os.path.join(input.file[0], "AllFractions.tsv")
        df = pd.read_csv(df_file, sep="\t")
        df = df[df["category"] != "Human Diseases"]
        df["ONTOLOGY"] = "KEGG"
        df["Description"] = df["Description"].str.replace(" - Mus musculus (house mouse)", "")
        df = df.sort_values(by="p.adjust", ascending=False)

        fig = enrichment_plot_from_cp_table(df, mode="bar", colors=(DEFAULT_COLORS[1], DEFAULT_COLORS[0]))
        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            font=config["fonts"]["default"],
            width=config["width"],
            height = config["natureMouseFig"]["A"]["height"]

        )
        fig.write_image(output.svg)

rule plotExampleDistributions:
    input:
        json = expand(rules.AnalyzeNatureWithRAPDOR.output.json, experiment="egf_liver")
    output:
        svg = "Pipeline/Paper/Subfigures/NatureMouseDistributions.svg"
    run:
        from RAPDOR.datastructures import RAPDORData

        from RAPDOR.plots import plot_protein_distributions
        rapdor_data = RAPDORData.from_file(input.json[0])
        df = rapdor_data.df
        dist_config = config["natureMouseFig"]["B"]
        topids = dist_config["ids"]
        fig = plot_protein_distributions(
            topids, rapdordata=rapdor_data, colors=COLOR_SCHEMES["Dolphin"], mode="bar",
            plot_type="mixed", column_widths=[0.5, 0.5], horizontal_spacing=0.1, title_col="RAPDORid", vertical_spacing=dist_config["vertical_spacing"]
        )
        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=dist_config["height"])
        fig.update_layout(
            legend=dict(font=config["fonts"]["legend"], y=1.015),
            legend2=dict(font=config["fonts"]["legend"], y=dist_config["legend2_y"]),
        )
        fig.update_annotations(font=config["fonts"]["annotations"])
        fig.update_traces(
            error_y=dict(width=3)
        )
        fig.update_layout(
            margin=dict(l=60, b=120)
        )
        for annotation in fig.layout.annotations:
            if annotation.text == "Fraction":
                annotation.update(y=-0.22)
        print(fig.layout.annotations)
        fig.update_yaxes(nticks=2, col=2)
        fig.update_yaxes(nticks=3, col=1)
        fig.write_image(output.svg)

rule joinNaturePlot:
    input:
        a = rules.plotKEGGEnrichment.output.svg,
        b = rules.plotExampleDistributions.output.svg
    output:
        svg = "Pipeline/Paper/FigureMouse.svg"
    run:
        from svgutils.compose import Figure, Panel, SVG, Text

        b_y = config["natureMouseFig"]["A"]["height"]
        ges_y = b_y +  config["natureMouseFig"]["B"]["height"]
        f = Figure("624px",f"{ges_y}px",
            Panel(
                SVG(input.a),
                Text("A",2,15,size=config["multipanel_font_size"],weight='bold',font="Arial")

            ),
            Panel(
                SVG(input.b),
                Text("B",2,2,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ).move(0,b_y),
        )
        svg_string = f.tostr()
        svg_string = svg_string.decode().replace("encoding='ASCII'","encoding='utf-8'")
        with open(output.svg,"w") as handle:
            handle.write(svg_string)


rule joinQCPlot:
    input:
        a = rules.sampleCorrelation.output.corr,
        b = rules.pca.output.svg,
        c = rules.sampleCorrelation.output.histo
    output:
        svg = "Pipeline/Paper/Supplementary/Figures/FigureS2.svg"
    run:
        from svgutils.compose import Figure, Panel, SVG, Text

        c_y = config["QCPlot"]["A"]["height"]
        b_x = config["width"] // 2
        ges_y = c_y + config["QCPlot"]["C"]["height"]
        f = Figure("624px",f"{ges_y}px",
            Panel(
                SVG(input.a),
                Text("A",2,15,size=config["multipanel_font_size"],weight='bold',font="Arial")

            ),
            Panel(
                SVG(input.b),
                Text("B",2, 15,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ).move(b_x, 0),
            Panel(
                SVG(input.c),
                Text("C",2,2,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ).move(0,c_y),
        )
        svg_string = f.tostr()
        svg_string = svg_string.decode().replace("encoding='ASCII'","encoding='utf-8'")
        with open(output.svg,"w") as handle:
            handle.write(svg_string)



rule copyTableS2:
    input:
        s2 = "Data/svm.tsv",
        s3 = config["intensities"]
    output:
        s2 = "Pipeline/Paper/Supplementary/Tables/TableS2.xlsx",
        s3  = "Pipeline/Paper/Supplementary/Tables/TableS3.xlsx"
    run:
        import pandas as pd
        df = pd.read_csv(input.s2, sep="\t")
        df.to_excel(output.s2, index=False)
        df2 = pd.read_csv(input.s3, sep="\t")
        df2.to_excel(output.s3, index=False)

rule postProcessRapdorData:
    input:
        json=rules.run_RAPDOR.output.json,
        drop="Data/drop_columns.txt"
    output:
        file = "Pipeline/Paper/Supplementary/JSON/synechoRAPDORGradRFile.json",
    run:

        from RAPDOR.datastructures import RAPDORData
        with open(input.drop) as handle:
            drop = [line.rstrip() for line in handle]
        data = RAPDORData.from_file(input.json)
        data.pca()

        data.df = data.df.drop(drop, axis=1)
        data.to_json(output.file)



rule copySubfigures:
    input:
        s3 = expand(rules.plotTopHitDistributions.output, distribution=["S3"]),
        f6 = expand(rules.plotTopHitDistributions.output, distribution=["D1"]),
        wf2 = "Workflow2.svg",
        wf = "Workflow.svg",
        f1 = "Data/gradients_and_gels.svg",
        s2 = "Data/Regression.svg",
        s1 = "Data/northernblot.svg",
    output:
        s3 = "Pipeline/Paper/Supplementary/Figures/FigureS4.svg",
        f6 = "Pipeline/Paper/Figure7.svg",
        wf2 = "Pipeline/Paper/Figure10.svg",
        wf = "Pipeline/Paper/Figure2.svg",
        f1 =  "Pipeline/Paper/Figure1.svg",
        s2 = "Pipeline/Paper/Supplementary/Figures/FigureS3.svg",
        s1 = "Pipeline/Paper/Supplementary/Figures/FigureS1.svg",
    shell:
        """
        cp {input.s3[0]} {output.s3}
        cp {input.f6[0]} {output.f6}
        cp {input.wf2} {output.wf2}
        cp {input.wf} {output.wf}
        cp {input.f1} {output.f1}
        cp {input.s2} {output.s2}
        cp {input.s1} {output.s1}
        """

rule joinSourceData:
    input:
        f3a = rules.rplaDistribution.output.tsv,
        f3b = rules.plotMeanDistribution.output.tsv,
        f3c = rules.plotBarcodePlot.output.tsv,
        f5a = rules.plotFigureX.output.source_data,
        f5c = rules.plotComparisonExample.output.source_data,
        f6 = expand(rules.createBubblePlot.output.source_data, highlight="topHits")[0],
        f7 = expand(rules.plotTopHitDistributions.output.source_data, distribution=["D1"])[0],
        f9a = rules.plotEGFHeLa.output.source_data,
        f9b = rules.plotEGFHeLa.output.source_datab,
        s4 = expand(rules.plotTopHitDistributions.output.source_data,distribution=["S3"])[0],
        s5gradr = rules.plotANOSIMRDistribution.output.source_gradr,
        s5rdeep = rules.plotANOSIMRDistribution.output.source_rdeep,
        s5egf2min = rules.plotANOSIMRDistribution.output.source_egf2,
        s5egf8min = rules.plotANOSIMRDistribution.output.source_egf8,
        s5egf20min = rules.plotANOSIMRDistribution.output.source_egf20,
        s5egf90min = rules.plotANOSIMRDistribution.output.source_egf90,
        sf6ac = rules.plotFigureX.output.source_data2,
        sf7 = rules.meanRuntimeAndRAM.output.source_data,
        sf9 = rules.plotJSDMScoreCorrelation.output.source_data
    output:
        xlsx = "Pipeline/Paper/SourceData/SourceDataFigures.xlsx"
    run:
        import pandas as pd
        f3a = pd.read_csv(input.f3a, sep="\t")
        f3b = pd.read_csv(input.f3b, sep="\t")
        f3c = pd.read_csv(input.f3c, sep="\t")
        f5a = pd.read_csv(input.f5a, sep="\t")
        f5c = pd.read_csv(input.f5c, sep="\t")
        f6 = pd.read_csv(input.f6, sep="\t")
        f7 = pd.read_csv(input.f7, sep="\t")
        f9a = pd.read_csv(input.f9a, sep="\t")
        f9b = pd.read_csv(input.f9b, sep="\t")
        s4 = pd.read_csv(input.s4, sep="\t")
        s5gradr = pd.read_csv(input.s5gradr, sep="\t")
        s5rdeep = pd.read_csv(input.s5rdeep, sep="\t")
        s5egf2min = pd.read_csv(input.s5egf2min, sep="\t")
        s5egf8min = pd.read_csv(input.s5egf8min, sep="\t")
        s5egf20min = pd.read_csv(input.s5egf20min, sep="\t")
        s5egf90min = pd.read_csv(input.s5egf90min, sep="\t")
        sf6ac = pd.read_csv(input.sf6ac, sep="\t")
        sf7 = pd.read_csv(input.sf7, sep="\t")
        sf9 = pd.read_csv(input.sf9, sep="\t")
        with pd.ExcelWriter(output.xlsx, engine="openpyxl") as writer:
            f3a.to_excel(writer,sheet_name="Figure3A",index=False)
            f3b.to_excel(writer,sheet_name="Figure3B",index=False)
            f3c.to_excel(writer,sheet_name="Figure3C",index=False)
            f5a.to_excel(writer,sheet_name="Figure5A",index=False)
            f5c.to_excel(writer,sheet_name="Figure5C",index=False)
            f6.to_excel(writer,sheet_name="Figure6",index=False)
            f7.to_excel(writer,sheet_name="Figure7",index=False)
            f9a.to_excel(writer,sheet_name="Figure9A",index=False)
            f9b.to_excel(writer,sheet_name="Figure9B",index=False)
            s4.to_excel(writer,sheet_name="SupplementaryFigure4",index=False)
            s5gradr.to_excel(writer,sheet_name="SupplementaryFigure5GradR",index=False)
            s5rdeep.to_excel(writer,sheet_name="SupplementaryFigure5RDeeP",index=False)
            s5egf2min.to_excel(writer,sheet_name="SupplementaryFigure5EGF2min",index=False)
            s5egf8min.to_excel(writer,sheet_name="SupplementaryFigure5EGF8min",index=False)
            s5egf20min.to_excel(writer,sheet_name="SupplementaryFigure5EGF20min",index=False)
            s5egf90min.to_excel(writer,sheet_name="SupplementaryFigure5EGF90min",index=False)
            sf6ac.to_excel(writer,sheet_name="SupplementaryFigure6AandC",index=False)
            sf7.to_excel(writer,sheet_name="SupplementaryFigure7",index=False)
            sf9.to_excel(writer,sheet_name="SupplementaryFigure9",index=False)


rule svgsToPngs:
    input:
        expand("Pipeline/Paper/Figure{i}.svg", i=list(range(1, 8)) + list(range(9, 11)) )
    output:
        expand("Pipeline/Paper/Figure{i}.png", i=list(range(1, 8)) + list(range(9, 11)))
    shell:
        """
        for svg in {input}; do
            png="${{svg%.svg}}.png"
            inkscape "$svg" --export-type=png --export-filename="$png" --export-dpi=600
        done
        """

rule svgsToPngs_supplementary:
    input:
        expand("Pipeline/Paper/Supplementary/Figures/FigureS{i}.svg", i=[1, 2, 3, 4, 5, 6, 7, 9])
    output:
        expand("Pipeline/Paper/Supplementary/Figures/FigureS{i}.png", i=[1, 2, 3, 4, 5, 6, 7, 9])
    shell:
        """
        for svg in {input}; do
            png="${{svg%.svg}}.png"
            inkscape "$svg" --export-type=png --export-filename="$png" --export-dpi=600
        done
        """

rule zipSupplementaryFile1:
    input:
        stress = expand(rules.cpRedistSynecho.output.json, condition=["COLD", "HEAT", "DARK", "N", "Fe"]),
        synechorapdor = rules.postProcessRapdorData.output.file,
        hela = expand(rules.calcMobilityScore.output.json, experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
        rdeep = rules.runRAPDORonRDeeP.output.json
    output:
        file = "Pipeline/Paper/Supplementary/SupplementaryFile1.zip"
    shell:
        "zip -j {output.file} {input.synechorapdor} {input.stress} {input.hela} {input.rdeep}"
