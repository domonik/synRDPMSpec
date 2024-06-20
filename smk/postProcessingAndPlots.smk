import os.path
import pandas as pd
import numpy as np

include: "rapdor.smk"
include: "benchmarking.smk"

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

        svm_table = svm_table[["old_locus_tag", "SVM Score", "SVM RNA-binding"]]
        svm_table = svm_table.sort_values(by="SVM Score", ascending=False)
        svm_table = svm_table.drop_duplicates("old_locus_tag")

        rapdor_table = rapdor_table.merge(svm_table, on="old_locus_tag", how="left")


        rapdor_table = rapdor_table.sort_values(by="SVM Score",ascending=False)
        rapdor_table["Gene"] = rapdor_table["Gene"].fillna("None")
        rapdor_table["Gene"] = rapdor_table["Gene"].replace("None", pd.NA)

        rapdor_table["SVM Rank"] = np.arange(1,len(rapdor_table) + 1)
        rapdor_table["RDeeP significant"] = (rapdor_table["ctrl_peak_p_value"] <= 0.05) | (rapdor_table["rnase_peak_p_value"] <= 0.05 )

        rapdor_table = rapdor_table.sort_values(by="Rank", ascending=False)

        rapdor_table.to_csv(output.file ,sep="\t", index=False)

rule extract_tophits:
    input:
        tsv = rules.joinAnalysisMethods.output.file
    output:
        tsv = "Pipeline/Paper/Tables/overlappingProteins.tsv",
        tsv2 = "Pipeline/Paper/Tables/overlappingProteins_no_ribosomes.tsv",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")

        svm_candidates = [300, 643, 421]
        top200SVM = set(
            df[df["RAPDORid"].isin(svm_candidates)]["old_locus_tag"]
        )

        topRDPM = set(df[df["ANOSIM R"] >= config["anosimCutoff"]]["old_locus_tag"])
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
        all = topRDPM | top200SVM
        all = df[df["old_locus_tag"].isin(all)]
        all.to_csv(output.tsv,sep="\t",index=False)
        all = all[all["ribosomal protein"] == False]
        all.to_csv(output.tsv2,sep="\t",index=False)


rule createTable1andTableS1:
    input:
        tsv = rules.extract_tophits.output.tsv2,
        tsv2 = rules.joinAnalysisMethods.output.file
    output:
        tsv = "Pipeline/Paper/Table1.tsv",
        tsv2 = "Pipeline/Paper/TableS2JoinedAllMethods.tsv",
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
        df.to_csv(output.tsv2, index=False, sep="\t")



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
        svg="Pipeline/Paper/Subfigures/rplAFigure.svg"
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
        subdata2 = data.array[data[rpla][0]]
        row_heights = [0.5, 0.15, 0.15, 0.1, 0.1]
        fig = make_subplots(rows=5, cols=1, shared_xaxes=True, vertical_spacing=0.01, row_heights=row_heights)
        idesign = data.internal_design_matrix[data.internal_design_matrix["Replicate"] == "BR1"]
        fig1 = plot_distribution(subdata, data.internal_design_matrix, offset=data.state.kernel_size // 2, colors=COLOR_SCHEMES["Dolphin"], show_outliers=False)
        fig2 = plot_barcode_plot(subdata2, idesign, colors=COLOR_SCHEMES["Dolphin"], fractions=data.fractions, scale_max=False)
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
        for i_idx, trace in enumerate(fig2["data"]):
            v = row_heights[1] + row_heights[2]

            trace.colorbar.update(
                len=v,
                yref="paper",
                y=row_heights[0] - v * (i_idx // 2),
                yanchor="top",
                nticks=5,
                x=1. + 0.05 if i_idx % 2 else 1.,
                showticklabels=True,
                thickness=0.05,
                thicknessmode="fraction",
                ticklabelposition="inside",
                title=dict(text="iBAQ BR1" if i_idx == 0 else "-", font=dict(size=12, color="rgba(0,0,0,0)" if i_idx else None)),
                tickfont=config["fonts"]["legend"],
                outlinewidth=1,
                outlinecolor="black"

            )
            fig.add_trace(trace,row=2 + i_idx, col=1)

        for row in [4, 5]:
            fig.update_yaxes(range=[0, 1], row=row)
            cat = ["Control" if row == 4 else "RNase"]
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
        for row in range(2, 6):
            fig.update_yaxes(row=row,showgrid=False)
            fig.update_xaxes(row=row,showgrid=False)


        ctrl_western = Image.open(input.ctrl_img)
        rnase_western = Image.open(input.rnase_img)
        fig.add_layout_image(dict(
            source=ctrl_western,
            x=0,
            y=0.5,
            xref="x4 domain",
            yref="y4 domain",
            xanchor="left",
            yanchor="middle",

        )
        )
        fig.add_layout_image(dict(
            source=rnase_western,
            x=0,
            y=0.5,
            xref="x5 domain",
            yref="y5 domain",
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
        fig.update_layout(template=DEFAULT_TEMPLATE, height=config["F3"]["A"]["height"], width=config["width"])
        fig.update_yaxes(row=4, showgrid=False)

        fig.update_yaxes( row=5, showgrid=False)
        fig.update_xaxes(title = "Fraction", row=5, range=[0.5, 20.5])
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
        svg = "Pipeline/Paper/Distribution{distribution}.svg"
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
        topids = df[df["old_locus_tag"].isin(dist_config["locus_tags"])]["RAPDORid"]
        fig = plot_protein_distributions(
            topids, rapdordata=rapdor_data, colors=COLOR_SCHEMES["Dolphin"],
            plot_type="zoomed", column_widths=[0.7, 0.3], horizontal_spacing=0.075, title_col="Protein name", vertical_spacing=dist_config["vertical_spacing"]
        )
        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=dist_config["height"])
        fig.update_layout(
            legend=dict(font=config["fonts"]["legend"], y=1.015),
            legend2=dict(font=config["fonts"]["legend"], y=dist_config["legend2_y"]),
        )
        fig.update_annotations(font=config["fonts"]["annotations"])
        fig.update_traces(line=dict(width=2),
                marker=dict(size=3),)
        fig.update_yaxes(nticks=2, col=2)
        fig.update_yaxes(nticks=3, col=1)
        fig.write_image(output.svg)

rule plotMeanDistribution:
    input:
        file= rules.run_RAPDOR.output.json,
    output:
        svg="Pipeline/Paper/Subfigures/MeanDistributions.svg"
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
            width=config["width"],
            template=DEFAULT_TEMPLATE,
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
        x = list(range(df["Rank"].min(), df.Rank.max()+1))
        small_subunit = df[df["small_ribosomal"] == True]["RAPDORid"]
        large_subunit = df[df["large_ribosomal"] == True]["RAPDORid"]
        photosystem = df[df["photosystem"] == True]["RAPDORid"]
        colors = (COLOR_SCHEMES["Dolphin"][0], COLOR_SCHEMES["Dolphin"][1], COLOR_SCHEMES["Viking"][0])

        d = {"small subunit": small_subunit, "large subunit": large_subunit, "photosystem": photosystem}
        fig = rank_plot(d, data, colors, orientation="h", triangles="inside", tri_x=17, tri_y=0.2)

        fig.update_layout(
            width=config["width"], height=config["F3"]["C"]["height"],
            legend=dict(bgcolor = 'rgba(0,0,0,0)', )

        )





        #fig.show()
        fig.write_image(output.svg)
        fig.write_html(output.html)

rule plotAllVenns:
    input:
        files = expand(rules.plotVennDiagramm.output.json, set=["others", "ribosomes", "rna_binding"])
    output:
        svg = "Pipeline/Paper/Subfigures/VennDiagramm_all.svg"
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
                text=f"<b>{annos[idx-1]}</b>",
                xref= f"x{idx} domain" if idx > 1 else "x domain",
                xanchor="left",
                x=0,
                yref=f"y{idx} domain" if idx > 1 else "y domain",
                yanchor="top",
                y=1,
                font=dict(size=config["multipanel_font_size"]),
                showarrow = False


            )
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
        html = "Pipeline/Paper/benchmark.html",
        svg = "Pipeline/Paper/benchmark.svg"
    run:
        import plotly.graph_objs as go
        import plotly.express as px
        import pandas as pd
        import plotly.subplots as sp
        import string

        annos = string.ascii_uppercase
        df = pd.read_csv(input.all_benchmark, sep="\t")
        fig = sp.make_subplots(rows=1,cols=2, x_title="Replicates per group", horizontal_spacing=0.2)
        names = ["s", "max_uss"]
        for col in range(1, 3):
            name = names[col-1]

            fig.add_trace(
                go.Bar(
                    x=df[(df["Tool"] == "RDeeP")]["Replicates per group"],
                    y=df[(df["Tool"] == "RDeeP")][name],
                    error_y=dict(
                        type='data',# value of error bar given in data coordinates
                        array=df[(df["Tool"] == "RDeeP")][f"stdev_{name}"],
                        visible=True,
                        color="black"
                    ),
                    legendgroup="RDeeP",
                    legendgrouptitle=dict(text="RDeeP"),
                    name="RDeeP",
                    marker_color=COLOR_SCHEMES["Viking"][0],
                    marker_line=dict(width=2,color='black'),
                    showlegend=True if col == 1 else False,

                ),
                col=col, row=1
            )


            fig.add_trace(
                go.Bar(
                    x=df[(df["ANOSIM"] == False) & (df["Tool"] == "RAPDOR")]["Replicates per group"],
                    y=df[(df["ANOSIM"] == False) & (df["Tool"] == "RAPDOR")][name],
                    error_y=dict(
                        type='data',# value of error bar given in data coordinates
                        array=df[(df["ANOSIM"] == False) & (df["Tool"] == "RAPDOR")][f"stdev_{name}"],
                        visible=True,
                        color="black"
                    ),

                    legendgroup="RAPDOR",
                    legendgrouptitle=dict(text="RAPDOR"),
                    marker_color=COLOR_SCHEMES["Dolphin"][0],
                    marker_line=dict(width=2,color='black'),
                    showlegend=True if col == 1 else False,


                    name="Ranking"
                ),
                col=col, row=1
            )
            fig.add_trace(
                go.Bar(
                    x=df[(df["ANOSIM"] == True) & (df["Tool"] == "RAPDOR")]["Replicates per group"],
                    y=df[(df["ANOSIM"] == True) & (df["Tool"] == "RAPDOR")][name],
                    error_y=dict(
                        type='data',# value of error bar given in data coordinates
                        array=df[(df["ANOSIM"] == True) & (df["Tool"] == "RAPDOR")][f"stdev_{name}"],
                        visible=True,
                        color="black"
                    ),
                    legendgroup="RAPDOR",
                    legendgrouptitle=dict(text="RAPDOR"),
                    marker_color=COLOR_SCHEMES["Dolphin"][1],
                    marker_line=dict(width=2,color='black'),
                    showlegend=True if col == 1 else False,


                    name="ANOSIM"
                ),
                col=col,row=1

            )
            fig.add_annotation(
                text=f"<b>{annos[col - 1]}</b>",
                xref=f"x{col} domain" if col > 1 else "x domain",
                xanchor="left",
                x=-0.4,
                yref=f"y{col} domain" if col > 1 else "y domain",
                yanchor="top",
                y=1.1,
                font=dict(size=config["multipanel_font_size"]),
                showarrow=False

            )

        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=300,
        )
        fig.update_yaxes(type="log", title=dict(text="Average Runtime [s]"), col=1)
        fig.update_yaxes(title=dict(text="Memory [Mb]"), col=2)

        fig.update_xaxes(type='linear')
        fig.write_image(output.svg)
        fig.write_html(output.html)


rule createBubblePlot:
    input:
        file=rules.run_RAPDOR.output.json,
        joined=rules.extract_tophits.output.tsv2,
    output:
        svg = "Pipeline/Paper/Subfigures/DimensionReduction_ids{highlight}.svg",
        html = "Pipeline/plots/DimensionReduction_ids{highlight}.html",
        json = "Pipeline/plots/DimensionReduction_ids{highlight}.json"
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
        data.df['Gene'] = data.df.apply(lambda row: row['old_locus_tag'] if pd.isnull(row['Gene']) or row['Gene'] == '' else row['Gene'], axis=1)
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
                annotation.update(showarrow=True,ax=-12.5, ay=-.2, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Tsf":
                annotation.update(showarrow=True,ay=.1,ax=annotation.x,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "FtsH2":
                annotation.update(showarrow=True,ay=-.1,ax=-8,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "CysC":
                annotation.update(showarrow=True,ay=.075,ax=1,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RpoB":
                annotation.update(showarrow=True,ay=annotation.y,ax=0, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "RbpA":
                annotation.update(showarrow=True,ay=-.35,ax=-7.5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0147":
                annotation.update(showarrow=True,ay=-.4,ax=-5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr1143":
                annotation.update(showarrow=True,ay=-.45,ax=1, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Hpf":
                annotation.update(showarrow=True,ay=.3,ax=-16, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1388":
                annotation.update(showarrow=True,ax=-14,ay=.4, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "ChlI":
                annotation.update(showarrow=True,ax=-9 ,ay=.5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll0921":
                annotation.update(showarrow=True,ax=annotation.x,ay=.5, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Pgm":
                annotation.update(showarrow=True,ax=-3.5, ay=.57, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll1967":
                annotation.update(showarrow=True,ax=0,ay=.55, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0782":
                annotation.update(showarrow=True,ax=0.5,ay=.4, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "QueF":
                annotation.update(showarrow=True,ax=1.3,ay=.32, axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "SecE":
                annotation.update(showarrow=True,ax=-16,ay=-.05,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Ffh":
                annotation.update(showarrow=True,ax=14,ay=-.35,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll0284":
                annotation.update(showarrow=True,ax=10,ay=.5,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "AroB":
                annotation.update(showarrow=True,ax=5,ay=-.375,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Slr0678":
                annotation.update(showarrow=True,ax=5.5,ay=.43,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)
            elif annotation.text == "Sll7087":
                annotation.update(showarrow=True,ax=-7,ay=.56,axref=annotation.xref,ayref=annotation.yref,arrowcolor='black',)


        fig.write_image(output.svg)
        fig.write_html(output.html)
        fig.write_json(output.json)


rule createFigure4:
    input:
        venns = rules.plotAllVenns.output.svg,
        bubble = expand(rules.createBubblePlot.output.svg, highlight="topHits")
    output:
        svg = "Pipeline/Paper/Figure4.svg"
    run:
        from svgutils.compose import Figure, Panel, SVG, Text

        f = Figure("624px","650px",
            Panel(
                SVG(input.venns),
            ),
            Panel(
                SVG(input.bubble[0]),
                Text("D", 10, 30,size=config["multipanel_font_size"],weight='bold', font="Arial")
            ).move(0, 250)
        )
        svg_string = f.tostr()

        svg_string = svg_string.decode().replace("encoding='ASCII'", "encoding='utf-8'")
        with open(output.svg, "w") as handle:
            handle.write(svg_string)


rule plotConditionedSynechochoColdRibo:
    input:
        json = expand(rules.runOnSynData.output.json, condition="COLD")
    output:
        svg = "Pipeline/Paper/FigureS2.svg",
        svg2 = "Pipeline/Paper/FigureS3.svg",
        html2 = "Pipeline/Paper/FigureS3.html",
        tsv = "Pipeline/Paper/TableforFigureS4.tsv",
        json = "Pipeline/Paper/RAPDORforFigureS4.json",
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
        d = {"ribosomal & membrane": r_to_membrane, "Wang et al. (2023)": identified_inpub}
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
        fig.update_layout(width=config["width"])
        fig.update_layout(height=config["height"] / 2)

        fig.write_image(output.svg2)
        fig.write_html(output.html2)
        data.df["mentioned in wang"] = data.df["Gene"].isin(names)
        data.export_csv(output.tsv, sep="\t")
        data.to_json(output.json)


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
        file = "Pipeline/Paper/TableS1BasedOnKArray.tsv",
        file2 = "Pipeline/Paper/TableS1BasedOnPeaks.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        from scipy.stats import pearsonr
        import pandas as pd
        import plotly.graph_objs as go

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
        df.to_csv(output.file, sep="\t")
        df.loc['Sum'] = df.sum(axis=0)
        df2.to_csv(output.file2, sep="\t")
        df2.loc['Sum'] = df2.sum(axis=0)

rule calcANOSIMDistribution:
    input:
        json=rules.run_RAPDOR.output.json,
        json_liver=expand(rules.AnalyzeNatureWithRAPDOR.output.json,experiment="egf_liver")
    output:
        gradr_json="Pipeline/analyzed/ANOSIMDistributionGradR.json",
        liver_json="Pipeline/analyzed/ANOSIMDistributionLiver.json",
    threads: 10
    run:
        from RAPDOR.datastructures import RAPDORData
        gradr_data = input.json
        liver_data = input.json_liver[0]
        for plot_id, data in enumerate([(gradr_data, 3, output.gradr_json), (liver_data, 4, output.liver_json), ]):
            data, samples, outfile = data
            data = RAPDORData.from_file(data)
            data: RAPDORData
            print("starting")
            data.calc_anosim_p_value(999,threads=threads,mode="global")
            distribution = data._anosim_distribution
            print(plot_id, "done")
            with open(outfile, "wb") as handle:
                np.save(handle, distribution)
            data._anosim_distribution = None
            print("saved")



rule plotANOSIMRDistribution:
    input:
        np=rules.calcANOSIMDistribution.output.gradr_json,
        np_liver=rules.calcANOSIMDistribution.output.liver_json,
        json=rules.run_RAPDOR.output.json,
        json_liver=expand(rules.AnalyzeNatureWithRAPDOR.output.json,experiment="egf_liver")
    output:
        svg = "Pipeline/Paper/ANOSIMDistribution.svg",
    threads: 1
    run:
        from RAPDOR.datastructures import RAPDORData
        import plotly.graph_objs as go
        from plotly.subplots import make_subplots
        from scipy.stats import spearmanr
        gradr_data = input.json
        liver_data = input.json_liver[0]
        fig = make_subplots(rows=2)
        for plot_id, data in enumerate([(gradr_data, 3, input.np), (liver_data, 4, input.np_liver)]):
            data, samples, array = data
            data = RAPDORData.from_file(data)
            print("data loaded")
            data: RAPDORData
            with open(array, "rb") as handle:
                distribution = np.load(handle)
            print("distribution")
            indices = data.df[data.df["min replicates per group"] == samples]["id"].to_numpy()
            print(data.df["min replicates per group"].min())
            distribution = distribution[:, indices].flatten()
            distribution = data.df["ANOSIM R"]
            y, x = np.histogram(
                distribution,
                bins=np.linspace(-1, 1, int(np.floor(2/0.05)))
            )
            y = y / np.sum(y)

            fig.add_trace(go.Bar(
                x=x,
                y=y,
                marker=dict(
                    line=dict(color="black", width=1),
                    color=DEFAULT_COLORS[0],

                ),
            ), row=plot_id+1, col=1)
        fig.update_xaxes(title="ANOSIM R")
        fig.update_yaxes(title="probability")
        fig.update_layout(template=DEFAULT_TEMPLATE)
        fig.update_layout(width=config["width"])
        fig.update_layout(height=config["height"] / 2)
        fig.update_layout(font=config["fonts"]["default"])
        fig.write_image(output.svg)


rule theoreticalDistinctRValues:
    input:
        json=rules.run_RAPDOR.output.json,
    output:
        tsv = "Pipeline/Paper/NrDistinctRValues.tsv"
    run:
        from RAPDOR.datastructures import RAPDORData
        import numpy as np
        import pandas as pd
        data = RAPDORData.from_file(input.json)

        def distinct_perms(n):
            return np.math.factorial(2 * n) / (2 * np.square(np.math.factorial(n)))

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
        svg = "Pipeline/Paper/artificialSoftArgMaxExample.svg"
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
    output:
        enriched = directory("Pipeline/NatureMouse/GO/{experiment}GOEnrichment/",)
    script: "../Rscripts/AnalyzeGoTermsMouse.R"


rule KEGGEnrichmentMouseNature:
    input:
        file = rules.AnalyzeNatureWithRAPDOR.output.tsv
    conda: "../envs/ClusterProfiler.yml"
    output:
        enriched = directory("Pipeline/NatureMouse/KEGG/{experiment}GOEnrichment/",)
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


rule postProcessRapdorData:
    input:
        json=rules.run_RAPDOR.output.json,
        drop="Data/drop_columns.txt"
    output:
        file = "Pipeline/Paper/synechoRAPDORFile.json",
    run:
        from RAPDOR.datastructures import RAPDORData
        with open(input.drop) as handle:
            drop = [line.rstrip() for line in handle]

        data = RAPDORData.from_file(input.json)
        data.df = data.df.drop(drop, axis=1)
        data.to_json(output.file)

