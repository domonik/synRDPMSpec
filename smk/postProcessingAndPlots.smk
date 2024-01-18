
include: "rapdor.smk"
include: "benchmarking.smk"

from RAPDOR.plots import DEFAULT_TEMPLATE, COLOR_SCHEMES


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

rule extract_overlapping_proteins:
    input:
        tsv = rules.joinAnalysisMethods.output.file
    output:
        tsv = "Pipeline/Paper/Tables/overlappingProteins.tsv",
        tsv2 = "Pipeline/Paper/Tables/overlappingProteins_no_ribosomes.tsv",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")
        top200SVM = set(
            df[df["SVM RNA-binding"] == True]["old_locus_tag"]
        )
        top200RDPM = set(df[df["Rank"] <= config["nr_proteins"]]["old_locus_tag"])
        df = df[
            [
                "RAPDORid",
                "old_locus_tag",
                "Gene",
                "Locus tag",
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
        all = top200RDPM.intersection(top200SVM)
        all = df[df["old_locus_tag"].isin(all)]
        all.to_csv(output.tsv,sep="\t",index=False)
        all = all[all["ribosomal protein"] == False]
        all.to_csv(output.tsv2,sep="\t",index=False)


rule createTable1:
    input:
        tsv = rules.extract_overlapping_proteins.output.tsv2
    output:
        tsv = "Pipeline/Paper/Table1.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")
        df = df[["old_locus_tag", "Gene", "Rank", "Mean Distance"]]
        df = df.rename({"old_locus_tag": "locus tag"}, axis=1)
        df = df.sort_values(by="Rank")
        df.loc[df.Gene.str.contains("sll|slr", na=False), "Gene"] = ""
        df.to_csv(output.tsv, sep="\t", index=False, float_format='%.2f')


rule plotVennDiagramm:
        input:
            tsv = rules.joinAnalysisMethods.output.file
        output:
            html = "Pipeline/plots/VennDiagramm_set{set}.html",
            svg = "Pipeline/plots/VennDiagramm_set{set}.svg",
            tsv = "Pipeline/postProcessing/topCandidates_set{set}.tsv",
            json = "Pipeline/postProcessing/VennDiagramm_set{set}.json",
        run:
            NUMBER_PROTEINS = config["nr_proteins"]
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
            top200RDPM = set(result_df[result_df["Rank"] <= NUMBER_PROTEINS]["old_locus_tag"])

            topRDeep = set(
                result_df[(result_df["rnase_peak_p_value"] <= 0.05) | (result_df["ctrl_peak_p_value"] <= 0.05)]["old_locus_tag"]
            )

            all = topRDeep.union(top200RDPM).union(top200SVM)
            all = result_df[result_df["old_locus_tag"].isin(all)]
            all = all[["RAPDORid", "old_locus_tag", "Gene",  "Mean Distance", "relative fraction shift",  "ANOSIM R", "RDeeP significant", "SVM RNA-binding", "Rank", "ribosomal protein"]]
            all.to_csv(output.tsv, sep="\t", index=False)
            colors = list(COLOR_SCHEMES["Dolphin"]) + [COLOR_SCHEMES["Viking"][0]]
            colors = [_color_to_calpha(col, alpha=0.6) for col in colors]
            fig = venn_to_plotly(L_sets=(top200RDPM, top200SVM, topRDeep), L_labels=("RAPDOR", "TripepSVM", "RDeeP"), L_color=colors)
            fig.update_layout(template=DEFAULT_TEMPLATE)
            fig.update_layout(
                font=dict(size=18),
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
                showarrow=False
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


rule plotDistribution:
    input:
        file = rules.run_RAPDOR.output.json,
        joined = rules.joinAnalysisMethods.output.file
    output:
        svg = "Pipeline/plots/Distribution_ids{distributionids}.svg"
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        import math
        from RAPDOR.plots import plot_protein_distributions, COLOR_SCHEMES, DEFAULT_TEMPLATE
        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        ids = config["distributions"][wildcards.distributionids]
        if isinstance(ids, list):
            rapdorids = ids
            tcol = "Gene"
        else:
            df = pd.read_csv(input.joined, sep="\t")
            for key, value in ids.items():
                df = df[df[key] == value]
            tcol = "old_locus_tag"
            rapdorids = df["RAPDORid"]
        rows = math.ceil(len(rapdorids)/2)
        fig = plot_protein_distributions(rapdorids, rapdordata=data, colors=COLOR_SCHEMES["Dolphin"], title_col=tcol, cols=2, rows=rows, vertical_spacing=0.1 * (1/rows), horizontal_spacing=0.075)
        fig.update_layout(template=DEFAULT_TEMPLATE, legend1=dict(y=1.02))
        fig.update_yaxes(mirror=True, tickfont=dict(size=12))
        fig.update_xaxes(mirror=True, dtick=1, tickfont=dict(size=12))
        fig.update_layout(
            font=dict(size=16),
        )
        fig.update_layout(width=1200, height=150 * rows,
            legend2=dict(
                   y=1 + 0.4 * (1/rows)
            )
        )

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
        from RAPDOR.plots import multi_means_and_histo, COLOR_SCHEMES, DEFAULT_TEMPLATE

        with open(input.file) as handle:
            data = RAPDORData.from_json(handle.read())
        ids = data.df[data.df["small_ribosomal"] == True]["RAPDORid"]
        ids2 = data.df[data.df["large_ribosomal"] == True]["RAPDORid"]
        ids3 = data.df[data.df["photosystem"] == True]["RAPDORid"]
        print(ids, ids2, ids3)
        d = {"small subunit": ids, "large subunit": ids2, "photosystem": ids3}
        fig = multi_means_and_histo(d, data, colors=COLOR_SCHEMES["Dolphin"] + COLOR_SCHEMES["Viking"])
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


rule plotBarcodePlot:
    input:
        file = rules.joinAnalysisMethods.output.file
    output:
        svg = "Pipeline/Paper/Subfigures/ribosomalIdentifier.svg"
    run:
        import plotly.graph_objs as go
        import pandas as pd
        import numpy as np
        df = pd.read_csv(input.file, sep="\t")
        #df = df[df["ribosomal protein"] == True]
        df = df.sort_values(by="Rank", ascending=True)
        df = df[~df["ANOSIM R"].isna()]
        x = list(range(df["Rank"].min(), df.Rank.max()+1))
        small_subunit = (df["small_ribosomal"] == True).to_numpy().astype(float)
        small_subunit[small_subunit == 0] = np.nan
        large_subunit = (df["large_ribosomal"] == True).to_numpy().astype(float)
        large_subunit[large_subunit == 0] = np.nan
        photosystem = (df["photosystem"] == True).to_numpy().astype(float)
        photosystem[photosystem == 0] = np.nan
        fig = go.Figure()
        colors = (COLOR_SCHEMES["Dolphin"][0], COLOR_SCHEMES["Dolphin"][1], COLOR_SCHEMES["Viking"][0])
        triangles = []
        tri_x = 25
        tri_y = 0.1

        for idx, (key, data) in enumerate({"small subunit": small_subunit, "large subunit": large_subunit, "photosystem": photosystem}.items()):

            fig.add_trace(
                go.Bar(
                    x=x,
                    y=data,
                    marker_color=colors[idx],
                    marker_line=dict(width=2,color=colors[idx]),
                    name=key,
                    hovertext=df["Gene"],

                )
            )
            tri = np.nanmedian((data * np.asarray(x)))
            triangles.append(
                go.Scatter(
                    x=[tri-tri_x, tri, tri+tri_x, tri-tri_x],
                    y=[0, tri_y, 0, 0],
                    mode="lines",
                    fill="toself",
                    fillcolor=colors[idx],
                    marker_color=colors[idx],
                    showlegend=False
                )
            )
        fig.add_traces(triangles)
        fig.update_layout(barmode="overlay", bargap=0, template=DEFAULT_TEMPLATE, width=config["width"], height=config["height"] / 2)
        fig.update_layout(margin=dict(r=10, l=10, t=0, b=0),
            legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1, xanchor="right"))
        fig.update_xaxes(range=(x[0], x[-1]), title="Rank", mirror=True, showline=True, linecolor="black", showgrid=False)
        fig.update_yaxes(range=(0, 1), mirror=True, showline=True, showticklabels=False, linecolor="black", showgrid=False)
        #fig.show()
        fig.write_image(output.svg)

rule plotAllVenns:
    input:
        files = expand(rules.plotVennDiagramm.output.json, set=["others", "ribosomes", "rna_binding"])
    output:
        svg = "Pipeline/plots/VennDiagramm_all.svg"
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

        print(multi_fig["layout"]["xaxis"])
        print(multi_fig["layout"]["xaxis2"])
        print(multi_fig["layout"]["yaxis"])
        print(multi_fig["layout"]["yaxis2"])
        multi_fig.update_layout(template=DEFAULT_TEMPLATE)
        multi_fig.update_layout(font=dict(size=config["font_size"]))
        multi_fig.update_layout(
            legend={'itemsizing': 'trace',"font": dict(size=config["font_size"]), "orientation": "h", "yanchor": "bottom", "y": 1.01},
            margin=dict(b=0,l=10,pad=0,r=10,t=40),
            width=config["width"], height=250
        )

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

        fig.layout.annotations[0].update(font=dict(size=config["font_size"]))
        fig.update_layout(template=DEFAULT_TEMPLATE, width=config["width"], height=300,
            margin=dict(r=10, l=70, t=30, b=50)
        )
        fig.update_yaxes(type="log", title=dict(text="Average Runtime [s]"), col=1)
        fig.update_yaxes(title=dict(text="Memory [Mb]"), col=2)
        fig.update_layout(font=dict(size=config["font_size"]))

        fig.update_xaxes(type='linear')
        fig.write_image(output.svg)
        fig.write_html(output.html)


rule createBubblePlot:
    input:
        file=rules.run_RAPDOR.output.json,
        joined=rules.extract_overlapping_proteins.output.tsv2,
    output:
        svg = "Pipeline/plots/DimensionReduction_ids{highlight}.svg",
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
        else:
            raise ValueError("Not supported")

        fig = plot_dimension_reduction(
            rapdordata=data,
            colors=COLOR_SCHEMES["Dolphin"],
            highlight=ids,
            title_col="Gene",
            legend_spread=0.1,
            legend_start=0.2
        )
        fig.update_annotations(font=dict(size=config["font_size"]))
        for annotation in fig.layout.annotations:
            print(annotation.text)
            if annotation.text == "Mean Distance":
                annotation.update(font=dict(size=config["font_size"] + 2))
        fig.layout.update(
            template=DEFAULT_TEMPLATE,
            width=config["width"],
            height=config["height"],
            legend=dict(visible=False),
            margin=dict(t=20, r=10, l=10, b=10),
            font=dict(size=config["font_size"])
        )

        fig.write_image(output.svg)
        fig.write_html(output.html)
        fig.write_json(output.json)


rule createFigure4:
    input:
        venns = rules.plotAllVenns.output.svg,
        bubble = expand(rules.createBubblePlot.output.svg, highlight="overlapping")
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
                Text("D", 10, 30,size=config["multipanel_font_size"],weight='bold')
            ).move(0, 250)
        )
        f.save(output.svg)


rule plotConditionedSynechochoColdRibo:
    input:
        json = expand(rules.runOnSynData.output.json, condition="COLD")
    output:
        svg = "Pipeline/Paper/FigureS2.svg",
        svg2 = "Pipeline/Paper/FigureS3.svg",
        html2 = "Pipeline/Paper/FigureS3.html",

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
        fig.update_layout(margin=dict(r=10, l=10, t=0, b=0), height=config["height"] / 2)

        fig.write_image(output.svg2)
        fig.write_html(output.html2)