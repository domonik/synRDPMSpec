
include: "rdpmspecidentifier.smk"

from RDPMSpecIdentifier.plots import DEFAULT_TEMPLATE


rule joinAnalysisMethods:
    input:
        svm = "Data/svm.tsv",
        rdpmspec = rules.run_rdpmspecidentifier.output.tsv,
        rdeep = rules.runRDeep.output.outfile,
        go_binders = rules.extractGORNABinding.output.file
    output:
        file = "Pipeline/postProcessing/joinedTable.tsv"
    run:
        import pandas as pd
        import numpy as np
        svm_table = pd.read_csv(input.svm, sep="\t",  decimal=",")
        rdpmspec_table = pd.read_csv(input.rdpmspec, sep="\t")
        rdeep_table = pd.read_csv(input.rdeep, sep=";")
        go_table = pd.read_csv(input.go_binders, sep="\t")

        indices = list(rdeep_table.columns[0:1]) + list(rdeep_table.columns[-14:])
        rdeep_table = rdeep_table[indices]
        rdeep_table = rdeep_table[rdeep_table["nb_shift"] > 0]
        rdeep_table = rdeep_table.sort_values(by=["rnase_peak_p_value", "ctrl_peak_p_value"])
        rdeep_table = rdeep_table.drop_duplicates("protein_name")
        rdeep_table = rdeep_table.rename({"protein_name": "RDPMSpecID"}, axis=1)

        rdpmspec_table = rdpmspec_table[["RDPMSpecID"] + [col for col in rdpmspec_table.columns if col != "RDPMSpecID"]]

        rdpmspec_table = rdpmspec_table.merge(rdeep_table, on="RDPMSpecID", how="left")

        svm_table = svm_table.rename({"Gene": "old_locus_tag", "Score": "SVM Score"}, axis=1)
        svm_table["SVM RNA-binding"] = svm_table["Prediction"] == "RNA-binding protein"

        svm_table = svm_table[["old_locus_tag", "SVM Score", "SVM RNA-binding"]]
        svm_table = svm_table.sort_values(by="SVM Score", ascending=False)
        svm_table = svm_table.drop_duplicates("old_locus_tag")

        rdpmspec_table = rdpmspec_table.merge(svm_table, on="old_locus_tag", how="left")

        rdpmspec_table = rdpmspec_table.merge(go_table, on="old_locus_tag", how="left")

        rdpmspec_table = rdpmspec_table.sort_values(by="SVM Score",ascending=False)
        rdpmspec_table["Gene"] = rdpmspec_table["Gene"].fillna("None")
        rdpmspec_table["ribosomal protein"] = ((rdpmspec_table["Gene"].str.contains('rpl|rps|Rpl|Rps', case=False)) | (rdpmspec_table['ProteinFunction'].str.contains('ribosomal protein', case=False)))
        rdpmspec_table["Gene"] = rdpmspec_table["Gene"].replace("None", pd.NA)

        rdpmspec_table["SVM Rank"] = np.arange(1,len(rdpmspec_table) + 1)
        rdpmspec_table["RDeeP significant"] = (rdpmspec_table["ctrl_peak_p_value"] <= 0.05) | (rdpmspec_table["rnase_peak_p_value"] <= 0.05 )

        rdpmspec_table = rdpmspec_table.sort_values(by="Rank", ascending=False)

        rdpmspec_table.to_csv(output.file ,sep="\t", index=False)


rule plotVennDiagramm:
        input:
            tsv = rules.joinAnalysisMethods.output.file
        output:
            html = "Pipeline/plots/VennDiagramm_set{set}.html",
            svg = "Pipeline/plots/VennDiagramm_set{set}.svg",
            tsv = "Pipeline/postProcessing/topCandidates_set{set}.tsv",
            json = "Pipeline/postProcessing/VennDiagramm_set{set}.json",
        run:
            NUMBER_PROTEINS = 200
            from RDPMSpecIdentifier.plots import COLOR_SCHEMES, _color_to_calpha
            import pandas as pd
            import numpy as np
            from pyfunctions.plotlyVenn import venn_to_plotly
            result_df = pd.read_csv(input.tsv, sep="\t")
            if wildcards.set == "others":
                result_df = result_df[~result_df["ribosomal protein"]]
            elif wildcards.set == "ribosomes":
                result_df = result_df[result_df["ribosomal protein"]]
            elif wildcards.set == "rna_binding":
                #result_df = result_df[~result_df["ribosomal protein"]]
                result_df = result_df[result_df["GO RNA-binding"] == True]

            top200SVM = set(
                result_df[result_df["SVM RNA-binding"] == True]["old_locus_tag"]
            )
            top200RDPM = set(result_df[result_df["Rank"] <= NUMBER_PROTEINS]["old_locus_tag"])

            topRDeep = set(
                result_df[(result_df["rnase_peak_p_value"] <= 0.05) | (result_df["ctrl_peak_p_value"] <= 0.05)]["old_locus_tag"]
            )

            print(topRDeep & top200RDPM & top200SVM)
            all = topRDeep.union(top200RDPM).union(top200SVM)
            all = result_df[result_df["old_locus_tag"].isin(all)]
            all = all[["RDPMSpecID", "old_locus_tag", "Gene",  "Mean Distance", "relative fraction shift",  "ANOSIM R", "RDeeP significant", "SVM RNA-binding", "Rank", "ribosomal protein"]]
            all.to_csv(output.tsv, sep="\t", index=False)
            colors = list(COLOR_SCHEMES["Dolphin"]) + [COLOR_SCHEMES["Viking"][0]]
            colors = [_color_to_calpha(col, alpha=0.6) for col in colors]
            fig = venn_to_plotly(L_sets=(top200RDPM, top200SVM, topRDeep), L_labels=("RDPMSpecIdentifier", "TripepSVM", "RDeeP"), L_color=colors)
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


rule plotDistribution:
    input:
        file = rules.run_rdpmspecidentifier.output.json,
        joined = rules.joinAnalysisMethods.output.file
    output:
        svg = "Pipeline/plots/Distribution_ids{distributionids}.svg"
    run:
        from RDPMSpecIdentifier.datastructures import RDPMSpecData
        import pandas as pd
        import math
        from RDPMSpecIdentifier.plots import plot_protein_distributions, COLOR_SCHEMES, DEFAULT_TEMPLATE
        with open(input.file) as handle:
            data = RDPMSpecData.from_json(handle.read())
        ids = config["distributions"][wildcards.distributionids]
        if isinstance(ids, list):
            rdpmsids = ids
            tcol = "Gene"
        else:
            df = pd.read_csv(input.joined, sep="\t")
            rdpmsids = df[df[ids] == True]["RDPMSpecID"]
            tcol = "old_locus_tag"
        rows = math.ceil(len(rdpmsids)/2)
        fig = plot_protein_distributions(rdpmsids, rdpmsdata=data, colors=COLOR_SCHEMES["Dolphin"], title_col=tcol, cols=2, rows=rows, vertical_spacing=0.01, horizontal_spacing=0.075)
        fig.update_layout(template=DEFAULT_TEMPLATE, legend2=dict(y=1.11), legend1=dict(y=1.02))
        fig.update_yaxes(mirror=True, tickfont=dict(size=12))
        fig.update_xaxes(mirror=True, dtick=1, tickfont=dict(size=12))
        fig.update_layout(
            font=dict(size=16),
        )
        fig.update_layout(width=1200, height=150 * rows)

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
                font=dict(size=22),
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
                font=dict(size=16),
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
        multi_fig.update_layout(font=dict(size=16))
        multi_fig.update_layout(
            legend={'itemsizing': 'trace',"font": dict(size=18), "orientation": "h", "yanchor": "bottom", "y": 1.01},
            margin=dict(b=0,l=10,pad=0,r=10,t=40),
            width=800, height=400
        )

        multi_fig.show()
        multi_fig.write_image(output.svg)





