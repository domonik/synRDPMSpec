import pandas as pd


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




rule downloadString:
    output:
        string = "Pipeline/STRING/stringdata.tsv",
        string_info = "Pipeline/STRING/string_info.tsv",
    shell: """
    wget -O {output.string} https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz
    wget -O {output.string_info} https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
    """





rule processHumanGOTerms:
    input:
        go_terms = "Data/GOAnnoHuman.tsv"
    output:
        go_tsv = "Pipeline/HumanGO/expandedGO.tsv",
        symbols = "Pipeline/HumanGO/GOfakeSymbols.tsv",
    run:
        import pandas as pd
        import numpy as np
        df = pd.read_csv(input.go_terms, sep="\t")
        df = df[["Entry Name", "Gene Ontology IDs"]]
        df = df.rename({"Entry Name": "old_locus_tag"}, axis=1)
        symbols = df[["old_locus_tag"]]
        symbols = symbols.drop_duplicates()
        symbols["random"] = np.arange(0, len(symbols))
        df["GOTerm"] = df["Gene Ontology IDs"].str.split("; ")
        df = df.explode("GOTerm")[["old_locus_tag", "GOTerm"]]
        df = df.dropna()
        df.to_csv(output.go_tsv, sep="\t", index=False)
        symbols.to_csv(output.symbols, sep="\t", index=False)

rule expandHumanGOTerms:
    input:
        go_table = rules.processHumanGOTerms.output.go_tsv
    output:
        all_terms = "Pipeline/HumanGO/expandedTerms.tsv"
    conda: "../envs/annotation_forge.yml"
    threads: 8
    script:
        "../Rscripts/getFullGOTerms.R"

rule processRNABinding:
    input:
        string = rules.downloadString.output.string,
        string_info = rules.downloadString.output.string_info,
        go = rules.expandHumanGOTerms.output.all_terms
    output:
        string = "Pipeline/RAPDORonRDeeP/stringProteinsWithRNABinding.tsv",
    run:
        binding = pd.read_csv(input.go,sep="\t")
        binding = binding[binding["GOTerm"] == "GO:0003723"]
        binding = set(binding["old_locus_tag"].unique().tolist()).union(set(binding["old_locus_tag"].str.split("_HUMAN").str[0].unique().tolist()))

        string_info = pd.read_csv(input.string_info, sep="\t", compression="gzip")
        string = pd.read_csv(input.string, sep=" ", compression="gzip")
        print(string.columns)

        string = string.merge(string_info, left_on='protein1', right_on='#string_protein_id').drop([ '#string_protein_id', 'protein1'] ,axis=1)
        string.rename(columns={'preferred_name': 'protein1_name'},inplace=True)
        print(string.columns)
        string = string.merge(string_info,left_on='protein2', right_on='#string_protein_id').drop(['protein2', '#string_protein_id'],axis=1)
        string.rename(columns={'preferred_name': 'protein2_name'},inplace=True)
        string = string[(string["protein1_name"].isin(binding)) | (string["protein2_name"].isin(binding))]
        string.to_csv(output.string, index=False, sep="\t")


rule fixDataLayout:
    input:
        tsv = "Data/proteinGroups.txt",
        go = rules.expandHumanGOTerms.output.all_terms
    output:
        intensities = "Pipeline/RAPDORonRDeeP/intensities.tsv",
        design = "Pipeline/RAPDORonRDeeP/design.tsv",
        rdeep_input="Pipeline/RDeePOriginal/input.csv"
    run:
        df = pd.read_csv(input.tsv, sep="\t")
        df["Protein"] = df["Fasta headers"].str.split(" ").str[0].str.split("|").str[-1]
        pattern = r'^Reporter intensity corrected [1-6] (F[1-9]|1[0-9]|2[0-5])$'

        data_cols = [col for col in df.columns if "Reporter intensity corrected " in col]
        df = df[["Protein", "Protein IDs"] + data_cols]
        df.loc[df["Protein"].isna(), "Protein"] = df.loc[df["Protein"].isna(), "Protein IDs"]
        data = {
            "Name": [],
            "Treatment": [],
            "Fraction": [],
            "Replicate":  []
        }
        repmap = {
            1: 1,
            2: 1,
            3: 2,
            4: 2,
            5: 3,
            6: 3
        }
        rdeep_df = df.rename({"Protein": "name"}, axis=1).drop("Protein IDs", axis=1)
        for col in df[data_cols]:
            frac = int(col.split("F")[-1])
            rep = col.split("Reporter intensity corrected ")[-1].split(" ")[0]
            treatment = "CTRL" if int(rep) in [1, 5, 3] else "RNase"
            rep = repmap[int(rep)]
            rdeep_name = f"{treatment.lower()}{rep}-{str(frac).zfill(2)}"
            rdeep_df = rdeep_df.rename({col: rdeep_name}, axis=1)

            data["Name"].append(col)
            data["Treatment"].append(treatment)
            data["Fraction"].append(frac)
            data["Replicate"].append(rep)
        design = pd.DataFrame(data)
        design.to_csv(output.design, index=False, sep="\t")

        binding = pd.read_csv(input.go,sep="\t")
        binding = binding[binding["GOTerm"] == "GO:0003723"]
        rna_binding = binding["old_locus_tag"].unique()
        df["RNA binding"] = df["Protein"].isin(rna_binding)
        df.to_csv(output.intensities, index=False, sep="\t")


        def sort_key(col):
            print(col)
            if col == 'name':
                return (0, 0, 0)  # Ensure 'name' is first
            parts = col.split('-')
            group = 1 if 'ctrl' in col else 2  # Sort 'ctrl' before 'rnase'
            sample = int(parts[0][-1])  # Extract sample number (e.g., 1, 2, 3)
            fraction = int(parts[1])  # Extract replicate number (e.g., 01, 02)
            return (fraction, sample,  group)

        print(rdeep_df.columns)
        # Sort columns
        sorted_columns = sorted(rdeep_df.columns, key=sort_key)
        rdeep_df = rdeep_df[sorted_columns]
        rdeep_df = rdeep_df.drop_duplicates("name")
        rdeep_df.to_csv(output.rdeep_input, sep=";", index=False)


rule originalRDeeP:
    input:
        file = rules.fixDataLayout.output.rdeep_input,
    output:
        outfile = "Pipeline/RAPDORonRDeeP/OriginalRDeePAnalysis.csv",
        normalized_counts = "Pipeline/RAPDORonRDeeP/OriginalRDeePNormalizedCounts.tsv",
    conda: "../envs/rdeep.yml"
    script: "../Rscripts/OriginalRDeeP.R"


# rule postProcessNormCounts:
#     input:
#         tsv = rules.originalRDeeP.output.normalized_counts
#     output:
#         intensities = "",
#         design =


rule runRAPDORonRDeeP:
    input:
        intensities = rules.fixDataLayout.output.intensities,
        design = rules.fixDataLayout.output.design,
        rdeeP_data = rules.originalRDeeP.output.outfile
    output:
        json = "Pipeline/RAPDORonRDeeP/RDeePRAPDOR.json",
        tsv = "Pipeline/RAPDORonRDeeP/RDeePRAPDOR.tsv",
        RDistribution = "Pipeline/RAPDORonRDeeP/RDeePRAPDORDistribution.npy",
    threads: 8
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd

        df = pd.read_csv(input.intensities,sep="\t")
        rdeep_df = pd.read_csv(input.rdeeP_data, sep=" ")
        rdeep_filter = ~(
                (rdeep_df["ctrl_peak_amount_loss"] <= 0) |
                (rdeep_df["rnase_peak_amount_gain"] <= 0) |
                (rdeep_df["dist"].abs() <= 1)
        )
        rdeep_sig = rdeep_df[rdeep_filter & (rdeep_df["rnase_peak_p_value"] < 0.05) & (rdeep_df["ctrl_peak_p_value"] < 0.05)]["protein_name"]
        df["RDeePSignificant"] = df["Protein"].isin(rdeep_sig)
        design = pd.read_csv(input.design,sep="\t")
        rbpmdata = RAPDORData(df=df, design=design)
        rbpmdata.normalize_and_get_distances(method="Jensen-Shannon-Distance", kernel=3)
        rbpmdata.calc_all_scores()
        rbpmdata.rank_table(["ANOSIM R", "Mean Distance"],ascending=(False, False))
        _, distribution = rbpmdata.calc_anosim_p_value(permutations=-1, mode="global", threads=threads)
        indices = rbpmdata.df[rbpmdata.df["min replicates per group"] == 3]["id"].to_numpy()
        distribution = distribution[:, indices].flatten()
        with open(output.RDistribution, "wb") as handle:
            np.save(handle, distribution)
        del distribution

        rbpmdata.to_json(output.json)
        rbpmdata.extra_df.to_csv(output.tsv, sep="\t", index=False)



rule sortAndRankRDeep:
    input:
        file = rules.originalRDeeP.output.outfile,
        rna_binding = rules.expandHumanGOTerms.output.all_terms,
        rapdor = rules.runRAPDORonRDeeP.output.tsv,
        string = rules.processRNABinding.output.string
    output:
        file = "Pipeline/RAPDORonRDeeP/RDeePRAPDORJoined.tsv",
    run:
        import numpy as np
        df = pd.read_csv(input.file, sep=" ")
        binding_df = pd.read_csv(input.rna_binding, sep="\t")
        rapdor_df = pd.read_csv(input.rapdor, sep="\t")
        string = pd.read_csv(input.string, sep="\t")
        binding = binding_df[binding_df["GOTerm"] == "GO:0003723"]
        na_binding = binding_df[binding_df["GOTerm"] == "GO:0003676"]
        print(df)
        df["maxpval"] = df[["rnase_peak_p_value", "ctrl_peak_p_value"]].max(axis=1)
        df["end"] = (df["ctrl_peak_amount_loss"] <= 0) | (df["rnase_peak_amount_gain"] <= 0) | (df["dist"].abs() <= 1)
        df = df.sort_values(["end", "maxpval", "rnase_peak_amount_gain", "ctrl_peak_amount_loss",  "right_shift"], ascending=[True, True, False, False,  True])
        df = df.drop_duplicates(subset=["protein_name"])
        df["RDeepRank"] = np.arange(0, df.shape[0])
        print(rapdor_df)
        df = pd.merge(rapdor_df, df, how="left", left_on="RAPDORid", right_on="protein_name")

        rna_binding = binding["old_locus_tag"].unique()
        na_binding = na_binding["old_locus_tag"].unique()
        df["RNA binding"] = df["RAPDORid"].isin(rna_binding)
        df["NA binding"] = df["RAPDORid"].isin(na_binding)
        df["RNA binding complex"] = (
                df["RAPDORid"].isin(string["protein1_name"]) |
                df["RAPDORid"].isin(string["protein2_name"]) |
                df["RAPDORid"].str.replace("_HUMAN", "").isin(string["protein1_name"]) |
                df["RAPDORid"].str.replace("_HUMAN", "").isin(string["protein2_name"])
        )
        df["RNA binding or complex"] = df["RNA binding"] | df["RNA binding complex"]
        df.to_csv(output.file, sep="\t", index=False)



rule plotRDeePDataVennDiagram:
    input:
        rapdorfile = rules.sortAndRankRDeep.output.file
    output:
        svg = "Pipeline/RAPDORonRDeeP/Plots/VennDiagramAll.svg",
        tsv = "Pipeline/RAPDORonRDeeP/Summary/Performance.tsv",
    run:
        from pyfunctions.plotlyVenn import venn_to_plotly
        from RAPDOR.plots import COLOR_SCHEMES
        import plotly.graph_objs as go
        import plotly.subplots as sp
        import string
        df = pd.read_csv(input.rapdorfile, sep="\t")
        rdeep = set(df[df["RDeePSignificant"]]["RAPDORid"].tolist())
        rdeep_correct = set(df[df["RDeePSignificant"] & df["RNA binding or complex"]]["RAPDORid"].tolist())
        rapdor = set(df[(df["ANOSIM R"] >= 0.8)]["RAPDORid"].tolist())
        rapdor_correct = set(df[(df["ANOSIM R"] >= 0.8) & df["RNA binding or complex"]]["RAPDORid"].tolist())
        colors = COLOR_SCHEMES["Dolphin"]

        rdeep_tp = len(rdeep_correct)
        rapdor_tp = len(rapdor_correct)
        ppv_rdeep = rdeep_tp / len(rdeep)
        ppv_rapdor = rapdor_tp / len(rapdor)
        data = {
            "Tool": ["RDeeP", "RAPDOR"],
            "True positives": [rdeep_tp, rapdor_tp],
            "False positives": [len(rdeep)- rdeep_tp, len(rapdor)-rapdor_tp],
            "Positive predictive value": [ppv_rdeep, ppv_rapdor]
        }
        df = pd.DataFrame(data)
        df.to_csv(output.tsv, sep="\t", index=False)

        fig1 = venn_to_plotly(L_sets=(rapdor, rdeep), L_labels=("RAPDOR", "RDeeP"), L_color=colors)
        fig2 = venn_to_plotly(L_sets=(rapdor_correct, rdeep_correct), L_labels=("RAPDOR", "RDeeP"), L_color=colors)
        multi_fig = sp.make_subplots(rows=1,cols=2)
        annos = string.ascii_uppercase

        for idx, (name, fig) in enumerate([("All", fig1), ("RNA dependent proteins", fig2)], 1):
            fig.update_layout(
                xaxis=dict(showgrid=False,zeroline=False,showline=False),
                yaxis=dict(showgrid=False,zeroline=False,showline=False)
            )
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
            for annotation in fig['layout']['annotations'][2:]:
                if annotation["xref"] == "x":
                    annotation["xref"] = f"x{idx}" if idx > 1 else "x"
                    annotation["yref"] = f"y{idx}" if idx > 1 else "y"
                else:
                    annotation["xref"] = f"x{idx} domain" if idx > 1 else "x domain"
                    annotation["yref"] = f"y{idx} domain" if idx > 1 else "y domain"
                multi_fig.add_annotation(annotation)

            multi_fig.add_annotation(
                text=f"{name}",
                xref=f"x{idx} domain" if idx > 1 else "x domain",
                xanchor="center",
                x=0.5,
                yref=f"y{idx} domain" if idx > 1 else "y domain",
                yanchor="bottom",
                y=0,
                font=dict(size=12),
                showarrow=False

            )

            multi_fig.update_xaxes(fig["layout"]["xaxis"],row=1,col=idx)
            multi_fig.update_yaxes(fig["layout"]["yaxis"],row=1,col=idx)
            multi_fig.update_yaxes(scaleanchor="x" if idx == 1 else f"x{idx}",scaleratio=1,col=idx)

        multi_fig.update_layout(template=DEFAULT_TEMPLATE)
        multi_fig.update_layout(
            legend={'itemsizing': 'trace', "orientation": "h", "yanchor": "bottom", "y": 1.01},
            width=config["width"],height=300,font=config["fonts"]["legend"]
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



rule plotRDeePRHistogram:
    input:
        rapdor = rules.runRAPDORonRDeeP.output.json
    output:
        html = "Pipeline/RAPDORonRDeeP/Plots/AnosimRDistribution.html"
    run:
        from RAPDOR.datastructures import RAPDORData
        import plotly.graph_objs as go
        rapdor_data = RAPDORData.from_file(input.rapdor)
        original_dist = rapdor_data.df["ANOSIM R"]
        fig = go.Figure()
        fig.add_trace(
            go.Histogram(
                x=original_dist,
            )
        )
        fig.write_html(output.html)




rule plotAUROC:
    input:
        tsv = rules.sortAndRankRDeep.output.file
    output:
        html = "Pipeline/RAPDORonRDeeP/Plots/AUROC.html",
        svg = "Pipeline/RAPDORonRDeeP/Plots/AUROC.svg",
    run:
        import plotly.graph_objs as go
        from plotly.express.colors import DEFAULT_PLOTLY_COLORS
        from sklearn.metrics import auc
        from RAPDOR.plots import COLOR_SCHEMES

        fig = go.Figure()
        df = pd.read_csv(input.tsv,sep="\t")
        data = (
            ("Rank", "RAPDOR", "RNA binding or complex"),
            ("RDeepRank", "RDeeP", "RNA binding or complex"),
        )
        colors = COLOR_SCHEMES["Dolphin"]
        for idx, (sort_col, name, term) in enumerate(data):
            df = df.sort_values(sort_col)
            iidx = int((np.nanmin(df[(df["maxpval"] > 0.05)][sort_col]) if name == "RDeeP" else df[(df["ANOSIM R"] == 1)][sort_col].max()) -1)
            print(name, iidx)
            p = np.sum(df[term])
            n = np.sum(~(df[term].astype(bool)))

            tp = np.cumsum(df[term])
            print(name, p)
            print(name, n)
            fp = np.cumsum(~(df[term].astype(bool)))
            tpr = tp / p
            fpr = fp / n
            a = auc(fpr,tpr)
            display_name = name + term
            print(fpr.iloc[iidx:iidx+10])

            fig.add_trace(go.Scatter(
                x=fpr,
                y=tpr,
                mode="lines",
                name=name + f" AUC: {a:.3f}",
                fill="tozeroy",
                line=dict(color=colors[idx]),
                legendgroup="AUC",
                legendgrouptitle=dict(text="AUC")


            ))
            fig.add_trace(
                go.Scatter(
                    x=[fpr.iloc[iidx]],
                    y=[tpr.iloc[iidx]],
                    marker=dict(size=10,color=colors[idx]),
                    mode="markers",
                    name=f"{name}: R=1" if name == "RAPDOR" else f"{name}: p-Value < 0.05",
                    legendgroup="Default cutoffs",
                    legendgrouptitle=dict(text="Default cutoffs")
                )
            )

        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            legend=dict(
                x=0,
                y=1,
                xref="paper",
                yref="paper",
                xanchor="left",
                yanchor="top"

            ),
            yaxis=dict(range=[0, 1], title=dict(text="True positive rate")),
            xaxis=dict(range=[0, 1], title=dict(text="False positive rate")),
            width=config["width"],
            height=400,
        )
        fig.write_html(output.html)
        fig.write_image(output.svg)
