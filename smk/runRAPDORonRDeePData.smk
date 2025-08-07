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
        df = df[["Entry", "Gene Ontology IDs"]]
        df = df.rename({"Entry": "old_locus_tag"}, axis=1)
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
        go = rules.expandHumanGOTerms.output.all_terms,
        uniprot = "Data/GOAnnoHuman.tsv"
    output:
        mapping = "Pipeline/RAPDORonRDeeP/ProteinsWithRNABinding.tsv",
    run:
        binding = pd.read_csv(input.go,sep="\t")
        binding = binding[binding["GOTerm"] == "GO:0003723"]
        binding = set(binding["old_locus_tag"].unique().tolist())
        uniprot = pd.read_csv(input.uniprot,sep="\t")
        string = pd.read_csv(input.string, sep=" ", compression="gzip")
        string = string[string["combined_score"] >= 700]
        mapping = uniprot[["Entry", "STRING"]]
        mapping["STRING"] = mapping["STRING"].str.split(";").str[0]
        string = string.merge(mapping, left_on='protein1', right_on='STRING', )
        string = string.merge(mapping, left_on='protein2', right_on='STRING', suffixes=("_protein1", "_protein2"))

        mapping["RNA binding"] = mapping["Entry"].isin(binding)


        string = string[(string["Entry_protein1"].isin(binding)) | (string["Entry_protein2"].isin(binding))]
        binding_complex = set(string["Entry_protein1"]).union(set(string["Entry_protein2"]))
        mapping["RNA binding complex"] = mapping["Entry"].isin(binding_complex)
        mapping["RNA binding or complex"] = mapping["RNA binding complex"] | mapping["RNA binding"]
        mapping.to_csv(output.mapping, index=False, sep="\t")


rule fixDataLayout:
    input:
        tsv = "Data/proteinGroups.txt",
    output:
        intensities = "Pipeline/RAPDORonRDeeP/intensities.tsv",
        design = "Pipeline/RAPDORonRDeeP/design.tsv",
        rdeep_input="Pipeline/RDeePOriginal/input.csv"
    run:
        df = pd.read_csv(input.tsv, sep="\t")
        df["Protein"] = df["Fasta headers"].str.split(" ").str[0].str.split("|").str[-1]
        pattern = r'^Reporter intensity corrected [1-6] (F[1-9]|1[0-9]|2[0-5])$'

        data_cols = [col for col in df.columns if "Reporter intensity corrected " in col]
        df.loc[df["Protein"].isna(), "Protein"] = df.loc[df["Protein"].isna(), "Protein IDs"]
        df = df[df["Razor + unique peptides"] >= 2]
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
        rnase_cols = []
        ctrl_cols = []
        for col in df[data_cols]:
            frac = int(col.split("F")[-1])
            rep = col.split("Reporter intensity corrected ")[-1].split(" ")[0]
            treatment = "CTRL" if int(rep) in [1, 5, 3] else "RNase"
            rep = repmap[int(rep)]
            rdeep_name = f"{treatment.lower()}{rep}-{str(frac).zfill(2)}"
            df = df.rename({col: rdeep_name}, axis=1)
            if treatment == "CTRL":
                ctrl_cols.append(rdeep_name)
            else:
                rnase_cols.append(rdeep_name)
            data["Name"].append(rdeep_name)
            data["Treatment"].append(treatment)
            data["Fraction"].append(frac)
            data["Replicate"].append(rep)
        design = pd.DataFrame(data)
        design.to_csv(output.design, index=False, sep="\t")

        df = df[["Protein", "Protein IDs", "Majority protein IDs"] + data["Name"]]
        df = df[(df[rnase_cols].sum(axis=1) != 0) & (df[ctrl_cols].sum(axis=1) != 0)]
        rdeep_df = df.rename({"Protein": "name"}, axis=1).drop(["Protein IDs", "Majority protein IDs"], axis=1)

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
        outfile2 = "Pipeline/RAPDORonRDeeP/OriginalRDeePAnalysisCTRLTable.csv",
        outfile3 = "Pipeline/RAPDORonRDeeP/OriginalRDeePAnalysisRNaseTable.csv",
        normalized_counts = "Pipeline/RAPDORonRDeeP/OriginalRDeePNormalizedCounts.tsv",
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/RDeePOriginalBenchmark.csv", config["benchmark_repeats"])
    conda: "../envs/rdeep.yml"
    script: "../Rscripts/OriginalRDeeP.R"


rule postProcessNormCounts:
    input:
        tsv = rules.originalRDeeP.output.normalized_counts,
        rapdor = rules.fixDataLayout.output.intensities,
        mapping = rules.processRNABinding.output.mapping
    output:
        intensities = "Pipeline/RAPDORonRDeeP/IntensitiesNorm.tsv",
    run:
        df = pd.read_csv(input.tsv, sep="\t")
        df.columns = df.columns.str.replace(".", "-")
        rapdor_df = pd.read_csv(input.rapdor, sep="\t")
        mapping = pd.read_csv(input.mapping, sep="\t")
        rapdor_df = rapdor_df.merge(df, left_on="Protein", right_index=True, suffixes=("_delete", ""))
        rapdor_df = rapdor_df[[col for col in rapdor_df.columns if  "_delete" not in col]]
        rapdor_df["Protein ID"] = rapdor_df["Majority protein IDs"].str.split(";").str[0]
        rapdor_df = rapdor_df.merge(mapping, left_on="Protein ID", right_on="Entry")
        print(rapdor_df.head())
        print(rapdor_df.columns)
        rapdor_df.to_csv(output.intensities, sep="\t", index=False)



rule runRAPDORonRDeeP:
    input:
        intensities = rules.postProcessNormCounts.output.intensities,
        design = rules.fixDataLayout.output.design,
        rdeeP_data = rules.originalRDeeP.output.outfile
    output:
        json = "Pipeline/Paper/Supplementary/JSON/RDeePRAPDOR.json",
        tsv = "Pipeline/RAPDORonRDeeP/RDeePRAPDOR.tsv",
        RDistribution = "Pipeline/RAPDORonRDeeP/RDeePRAPDORDistribution.npy",
    threads: 999
    benchmark:
        repeat("Pipeline/benchmarks/RAPDORonRDeePBenchmark.csv", config["benchmark_repeats"])
    run:
        from RAPDOR.datastructures import RAPDORData
        import pandas as pd
        from multiprocessing import set_start_method
        set_start_method("spawn", force=True)
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
        _, distribution = rbpmdata.calc_anosim_p_value(permutations=-1, mode="global", threads=1)
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
        rapdor_json = rules.runRAPDORonRDeeP.output.json,
    output:
        file = "Pipeline/RAPDORonRDeeP/RDeePRAPDORJoined.tsv",
        json = "Pipeline/RAPDORonRDeeP/RDeePRAPDORJoined.json"
    run:
        import numpy as np
        from RAPDOR.datastructures import RAPDORData
        df = pd.read_csv(input.file, sep=" ")
        binding_df = pd.read_csv(input.rna_binding, sep="\t")
        rapdor_df = pd.read_csv(input.rapdor, sep="\t")
        print(df)
        df["maxpval"] = df[["rnase_peak_p_value", "ctrl_peak_p_value"]].max(axis=1)
        df["end"] = (df["ctrl_peak_amount_loss"] <= 0) | (df["rnase_peak_amount_gain"] <= 0) | (df["dist"].abs() <= 1)
        df = df.sort_values(["end", "maxpval", "rnase_peak_amount_gain", "ctrl_peak_amount_loss",  "right_shift"], ascending=[True, True, False, False,  True])
        df = df.drop_duplicates(subset=["protein_name"])
        df["RDeepRank"] = np.arange(0, df.shape[0])
        print(rapdor_df)
        df = pd.merge(rapdor_df, df, how="left", left_on="RAPDORid", right_on="protein_name")
        df.to_csv(output.file, sep="\t", index=False)
        add_cols = df[["RAPDORid","maxpval", "RDeepRank", "pb_fit", "rnase_peak", "ctrl_peak"]]
        data = RAPDORData.from_file(input.rapdor_json)
        data.df = data.df.merge(add_cols, on="RAPDORid")
        data.to_json(output.json)


rule plotComparisonExample:
    input:
        rapdor = rules.sortAndRankRDeep.output.json
    output:
        svg = "Pipeline/RAPDORonRDeeP/Plots/ComparisonExample.svg",
        html = "Pipeline/RAPDORonRDeeP/Plots/ComparisonExample.html",
        source_data = "Pipeline/Paper/SourceData/F5C.tsv",

    run:
        from RAPDOR.datastructures import RAPDORData
        from RAPDOR.plots import plot_protein_distributions, COLOR_SCHEMES
        colors = COLOR_SCHEMES["Dolphin"]
        data = RAPDORData.from_file(input.rapdor)
        data.df["Protein"] = data.df["RAPDORid"].str.split("_HUMAN").str[0]
        ids = ["NHP2_HUMAN", "EF2_HUMAN", "TIM50_HUMAN"]
        source_data = data.df.loc[data.df["RAPDORid"].isin(ids), ["Protein"]]
        indices = source_data.index

        source_data = source_data.loc[source_data.index.repeat(6)].reset_index(drop=True)
        source_data["type"] = (["Control"] * 3 + ["RNase"] * 3)  * 3

        subdata = data.norm_array[indices]
        rsubdata = subdata.reshape(3 * 6,-1)

        source_data.loc[:, [f"rel. Intensity Fraction {x+2}" for x in range(rsubdata.shape[1])]] = rsubdata
        source_data.to_csv(output.source_data, sep="\t", index=False)

        fig = plot_protein_distributions(
            ids,
            data,
            colors=COLOR_SCHEMES["Dolphin"],
            title_col="Protein",
            vertical_spacing=0.05
        )
        for idx, name in enumerate(ids):
            sdf = data.df[data.df["RAPDORid"] == name]
            if name != "TIM50_HUMAN":
                fig.add_vline(
                    x=sdf["rnase_peak"].iloc[0],
                    line=dict(color=colors[1]),
                    row=idx+1,
                    col=1
                )
            fig.add_vline(
                x=sdf["ctrl_peak"].iloc[0],
                line=dict(color=colors[0]),
                row=idx+1,
                col=1
            )
        fig.update_layout(
            template=DEFAULT_TEMPLATE,
            width=config["width"],
            height=config["FR1"]["C"]["height"],
        )
        fig.update_annotations(font=config["fonts"]["annotations"])
        fig.update_legends(
            font=config["fonts"]["legend"]
        )
        fig.write_image(
            output.svg
        )
        fig.write_html(output.html)


rule plotRDeePDataVennDiagram:
    input:
        rapdorfile = rules.sortAndRankRDeep.output.file
    output:
        json_all = "Pipeline/RAPDORonRDeeP/Plots/VennDiagramAll.json",
        json_rbp = "Pipeline/RAPDORonRDeeP/Plots/VennDiagramRBP.json",
        json_rdp = "Pipeline/RAPDORonRDeeP/Plots/VennDiagramRDP.json",
    run:
        from pyfunctions.plotlyVenn import venn_to_plotly
        from RAPDOR.plots import COLOR_SCHEMES
        import plotly.graph_objs as go
        import plotly.subplots as sp
        import string
        import plotly.io
        df = pd.read_csv(input.rapdorfile, sep="\t")
        rdeep = set(df[df["RDeePSignificant"]]["RAPDORid"].tolist())
        rdeep_rdp = set(df[df["RDeePSignificant"] & df["RNA binding or complex"]]["RAPDORid"].tolist())
        rdeep_rbp = set(df[df["RDeePSignificant"] & df["RNA binding"]]["RAPDORid"].tolist())
        rapdor = set(df[(df["ANOSIM R"] >= 1)]["RAPDORid"].tolist())
        rapdor_rdp = set(df[(df["ANOSIM R"] >= 1) & df["RNA binding or complex"]]["RAPDORid"].tolist())
        rapdor_rbp = set(df[(df["ANOSIM R"] >= 1) & df["RNA binding"]]["RAPDORid"].tolist())
        colors = COLOR_SCHEMES["Dolphin"]

        fig1 = venn_to_plotly(L_sets=(rapdor, rdeep), L_labels=("RAPDOR", "RDeeP"), L_color=colors)
        plotly.io.write_json(fig1, output.json_all)
        fig2 = venn_to_plotly(L_sets=(rapdor_rbp, rdeep_rbp), L_labels=("RAPDOR", "RDeeP"), L_color=colors)
        plotly.io.write_json(fig2, output.json_rbp)

        fig3 = venn_to_plotly(L_sets=(rapdor_rdp, rdeep_rdp), L_labels=("RAPDOR", "RDeeP"), L_color=colors)
        plotly.io.write_json(fig3, output.json_rdp)
        #
        # multi_fig = sp.make_subplots(rows=1,cols=3, horizontal_spacing=0.01)
        # annos = string.ascii_uppercase
        #
        # for idx, (name, fig) in enumerate([("All", fig1), ("RNA binding proteins", fig2), ("RNA dependent proteins", fig3)], 1):
        #     fig.update_layout(
        #         xaxis=dict(showgrid=False,zeroline=False,showline=False),
        #         yaxis=dict(showgrid=False,zeroline=False,showline=False)
        #     )
        #     for iidx, shape in enumerate(fig['layout']['shapes']):
        #         shape["xref"] = f"x{idx}" if idx > 1 else "x"
        #         shape["yref"] = f"y{idx}" if idx > 1 else "y"
        #         multi_fig.add_shape(shape)
        #         color = shape["fillcolor"]
        #         if idx == 1 and color != "rgba(255,0,0,0)":
        #             multi_fig.add_trace(
        #                 go.Scatter(
        #                     x=[None],
        #                     y=[None],
        #                     mode="markers",
        #                     name=fig['layout']['annotations'][iidx]["text"],
        #                     marker=dict(
        #                         color=color,
        #                         size=900,
        #                         # line=dict(
        #                         #     color="black",
        #                         #     width=3
        #                         # ),
        #                     ),
        #                     showlegend=True
        #
        #                 )
        #             )
        #     for annotation in fig['layout']['annotations'][2:]:
        #         if annotation["xref"] == "x":
        #             annotation["xref"] = f"x{idx}" if idx > 1 else "x"
        #             annotation["yref"] = f"y{idx}" if idx > 1 else "y"
        #         else:
        #             annotation["xref"] = f"x{idx} domain" if idx > 1 else "x domain"
        #             annotation["yref"] = f"y{idx} domain" if idx > 1 else "y domain"
        #         multi_fig.add_annotation(annotation)
        #
        #     multi_fig.add_annotation(
        #         text=f"{name}",
        #         xref=f"x{idx} domain" if idx > 1 else "x domain",
        #         xanchor="center",
        #         x=0.5,
        #         yref=f"y{idx} domain" if idx > 1 else "y domain",
        #         yanchor="bottom",
        #         y=0,
        #         font=dict(size=12),
        #         showarrow=False
        #
        #     )
        #
        #     multi_fig.update_xaxes(fig["layout"]["xaxis"],row=1,col=idx)
        #     multi_fig.update_yaxes(fig["layout"]["yaxis"],row=1,col=idx)
        #     multi_fig.update_yaxes(scaleanchor="x" if idx == 1 else f"x{idx}",scaleratio=1,col=idx)
        #
        # multi_fig.update_layout(template=DEFAULT_TEMPLATE)
        # multi_fig.update_layout(
        #     legend={'itemsizing': 'trace', "orientation": "h", "yanchor": "bottom", "y": 1.01},
        #     width=config["width"],height=300,font=config["fonts"]["legend"],
        #     margin=dict(r=20, l=20, b=20, t=20)
        # )
        # for annotation in multi_fig.layout.annotations:
        #     if annotation.text.startswith("<b>"):
        #         pass
        #     elif annotation.text.startswith("n:"):
        #         annotation.update(font=config["fonts"]["annotations"])
        #     else:
        #         annotation.update(font=config["fonts"]["annotations"], xanchor="center", yanchor="middle")
        #         if annotation.text == "0":
        #             if annotation.xref == "x":
        #                 annotation.update(showarrow=True, ay=-0.5, ax=0.5, axref=annotation.xref, ayref=annotation.yref, arrowcolor='black')
        #         if annotation.text == "1":
        #             if annotation.xref == "x":
        #                 annotation.update(showarrow=True, ay=-0.3, ax=0.7, axref=annotation.xref, ayref=annotation.yref, arrowcolor='black')
        #         if annotation.text == "4":
        #             annotation.update(x=annotation.x - 0.1, y=annotation.y - 0.1)
        #         if annotation.text == "8":
        #             annotation.update(x=annotation.x + 0.05)
        #         if annotation.text == "12":
        #             annotation.update(y=annotation.y + 0.05)
        # multi_fig.write_image(output.svg)



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



rule plotFigureX:
    input:
        tsv = rules.sortAndRankRDeep.output.file,
        json_all = rules.plotRDeePDataVennDiagram.output.json_all,
        json_rbp = rules.plotRDeePDataVennDiagram.output.json_rbp,
        json_rdp = rules.plotRDeePDataVennDiagram.output.json_rdp,
    output:
        df = "Pipeline/Paper/Supplementary/Tables/SupplementaryTableX.tsv",
        source_data = "Pipeline/Paper/SourceData/F5A.tsv",
        figurex = "Pipeline/Paper/Subfigures/FigureX.svg",
        supplementary_figurex = "Pipeline/Paper/Supplementary/Figures/FigureS6.svg",

    run:
        import plotly.graph_objs as go
        import plotly.io
        from plotly.express.colors import DEFAULT_PLOTLY_COLORS
        from sklearn.metrics import auc
        from RAPDOR.plots import COLOR_SCHEMES
        from plotly.subplots import make_subplots

        multi_fig = make_subplots(rows=1, cols=3, horizontal_spacing= 0)
        fig2 = make_subplots(
            rows=2, cols=3, horizontal_spacing= 0.1, vertical_spacing=0.25,
            row_heights=[0.7, 0.3],
            specs=[
                [{}, {}, {}],  # Row 1: Three separate subplots
                [{"colspan": 3, "type": "domain"}, None, None]  # Row 2: One table spanning all 3 columns
            ],
        )

        fig2.add_annotation(
            text="<b>A</b>",
            font=dict(size=config["multipanel_font_size"],family="Arial"),
            x=0 - 0.1,
            y=.99,
            xref="paper",
            yref="paper",
            xanchor="left",
            yanchor="bottom",
            showarrow=False
        )
        fig2.add_annotation(
            text="<b>B</b>",
            font=dict(size=config["multipanel_font_size"],family="Arial"),
            x=2 / 3 + 0.02,
            y=.99,
            xref="paper",
            yref="paper",
            xanchor="left",
            yanchor="bottom",
            showarrow=False

        )
        fig2.add_annotation(
            text="<b>C</b>",
            font=dict(size=config["multipanel_font_size"],family="Arial"),
            x=0 - 0.1,
            y=.25,
            xref="paper",
            yref="paper",
            xanchor="left",
            yanchor="bottom",
            showarrow=False

        )
        df = pd.read_csv(input.tsv,sep="\t")
        data = (
            ("Rank", "RAPDOR", "RNA binding or complex", COLOR_SCHEMES["Dolphin"][0]),
            ("RDeepRank", "RDeeP", "RNA binding or complex", COLOR_SCHEMES["Dolphin"][1]),
            ("Rank", "RAPDOR", "RNA binding", COLOR_SCHEMES["Dolphin"][0]),
            ("RDeepRank", "RDeeP", "RNA binding", COLOR_SCHEMES["Dolphin"][1]),
        )
        colors = list(COLOR_SCHEMES["Dolphin"]) + list(COLOR_SCHEMES["Viking"])
        cutoffs = []

        out_df = pd.DataFrame(
            {
                ("RBP", "AUROC"): [pd.NA, pd.NA],
                ("RBP", "AUPRC"): [pd.NA, pd.NA],
                ("RDP", "AUROC"): [pd.NA, pd.NA],
                ("RDP", "AUPRC"): [pd.NA, pd.NA],
            },
            index=["RAPDOR", "RDeeP"]
        )
        source_data = pd.DataFrame()
        #out_df["AUROC"] = pd.NA
        #out_df["AUPRC"] = pd.NA
        for idx, (sort_col, name, term, color) in enumerate(data):
            df = df.sort_values(sort_col)
            #iidx = int((np.nanmin(df[(df["maxpval"] > 0.05)][sort_col]) if name == "RDeeP" else df[(df["ANOSIM R"] == 1)][sort_col].max()) -1)
            #print(name, iidx)
            p = np.sum(df[term])
            n = np.sum(~(df[term].astype(bool)))

            tp = np.cumsum(df[term])
            print(name, p)
            print(name, n)
            fp = np.cumsum(~(df[term].astype(bool)))
            tpr = tp / p
            fpr = fp / n
            precision = tp / (tp + fp)
            auroc = auc(fpr,tpr)
            auprc = auc(tpr, precision)
            display_name = name + term
            lgroup = "RBP" if term == "RNA binding" else "RDP"
            out_df.loc[name, (lgroup, "AUROC")] = auroc
            out_df.loc[name, (lgroup, "AUPRC")] = auprc
            if lgroup == "RBP":
                source_data[f"{name} - False positive rate"] = fpr
                source_data[f"{name} - True positive rate"] = tpr
                multi_fig.add_trace(go.Scatter(x=[None],y=[None],mode="markers",marker=dict(color=color),name=name, showlegend=True))

                multi_fig.add_trace(go.Scatter(
                    x=fpr,
                    y=tpr,
                    mode="lines",
                    name=name + f" AUC: {auroc:.3f}",
                    #fill="tozeroy",
                    line=dict(color=color, dash=None if term == "RNA binding" else "dash"),
                    showlegend=False,

                    hovertext=df[sort_col]


                ))
                multi_fig.add_annotation(
                    text=f"AUC: {auroc:.3f}",
                    font=dict(color=color),
                    x=0.45 if name == "RAPDOR" else 0.35,
                    y=0.8 if name == "RAPDOR" else 0.5,
                    xref="x domain",
                    yref="y domain",
                    xanchor="right" if name == "RAPDOR" else "left",
                    yanchor="middle",
                    showarrow=False

                )
            if lgroup == "RDP":
                fig2.add_trace(
                    go.Scatter(
                        x=fpr,
                        y=tpr,
                        mode="lines",
                        name=name,
                        # fill="tozeroy",
                        line=dict(color=color,dash=None if term == "RNA binding" else "dash"),
                        showlegend=False,
                        hovertext=df[sort_col],

                    ), col=1, row=1
                )
            fig2.add_trace(
                go.Scatter(
                    x=tpr,
                    y=precision,
                    mode="lines",
                    name=name,
                    # fill="tozeroy",
                    line=dict(color=color, dash=None if term == "RNA binding" else "dash"),
                    showlegend=False,

                    hovertext=df[sort_col],

                ),col=2,row=1
            )
        fig2.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=dict(color="black", dash=None), name="RBP"))
        fig2.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=dict(color="black", dash="dash"), name="RDP"))
        fig2.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(color=colors[0]), name="RAPDOR"))
        fig2.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(color=colors[1]), name="RDeeP"))

        all_fig = plotly.io.read_json(input.json_all)
        rbp_fig = plotly.io.read_json(input.json_rbp)
        rdp_fig = plotly.io.read_json(input.json_rdp)
        source_data.to_csv(output.source_data, sep="\t", index=False)

        for idx, (name, fig) in enumerate([("All", all_fig), ("RBP", rbp_fig), ("RDP", rdp_fig)], 1):
            if name == "RDP":
                ref_fig = fig2
                col=3
            else:
                ref_fig = multi_fig
                col = 2 if name == "All" else 3
            fig.update_layout(
                xaxis=dict(showgrid=False,zeroline=False,showline=False),
                yaxis=dict(showgrid=False,zeroline=False,showline=False)
            )
            for iidx, shape in enumerate(fig['layout']['shapes']):
                shape["xref"] = f"x{col}"
                shape["yref"] = f"y{col}"
                ref_fig.add_shape(shape)
                color = shape["fillcolor"]

            for annotation in fig['layout']['annotations'][2:]:
                if annotation["xref"] == "x":
                    annotation["xref"] = f"x{col}"
                    annotation["yref"] = f"y{col}"
                else:
                    annotation["xref"] = f"x{col} domain"
                    annotation["yref"] = f"y{col} domain"
                ref_fig.add_annotation(annotation)

            ref_fig.add_annotation(
                text=f"{name}",
                xref=f"x{col} domain",
                xanchor="center",
                x=0.5,
                yref=f"y{col} domain",
                yanchor="bottom",
                y=0,
                font=dict(size=12),
                showarrow=False

            )

            ref_fig.update_xaxes(fig["layout"]["xaxis"],row=1,col=col)
            ref_fig.update_yaxes(fig["layout"]["yaxis"],row=1,col=col)
            ref_fig.update_yaxes(scaleanchor=f"x{col}",scaleratio=1,col=col)
        table = go.Table(
            header=dict(values=["Tool"] + ["-".join(i) for i in out_df.columns]),
            cells=dict(values=[out_df.index.tolist()] + [out_df[col].astype(float).round(2) for col in out_df.columns])
        )
        fig2.add_trace(
            table, row=2, col=1
        )


        multi_fig.update_layout(
            template=DEFAULT_TEMPLATE,
            legend=dict(
                x=0.5,
                y=0,
                xref="container",
                yref="container",
                xanchor="center",
                yanchor="bottom",
                orientation="h"

            ),
            yaxis=dict(range=[0, 1], title=dict(text="True positive rate")),
            xaxis=dict(range=[0, 1], title=dict(text="False positive rate")),
            width=config["width"],
            margin=config["margin"],
            height=config["FR1"]["AB"]["height"],
        )
        title_standoff = 3
        fig2.update_layout(
            template=DEFAULT_TEMPLATE,
            legend=dict(orientation="h", y=0.99, yref="container", yanchor="top"),
            yaxis=dict(range=[0, 1], title=dict(text="True positive rate"), title_standoff = title_standoff),
            xaxis=dict(range=[0, 1], title=dict(text="False positive rate"),),
            width=config["width"],
            margin=config["margin"],
            height=400,
        )
        fig2.update_xaxes(title="True positive rate", col=2, row=1)
        fig2.update_yaxes(title="Positive predictive value ", title_standoff = title_standoff, col=2, row=1)
        fig2.write_image(output.supplementary_figurex)
        multi_fig.write_image(output.figurex)
        out_df.to_csv(output.df, sep="\t")


rule plotFigureR1:
    input:
        ab = rules.plotFigureX.output.figurex,
        c = rules.plotComparisonExample.output.svg
    output:
        svg = "Pipeline/Paper/Figure5.svg"
    run:
        from svgutils.compose import Figure, Panel, SVG, Text

        c_y = config["FR1"]["C"]["height"]
        ges_y = c_y + config["FR1"]["AB"]["height"]
        letter_b_x = config["width"] // 3 + 50
        f = Figure("624px",f"{ges_y}px",
            Panel(
                SVG(input.ab),
                Text("A",2,15,size=config["multipanel_font_size"],weight='bold',font="Arial"),
                Text("B",letter_b_x,15,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ),
            Panel(
                SVG(input.c),
                Text("C",2,2,size=config["multipanel_font_size"],weight='bold',font="Arial")
            ).move(0,c_y),
        )
        svg_string = f.tostr()
        svg_string = svg_string.decode().replace("encoding='ASCII'","encoding='utf-8'")
        with open(output.svg,"w") as handle:
            handle.write(svg_string)

