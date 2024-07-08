import os

ORGANISMID = f"{config['genus']}.{config['species']}_taxid{config['tax_id']}"


rule multiplyData:
    input:
        intensities = config["intensities"],
        design = config["design"]
    output:
        multipleIntensities = "Pipeline/BenchmarkData/synIntensitiesFake.tsv",
        multipleDesign = "Pipeline/BenchmarkData/synDesignFake.tsv",
    run:
        import pandas as pd
        intensities = pd.read_csv(input.intensities, sep="\t")
        design = pd.read_csv(input.design, sep="\t")

        names = []
        rnases = []
        fractions = []
        replicates = []
        newcols = {"index": intensities.iloc[:, 0]}
        for x in intensities.columns:
            if x.startswith("ctrl") or x.startswith("rnase"):
                d = x.split("-")
                print(d)
                frac = int(d[1])
                d = d[0]

                if d.startswith("ctrl"):
                    RNase = False
                else:
                    RNase = True
                for nr in range(3):
                    br = d[-1] + f".{nr}"
                    fakedata = intensities[x].sample(frac=1).values
                    name = x + f".{nr}"

                    names.append(name)
                    rnases.append(RNase)
                    fractions.append(frac)
                    replicates.append(br)
                    newcols[name] = fakedata

        design = {"Name": names, "Treatment": rnases, "Fraction": fractions, "Replicate": replicates}
        df = pd.DataFrame(design)
        df.to_csv(output.multipleDesign, sep="\t", index=False)
        intensities = pd.DataFrame(newcols)
        intensities.to_csv(output.multipleIntensities, sep="\t", index=False)



rule processGOTerms:
    input:
        go_terms = "Data/GOAnno.tsv"
    output:
        go_tsv = "Pipeline/GO/expandedGO.tsv",
        symbols = "Pipeline/GO/GOfakeSymbols.tsv",
    run:
        import pandas as pd
        import numpy as np
        df = pd.read_csv(input.go_terms, sep="\t")
        df = df[["Gene Names (ordered locus)", "Gene Ontology IDs"]]
        df = df.rename({"Gene Names (ordered locus)": "old_locus_tag"}, axis=1)
        symbols = df[["old_locus_tag"]]
        symbols = symbols.drop_duplicates()
        symbols["random"] = np.arange(0, len(symbols))


        df["GOTerm"] = df["Gene Ontology IDs"].str.split("; ")
        df = df.explode("GOTerm")[["old_locus_tag", "GOTerm"]]
        df = df.dropna()
        df.to_csv(output.go_tsv, sep="\t", index=False)
        symbols.to_csv(output.symbols, sep="\t", index=False)

rule tmpbuildGOTable:
    input:
        go_table = rules.processGOTerms.output.go_tsv
    output:
        all_terms = "Pipeline/GO/expandedTerms.tsv"
    conda: "../envs/annotation_forge.yml"
    script:
        "../Rscripts/getFullGOTerms.R"


rule generateOrgDB:
    input:
        symbols = rules.processGOTerms.output.symbols,
        go_terms = rules.processGOTerms.output.go_tsv
    params:
        species=config["species"],
        genus=config["genus"]
    conda:
        "../envs/annotation_forge.yml"
    output:
        annotation_db = directory(
            os.path.join("Pipeline/AnnotationDBs", ORGANISMID)
        ),
        finished_file = temporary(
            os.path.join( "Pipeline/AnnotationDBs/finished_" + ORGANISMID)
        )
    script:
        "../Rscripts/createAnnotation.R"



rule extractGORNABinding:
    input:
        file = rules.tmpbuildGOTable.output.all_terms
    output:
        file = "Pipeline/GO/RNABindingGO.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.file, sep="\t")

        binders = df[df["GOTerm"] == "GO:0003723"]
        df["GO RNA-binding"] = df["old_locus_tag"].isin(binders["old_locus_tag"])
        df = df.drop_duplicates("old_locus_tag")
        df = df[["old_locus_tag", "GO RNA-binding"]]
        df.to_csv(output.file, sep="\t", index=False)


rule prepareSalmonellaData:
    input:
        file = "Data/salmonella_table.tsv"
    output:
        file = "Pipeline/Salmonella/preparedIntensities.tsv",
        file2 = "Pipeline/Salmonella/preparedDesing.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.file, sep="\t")
        df = df[["UniprotID", "gene"] + list(df.columns[df.columns.str.contains("norm_to_spike-in fraction")])]
        df = df.rename({"norm_to_spike-in fraction PR": "norm_to_spike-in fraction 21R"}, axis=1)
        df = df.rename({"norm_to_spike-in fraction P": "norm_to_spike-in fraction 21"}, axis=1)
        df = df[~df["UniprotID"].str.contains("_HUMAN_")]

        design = {"Name": [], "Treatment": [], "Fraction": [], "Replicate": []}
        for orig_col in df.columns[2:]:
            ncol = orig_col
            if ncol[-1] == "R":
                rnase = "RNase"
                ncol = ncol[0:-1]
            else:
                rnase = "Control"
            fraction = int(ncol.split(" ")[-1])
            rep = 1
            design["Name"].append(orig_col)
            design["Treatment"].append(rnase)
            design["Fraction"].append(fraction)
            design["Replicate"].append(rep)
        design = pd.DataFrame(design)
        df.to_csv(output.file, sep="\t", index=False)
        design.to_csv(output.file2, sep="\t", index=False)


rule prepareinitialData:
    input:
        file = config["intensities"],
        design = config["design"],
        uniprot = "Data/GOAnno.tsv",
        go_binders = rules.extractGORNABinding.output.file
    output:
        sanitized_df = "Pipeline/sanitized_df.tsv",
        design = "Pipeline/sanitized_design.tsv"
    run:
        import pandas as pd
        go_table = pd.read_csv(input.go_binders, sep="\t")
        design = pd.read_csv(input.design, sep="\t")

        rapdor_table = pd.read_csv(input.file, sep="\t")
        uniprot_anno = pd.read_csv(input.uniprot, sep="\t")
        uniprot_anno = uniprot_anno.rename({"Gene Names (ordered locus)": "old_locus_tag"}, axis=1)
        uniprot_anno = uniprot_anno[["Gene Names", "old_locus_tag"]]
        uniprot_anno["Gene Names"] = uniprot_anno["Gene Names"].str.split(" ").str[0]
        uniprot_anno = uniprot_anno[~uniprot_anno["old_locus_tag"].isna()]
        rapdor_table = rapdor_table.merge(uniprot_anno, on="old_locus_tag", how="left")

        rapdor_table['Gene'] = rapdor_table['Gene Names'].fillna(rapdor_table['Gene'])
        def capitalize_string(s):
            if pd.isna(s):
                return s
            else:
                return s[0].upper() + s[1:3] + s[3:].upper()
        rapdor_table["Protein name"] = rapdor_table["Gene"].apply(capitalize_string)
        rapdor_table["small_ribosomal"] = rapdor_table["Gene"].str.contains(r'\b(rps|Rps)(?![1-9][A-Za-z])',case=False)
        rapdor_table.loc[rapdor_table["Gene"].str.contains(r'rps1A|rps1b') == True, "small_ribosomal"] = False
        rapdor_table["large_ribosomal"] = rapdor_table["Gene"].str.contains(r'\b(rpl|Rpl)', case=False)
        rapdor_table["photosystem"] = rapdor_table["Gene"].str.contains(r'\b(psb|psa)', case=False)
        poly_tags = r"|".join([
            "ssl2982", "sll1789", "sll1689", "slr0653", "sll0722", "sll0184", "slr1564", "slr1545", "slr1861", "ssr1600", "slr1859"
        ])
        rapdor_table["RNA ploymerase"] = (rapdor_table["Gene"].str.contains(r'\b(Rpo)', case=False) | rapdor_table["old_locus_tag"].str.contains(poly_tags, case=False))
        rapdor_table["ribosomal protein"] = (rapdor_table["large_ribosomal"] | rapdor_table["small_ribosomal"])
        rapdor_table = rapdor_table.merge(go_table, on="old_locus_tag", how="left")
        rapdor_table["ribosomal protein"] = ((rapdor_table["Gene"].str.contains('rpl|rps|Rpl|Rps', case=False)) | (rapdor_table['ProteinFunction'].str.contains('ribosomal protein', case=False)))
        rapdor_table.loc[rapdor_table["Gene"].str.contains(r'rps1A|rps1b') == True, "ribosomal protein"] = False
        # rapdor_table = rapdor_table.drop([f"rnase{idx}-20" for idx in range(1, 4)], axis=1)
        # rapdor_table = rapdor_table.drop([f"ctrl{idx}-20" for idx in range(1, 4)], axis=1)
        # design = design[~design["Name"].isin([f"ctrl{idx}-20" for idx in range(1, 4)])]
        # design = design[~design["Name"].isin([f"rnase{idx}-20" for idx in range(1, 4)])]

        rapdor_table = rapdor_table[~rapdor_table["id"].isin([1714, 1743, 0])]

        rapdor_table.to_csv(output.sanitized_df, sep="\t", index=False)
        design.to_csv(output.design, sep="\t", index=False)



rule prepareSynConditionedMS:
    input:
        membrane = "Data/MembraneProteins.csv",
        soluble = "Data/SolubleProteins.csv"
    output:
        intensities = "Pipeline/ConditionedSynechocystis/csv/{condition}_intensities.tsv",
        design = "Pipeline/ConditionedSynechocystis/csv/{condition}_design.tsv",
    run:
        import pandas as pd
        def multi_delimiter_split(string, delimiters):
            for delimiter in delimiters:
                string = " ".join(string.split(delimiter))

            result = string.split()
            return result
        membrane = pd.read_csv(input.membrane, sep="\t")
        membrane = membrane.rename(columns=lambda x: x.replace('Fe', 'Fe-'))


        membrane = membrane[["Accession", "Gene_symbol"] + list(membrane.columns[membrane.columns.str.contains(f"Ctrl1_|{wildcards.condition}-|{wildcards.condition}_")])]
        membrane.columns = ["Accession", "Gene_symbol"] + [f"{col}_Membrane" for col in membrane.columns[2:]]
        membrane = membrane.dropna()

        soluble = pd.read_csv(input.soluble, sep="\t")
        soluble = soluble[["Accession", "Gene_symbol"] + list(soluble.columns[soluble.columns.str.contains(f"Ctrl1_|{wildcards.condition}-|{wildcards.condition}_")])]
        soluble.columns = ["Accession", "Gene_symbol"] + [f"{col}_Soluble" for col in soluble.columns[2:]]
        if wildcards.condition == "HEAT":
            soluble = soluble.drop("HEAT_3_Soluble", axis=1)
        soluble = soluble.dropna()


        df = soluble.merge(membrane, on=["Accession", "Gene_symbol"])
        design = {"Name": [], "Treatment": [], "Fraction": [], "Replicate": []}
        for ncol in df.columns[2:]:
            rnase, replicate, fraction = multi_delimiter_split(ncol, ["-", "_"])

            design["Treatment"].append(wildcards.condition if rnase == wildcards.condition else "Control")
            design["Replicate"].append(replicate)
            design["Fraction"].append(fraction)
            design["Name"].append(ncol)
        design = pd.DataFrame(design)
        df = df.rename({"Gene_symbol": "Gene"}, axis=1)
        df.to_csv(output.intensities,sep="\t",index=False)
        design.to_csv(output.design ,sep="\t",index=False)

rule downloadNatureData:
    output:
        muscle="Pipeline/NatureMouse/RawData/{experiment}.xlsx",
    params:
        link = lambda wildcards: config["Mouse"][wildcards.experiment]["link"]
    shell:
        """
        wget {params.link} -O {output.muscle}
        """

def split_muscle(df):
    names, treatments, replicates, fractions = ([] for _ in range(4))
    for col in df.columns[1:]:
        ncol = col.split("_")
        names.append(col)
        replicates.append(ncol[-2])
        fractions.append(ncol[-1])
        treatments.append("Contracted" if len(ncol) == 4 else "Control")
        df[col][df[col] == 0] = np.nan
    return names, treatments, replicates, fractions, df

def split_kidney(df):
    names, treatments, replicates, fractions = ([] for _ in range(4))
    for col in df.columns[1:]:
        ncol = col.split("_")
        names.append(col)
        replicates.append(ncol[2])
        fractions.append(ncol[1])
        treatments.append(ncol[0])
        df[col][df[col] == 0] = np.nan
    return names, treatments, replicates, fractions, df

def split_liver(df):
    df = df.rename({"CTRL_FR1_REP3_01_20201002195754": "CTRL_FR1_REP4_01"}, axis=1)
    return *split_kidney(df),

def split_min(min):
    def split_8min(df):
        names, treatments, replicates, fractions = ([] for _ in range(4))
        df = df.loc[:, df.columns.str.contains(f"{min}|Gene.name|CTRL")]
        for col in df.columns[1:]:
            ncol = col.split("_")
            names.append(col)
            replicates.append(ncol[-1])
            fractions.append(ncol[-2])
            treatments.append("EGF" if ncol[1] == f"{min}" else "CTRL")
            df[col][df[col] == 0] = np.nan
        return names, treatments, replicates, fractions, df
    return split_8min



SPLITFCTS = {
    "muscle": split_muscle,
    "egf_kidney": split_kidney,
    "egf_liver": split_liver,
    "egf_8min": split_min("8min"),
    "egf_2min": split_min("2min"),
    "egf_20min": split_min("20min"),
    "egf_90min": split_min("90min"),
}

RENAMEFRCT = {
    "FR1": "Cytosolic (FR1)",
    "FR2": "Cytosolic (FR2)",
    "FR3": "Membrane-bound organelles (FR3)",
    "FR4": "Membrane-bound organelles (FR4)",
    "FR5": "Nucleus & Nucleolus (FR5)",
    "FR6": "Nucleus & Nucleolus (FR6)",
}

rule AnalyzeNatureWithRAPDOR:
    input:
        xlsx = rules.downloadNatureData.output.muscle
    output:
        design = "Pipeline/NatureMouse/prepared/{experiment}_design.tsv",
        df = "Pipeline/NatureMouse/prepared/{experiment}_intensities.tsv",
        tsv= "Pipeline/NatureMouse/analyzed/{experiment}_intensities_analyzed.tsv",
        json = "Pipeline/NatureMouse/analyzed/{experiment}_intensities.json"
    run:
        import pandas as pd
        from RAPDOR.datastructures import RAPDORData
        cfg = config["Mouse"][wildcards.experiment]
        df = pd.read_excel(input.xlsx, sheet_name=cfg["data_sheet"])
        addinfo = pd.read_excel(input.xlsx, sheet_name=cfg["col_sheet"])
        addinfo = addinfo[["PG.Genes", "PG.ProteinDescriptions", "PG.ProteinGroups"]]
        addinfo["PG.ProteinGroups"] = addinfo["PG.ProteinGroups"].str.split(";").str[0]
        split_fct = SPLITFCTS[wildcards.experiment]
        names, treatments, replicates, fractions, df = split_fct(df)
        fractions = [RENAMEFRCT[frac] for frac in fractions]
        design = pd.DataFrame(
            {"Name": names, "Treatment": treatments, "Fraction": fractions, "Replicate": replicates}
        )
        design.to_csv(output.design, sep="\t", index=False)
        if "min" in wildcards.experiment and "egf" in wildcards.experiment:
            addinfo["PG.Genes"] = addinfo["PG.Genes"].str.split(";").str[0]
            df = df.merge(addinfo, left_on="Gene.names", how="left", right_on="PG.Genes")

        else:
            df = df.merge(addinfo, on="PG.Genes", how="left")

        print(design)
        print(df)

        df.to_csv(output.df, sep="\t", index=False)
        data = RAPDORData(df = df, design=design, logbase=2)
        data.normalize_and_get_distances(method="Jensen-Shannon-Distance")

        data.calc_all_scores()
        data.rank_table(["ANOSIM R", "Mean Distance"], ascending=(False, False))
        mean1 = np.nansum(np.nanmean(data.array[:, data.indices[0]],axis=-2),axis=-1)
        mean2 = np.nansum(np.nanmean(data.array[:, data.indices[1]],axis=-2),axis=-1)

        data.df["rawDiff"] = np.abs(mean1 - mean2)

        data.to_json(output.json)
        data.export_csv(output.tsv ,sep="\t")

