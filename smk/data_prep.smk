import os

ORGANISMID = f"{config['genus']}.{config['species']}_taxid{config['tax_id']}"


rule multiplyData:
    input:
        intensities = "Data/synIntensities.tsv",
        design = "Data/synDesign.tsv"
    output:
        multipleIntensities = "Pipeline/BenchmarkData/synIntensitiesFake.tsv",
        multipleDesign = "Pipeline/BenchmarkData/synDesignFake.tsv",
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
            col = orig_col
            if col[-1] == "R":
                rnase = True
                col = col[0:-1]
            else:
                rnase = False
            fraction = int(col.split(" ")[-1])
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
        file = "Data/synIntensities.tsv",
        uniprot = "Data/GOAnno.tsv",
        go_binders = rules.extractGORNABinding.output.file
    output:
        sanitized_df = "Pipeline/sanitized_df.tsv"
    run:
        import pandas as pd
        go_table = pd.read_csv(input.go_binders, sep="\t")

        rapdor_table = pd.read_csv(input.file, sep="\t")
        uniprot_anno = pd.read_csv(input.uniprot, sep="\t")
        uniprot_anno = uniprot_anno.rename({"Gene Names (ordered locus)": "old_locus_tag"}, axis=1)
        uniprot_anno = uniprot_anno[["Gene Names", "old_locus_tag"]]
        uniprot_anno["Gene Names"] = uniprot_anno["Gene Names"].str.split(" ").str[0]
        uniprot_anno = uniprot_anno[~uniprot_anno["old_locus_tag"].isna()]
        rapdor_table = rapdor_table.merge(uniprot_anno, on="old_locus_tag", how="left")

        rapdor_table['Gene'] = rapdor_table['Gene Names'].fillna(rapdor_table['Gene'])
        rapdor_table["small_ribosomal"] = rapdor_table["Gene"].str.contains(r'\b(rps|Rps)(?![1-9][A-Za-z])',case=False)
        rapdor_table["large_ribosomal"] = rapdor_table["Gene"].str.contains(r'\b(rpl|Rpl)', case=False)
        rapdor_table["photosystem"] = rapdor_table["Gene"].str.contains(r'\b(psb|psa)', case=False)
        rapdor_table["ribosomal protein"] = (rapdor_table["large_ribosomal"] | rapdor_table["small_ribosomal"])
        rapdor_table = rapdor_table.merge(go_table, on="old_locus_tag", how="left")

        rapdor_table["ribosomal protein"] = ((rapdor_table["Gene"].str.contains('rpl|rps|Rpl|Rps', case=False)) | (rapdor_table['ProteinFunction'].str.contains('ribosomal protein', case=False)))
        rapdor_table.to_csv(output.sanitized_df, sep="\t", index=False)



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
        for col in df.columns[2:]:
            rnase, replicate, fraction = multi_delimiter_split(col, ["-", "_"])

            design["Treatment"].append(wildcards.condition if rnase == wildcards.condition else "Control")
            design["Replicate"].append(replicate)
            design["Fraction"].append(fraction)
            design["Name"].append(col)
        design = pd.DataFrame(design)
        df = df.rename({"Gene_symbol": "Gene"}, axis=1)
        df.to_csv(output.intensities,sep="\t",index=False)
        design.to_csv(output.design ,sep="\t",index=False)

