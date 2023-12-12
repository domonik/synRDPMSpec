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

        design = {"Name": names, "RNase": rnases, "Fraction": fractions, "Replicate": replicates}
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
