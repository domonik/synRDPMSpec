
include: "smk/rdeep.smk"
include: "smk/rapdor.smk"
include: "smk/data_prep.smk"
include: "smk/postProcessingAndPlots.smk"

rule all:
    input:
        rdeep =  rules.runRDeep.output,
        rdpmspec = rules.plotMeanDistribution.output,
        p = rules.plotBarcodePlot.output,

        #joined = expand(rules.plotVennDiagramm.output, ribo=[False, True, "only"]),
        #allvenns = rules.plotAllVenns.output,
        #distribution = expand(rules.plotDistribution.output, distributionids=config["distributions"].keys()),
        bubble_plot = expand(rules.createBubblePlot.output, highlight=["overlapping"]),
        #rna_binding = rules.extractGORNABinding.output,
        #enrich = rules.GOEnrichment.output,
        #benchmark = rules.plotRuntime.output,
        #salmonella = rules.runIdentifierOnSalmonella.output,
        syn_cond = expand(rules.runOnSynData.output, condition=["COLD", "HEAT", "DARK", "N", "Fe"]),
        overlapping_data = rules.extract_tophits.output,
        fig4 = rules.createFigure4.output,
        table1 = rules.createTable1andTableS1.output,
        figs2 = rules.plotConditionedSynechochoColdRibo.output,
        rdeep_rapdor = rules.runRAPDORonRDeePNorm.output,
        rpla = rules.combineFigure3.output,
        fig = rules.plotSpearmans.output,
        ts1 = rules.detectedProteinTable.output







