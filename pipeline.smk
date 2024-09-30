
include: "smk/rdeep.smk"
include: "smk/rapdor.smk"
include: "smk/data_prep.smk"
include: "smk/postProcessingAndPlots.smk"

rule all:
    input:
        rdeep =  rules.runRDeep.output,
        rdpmspec = rules.plotMeanDistribution.output,
        p = rules.plotBarcodePlot.output,
        df = rules.prepareinitialData.output,

        #joined = expand(rules.plotVennDiagramm.output, ribo=[False, True, "only"]),
        #allvenns = rules.plotAllVenns.output,
        #distribution = expand(rules.plotDistribution.output, distributionids=config["distributions"].keys()),
        bubble_plot = expand(rules.createBubblePlot.output, highlight=["overlapping"]),
        #rna_binding = rules.extractGORNABinding.output,
        #enrich = rules.GOEnrichment.output,
        benchmark = rules.plotRuntime.output,
        salmonella = rules.runIdentifierOnSalmonella.output,
        syn_cond = expand(rules.runOnSynData.output, condition=["COLD", "HEAT", "DARK", "N", "Fe"]),
        overlapping_data = rules.extract_tophits.output,
        fig5 = rules.createFigure5.output,
        fig4 = rules.plotAllVenns.output,
        table1 = rules.createTable1andTableS1.output,
        figs2 = rules.plotConditionedSynechochoColdRibo.output,
        rdeep_rapdor = rules.runRAPDORonRDeePNorm.output,
        rpla = rules.combineFigure3.output,
        fig = rules.plotSpearmans.output,
        ts1 = rules.detectedProteinTable.output,
        tophits = expand(rules.plotTopHitDistributions.output, distribution=config["distributions"].keys()),
        copySFigures = rules.copySubfigures.output,
        #fit = rules.fitMultiGaussian.output
        final_data = rules.postProcessRapdorData.output,
        histo = rules.plotANOSIMRDistribution.output,
        qq = rules.produceArtificialSoftArgMaxExample.output,
        distinct = rules.theoreticalDistinctRValues.output,
        #naturego = expand(rules.GOTermEnrichmentMouseNature.output, experiment=config["Mouse"].keys()),
        #naturekegg = expand(rules.KEGGEnrichmentMouseNature.output, experiment=config["Mouse"].keys()),
        naturrap = expand(rules.AnalyzeNatureWithRAPDOR.output, experiment=config["Mouse"].keys()),
        HelaPlot = rules.joinHeLaPlot.output,
        mobility = expand(rules.calcMobilityScore.output, experiment=[f"egf_{x}min" for x in (2, 8, 20, 90)]),
        #limma = expand(rules.collectLimmaResults.output, experiment=[f"egf_{x}min" for x in (2,)]),
        limma = expand(rules.extractNatureForLimma.output, experiment=[f"egf_{x}min" for x in (2,)]),
        sup = rules.copyTableS2.output,
        file1 = rules.zipSupplementaryFile1.output








