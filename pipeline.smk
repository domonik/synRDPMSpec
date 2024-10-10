
include: "smk/rdeep.smk"
include: "smk/rapdor.smk"
include: "smk/data_prep.smk"
include: "smk/postProcessingAndPlots.smk"

rule all:
    input:
        rdeep =  rules.runRDeep.output,
        rdpmspec = rules.plotMeanDistribution.output,
        #p = rules.plotBarcodePlot.output,
        #df = rules.prepareinitialData.output,
        #bubble_plot = expand(rules.createBubblePlot.output, highlight=["overlapping"]),
        benchmark = rules.plotRuntime.output,
        #salmonella = rules.runIdentifierOnSalmonella.output,
        #syn_cond = expand(rules.runOnSynData.output, condition=["COLD", "HEAT", "DARK", "N", "Fe"]),
        #overlapping_data = rules.extract_tophits.output,
        fig5 = rules.createFigure5.output,
        fig4 = rules.plotAllVenns.output,
        table1 = rules.createTable1andTableS1.output,
        figs2 = rules.plotConditionedSynechochoColdRibo.output,
        #rdeep_rapdor = rules.runRAPDORonRDeePNorm.output,
        rpla = rules.combineFigure3.output,
        #fig = rules.plotSpearmans.output,
        ts1 = rules.detectedProteinTable.output,
        #tophits = expand(rules.plotTopHitDistributions.output, distribution=config["distributions"].keys()),
        copySFigures = rules.copySubfigures.output,
        final_data = rules.postProcessRapdorData.output,
        histo = rules.plotANOSIMRDistribution.output,
        #qq = rules.produceArtificialSoftArgMaxExample.output,
        #distinct = rules.theoreticalDistinctRValues.output,
        #naturrap = expand(rules.AnalyzeNatureWithRAPDOR.output, experiment=config["Mouse"].keys()),
        HelaPlot = rules.joinHeLaPlot.output,
        sup = rules.copyTableS2.output,
        file1 = rules.zipSupplementaryFile1.output








