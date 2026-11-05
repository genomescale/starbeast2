open module starbeast2 {
    requires beast.pkgmgmt;
    requires beast.base;
    requires static beast.fx;
    requires static javafx.controls;
    requires sampled.ancestors;
    requires morph.models;
    requires com.google.common;
    requires commons.math3;
    requires org.apache.commons.statistics.distribution;
    requires org.apache.commons.numbers.gamma;
    requires java.desktop;

    exports starbeast2;
    exports starbeast2.aimannotator;
    exports starbeast2.aimoperator;
    exports starbeast2.app.beauti;
    exports starbeast2.math.distributions;
    exports starbeast2.utils;

    provides beast.base.core.BEASTInterface with
        starbeast2.app.beauti.StarBeastAlignmentProvider,
        starbeast2.app.beauti.StarBeastMorphModelAlignmentProvider,
        starbeast2.math.distributions.AltInverseGamma,
        starbeast2.AllEqual,
        starbeast2.ConstantPopulations,
        starbeast2.ConstantWithGeneFlow,
        starbeast2.CoordinatedExchange,
        starbeast2.CoordinatedExponential,
        starbeast2.CoordinatedUniform,
        starbeast2.DiscreteRateCycle,
        starbeast2.DiscreteRateUniform,
        starbeast2.DummyModel,
        starbeast2.GeneTree,
        starbeast2.GeneTreeWithMigration,
        starbeast2.LinearWithConstantRoot,
        starbeast2.MinimalBranchLength,
        starbeast2.MultispeciesCoalescent,
        starbeast2.NetworkRateExchange,
        starbeast2.NodeReheight2,
        starbeast2.Overlap,
        starbeast2.PassthroughModel,
        starbeast2.RandomLocalRates,
        starbeast2.RealCycle,
        starbeast2.SAMau1999,
        starbeast2.STDirectionalLogger,
        starbeast2.SpeciesTree,
        starbeast2.SpeciesTreeLogger,
        starbeast2.SpeciesTreeLoggerWithGeneFlow,
        starbeast2.SpeciesTreeParser,
        starbeast2.StarBeastClock,
        starbeast2.StarBeastInitializer,
        starbeast2.StarBeastTaxonSet,
        starbeast2.TreeLengthLogger,
        starbeast2.UncorrelatedRates,
        starbeast2.UniformPopulations,
        starbeast2.aimoperator.ExchangeAndSwap,
        starbeast2.aimoperator.RankingAwareOperator,
        starbeast2.aimoperator.SubtreeSlideAndSwap,
        starbeast2.aimoperator.UniformAndSwap,
        starbeast2.aimoperator.WilsonBaldingAndSwap,
        starbeast2.utils.BernoulliDistribution,
        starbeast2.utils.GeneTreeSimulator,
        starbeast2.utils.NodeHeightLogger,
        starbeast2.utils.SimulatedGeneTree,
        starbeast2.utils.SimulatedGeneTreeLogger,
        starbeast2.utils.TreeTopologyDistGenerator,
        starbeast2.utils.TreeTopologyDistLogger,
        starbeast2.utils.TreeTraceAnalysisWithError;

    provides beastfx.app.inputeditor.InputEditor with
        starbeast2.app.beauti.StarBeastTipDatesInputEditor;
}
