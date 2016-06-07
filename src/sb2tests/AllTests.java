package sb2tests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantIOTest.class,
    ConstantPopulationTest.class,
    ContinuousRatesTest.class,
    UncorrelatedRatesTest.class,
    StarbeastClockTest.class,
    LinearPopulationTest.class,
    LinearWithConstantRootTest.class,
    IncompatibleTreeTest.class,
    SmallCoordinatedExchangeTest.class,
    BigCoordinatedExchangeTest.class,
    WideExchangeTest.class,
    MissingDataCoordinatedExchange.class,
    MissingDataConstantIO.class,
})
public class AllTests {

}
