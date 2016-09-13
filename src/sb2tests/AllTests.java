package sb2tests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantIOTest.class,
    ConstantPopulationTest.class,
    UncorrelatedRatesTest.class,
    StarbeastClockTest.class,
    LinearWithConstantRootTest.class,
    IncompatibleTreeTest.class,
    SmallCoordinatedExchangeTest.class,
    BigCoordinatedExchangeTest.class,
    WideExchangeTest.class,
    MissingDataCoordinatedExchange.class,
    MissingDataConstantIO.class,
    TreeLengthLoggerTest.class,
})
public class AllTests {

}
