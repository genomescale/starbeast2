package snetworktests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantIOTest.class,
    ConstantPopulationTest.class,
    ContinuousRatesTest.class,
    DiscreteRatesTest.class,
    StarbeastClockTest.class,
    IncompatibleTreeTest.class,
    SmallCoordinatedExchangeTest.class,
    BigCoordinatedExchangeTest.class,
    MissingDataCoordinatedExchange.class,
    MissingDataConstantIO.class,
})
public class AllTests {

}
