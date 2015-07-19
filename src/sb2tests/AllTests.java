package sb2tests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantIOTest.class,
    ConstantPopulationTest.class,
    LinearPopulationTest.class,
    LinearWithConstantRootTest.class,
    IncompatibleTreeTest.class,
    CoordinatedExchangeTest.class
})
public class AllTests {

}
