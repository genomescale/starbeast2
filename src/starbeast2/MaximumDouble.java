package starbeast2;

// store and return a single double value
// value if never set() is negative infinity
// if set() is called multiple times, the largest value will be stored
final class MaximumDouble {
    private double storedDouble;

    public MaximumDouble() {
        storedDouble = Double.NEGATIVE_INFINITY;
    }

    public void set(double inputDouble) {
        if (inputDouble > storedDouble) {
            storedDouble = inputDouble;
        }
    }

    public double get() {
        return storedDouble;
    }
}
