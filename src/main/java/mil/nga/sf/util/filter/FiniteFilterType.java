package mil.nga.sf.util.filter;

/**
 * Finite Filter Type, including finite values and optionally one of either
 * infinite or NaN values
 * 
 * @author osbornb
 * @since 2.0.3
 */
public enum FiniteFilterType {

	/**
	 * Accept only finite values
	 */
	FINITE,

	/**
	 * Accept finite and infinite values
	 */
	FINITE_AND_INFINITE,

	/**
	 * Accept finite and Not a Number values
	 */
	FINITE_AND_NAN;

}
