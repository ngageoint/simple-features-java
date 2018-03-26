package mil.nga.sf;

/**
 * The base type for all 1-dimensional geometry types. A 1-dimensional geometry
 * is a geometry that has a length, but no area. A curve is considered simple if
 * it does not intersect itself (except at the start and end point). A curve is
 * considered closed its start and end point are coincident. A simple, closed
 * curve is called a ring.
 * 
 * @author osbornb
 */
public abstract class Curve extends Geometry {

	/**
	 * Constructor
	 * 
	 * @param type
	 *            geometry type
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 */
	protected Curve(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the start Point of this Curve
	 * 
	 * @return start point
	 */
	public abstract Point startPoint();

	/**
	 * Get the end Point of this Curve
	 * 
	 * @return end point
	 */
	public abstract Point endPoint();

	/**
	 * Determine if this Curve is closed (start point = end point)
	 * 
	 * @return true if closed
	 */
	public boolean isClosed() {
		return !isEmpty() && startPoint().equals(endPoint());
	}

	/**
	 * Determine if this Curve is a ring (closed and simple)
	 * 
	 * @return true if a ring
	 */
	public boolean isRing() {
		return isClosed() && isSimple();
	}

}
