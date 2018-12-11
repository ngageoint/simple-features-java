package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * Circular String, Curve sub type
 * 
 * @author osbornb
 */
public class CircularString extends LineString {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Constructor
	 */
	public CircularString() {
		this(false, false);
	}

	/**
	 * Constructor
	 * 
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 */
	public CircularString(boolean hasZ, boolean hasM) {
		super(GeometryType.CIRCULARSTRING, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param points
	 *            list of points
	 */
	public CircularString(List<Point> points) {
		this(GeometryUtils.hasZ(points), GeometryUtils.hasM(points));
		setPoints(points);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param circularString
	 *            circular string to copy
	 */
	public CircularString(CircularString circularString) {
		this(circularString.hasZ(), circularString.hasM());
		for (Point point : circularString.getPoints()) {
			addPoint((Point) point.copy());
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new CircularString(this);
	}

}
