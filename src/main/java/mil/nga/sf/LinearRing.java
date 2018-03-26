package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.SFException;

/**
 * A LineString that is both closed and simple.
 * 
 * @author osbornb
 */
public class LinearRing extends LineString {

	/**
	 * Constructor
	 */
	public LinearRing() {
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
	public LinearRing(boolean hasZ, boolean hasM) {
		super(hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param points
	 *            list of points
	 */
	public LinearRing(List<Point> points) {
		this(GeometryUtils.hasZ(points), GeometryUtils.hasM(points));
		setPoints(points);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param linearRing
	 *            linear ring to copy
	 */
	public LinearRing(LinearRing linearRing) {
		this(linearRing.hasZ(), linearRing.hasM());
		for (Point point : linearRing.getPoints()) {
			addPoint((Point) point.copy());
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setPoints(List<Point> points) {
		super.setPoints(points);
		if (!isEmpty()) {
			if (!isClosed()) {
				addPoint(points.get(0));
			}
			if (numPoints() < 4) {
				throw new SFException(
						"A closed linear ring must have at least four points.");
			}
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new LinearRing(this);
	}

}
