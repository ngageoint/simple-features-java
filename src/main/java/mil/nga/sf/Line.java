package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.SFException;

/**
 * A LineString with exactly 2 Points.
 * 
 * @author osbornb
 */
public class Line extends LineString {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Constructor
	 */
	public Line() {
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
	public Line(boolean hasZ, boolean hasM) {
		super(hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param points
	 *            list of points
	 */
	public Line(List<Point> points) {
		this(GeometryUtils.hasZ(points), GeometryUtils.hasM(points));
		setPoints(points);
	}

	/**
	 * Constructor
	 * 
	 * @param point1
	 *            first point
	 * @param point2
	 *            second point
	 * @since 2.2.0
	 */
	public Line(Point point1, Point point2) {
		this(point1.hasZ() || point2.hasZ(), point1.hasM() || point2.hasM());
		addPoint(point1);
		addPoint(point2);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param line
	 *            line to copy
	 */
	public Line(Line line) {
		this(line.hasZ(), line.hasM());
		for (Point point : line.getPoints()) {
			addPoint((Point) point.copy());
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setPoints(List<Point> points) {
		super.setPoints(points);
		if (numPoints() != 2) {
			throw new SFException("A line must have exactly 2 points.");
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new Line(this);
	}

}
