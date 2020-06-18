package mil.nga.sf;

import java.util.ArrayList;
import java.util.List;

import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.sweep.ShamosHoey;

/**
 * A Curve that connects two or more points in space.
 * 
 * @author osbornb
 */
public class LineString extends Curve {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * List of points
	 */
	private List<Point> points = new ArrayList<Point>();

	/**
	 * Constructor
	 */
	public LineString() {
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
	public LineString(boolean hasZ, boolean hasM) {
		super(GeometryType.LINESTRING, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param points
	 *            list of points
	 */
	public LineString(List<Point> points) {
		this(GeometryUtils.hasZ(points), GeometryUtils.hasM(points));
		setPoints(points);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param lineString
	 *            line string to copy
	 */
	public LineString(LineString lineString) {
		this(lineString.hasZ(), lineString.hasM());
		for (Point point : lineString.getPoints()) {
			addPoint((Point) point.copy());
		}
	}

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
	protected LineString(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the points
	 * 
	 * @return points
	 */
	public List<Point> getPoints() {
		return points;
	}

	/**
	 * Set the points
	 * 
	 * @param points
	 *            points
	 */
	public void setPoints(List<Point> points) {
		this.points = points;
	}

	/**
	 * Add a point
	 * 
	 * @param point
	 *            point
	 */
	public void addPoint(Point point) {
		points.add(point);
		updateZM(point);
	}

	/**
	 * Add points
	 * 
	 * @param points
	 *            points
	 */
	public void addPoints(List<Point> points) {
		for (Point point : points) {
			addPoint(point);
		}
	}

	/**
	 * Get the number of points
	 * 
	 * @return number of points
	 */
	public int numPoints() {
		return points.size();
	}

	/**
	 * Returns the Nth point
	 * 
	 * @param n
	 *            nth point to return
	 * @return point
	 */
	public Point getPoint(int n) {
		return points.get(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Point startPoint() {
		Point startPoint = null;
		if (!isEmpty()) {
			startPoint = points.get(0);
		}
		return startPoint;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Point endPoint() {
		Point endPoint = null;
		if (!isEmpty()) {
			endPoint = points.get(points.size() - 1);
		}
		return endPoint;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		return ShamosHoey.simplePolygonPoints(points);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new LineString(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return points.isEmpty();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((points == null) ? 0 : points.hashCode());
		return result;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		LineString other = (LineString) obj;
		if (points == null) {
			if (other.points != null)
				return false;
		} else if (!points.equals(other.points))
			return false;
		return true;
	}

}
