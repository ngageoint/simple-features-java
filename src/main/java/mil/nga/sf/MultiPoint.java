package mil.nga.sf;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import mil.nga.sf.util.GeometryUtils;

/**
 * A restricted form of GeometryCollection where each Geometry in the collection
 * must be of type Point.
 * 
 * @author osbornb
 */
public class MultiPoint extends GeometryCollection<Point> {

	/**
	 * Constructor
	 */
	public MultiPoint() {
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
	public MultiPoint(boolean hasZ, boolean hasM) {
		super(GeometryType.MULTIPOINT, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param points
	 *            list of points
	 */
	public MultiPoint(List<Point> points) {
		this(GeometryUtils.hasZ(points), GeometryUtils.hasM(points));
		setPoints(points);
	}

	/**
	 * Constructor
	 * 
	 * @param point
	 *            point
	 */
	public MultiPoint(Point point) {
		this(point.hasZ(), point.hasM());
		addPoint(point);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param multiPoint
	 *            multi point to copy
	 */
	public MultiPoint(MultiPoint multiPoint) {
		this(multiPoint.hasZ(), multiPoint.hasM());
		for (Point point : multiPoint.getPoints()) {
			addPoint((Point) point.copy());
		}
	}

	/**
	 * Get the points
	 * 
	 * @return points
	 */
	public List<Point> getPoints() {
		return getGeometries();
	}

	/**
	 * Set the points
	 * 
	 * @param points
	 *            points
	 */
	public void setPoints(List<Point> points) {
		setGeometries(points);
	}

	/**
	 * Add a point
	 * 
	 * @param point
	 *            point
	 */
	public void addPoint(Point point) {
		addGeometry(point);
	}

	/**
	 * Add points
	 * 
	 * @param points
	 *            points
	 */
	public void addPoints(List<Point> points) {
		addGeometries(points);
	}

	/**
	 * Get the number of points
	 * 
	 * @return number of points
	 */
	public int numPoints() {
		return numGeometries();
	}

	/**
	 * Returns the Nth point
	 * 
	 * @param n
	 *            nth point to return
	 * @return point
	 */
	public Point getPoint(int n) {
		return getGeometry(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new MultiPoint(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		Set<Point> points = new HashSet<>(getPoints());
		return points.size() == numPoints();
	}

}
