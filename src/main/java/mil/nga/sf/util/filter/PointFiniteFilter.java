package mil.nga.sf.util.filter;

import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryType;
import mil.nga.sf.Point;
import mil.nga.sf.util.SFException;

/**
 * Point filter for finite checks on x and y properties, optionally filter on z
 * and m properties and non finite values (NaN or infinity)
 * 
 * @author osbornb
 * @since 2.0.3
 */
public class PointFiniteFilter implements GeometryFilter {

	/**
	 * Finite Filter type
	 */
	private FiniteFilterType type = FiniteFilterType.FINITE;

	/**
	 * Include z values in filtering
	 */
	boolean filterZ = false;

	/**
	 * Include m values in filtering
	 */
	boolean filterM = false;

	/**
	 * Default Constructor, filter on x and y, allowing only finite values
	 */
	public PointFiniteFilter() {

	}

	/**
	 * Constructor, filter on x and y
	 * 
	 * @param type
	 *            finite filter type
	 */
	public PointFiniteFilter(FiniteFilterType type) {
		setType(type);
	}

	/**
	 * Constructor, filter on x, y, and optionally z
	 * 
	 * @param type
	 *            finite filter type
	 * @param filterZ
	 *            filter z values mode
	 */
	public PointFiniteFilter(FiniteFilterType type, boolean filterZ) {
		setType(type);
		setFilterZ(filterZ);
	}

	/**
	 * Constructor, filter on x, y, and optionally z and m
	 * 
	 * @param type
	 *            finite filter type
	 * @param filterZ
	 *            filter z values mode
	 * @param filterM
	 *            filter m values mode
	 */
	public PointFiniteFilter(FiniteFilterType type, boolean filterZ,
			boolean filterM) {
		setType(type);
		setFilterZ(filterZ);
		setFilterM(filterM);
	}

	/**
	 * Constructor, filter on x, y, and optionally z
	 * 
	 * @param filterZ
	 *            filter z values mode
	 */
	public PointFiniteFilter(boolean filterZ) {
		setFilterZ(filterZ);
	}

	/**
	 * Constructor, filter on x, y, and optionally z and m
	 * 
	 * @param filterZ
	 *            filter z values mode
	 * @param filterM
	 *            filter m values mode
	 */
	public PointFiniteFilter(boolean filterZ, boolean filterM) {
		setFilterZ(filterZ);
		setFilterM(filterM);
	}

	/**
	 * Get the finite filter type
	 * 
	 * @return finite filter type
	 */
	public FiniteFilterType getType() {
		return type;
	}

	/**
	 * Set the finite filter type, null defaults to
	 * {@link FiniteFilterType#FINITE}
	 * 
	 * @param type
	 *            finite filter type
	 */
	public void setType(FiniteFilterType type) {
		this.type = type;
	}

	/**
	 * Is filtering for z values enabled?
	 * 
	 * @return true if z filtering
	 */
	public boolean isFilterZ() {
		return filterZ;
	}

	/**
	 * Set the z value filtering mode
	 * 
	 * @param filterZ
	 *            true to z filter
	 */
	public void setFilterZ(boolean filterZ) {
		this.filterZ = filterZ;
	}

	/**
	 * Is filtering for m values enabled?
	 * 
	 * @return true if m filtering
	 */
	public boolean isFilterM() {
		return filterM;
	}

	/**
	 * Set the m value filtering mode
	 * 
	 * @param filterM
	 *            true to m filter
	 */
	public void setFilterM(boolean filterM) {
		this.filterM = filterM;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean filter(GeometryType containingType, Geometry geometry) {
		return geometry.getGeometryType() != GeometryType.POINT
				|| !(geometry instanceof Point) || filter((Point) geometry);
	}

	/**
	 * Filter the point
	 * 
	 * @param point
	 *            point
	 * @return true if passes filter and point should be included
	 */
	private boolean filter(Point point) {
		return filter(point.getX()) && filter(point.getY()) && filterZ(point)
				&& filterM(point);
	}

	/**
	 * Filter the double value
	 * 
	 * @param value
	 *            double value
	 * @return
	 */
	private boolean filter(double value) {
		boolean passes;
		switch (type) {
		case FINITE:
			passes = Double.isFinite(value);
			break;
		case FINITE_AND_INFINITE:
			passes = !Double.isNaN(value);
			break;
		case FINITE_AND_NAN:
			passes = !Double.isInfinite(value);
			break;
		default:
			throw new SFException("Unsupported filter type: " + type);
		}
		return passes;
	}

	/**
	 * Filter the Z value
	 * 
	 * @param point
	 *            point
	 * @return true if passes
	 */
	private boolean filterZ(Point point) {
		return !filterZ || !point.hasZ() || filter(point.getZ());
	}

	/**
	 * Filter the M value
	 * 
	 * @param point
	 *            point
	 * @return true if passes
	 */
	private boolean filterM(Point point) {
		return !filterM || !point.hasM() || filter(point.getM());
	}

	/**
	 * Filter the Double value
	 * 
	 * @param value
	 *            Double value
	 * @return true if passes
	 */
	private boolean filter(Double value) {
		return value == null || filter(value.doubleValue());
	}

}
