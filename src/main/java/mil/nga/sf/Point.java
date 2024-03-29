package mil.nga.sf;

/**
 * A single location in space. Each point has an X and Y coordinate. A point MAY
 * optionally also have a Z and/or an M value.
 * 
 * @author osbornb
 */
public class Point extends Geometry {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * X coordinate
	 */
	private double x;

	/**
	 * Y coordinate
	 */
	private double y;

	/**
	 * Z coordinate
	 */
	private Double z;

	/**
	 * M value
	 */
	private Double m;

	/**
	 * Constructor
	 */
	public Point() {
		this(0.0, 0.0);
	}

	/**
	 * Constructor
	 * 
	 * @param x
	 *            x coordinate
	 * @param y
	 *            y coordinate
	 */
	public Point(double x, double y) {
		this(false, false, x, y);
	}

	/**
	 * Constructor
	 * 
	 * @param x
	 *            x coordinate
	 * @param y
	 *            y coordinate
	 * @param z
	 *            z coordinate
	 */
	public Point(double x, double y, Double z) {
		this(x, y, z, null);
	}

	/**
	 * Constructor
	 * 
	 * @param x
	 *            x coordinate
	 * @param y
	 *            y coordinate
	 * @param z
	 *            z coordinate
	 * @param m
	 *            m coordinate
	 */
	public Point(double x, double y, Double z, Double m) {
		super(GeometryType.POINT, z != null, m != null);
		this.x = x;
		this.y = y;
		this.z = z;
		this.m = m;
	}

	/**
	 * Constructor
	 * 
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 * @param x
	 *            x coordinate
	 * @param y
	 *            y coordinate
	 */
	public Point(boolean hasZ, boolean hasM, double x, double y) {
		super(GeometryType.POINT, hasZ, hasM);
		this.x = x;
		this.y = y;
	}

	/**
	 * Copy Constructor
	 * 
	 * @param point
	 *            point to copy
	 */
	public Point(Point point) {
		this(point.hasZ(), point.hasM(), point.getX(), point.getY());
		setZ(point.getZ());
		setM(point.getM());
	}

	/**
	 * Get x
	 * 
	 * @return x
	 */
	public double getX() {
		return x;
	}

	/**
	 * Set x
	 * 
	 * @param x
	 *            x coordinate
	 */
	public void setX(double x) {
		this.x = x;
	}

	/**
	 * Get y
	 * 
	 * @return y
	 */
	public double getY() {
		return y;
	}

	/**
	 * Set y
	 * 
	 * @param y
	 *            y coordinate
	 */
	public void setY(double y) {
		this.y = y;
	}

	/**
	 * Get z
	 * 
	 * @return z
	 */
	public Double getZ() {
		return z;
	}

	/**
	 * Set z
	 * 
	 * @param z
	 *            z coordinate
	 */
	public void setZ(Double z) {
		this.z = z;
		setHasZ(z != null);
	}

	/**
	 * Get m
	 * 
	 * @return m
	 */
	public Double getM() {
		return m;
	}

	/**
	 * Set m
	 * 
	 * @param m
	 *            m coordinate
	 */
	public void setM(Double m) {
		this.m = m;
		setHasM(m != null);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new Point(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return false;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		return true;
	}

	/**
	 * Indicates if x values are equal
	 * 
	 * @param point
	 *            point to compare
	 * @return true if x is equal
	 * @since 2.2.1
	 */
	public boolean equalsX(Point point) {
		return Double.doubleToLongBits(x) == Double.doubleToLongBits(point.x);
	}

	/**
	 * Indicates if y values are equal
	 * 
	 * @param point
	 *            point to compare
	 * @return true if y is equal
	 * @since 2.2.1
	 */
	public boolean equalsY(Point point) {
		return Double.doubleToLongBits(y) == Double.doubleToLongBits(point.y);
	}

	/**
	 * Indicates if x and y values are equal
	 * 
	 * @param point
	 *            point to compare
	 * @return true if x and y are equal
	 * @since 2.2.1
	 */
	public boolean equalsXY(Point point) {
		return equalsX(point) && equalsY(point);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((m == null) ? 0 : m.hashCode());
		long temp;
		temp = Double.doubleToLongBits(x);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(y);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((z == null) ? 0 : z.hashCode());
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
		Point other = (Point) obj;
		if (m == null) {
			if (other.m != null)
				return false;
		} else if (!m.equals(other.m))
			return false;
		if (!equalsXY(other))
			return false;
		if (z == null) {
			if (other.z != null)
				return false;
		} else if (!z.equals(other.z))
			return false;
		return true;
	}

}
