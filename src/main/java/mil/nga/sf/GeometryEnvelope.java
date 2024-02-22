package mil.nga.sf;

import java.io.Serializable;

import mil.nga.sf.util.GeometryEnvelopeBuilder;

/**
 * Geometry envelope
 * 
 * @author osbornb
 */
public class GeometryEnvelope implements Serializable {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Min X
	 */
	private double minX;

	/**
	 * Max X
	 */
	private double maxX;

	/**
	 * Min Y
	 */
	private double minY;

	/**
	 * Max Y
	 */
	private double maxY;

	/**
	 * True if has z coordinates
	 */
	private boolean hasZ = false;

	/**
	 * Min Z
	 */
	private Double minZ;

	/**
	 * Max Z
	 */
	private Double maxZ;

	/**
	 * True if has M measurements
	 */
	private boolean hasM = false;

	/**
	 * Min M
	 */
	private Double minM;

	/**
	 * Max M
	 */
	private Double maxM;

	/**
	 * Constructor
	 */
	public GeometryEnvelope() {

	}

	/**
	 * Constructor
	 * 
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 */
	public GeometryEnvelope(boolean hasZ, boolean hasM) {
		this.hasZ = hasZ;
		this.hasM = hasM;
	}

	/**
	 * Constructor
	 * 
	 * @param minX
	 *            min x
	 * @param minY
	 *            min y
	 * @param maxX
	 *            max x
	 * @param maxY
	 *            max y
	 */
	public GeometryEnvelope(double minX, double minY, double maxX,
			double maxY) {
		this(minX, minY, null, null, maxX, maxY, null, null);
	}

	/**
	 * Constructor
	 * 
	 * @param minX
	 *            min x
	 * @param minY
	 *            min y
	 * @param minZ
	 *            min z
	 * @param maxX
	 *            max x
	 * @param maxY
	 *            max y
	 * @param maxZ
	 *            max z
	 */
	public GeometryEnvelope(double minX, double minY, Double minZ, double maxX,
			double maxY, Double maxZ) {
		this(minX, minY, minZ, null, maxX, maxY, maxZ, null);
	}

	/**
	 * Constructor
	 * 
	 * @param minX
	 *            min x
	 * @param minY
	 *            min y
	 * @param minZ
	 *            min z
	 * @param minM
	 *            min m
	 * @param maxX
	 *            max x
	 * @param maxY
	 *            max y
	 * @param maxZ
	 *            max z
	 * @param maxM
	 *            max m
	 */
	public GeometryEnvelope(double minX, double minY, Double minZ, Double minM,
			double maxX, double maxY, Double maxZ, Double maxM) {
		this.minX = minX;
		this.minY = minY;
		this.minZ = minZ;
		this.minM = minM;
		this.maxX = maxX;
		this.maxY = maxY;
		this.maxZ = maxZ;
		this.maxM = maxM;
		this.hasZ = minZ != null || maxZ != null;
		this.hasM = minM != null || maxM != null;
	}

	/**
	 * Copy Constructor
	 * 
	 * @param envelope
	 *            Envelope to copy
	 */
	public GeometryEnvelope(GeometryEnvelope envelope) {
		this.minX = envelope.minX;
		this.maxX = envelope.maxX;
		this.minY = envelope.minY;
		this.maxY = envelope.maxY;
		this.hasZ = envelope.hasZ;
		this.minZ = envelope.minZ;
		this.maxZ = envelope.maxZ;
		this.hasM = envelope.hasM;
		this.minM = envelope.minM;
		this.maxM = envelope.maxM;
	}

	/**
	 * True if has Z coordinates
	 * 
	 * @return has z
	 */
	public boolean hasZ() {
		return hasZ;
	}

	/**
	 * True if has Z coordinates
	 * 
	 * @return has z
	 * @see #hasZ()
	 */
	public boolean is3D() {
		return hasZ();
	}

	/**
	 * True if has M measurements
	 * 
	 * @return has m
	 */
	public boolean hasM() {
		return hasM;
	}

	/**
	 * True if has M measurements
	 * 
	 * @return has m
	 * @see #hasM()
	 */
	public boolean isMeasured() {
		return hasM();
	}

	/**
	 * Get min x
	 * 
	 * @return min x
	 */
	public double getMinX() {
		return minX;
	}

	/**
	 * Set min x
	 * 
	 * @param minX
	 *            min x
	 */
	public void setMinX(double minX) {
		this.minX = minX;
	}

	/**
	 * Get max x
	 * 
	 * @return max x
	 */
	public double getMaxX() {
		return maxX;
	}

	/**
	 * Set max x
	 * 
	 * @param maxX
	 *            max x
	 */
	public void setMaxX(double maxX) {
		this.maxX = maxX;
	}

	/**
	 * Get min y
	 * 
	 * @return min y
	 */
	public double getMinY() {
		return minY;
	}

	/**
	 * Set min y
	 * 
	 * @param minY
	 *            min y
	 */
	public void setMinY(double minY) {
		this.minY = minY;
	}

	/**
	 * Get max y
	 * 
	 * @return max y
	 */
	public double getMaxY() {
		return maxY;
	}

	/**
	 * Set max y
	 * 
	 * @param maxY
	 *            max y
	 */
	public void setMaxY(double maxY) {
		this.maxY = maxY;
	}

	/**
	 * Has z coordinates
	 * 
	 * @return true if has z coordinates
	 */
	public boolean isHasZ() {
		return hasZ;
	}

	/**
	 * Set has z coordinates
	 * 
	 * @param hasZ
	 *            has z
	 */
	public void setHasZ(boolean hasZ) {
		this.hasZ = hasZ;
	}

	/**
	 * Get min z
	 * 
	 * @return min z
	 */
	public Double getMinZ() {
		return minZ;
	}

	/**
	 * Set min z
	 * 
	 * @param minZ
	 *            min z
	 */
	public void setMinZ(Double minZ) {
		this.minZ = minZ;
	}

	/**
	 * Get max z
	 * 
	 * @return max z
	 */
	public Double getMaxZ() {
		return maxZ;
	}

	/**
	 * Set max z
	 * 
	 * @param maxZ
	 *            max z
	 */
	public void setMaxZ(Double maxZ) {
		this.maxZ = maxZ;
	}

	/**
	 * Has m coordinates
	 * 
	 * @return true if has m coordinates
	 */
	public boolean isHasM() {
		return hasM;
	}

	/**
	 * Set has m coordinates
	 * 
	 * @param hasM
	 *            has m
	 */
	public void setHasM(boolean hasM) {
		this.hasM = hasM;
	}

	/**
	 * Get min m
	 * 
	 * @return min m
	 */
	public Double getMinM() {
		return minM;
	}

	/**
	 * Set min m
	 * 
	 * @param minM
	 *            min m
	 */
	public void setMinM(Double minM) {
		this.minM = minM;
	}

	/**
	 * Get max m
	 * 
	 * @return max m
	 */
	public Double getMaxM() {
		return maxM;
	}

	/**
	 * Set max m
	 * 
	 * @param maxM
	 *            max m
	 */
	public void setMaxM(Double maxM) {
		this.maxM = maxM;
	}

	/**
	 * Get the x range
	 * 
	 * @return x range
	 * @since 2.0.5
	 */
	public double getXRange() {
		return maxX - minX;
	}

	/**
	 * Get the y range
	 * 
	 * @return y range
	 * @since 2.0.5
	 */
	public double getYRange() {
		return maxY - minY;
	}

	/**
	 * Get the z range
	 * 
	 * @return z range
	 * @since 2.0.5
	 */
	public Double getZRange() {
		Double range = null;
		if (minZ != null && maxZ != null) {
			range = maxZ - minZ;
		}
		return range;
	}

	/**
	 * Get the m range
	 * 
	 * @return m range
	 * @since 2.0.5
	 */
	public Double getMRange() {
		Double range = null;
		if (minM != null && maxM != null) {
			range = maxM - minM;
		}
		return range;
	}

	/**
	 * Determine if the envelope is of a single point
	 * 
	 * @return true if a single point bounds
	 * @since 2.0.5
	 */
	public boolean isPoint() {
		return Double.compare(minX, maxX) == 0
				&& Double.compare(minY, maxY) == 0;
	}

	/**
	 * Get the top left point
	 * 
	 * @return top left point
	 * @since 2.2.0
	 */
	public Point getTopLeft() {
		return new Point(minX, maxY);
	}

	/**
	 * Get the bottom left point
	 * 
	 * @return bottom left point
	 * @since 2.2.0
	 */
	public Point getBottomLeft() {
		return new Point(minX, minY);
	}

	/**
	 * Get the bottom right point
	 * 
	 * @return bottom right point
	 * @since 2.2.0
	 */
	public Point getBottomRight() {
		return new Point(maxX, minY);
	}

	/**
	 * Get the top right point
	 * 
	 * @return top right point
	 * @since 2.2.0
	 */
	public Point getTopRight() {
		return new Point(maxX, maxY);
	}

	/**
	 * Get the left mid point
	 * 
	 * @return left mid point
	 * @since 2.2.2
	 */
	public Point getLeftMid() {
		return new Point(minX, getMidY());
	}

	/**
	 * Get the bottom mid point
	 * 
	 * @return bottom mid point
	 * @since 2.2.2
	 */
	public Point getBottomMid() {
		return new Point(getMidX(), minY);
	}

	/**
	 * Get the right mid point
	 * 
	 * @return right mid point
	 * @since 2.2.2
	 */
	public Point getRightMid() {
		return new Point(maxX, getMidY());
	}

	/**
	 * Get the top mid point
	 * 
	 * @return top mid point
	 * @since 2.2.2
	 */
	public Point getTopMid() {
		return new Point(getMidX(), maxY);
	}

	/**
	 * Get the left line
	 * 
	 * @return left line
	 * @since 2.2.0
	 */
	public Line getLeft() {
		return new Line(getTopLeft(), getBottomLeft());
	}

	/**
	 * Get the bottom line
	 * 
	 * @return bottom line
	 * @since 2.2.0
	 */
	public Line getBottom() {
		return new Line(getBottomLeft(), getBottomRight());
	}

	/**
	 * Get the right line
	 * 
	 * @return right line
	 * @since 2.2.0
	 */
	public Line getRight() {
		return new Line(getBottomRight(), getTopRight());
	}

	/**
	 * Get the top line
	 * 
	 * @return top line
	 * @since 2.2.0
	 */
	public Line getTop() {
		return new Line(getTopRight(), getTopLeft());
	}

	/**
	 * Get the envelope mid x
	 * 
	 * @return mid x
	 * @since 2.1.0
	 */
	public double getMidX() {
		return (minX + maxX) / 2.0;
	}

	/**
	 * Get the envelope mid y
	 * 
	 * @return mid y
	 * @since 2.1.0
	 */
	public double getMidY() {
		return (minY + maxY) / 2.0;
	}

	/**
	 * Get the envelope centroid point
	 * 
	 * @return centroid point
	 * @since 2.0.5
	 */
	public Point getCentroid() {
		return new Point(getMidX(), getMidY());
	}

	/**
	 * Determine if the envelope is empty
	 * 
	 * @return true if empty
	 * @since 2.1.0
	 */
	public boolean isEmpty() {
		return getXRange() <= 0.0 || getYRange() <= 0.0;
	}

	/**
	 * Determine if intersects with the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @return true if intersects
	 * @since 2.0.1
	 */
	public boolean intersects(GeometryEnvelope envelope) {
		return overlap(envelope) != null;
	}

	/**
	 * Determine if intersects with the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @param allowEmpty
	 *            allow empty ranges when determining intersection
	 *
	 * @return true if intersects
	 * @since 2.0.1
	 */
	public boolean intersects(GeometryEnvelope envelope, boolean allowEmpty) {
		return overlap(envelope, allowEmpty) != null;
	}

	/**
	 * Get the overlapping geometry envelope with the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @return geometry envelope
	 */
	public GeometryEnvelope overlap(GeometryEnvelope envelope) {
		return overlap(envelope, false);
	}

	/**
	 * Get the overlapping geometry envelope with the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @param allowEmpty
	 *            allow empty ranges when determining intersection
	 * @return geometry envelope
	 * @since 2.0.1
	 */
	public GeometryEnvelope overlap(GeometryEnvelope envelope,
			boolean allowEmpty) {

		double minX = Math.max(getMinX(), envelope.getMinX());
		double maxX = Math.min(getMaxX(), envelope.getMaxX());
		double minY = Math.max(getMinY(), envelope.getMinY());
		double maxY = Math.min(getMaxY(), envelope.getMaxY());

		GeometryEnvelope overlap = null;

		if ((minX < maxX && minY < maxY)
				|| (allowEmpty && minX <= maxX && minY <= maxY)) {
			overlap = new GeometryEnvelope(minX, minY, maxX, maxY);
		}

		return overlap;
	}

	/**
	 * Get the union geometry envelope combined with the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @return geometry envelope
	 */
	public GeometryEnvelope union(GeometryEnvelope envelope) {

		double minX = Math.min(getMinX(), envelope.getMinX());
		double maxX = Math.max(getMaxX(), envelope.getMaxX());
		double minY = Math.min(getMinY(), envelope.getMinY());
		double maxY = Math.max(getMaxY(), envelope.getMaxY());

		GeometryEnvelope union = null;

		if (minX < maxX && minY < maxY) {
			union = new GeometryEnvelope(minX, minY, maxX, maxY);
		}

		return union;
	}

	/**
	 * Determine if contains the point
	 *
	 * @param point
	 *            point
	 * @return true if contains
	 * @since 2.1.0
	 */
	public boolean contains(Point point) {
		return contains(point, 0.0);
	}

	/**
	 * Determine if contains the point
	 *
	 * @param point
	 *            point
	 * @param epsilon
	 *            epsilon equality tolerance
	 * @return true if contains
	 * @since 2.2.0
	 */
	public boolean contains(Point point, double epsilon) {
		return contains(point.getX(), point.getY(), epsilon);
	}

	/**
	 * Determine if contains the coordinate
	 *
	 * @param x
	 *            x value
	 * @param y
	 *            y value
	 * @return true if contains
	 * @since 2.1.0
	 */
	public boolean contains(double x, double y) {
		return contains(x, y, 0.0);
	}

	/**
	 * Determine if contains the coordinate
	 *
	 * @param x
	 *            x value
	 * @param y
	 *            y value
	 * @param epsilon
	 *            epsilon equality tolerance
	 * @return true if contains
	 * @since 2.2.0
	 */
	public boolean contains(double x, double y, double epsilon) {
		return x >= getMinX() - epsilon && x <= getMaxX() + epsilon
				&& y >= getMinY() - epsilon && y <= getMaxY() + epsilon;
	}

	/**
	 * Determine if inclusively contains the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @return true if contains
	 * @since 2.0.1
	 */
	public boolean contains(GeometryEnvelope envelope) {
		return contains(envelope, 0.0);
	}

	/**
	 * Determine if inclusively contains the provided envelope
	 *
	 * @param envelope
	 *            geometry envelope
	 * @param epsilon
	 *            epsilon equality tolerance
	 * @return true if contains
	 * @since 2.2.0
	 */
	public boolean contains(GeometryEnvelope envelope, double epsilon) {
		return getMinX() - epsilon <= envelope.getMinX()
				&& getMaxX() + epsilon >= envelope.getMaxX()
				&& getMinY() - epsilon <= envelope.getMinY()
				&& getMaxY() + epsilon >= envelope.getMaxY();
	}

	/**
	 * Build a geometry representation of the geometry envelope
	 * 
	 * @return geometry, polygon or point
	 * @since 2.0.5
	 */
	public Geometry buildGeometry() {
		return GeometryEnvelopeBuilder.buildGeometry(this);
	}

	/**
	 * Copy the geometry envelope
	 * 
	 * @return geometry envelope copy
	 */
	public GeometryEnvelope copy() {
		return new GeometryEnvelope(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (hasM ? 1231 : 1237);
		result = prime * result + (hasZ ? 1231 : 1237);
		result = prime * result + ((maxM == null) ? 0 : maxM.hashCode());
		long temp;
		temp = Double.doubleToLongBits(maxX);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(maxY);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((maxZ == null) ? 0 : maxZ.hashCode());
		result = prime * result + ((minM == null) ? 0 : minM.hashCode());
		temp = Double.doubleToLongBits(minX);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(minY);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((minZ == null) ? 0 : minZ.hashCode());
		return result;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GeometryEnvelope other = (GeometryEnvelope) obj;
		if (hasM != other.hasM)
			return false;
		if (hasZ != other.hasZ)
			return false;
		if (maxM == null) {
			if (other.maxM != null)
				return false;
		} else if (!maxM.equals(other.maxM))
			return false;
		if (Double.doubleToLongBits(maxX) != Double
				.doubleToLongBits(other.maxX))
			return false;
		if (Double.doubleToLongBits(maxY) != Double
				.doubleToLongBits(other.maxY))
			return false;
		if (maxZ == null) {
			if (other.maxZ != null)
				return false;
		} else if (!maxZ.equals(other.maxZ))
			return false;
		if (minM == null) {
			if (other.minM != null)
				return false;
		} else if (!minM.equals(other.minM))
			return false;
		if (Double.doubleToLongBits(minX) != Double
				.doubleToLongBits(other.minX))
			return false;
		if (Double.doubleToLongBits(minY) != Double
				.doubleToLongBits(other.minY))
			return false;
		if (minZ == null) {
			if (other.minZ != null)
				return false;
		} else if (!minZ.equals(other.minZ))
			return false;
		return true;
	}

}
