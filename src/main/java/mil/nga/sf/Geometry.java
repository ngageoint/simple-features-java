package mil.nga.sf;

import java.io.Serializable;

import mil.nga.sf.util.GeometryEnvelopeBuilder;
import mil.nga.sf.util.GeometryUtils;

/**
 * The root of the geometry type hierarchy
 * 
 * @author osbornb
 */
public abstract class Geometry implements Serializable {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Geometry type
	 */
	private final GeometryType geometryType;

	/**
	 * Has z coordinates
	 */
	private boolean hasZ;

	/**
	 * Has m values
	 */
	private boolean hasM;

	/**
	 * Constructor
	 * 
	 * @param geometryType
	 *            geometry type
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 */
	protected Geometry(GeometryType geometryType, boolean hasZ, boolean hasM) {
		this.geometryType = geometryType;
		this.hasZ = hasZ;
		this.hasM = hasM;
	}

	/**
	 * Get the geometry type
	 * 
	 * @return geometry type
	 */
	public GeometryType getGeometryType() {
		return geometryType;
	}

	/**
	 * Does the geometry have z coordinates
	 * 
	 * @return true if has z coordinates
	 */
	public boolean hasZ() {
		return hasZ;
	}

	/**
	 * Set if the geometry has z coordinates
	 * 
	 * @param hasZ
	 *            true if has z coordinates
	 * @since 2.0.3
	 */
	public void setHasZ(boolean hasZ) {
		this.hasZ = hasZ;
	}

	/**
	 * Does the geometry have z coordinates
	 * 
	 * @return true if has z coordinates
	 * @see #hasZ()
	 */
	public boolean is3D() {
		return hasZ();
	}

	/**
	 * Does the geometry have m coordinates
	 * 
	 * @return true if has m coordinates
	 */
	public boolean hasM() {
		return hasM;
	}

	/**
	 * Set if the geometry has m coordinates
	 * 
	 * @param hasM
	 *            true if has m coordinates
	 * @since 2.0.3
	 */
	public void setHasM(boolean hasM) {
		this.hasM = hasM;
	}

	/**
	 * Update currently false hasZ and hasM values using the provided geometry
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.0.3
	 */
	protected void updateZM(Geometry geometry) {
		if (!hasZ()) {
			setHasZ(geometry.hasZ());
		}
		if (!hasM()) {
			setHasM(geometry.hasM());
		}
	}

	/**
	 * Does the geometry have m coordinates.
	 * 
	 * @return true if has m coordinates
	 * @see #hasM()
	 */
	public boolean isMeasured() {
		return hasM();
	}

	/**
	 * Get the minimum bounding box for this Geometry
	 * 
	 * @return geometry envelope
	 */
	public GeometryEnvelope getEnvelope() {
		return GeometryEnvelopeBuilder.buildEnvelope(this);
	}

	/**
	 * Get the inherent dimension (0, 1, or 2) for this Geometry
	 * 
	 * @return dimension
	 */
	public int getDimension() {
		return GeometryUtils.getDimension(this);
	}

	/**
	 * Get the mathematical centroid point of a 2 dimensional representation of
	 * the Geometry (balancing point of a 2d cutout of the geometry). Only the x
	 * and y coordinate of the resulting point are calculated and populated. The
	 * resulting {@link Point#getZ()} and {@link Point#getM()} methods will
	 * always return null.
	 * 
	 * @return centroid point
	 */
	public Point getCentroid() {
		return GeometryUtils.getCentroid(this);
	}

	/**
	 * Get the geographic centroid point of a 2 dimensional representation of
	 * the degree unit Geometry. Only the x and y coordinate of the resulting
	 * point are calculated and populated. The resulting {@link Point#getZ()}
	 * and {@link Point#getM()} methods will always return null.
	 * 
	 * @return centroid point
	 * @since 2.0.5
	 */
	public Point getDegreesCentroid() {
		return GeometryUtils.getDegreesCentroid(this);
	}

	/**
	 * Copy the geometry
	 * 
	 * @return geometry copy
	 */
	public abstract Geometry copy();

	/**
	 * Is the Geometry empty
	 * 
	 * @return true if empty
	 */
	public abstract boolean isEmpty();

	/**
	 * Determine if this Geometry has no anomalous geometric points, such as
	 * self intersection or self tangency
	 * 
	 * @return true if simple
	 */
	public abstract boolean isSimple();

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((geometryType == null) ? 0 : geometryType.hashCode());
		result = prime * result + (hasM ? 1231 : 1237);
		result = prime * result + (hasZ ? 1231 : 1237);
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
		Geometry other = (Geometry) obj;
		if (geometryType != other.geometryType)
			return false;
		if (hasM != other.hasM)
			return false;
		if (hasZ != other.hasZ)
			return false;
		return true;
	}

}
