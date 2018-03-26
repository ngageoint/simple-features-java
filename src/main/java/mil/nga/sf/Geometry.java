package mil.nga.sf;

import mil.nga.sf.util.GeometryEnvelopeBuilder;
import mil.nga.sf.util.GeometryUtils;

/**
 * The root of the geometry type hierarchy
 * 
 * @author osbornb
 */
public abstract class Geometry {

	/**
	 * Geometry type
	 */
	private final GeometryType geometryType;

	/**
	 * Has z coordinates
	 */
	private final boolean hasZ;

	/**
	 * Has m values
	 */
	private final boolean hasM;

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
	 * Get the mathematical centroid for this Geometry as a Point
	 * 
	 * @return centroid point
	 */
	public Point getCentroid() {
		return GeometryUtils.getCentroid(this);
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

}
