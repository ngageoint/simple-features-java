package mil.nga.sf;

import java.util.ArrayList;
import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * A collection of zero or more Geometry instances.
 * 
 * @author osbornb
 */
public class GeometryCollection<T extends Geometry> extends Geometry {

	/**
	 * List of geometries
	 */
	private List<T> geometries = new ArrayList<T>();

	/**
	 * Constructor
	 */
	public GeometryCollection() {
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
	public GeometryCollection(boolean hasZ, boolean hasM) {
		super(GeometryType.GEOMETRYCOLLECTION, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param geometries
	 *            list of geometries
	 */
	public GeometryCollection(List<T> geometries) {
		this(GeometryUtils.hasZ(geometries), GeometryUtils.hasM(geometries));
		setGeometries(geometries);
	}

	/**
	 * Constructor
	 * 
	 * @param geometry
	 *            geometry
	 */
	public GeometryCollection(T geometry) {
		this(geometry.hasZ(), geometry.hasM());
		addGeometry(geometry);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param geometryCollection
	 *            geometry collection to copy
	 */
	public GeometryCollection(GeometryCollection<T> geometryCollection) {
		this(geometryCollection.hasZ(), geometryCollection.hasM());
		for (T geometry : geometryCollection.getGeometries()) {
			@SuppressWarnings("unchecked")
			T geometryCopy = (T) geometry.copy();
			addGeometry(geometryCopy);
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
	protected GeometryCollection(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the list of geometries
	 * 
	 * @return geometries
	 */
	public List<T> getGeometries() {
		return geometries;
	}

	/**
	 * Set the geometries
	 * 
	 * @param geometries
	 *            geometries
	 */
	public void setGeometries(List<T> geometries) {
		this.geometries = geometries;
	}

	/**
	 * Add a geometry
	 * 
	 * @param geometry
	 *            geometry
	 */
	public void addGeometry(T geometry) {
		geometries.add(geometry);
	}

	/**
	 * Get the number of geometries in the collection
	 * 
	 * @return number of geometries
	 */
	public int numGeometries() {
		return geometries.size();
	}

	/**
	 * Returns the Nth geometry
	 * 
	 * @param n
	 *            nth geometry to return
	 * @return geometry
	 */
	public Geometry getGeometry(int n) {
		return geometries.get(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new GeometryCollection<T>(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return geometries.isEmpty();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		throw new UnsupportedOperationException(
				"Is Simple not implemented for " + getClass().getSimpleName());
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result
				+ ((geometries == null) ? 0 : geometries.hashCode());
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
		@SuppressWarnings("unchecked")
		GeometryCollection<T> other = (GeometryCollection<T>) obj;
		if (geometries == null) {
			if (other.geometries != null)
				return false;
		} else if (!geometries.equals(other.geometries))
			return false;
		return true;
	}

}
