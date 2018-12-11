package mil.nga.sf.extended;

import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryCollection;
import mil.nga.sf.GeometryType;
import mil.nga.sf.util.SFException;

/**
 * Extended Geometry Collection providing abstract geometry collection type
 * support
 * 
 * @author osbornb
 *
 * @param <T>
 *            geometry type
 */
public class ExtendedGeometryCollection<T extends Geometry> extends
		GeometryCollection<T> {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Extended geometry collection geometry type
	 */
	private GeometryType geometryType = GeometryType.GEOMETRYCOLLECTION;

	/**
	 * Constructor, wraps a geometry collection as extended
	 * 
	 * @param geometryCollection
	 *            geometry collection
	 */
	public ExtendedGeometryCollection(GeometryCollection<T> geometryCollection) {
		super(GeometryType.GEOMETRYCOLLECTION, geometryCollection.hasZ(),
				geometryCollection.hasM());
		setGeometries(geometryCollection.getGeometries());
		updateGeometryType();
	}

	/**
	 * Copy Constructor
	 * 
	 * @param extendedGeometryCollection
	 *            extended geometry collection to copy
	 */
	public ExtendedGeometryCollection(
			ExtendedGeometryCollection<T> extendedGeometryCollection) {
		super(GeometryType.GEOMETRYCOLLECTION, extendedGeometryCollection
				.hasZ(), extendedGeometryCollection.hasM());
		for (T geometry : extendedGeometryCollection.getGeometries()) {
			@SuppressWarnings("unchecked")
			T geometryCopy = (T) geometry.copy();
			addGeometry(geometryCopy);
		}
		geometryType = extendedGeometryCollection.getGeometryType();
	}

	/**
	 * Update the extended geometry type based upon the contained geometries
	 */
	public void updateGeometryType() {
		GeometryType geometryType = getCollectionType();
		switch (geometryType) {
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			break;
		case MULTIPOINT:
			geometryType = GeometryType.GEOMETRYCOLLECTION;
			break;
		case MULTILINESTRING:
			geometryType = GeometryType.MULTICURVE;
			break;
		case MULTIPOLYGON:
			geometryType = GeometryType.MULTISURFACE;
			break;
		default:
			throw new SFException(
					"Unsupported extended geometry collection geometry type: "
							+ geometryType);
		}
		this.geometryType = geometryType;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public GeometryType getGeometryType() {
		return geometryType;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new ExtendedGeometryCollection<T>(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result
				+ ((geometryType == null) ? 0 : geometryType.hashCode());
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
		ExtendedGeometryCollection<T> other = (ExtendedGeometryCollection<T>) obj;
		if (geometryType != other.geometryType)
			return false;
		return true;
	}

}
