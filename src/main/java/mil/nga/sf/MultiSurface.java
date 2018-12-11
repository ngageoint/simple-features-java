package mil.nga.sf;

import java.util.List;

/**
 * A restricted form of GeometryCollection where each Geometry in the collection
 * must be of type Surface.
 * 
 * @author osbornb
 * @param <T>
 *            surface type
 */
public abstract class MultiSurface<T extends Surface> extends
		GeometryCollection<T> {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

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
	protected MultiSurface(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the surfaces
	 * 
	 * @return surfaces
	 */
	public List<T> getSurfaces() {
		return getGeometries();
	}

	/**
	 * Set the surfaces
	 * 
	 * @param surfaces
	 *            surfaces
	 */
	public void setSurfaces(List<T> surfaces) {
		setGeometries(surfaces);
	}

	/**
	 * Add a surface
	 * 
	 * @param surface
	 *            surface
	 */
	public void addSurface(T surface) {
		addGeometry(surface);
	}

	/**
	 * Add surfaces
	 * 
	 * @param surfaces
	 *            surfaces
	 */
	public void addSurfaces(List<T> surfaces) {
		addGeometries(surfaces);
	}

	/**
	 * Get the number of surfaces
	 * 
	 * @return number of surfaces
	 */
	public int numSurfaces() {
		return numGeometries();
	}

	/**
	 * Returns the Nth surface
	 * 
	 * @param n
	 *            nth line surface to return
	 * @return surface
	 */
	public T getSurface(int n) {
		return getGeometry(n);
	}

}
