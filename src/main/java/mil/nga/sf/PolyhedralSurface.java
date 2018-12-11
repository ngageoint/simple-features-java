package mil.nga.sf;

import java.util.ArrayList;
import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * Contiguous collection of polygons which share common boundary segments.
 * 
 * @author osbornb
 */
public class PolyhedralSurface extends Surface {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * List of polygons
	 */
	private List<Polygon> polygons = new ArrayList<Polygon>();

	/**
	 * Constructor
	 */
	public PolyhedralSurface() {
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
	public PolyhedralSurface(boolean hasZ, boolean hasM) {
		super(GeometryType.POLYHEDRALSURFACE, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param polygons
	 *            list of polygons
	 */
	public PolyhedralSurface(List<Polygon> polygons) {
		this(GeometryUtils.hasZ(polygons), GeometryUtils.hasM(polygons));
		setPolygons(polygons);
	}

	/**
	 * Constructor
	 * 
	 * @param polygon
	 *            polygon
	 */
	public PolyhedralSurface(Polygon polygon) {
		this(polygon.hasZ(), polygon.hasM());
		addPolygon(polygon);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface to copy
	 */
	public PolyhedralSurface(PolyhedralSurface polyhedralSurface) {
		this(polyhedralSurface.hasZ(), polyhedralSurface.hasM());
		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			addPolygon((Polygon) polygon.copy());
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
	protected PolyhedralSurface(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get polygons
	 * 
	 * @return polygons
	 */
	public List<Polygon> getPolygons() {
		return polygons;
	}

	/**
	 * Get patches
	 * 
	 * @return patches
	 * @see #getPolygons()
	 */
	public List<Polygon> getPatches() {
		return getPolygons();
	}

	/**
	 * Set polygons
	 * 
	 * @param polygons
	 *            polygons
	 */
	public void setPolygons(List<Polygon> polygons) {
		this.polygons = polygons;
	}

	/**
	 * Set patches
	 * 
	 * @param patches
	 *            patches
	 * @see #setPolygons(List)
	 */
	public void setPatches(List<Polygon> patches) {
		setPolygons(patches);
	}

	/**
	 * Add polygon
	 * 
	 * @param polygon
	 *            polygon
	 */
	public void addPolygon(Polygon polygon) {
		polygons.add(polygon);
	}

	/**
	 * Add patch
	 * 
	 * @param patch
	 *            patch
	 * @see #addPolygon(Polygon)
	 */
	public void addPatch(Polygon patch) {
		addPolygon(patch);
	}

	/**
	 * Add polygons
	 * 
	 * @param polygons
	 *            polygons
	 */
	public void addPolygons(List<Polygon> polygons) {
		this.polygons.addAll(polygons);
	}

	/**
	 * Add patches
	 * 
	 * @param patches
	 *            patches
	 * @see #addPolygons(List)
	 */
	public void addPatches(List<Polygon> patches) {
		addPolygons(patches);
	}

	/**
	 * Get the number of polygons
	 * 
	 * @return number of polygons
	 */
	public int numPolygons() {
		return polygons.size();
	}

	/**
	 * Get the number of polygons
	 * 
	 * @return number of polygons
	 * @see #numPolygons()
	 */
	public int numPatches() {
		return numPolygons();
	}

	/**
	 * Get the Nth polygon
	 * 
	 * @param n
	 *            nth polygon to return
	 * @return polygon
	 */
	public Polygon getPolygon(int n) {
		return polygons.get(n);
	}

	/**
	 * Get the Nth polygon patch
	 * 
	 * @param n
	 *            nth polygon patch to return
	 * @return polygon patch
	 * @see #getPolygon(int)
	 */
	public Polygon getPatch(int n) {
		return getPolygon(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new PolyhedralSurface(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return polygons.isEmpty();
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
				+ ((polygons == null) ? 0 : polygons.hashCode());
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
		PolyhedralSurface other = (PolyhedralSurface) obj;
		if (polygons == null) {
			if (other.polygons != null)
				return false;
		} else if (!polygons.equals(other.polygons))
			return false;
		return true;
	}

}
