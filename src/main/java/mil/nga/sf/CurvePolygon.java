package mil.nga.sf;

import java.util.ArrayList;
import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * A planar surface defined by an exterior ring and zero or more interior ring.
 * Each ring is defined by a Curve instance.
 * 
 * @author osbornb
 */
public class CurvePolygon<T extends Curve> extends Surface {

	/**
	 * List of rings
	 */
	private List<T> rings = new ArrayList<T>();

	/**
	 * Constructor
	 */
	public CurvePolygon() {
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
	public CurvePolygon(boolean hasZ, boolean hasM) {
		super(GeometryType.CURVEPOLYGON, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param rings
	 *            list of rings
	 */
	public CurvePolygon(List<T> rings) {
		this(GeometryUtils.hasZ(rings), GeometryUtils.hasM(rings));
		setRings(rings);
	}

	/**
	 * Constructor
	 * 
	 * @param ring
	 *            ring
	 */
	public CurvePolygon(T ring) {
		this(ring.hasZ(), ring.hasM());
		addRing(ring);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param curvePolygon
	 *            curve polygon to copy
	 */
	public CurvePolygon(CurvePolygon<T> curvePolygon) {
		this(curvePolygon.hasZ(), curvePolygon.hasM());
		for (T ring : curvePolygon.getRings()) {
			@SuppressWarnings("unchecked")
			T ringCopy = (T) ring.copy();
			addRing(ringCopy);
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
	protected CurvePolygon(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the rings
	 * 
	 * @return rings
	 */
	public List<T> getRings() {
		return rings;
	}

	/**
	 * Set the rings
	 * 
	 * @param rings
	 *            rings
	 */
	public void setRings(List<T> rings) {
		this.rings = rings;
	}

	/**
	 * Add a ring
	 * 
	 * @param ring
	 *            ring
	 */
	public void addRing(T ring) {
		rings.add(ring);
	}

	/**
	 * Add rings
	 * 
	 * @param rings
	 *            rings
	 */
	public void addRings(List<T> rings) {
		this.rings.addAll(rings);
	}

	/**
	 * Get the number of rings including exterior and interior
	 * 
	 * @return number of rings
	 */
	public int numRings() {
		return rings.size();
	}

	/**
	 * Returns the Nth ring where the exterior ring is at 0, interior rings
	 * begin at 1
	 * 
	 * @param n
	 *            nth ring to return
	 * @return ring
	 */
	public T getRing(int n) {
		return rings.get(n);
	}

	/**
	 * Get the exterior ring
	 * 
	 * @return exterior ring
	 */
	public T getExteriorRing() {
		return rings.get(0);
	}

	/**
	 * Get the number of interior rings
	 * 
	 * @return number of interior rings
	 */
	public int numInteriorRings() {
		return rings.size() - 1;
	}

	/**
	 * Returns the Nth interior ring for this Polygon
	 * 
	 * @param n
	 *            interior ring number
	 * @return interior ring
	 */
	public T getInteriorRing(int n) {
		return rings.get(n + 1);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new CurvePolygon<T>(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return rings.isEmpty();
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
		result = prime * result + ((rings == null) ? 0 : rings.hashCode());
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
		CurvePolygon<T> other = (CurvePolygon<T>) obj;
		if (rings == null) {
			if (other.rings != null)
				return false;
		} else if (!rings.equals(other.rings))
			return false;
		return true;
	}

}
