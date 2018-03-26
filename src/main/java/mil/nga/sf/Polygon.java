package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.sweep.ShamosHoey;

/**
 * A restricted form of CurvePolygon where each ring is defined as a simple,
 * closed LineString.
 * 
 * @author osbornb
 */
public class Polygon extends CurvePolygon<LineString> {

	/**
	 * Constructor
	 */
	public Polygon() {
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
	public Polygon(boolean hasZ, boolean hasM) {
		super(GeometryType.POLYGON, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param rings
	 *            list of rings
	 */
	public Polygon(List<LineString> rings) {
		this(GeometryUtils.hasZ(rings), GeometryUtils.hasM(rings));
		setRings(rings);
	}

	/**
	 * Constructor
	 * 
	 * @param ring
	 *            ring
	 */
	public Polygon(LineString ring) {
		this(ring.hasZ(), ring.hasM());
		addRing(ring);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param polygon
	 *            polygon to copy
	 */
	public Polygon(Polygon polygon) {
		this(polygon.hasZ(), polygon.hasM());
		for (LineString ring : polygon.getRings()) {
			addRing((LineString) ring.copy());
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
	protected Polygon(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new Polygon(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		return ShamosHoey.simplePolygon(this);
	}

}
