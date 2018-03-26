package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * A tetrahedron (4 triangular faces), corner at the origin and each unit
 * coordinate digit.
 * 
 * @author osbornb
 */
public class TIN extends PolyhedralSurface {

	/**
	 * Constructor
	 */
	public TIN() {
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
	public TIN(boolean hasZ, boolean hasM) {
		super(GeometryType.TIN, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param polygons
	 *            list of polygons
	 */
	public TIN(List<Polygon> polygons) {
		this(GeometryUtils.hasZ(polygons), GeometryUtils.hasM(polygons));
		setPolygons(polygons);
	}

	/**
	 * Constructor
	 * 
	 * @param polygon
	 *            polygon
	 */
	public TIN(Polygon polygon) {
		this(polygon.hasZ(), polygon.hasM());
		addPolygon(polygon);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param tin
	 *            tin to copy
	 */
	public TIN(TIN tin) {
		this(tin.hasZ(), tin.hasM());
		for (Polygon polygon : tin.getPolygons()) {
			addPolygon((Polygon) polygon.copy());
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new TIN(this);
	}

}
