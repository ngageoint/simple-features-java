package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * Triangle
 * 
 * @author osbornb
 */
public class Triangle extends Polygon {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Constructor
	 */
	public Triangle() {
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
	public Triangle(boolean hasZ, boolean hasM) {
		super(GeometryType.TRIANGLE, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param rings
	 *            list of rings
	 */
	public Triangle(List<LineString> rings) {
		this(GeometryUtils.hasZ(rings), GeometryUtils.hasM(rings));
		setRings(rings);
	}

	/**
	 * Constructor
	 * 
	 * @param ring
	 *            ring
	 */
	public Triangle(LineString ring) {
		this(ring.hasZ(), ring.hasM());
		addRing(ring);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param triangle
	 *            triangle to copy
	 */
	public Triangle(Triangle triangle) {
		this(triangle.hasZ(), triangle.hasM());
		for (LineString ring : triangle.getRings()) {
			addRing((LineString) ring.copy());
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new Triangle(this);
	}

}
