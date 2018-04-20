package mil.nga.sf;

import java.util.List;

import mil.nga.sf.util.GeometryUtils;

/**
 * A restricted form of MultiCurve where each Curve in the collection must be of
 * type LineString.
 * 
 * @author osbornb
 */
public class MultiLineString extends MultiCurve<LineString> {

	/**
	 * Constructor
	 */
	public MultiLineString() {
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
	public MultiLineString(boolean hasZ, boolean hasM) {
		super(GeometryType.MULTILINESTRING, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param lineStrings
	 *            list of line strings
	 */
	public MultiLineString(List<LineString> lineStrings) {
		this(GeometryUtils.hasZ(lineStrings), GeometryUtils.hasM(lineStrings));
		setLineStrings(lineStrings);
	}

	/**
	 * Constructor
	 * 
	 * @param lineString
	 *            line string
	 */
	public MultiLineString(LineString lineString) {
		this(lineString.hasZ(), lineString.hasM());
		addLineString(lineString);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param multiLineString
	 *            multi line string to copy
	 */
	public MultiLineString(MultiLineString multiLineString) {
		this(multiLineString.hasZ(), multiLineString.hasM());
		for (LineString lineString : multiLineString.getLineStrings()) {
			addLineString((LineString) lineString.copy());
		}
	}

	/**
	 * Get the line strings
	 * 
	 * @return line strings
	 */
	public List<LineString> getLineStrings() {
		return getCurves();
	}

	/**
	 * Set the line strings
	 * 
	 * @param lineStrings
	 *            line strings
	 */
	public void setLineStrings(List<LineString> lineStrings) {
		setCurves(lineStrings);
	}

	/**
	 * Add a line string
	 * 
	 * @param lineString
	 *            line string
	 */
	public void addLineString(LineString lineString) {
		addCurve(lineString);
	}

	/**
	 * Add line strings
	 * 
	 * @param lineStrings
	 *            line strings
	 */
	public void addLineStrings(List<LineString> lineStrings) {
		addCurves(lineStrings);
	}

	/**
	 * Get the number of line strings
	 * 
	 * @return number of line strings
	 */
	public int numLineStrings() {
		return numCurves();
	}

	/**
	 * Returns the Nth line string
	 * 
	 * @param n
	 *            nth line string to return
	 * @return line string
	 */
	public LineString getLineString(int n) {
		return getCurve(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new MultiLineString(this);
	}

}
