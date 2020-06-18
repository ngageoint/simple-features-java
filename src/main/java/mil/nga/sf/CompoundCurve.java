package mil.nga.sf;

import java.util.ArrayList;
import java.util.List;

import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.sweep.ShamosHoey;

/**
 * Compound Curve, Curve sub type
 * 
 * @author osbornb
 */
public class CompoundCurve extends Curve {

	/**
	 * Serial Version UID
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * List of line strings
	 */
	private List<LineString> lineStrings = new ArrayList<LineString>();

	/**
	 * Constructor
	 */
	public CompoundCurve() {
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
	public CompoundCurve(boolean hasZ, boolean hasM) {
		super(GeometryType.COMPOUNDCURVE, hasZ, hasM);
	}

	/**
	 * Constructor
	 * 
	 * @param lineStrings
	 *            list of line strings
	 */
	public CompoundCurve(List<LineString> lineStrings) {
		this(GeometryUtils.hasZ(lineStrings), GeometryUtils.hasM(lineStrings));
		setLineStrings(lineStrings);
	}

	/**
	 * Constructor
	 * 
	 * @param lineString
	 *            line string
	 */
	public CompoundCurve(LineString lineString) {
		this(lineString.hasZ(), lineString.hasM());
		addLineString(lineString);
	}

	/**
	 * Copy Constructor
	 * 
	 * @param compoundCurve
	 *            compound Curve to copy
	 */
	public CompoundCurve(CompoundCurve compoundCurve) {
		this(compoundCurve.hasZ(), compoundCurve.hasM());
		for (LineString lineString : compoundCurve.getLineStrings()) {
			addLineString((LineString) lineString.copy());
		}
	}

	/**
	 * Get the line strings
	 * 
	 * @return line strings
	 */
	public List<LineString> getLineStrings() {
		return lineStrings;
	}

	/**
	 * Set the line strings
	 * 
	 * @param lineStrings
	 *            line strings
	 */
	public void setLineStrings(List<LineString> lineStrings) {
		this.lineStrings = lineStrings;
	}

	/**
	 * Add a line string
	 * 
	 * @param lineString
	 *            line string
	 */
	public void addLineString(LineString lineString) {
		lineStrings.add(lineString);
		updateZM(lineString);
	}

	/**
	 * Add line strings
	 * 
	 * @param lineStrings
	 *            line strings
	 */
	public void addLineStrings(List<LineString> lineStrings) {
		for (LineString lineString : lineStrings) {
			addLineString(lineString);
		}
	}

	/**
	 * Get the number of line strings
	 * 
	 * @return number of line strings
	 */
	public int numLineStrings() {
		return lineStrings.size();
	}

	/**
	 * Returns the Nth line string
	 * 
	 * @param n
	 *            nth line string to return
	 * @return line string
	 */
	public LineString getLineString(int n) {
		return lineStrings.get(n);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Point startPoint() {
		Point startPoint = null;
		if (!isEmpty()) {
			for (LineString lineString : lineStrings) {
				if (!lineString.isEmpty()) {
					startPoint = lineString.startPoint();
					break;
				}
			}
		}
		return startPoint;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Point endPoint() {
		Point endPoint = null;
		if (!isEmpty()) {
			for (int i = lineStrings.size() - 1; i >= 0; i--) {
				LineString lineString = lineStrings.get(i);
				if (!lineString.isEmpty()) {
					endPoint = lineString.endPoint();
					break;
				}
			}
		}
		return endPoint;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isSimple() {
		return ShamosHoey.simplePolygon(lineStrings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Geometry copy() {
		return new CompoundCurve(this);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean isEmpty() {
		return lineStrings.isEmpty();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result
				+ ((lineStrings == null) ? 0 : lineStrings.hashCode());
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
		CompoundCurve other = (CompoundCurve) obj;
		if (lineStrings == null) {
			if (other.lineStrings != null)
				return false;
		} else if (!lineStrings.equals(other.lineStrings))
			return false;
		return true;
	}

}
