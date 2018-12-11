package mil.nga.sf;

import java.util.List;

/**
 * A restricted form of GeometryCollection where each Geometry in the collection
 * must be of type Curve.
 * 
 * @author osbornb
 * @param <T>
 *            curve type
 */
public abstract class MultiCurve<T extends Curve> extends GeometryCollection<T> {

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
	protected MultiCurve(GeometryType type, boolean hasZ, boolean hasM) {
		super(type, hasZ, hasM);
	}

	/**
	 * Get the curves
	 * 
	 * @return curves
	 */
	public List<T> getCurves() {
		return getGeometries();
	}

	/**
	 * Set the curves
	 * 
	 * @param curves
	 *            curves
	 */
	public void setCurves(List<T> curves) {
		setGeometries(curves);
	}

	/**
	 * Add a curve
	 * 
	 * @param curve
	 *            curve
	 */
	public void addCurve(T curve) {
		addGeometry(curve);
	}

	/**
	 * Add curves
	 * 
	 * @param curves
	 *            curves
	 */
	public void addCurves(List<T> curves) {
		addGeometries(curves);
	}

	/**
	 * Get the number of curves
	 * 
	 * @return number of curves
	 */
	public int numCurves() {
		return numGeometries();
	}

	/**
	 * Returns the Nth curve
	 * 
	 * @param n
	 *            nth line curve to return
	 * @return curve
	 */
	public T getCurve(int n) {
		return getGeometry(n);
	}

	/**
	 * Determine if this Multi Curve is closed for each Curve (start point = end
	 * point)
	 * 
	 * @return true if closed
	 */
	public boolean isClosed() {
		boolean closed = true;
		for (Curve curve : getGeometries()) {
			if (!curve.isClosed()) {
				closed = false;
				break;
			}
		}
		return closed;
	}

}
