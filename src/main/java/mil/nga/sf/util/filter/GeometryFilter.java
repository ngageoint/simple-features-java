package mil.nga.sf.util.filter;

import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryType;

/**
 * Geometry Filter to filter included geometries and modify them during
 * construction
 * 
 * @author osbornb
 * @since 2.0.3
 */
public interface GeometryFilter {

	/**
	 * Filter the geometry
	 * 
	 * @param containingType
	 *            geometry type of the geometry containing this geometry
	 *            element, null if geometry is top level
	 * @param geometry
	 *            geometry, may be modified
	 * @return true if passes filter and geometry should be included
	 */
	public boolean filter(GeometryType containingType, Geometry geometry);

}
