package mil.nga.sf.util.centroid;

import mil.nga.sf.CompoundCurve;
import mil.nga.sf.Curve;
import mil.nga.sf.CurvePolygon;
import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryCollection;
import mil.nga.sf.GeometryType;
import mil.nga.sf.LineString;
import mil.nga.sf.Point;
import mil.nga.sf.Polygon;
import mil.nga.sf.PolyhedralSurface;
import mil.nga.sf.util.GeometryUtils;
import mil.nga.sf.util.SFException;

/**
 * Centroid calculations for geometries in degrees
 * 
 * @author osbornb
 * @since 2.0.5
 */
public class DegreesCentroid {

	/**
	 * Geometry
	 */
	private final Geometry geometry;

	/**
	 * Number of points
	 */
	private int points = 0;

	/**
	 * x sum
	 */
	private double x = 0;

	/**
	 * y sum
	 */
	private double y = 0;

	/**
	 * z sum
	 */
	private double z = 0;

	/**
	 * Get the degree geometry centroid
	 * 
	 * @param geometry
	 *            geometry
	 * @return centroid point
	 */
	public static Point getCentroid(Geometry geometry) {
		return new DegreesCentroid(geometry).getCentroid();
	}

	/**
	 * Constructor
	 * 
	 * @param geometry
	 *            geometry
	 */
	public DegreesCentroid(Geometry geometry) {
		this.geometry = geometry;
	}

	/**
	 * Get the centroid point
	 * 
	 * @return centroid point
	 */
	public Point getCentroid() {

		Point centroid = null;

		if (geometry.getGeometryType() == GeometryType.POINT) {
			centroid = (Point) geometry;
		} else {

			calculate(geometry);

			x = x / points;
			y = y / points;
			z = z / points;

			double centroidLongitude = Math.atan2(y, x);
			double centroidLatitude = Math.atan2(z, Math.sqrt(x * x + y * y));

			centroid = new Point(
					GeometryUtils.radiansToDegrees(centroidLongitude),
					GeometryUtils.radiansToDegrees(centroidLatitude));
		}

		return centroid;
	}

	/**
	 * Add to the centroid calculation for the Geometry
	 * 
	 * @param geometry
	 *            Geometry
	 */
	private void calculate(Geometry geometry) {

		GeometryType geometryType = geometry.getGeometryType();

		switch (geometryType) {

		case GEOMETRY:
			throw new SFException("Unexpected Geometry Type of "
					+ geometryType.name() + " which is abstract");
		case POINT:
			calculatePoint((Point) geometry);
			break;
		case LINESTRING:
		case CIRCULARSTRING:
			calculateLineString((LineString) geometry);
			break;
		case POLYGON:
		case TRIANGLE:
			calculatePolygon((Polygon) geometry);
			break;
		case GEOMETRYCOLLECTION:
		case MULTIPOINT:
		case MULTICURVE:
		case MULTILINESTRING:
		case MULTISURFACE:
		case MULTIPOLYGON:
			calculateGeometryCollection((GeometryCollection<?>) geometry);
			break;
		case COMPOUNDCURVE:
			calculateCompoundCurve((CompoundCurve) geometry);
			break;
		case CURVEPOLYGON:
			calculateCurvePolygon((CurvePolygon<?>) geometry);
			break;
		case CURVE:
			throw new SFException("Unexpected Geometry Type of "
					+ geometryType.name() + " which is abstract");
		case SURFACE:
			throw new SFException("Unexpected Geometry Type of "
					+ geometryType.name() + " which is abstract");
		case POLYHEDRALSURFACE:
		case TIN:
			calculatePolyhedralSurface((PolyhedralSurface) geometry);
			break;
		default:
			throw new SFException(
					"Geometry Type not supported: " + geometryType);
		}
	}

	/**
	 * Add to the centroid calculation for the Point
	 * 
	 * @param point
	 *            Point
	 */
	private void calculatePoint(Point point) {
		double latitude = GeometryUtils.degreesToRadians(point.getY());
		double longitude = GeometryUtils.degreesToRadians(point.getX());
		double cosLatitude = Math.cos(latitude);
		x += cosLatitude * Math.cos(longitude);
		y += cosLatitude * Math.sin(longitude);
		z += Math.sin(latitude);
		points++;
	}

	/**
	 * Add to the centroid calculation for the Line String
	 * 
	 * @param lineString
	 *            Line String
	 */
	private void calculateLineString(LineString lineString) {

		for (Point point : lineString.getPoints()) {
			calculatePoint(point);
		}

	}

	/**
	 * Add to the centroid calculation for the Polygon
	 * 
	 * @param polygon
	 *            Polygon
	 */
	private void calculatePolygon(Polygon polygon) {

		if (polygon.numRings() > 0) {
			LineString exteriorRing = polygon.getExteriorRing();
			int numPoints = exteriorRing.numPoints();
			if (GeometryUtils.closedPolygon(exteriorRing)) {
				numPoints--;
			}
			for (int i = 0; i < numPoints; i++) {
				calculatePoint(exteriorRing.getPoint(i));
			}
		}

	}

	/**
	 * Add to the centroid calculation for the Geometry Collection
	 * 
	 * @param geometryCollection
	 *            Geometry Collection
	 */
	private void calculateGeometryCollection(
			GeometryCollection<?> geometryCollection) {

		for (Geometry geometry : geometryCollection.getGeometries()) {
			calculate(geometry);
		}

	}

	/**
	 * Add to the centroid calculation for the Compound Curve
	 * 
	 * @param compoundCurve
	 *            Compound Curve
	 */
	private void calculateCompoundCurve(CompoundCurve compoundCurve) {

		for (LineString lineString : compoundCurve.getLineStrings()) {
			calculateLineString(lineString);
		}

	}

	/**
	 * Add to the centroid calculation for the Curve Polygon
	 * 
	 * @param curvePolygon
	 *            Curve Polygon
	 */
	private void calculateCurvePolygon(CurvePolygon<?> curvePolygon) {

		for (Curve ring : curvePolygon.getRings()) {
			calculate(ring);
		}

	}

	/**
	 * Add to the centroid calculation for the Polyhedral Surface
	 * 
	 * @param polyhedralSurface
	 *            Polyhedral Surface
	 */
	private void calculatePolyhedralSurface(
			PolyhedralSurface polyhedralSurface) {

		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			calculatePolygon(polygon);
		}

	}

}
