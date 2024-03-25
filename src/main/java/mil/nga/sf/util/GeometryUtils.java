package mil.nga.sf.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import mil.nga.sf.CircularString;
import mil.nga.sf.CompoundCurve;
import mil.nga.sf.Curve;
import mil.nga.sf.CurvePolygon;
import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryCollection;
import mil.nga.sf.GeometryEnvelope;
import mil.nga.sf.GeometryType;
import mil.nga.sf.Line;
import mil.nga.sf.LineString;
import mil.nga.sf.MultiLineString;
import mil.nga.sf.MultiPoint;
import mil.nga.sf.MultiPolygon;
import mil.nga.sf.Point;
import mil.nga.sf.Polygon;
import mil.nga.sf.PolyhedralSurface;
import mil.nga.sf.TIN;
import mil.nga.sf.Triangle;
import mil.nga.sf.util.centroid.CentroidCurve;
import mil.nga.sf.util.centroid.CentroidPoint;
import mil.nga.sf.util.centroid.CentroidSurface;
import mil.nga.sf.util.centroid.DegreesCentroid;

/**
 * Utilities for Geometry objects
 * 
 * @author osbornb
 * @since 1.0.3
 */
public class GeometryUtils {

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(GeometryUtils.class.getName());

	/**
	 * Get the dimension of the Geometry, 0 for points, 1 for curves, 2 for
	 * surfaces. If a collection, the largest dimension is returned.
	 * 
	 * @param geometry
	 *            geometry object
	 * @return dimension (0, 1, or 2)
	 */
	public static int getDimension(Geometry geometry) {

		int dimension = -1;

		GeometryType geometryType = geometry.getGeometryType();
		switch (geometryType) {
		case POINT:
		case MULTIPOINT:
			dimension = 0;
			break;
		case LINESTRING:
		case MULTILINESTRING:
		case CIRCULARSTRING:
		case COMPOUNDCURVE:
			dimension = 1;
			break;
		case POLYGON:
		case CURVEPOLYGON:
		case MULTIPOLYGON:
		case POLYHEDRALSURFACE:
		case TIN:
		case TRIANGLE:
			dimension = 2;
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			List<Geometry> geometries = geomCollection.getGeometries();
			for (Geometry subGeometry : geometries) {
				dimension = Math.max(dimension, getDimension(subGeometry));
			}
			break;
		default:
			throw new SFException("Unsupported Geometry Type: " + geometryType);
		}

		return dimension;
	}

	/**
	 * Get the Pythagorean theorem distance between two points
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return distance
	 */
	public static double distance(Point point1, Point point2) {
		double diffX = point1.getX() - point2.getX();
		double diffY = point1.getY() - point2.getY();
		return Math.sqrt(diffX * diffX + diffY * diffY);
	}

	/**
	 * Get the Pythagorean theorem distance between the line end points
	 * 
	 * @param line
	 *            line
	 * @return distance
	 * @since 2.2.0
	 */
	public static double distance(Line line) {
		return distance(line.startPoint(), line.endPoint());
	}

	/**
	 * Get the distance in meters between two points in degrees using the
	 * Haversine formula
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return distance in meters
	 * @since 2.2.2
	 */
	public static double distanceHaversine(Point point1, Point point2) {
		double lat1 = point1.getY();
		double lon1 = point1.getX();
		double lat2 = point2.getY();
		double lon2 = point2.getX();
		double diffLat = degreesToRadians(lat2 - lat1);
		double diffLon = degreesToRadians(lon2 - lon1);
		double a = Math.sin(diffLat / 2) * Math.sin(diffLat / 2)
				+ Math.cos(degreesToRadians(lat1))
						* Math.cos(degreesToRadians(lat2))
						* Math.sin(diffLon / 2) * Math.sin(diffLon / 2);
		double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
		return GeometryConstants.EARTH_RADIUS * c;
	}

	/**
	 * Get the distance in meters between the end points of a line in degrees
	 * using the Haversine formula
	 * 
	 * @param line
	 *            line
	 * @return distance in meters
	 * @since 2.2.2
	 */
	public static double distanceHaversine(Line line) {
		return distanceHaversine(line.startPoint(), line.endPoint());
	}

	/**
	 * Get the bearing heading in degrees between two points in degrees
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return bearing angle in degrees between 0 and 360
	 * @since 2.2.0
	 */
	public static double bearing(Point point1, Point point2) {
		double y1 = degreesToRadians(point1.getY());
		double y2 = degreesToRadians(point2.getY());
		double xDiff = degreesToRadians(point2.getX() - point1.getX());
		double y = Math.sin(xDiff) * Math.cos(y2);
		double x = Math.cos(y1) * Math.sin(y2)
				- Math.sin(y1) * Math.cos(y2) * Math.cos(xDiff);
		return (radiansToDegrees(Math.atan2(y, x)) + 360) % 360;
	}

	/**
	 * Get the geodesic midpoint in degrees between two points in degrees
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return geodesic midpoint in degrees
	 * @since 2.2.2
	 */
	public static Point geodesicMidpoint(Point point1, Point point2) {
		Point point1Radians = degreesToRadians(point1);
		Point point2Radians = degreesToRadians(point2);
		Point midpointRadians = geodesicMidpointRadians(point1Radians,
				point2Radians);
		return radiansToDegrees(midpointRadians);
	}

	/**
	 * Get the geodesic midpoint in radians between two points in radians
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return geodesic midpoint in radians
	 * @since 2.2.2
	 */
	public static Point geodesicMidpointRadians(Point point1, Point point2) {

		double xDiff = point2.getX() - point1.getX();
		double y1 = point1.getY();
		double y2 = point2.getY();
		double x1 = point1.getX();

		double bx = Math.cos(y2) * Math.cos(xDiff);
		double by = Math.cos(y2) * Math.sin(xDiff);

		double y = Math.atan2(Math.sin(y1) + Math.sin(y2),
				Math.sqrt((Math.cos(y1) + bx) * (Math.cos(y1) + bx) + by * by));
		double x = x1 + Math.atan2(by, Math.cos(y1) + bx);
		Point midpoint = new Point(x, y);

		return midpoint;
	}

	/**
	 * Get the bearing heading in degrees between line end points in degrees
	 * 
	 * @param line
	 *            line
	 * @return bearing angle in degrees between 0 inclusively and 360
	 *         exclusively
	 * @since 2.2.0
	 */
	public static double bearing(Line line) {
		return bearing(line.startPoint(), line.endPoint());
	}

	/**
	 * Determine if the bearing is in any north direction
	 * 
	 * @param bearing
	 *            bearing angle in degrees
	 * @return true if north bearing
	 * @since 2.2.0
	 */
	public static boolean isNorthBearing(double bearing) {
		bearing %= 360.0;
		return bearing < GeometryConstants.BEARING_EAST
				|| bearing > GeometryConstants.BEARING_WEST;
	}

	/**
	 * Determine if the bearing is in any east direction
	 * 
	 * @param bearing
	 *            bearing angle in degrees
	 * @return true if east bearing
	 * @since 2.2.0
	 */
	public static boolean isEastBearing(double bearing) {
		bearing %= 360.0;
		return bearing > GeometryConstants.BEARING_NORTH
				&& bearing < GeometryConstants.BEARING_SOUTH;
	}

	/**
	 * Determine if the bearing is in any south direction
	 * 
	 * @param bearing
	 *            bearing angle in degrees
	 * @return true if south bearing
	 * @since 2.2.0
	 */
	public static boolean isSouthBearing(double bearing) {
		bearing %= 360.0;
		return bearing > GeometryConstants.BEARING_EAST
				&& bearing < GeometryConstants.BEARING_WEST;
	}

	/**
	 * Determine if the bearing is in any west direction
	 * 
	 * @param bearing
	 *            bearing angle in degrees
	 * @return true if west bearing
	 * @since 2.2.0
	 */
	public static boolean isWestBearing(double bearing) {
		return (bearing % 360.0) > GeometryConstants.BEARING_SOUTH;
	}

	/**
	 * Convert degrees to radians
	 * 
	 * @param degrees
	 *            degrees
	 * @return radians
	 * @since 2.2.0
	 */
	public static double degreesToRadians(double degrees) {
		return degrees * GeometryConstants.DEGREES_TO_RADIANS;
	}

	/**
	 * Convert radians to degrees
	 * 
	 * @param radians
	 *            radians
	 * @return degrees
	 * @since 2.2.0
	 */
	public static double radiansToDegrees(double radians) {
		return radians * GeometryConstants.RADIANS_TO_DEGREES;
	}

	/**
	 * Convert point in degrees to radians
	 * 
	 * @param point
	 *            point in degrees
	 * @return point in radians
	 * @since 2.2.2
	 */
	public static Point degreesToRadians(Point point) {
		double x = degreesToRadians(point.getX());
		double y = degreesToRadians(point.getY());
		return new Point(x, y);
	}

	/**
	 * Convert point in radians to degrees
	 * 
	 * @param point
	 *            point in radians
	 * @return point in degrees
	 * @since 2.2.2
	 */
	public static Point radiansToDegrees(Point point) {
		double x = radiansToDegrees(point.getX());
		double y = radiansToDegrees(point.getY());
		return new Point(x, y);
	}

	/**
	 * Get the centroid point of a 2 dimensional representation of the Geometry
	 * (balancing point of a 2d cutout of the geometry). Only the x and y
	 * coordinate of the resulting point are calculated and populated. The
	 * resulting {@link Point#getZ()} and {@link Point#getM()} methods will
	 * always return null.
	 * 
	 * @param geometry
	 *            geometry object
	 * @return centroid point
	 */
	public static Point getCentroid(Geometry geometry) {
		Point centroid = null;
		int dimension = getDimension(geometry);
		switch (dimension) {
		case 0:
			CentroidPoint point = new CentroidPoint(geometry);
			centroid = point.getCentroid();
			break;
		case 1:
			CentroidCurve curve = new CentroidCurve(geometry);
			centroid = curve.getCentroid();
			break;
		case 2:
			CentroidSurface surface = new CentroidSurface(geometry);
			centroid = surface.getCentroid();
			break;
		}
		return centroid;
	}

	/**
	 * Get the geographic centroid point of a 2 dimensional representation of
	 * the degree unit Geometry. Only the x and y coordinate of the resulting
	 * point are calculated and populated. The resulting {@link Point#getZ()}
	 * and {@link Point#getM()} methods will always return null.
	 * 
	 * @param geometry
	 *            geometry object
	 * @return centroid point
	 * @since 2.0.5
	 */
	public static Point getDegreesCentroid(Geometry geometry) {
		return DegreesCentroid.getCentroid(geometry);
	}

	/**
	 * Minimize the WGS84 geometry using the shortest x distance between each
	 * connected set of points. Resulting x values will be in the range: -540.0
	 * &lt;= x &lt;= 540.0
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void minimizeWGS84(Geometry geometry) {
		minimize(geometry, GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH);
	}

	/**
	 * Minimize the Web Mercator geometry using the shortest x distance between
	 * each connected set of points. Resulting x values will be in the range:
	 * -60112525.028367732 &lt;= x &lt;= 60112525.028367732
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void minimizeWebMercator(Geometry geometry) {
		minimize(geometry, GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH);
	}

	/**
	 * Minimize the geometry using the shortest x distance between each
	 * connected set of points. The resulting geometry point x values will be in
	 * the range: (3 * min value &lt;= x &lt;= 3 * max value
	 *
	 * Example: For WGS84 provide a max x of
	 * {@link GeometryConstants#WGS84_HALF_WORLD_LON_WIDTH}. Resulting x values
	 * will be in the range: -540.0 &lt;= x &lt;= 540.0
	 *
	 * Example: For web mercator provide a world width of
	 * {@link GeometryConstants#WEB_MERCATOR_HALF_WORLD_WIDTH}. Resulting x
	 * values will be in the range: -60112525.028367732 &lt;= x &lt;=
	 * 60112525.028367732
	 *
	 * @param geometry
	 *            geometry
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	public static void minimizeGeometry(Geometry geometry, double maxX) {
		minimize(geometry, maxX);
	}

	/**
	 * Minimize the geometry using the shortest x distance between each
	 * connected set of points. The resulting geometry point x values will be in
	 * the range: (3 * min value &lt;= x &lt;= 3 * max value
	 *
	 * Example: For WGS84 provide a max x of
	 * {@link GeometryConstants#WGS84_HALF_WORLD_LON_WIDTH}. Resulting x values
	 * will be in the range: -540.0 &lt;= x &lt;= 540.0
	 *
	 * Example: For web mercator provide a world width of
	 * {@link GeometryConstants#WEB_MERCATOR_HALF_WORLD_WIDTH}. Resulting x
	 * values will be in the range: -60112525.028367732 &lt;= x &lt;=
	 * 60112525.028367732
	 *
	 * @param geometry
	 *            geometry
	 * @param maxX
	 *            max positive x value in the geometry projection
	 * @since 2.2.0
	 */
	public static void minimize(Geometry geometry, double maxX) {

		GeometryType geometryType = geometry.getGeometryType();
		switch (geometryType) {
		case LINESTRING:
			minimize((LineString) geometry, maxX);
			break;
		case POLYGON:
			minimize((Polygon) geometry, maxX);
			break;
		case MULTILINESTRING:
			minimize((MultiLineString) geometry, maxX);
			break;
		case MULTIPOLYGON:
			minimize((MultiPolygon) geometry, maxX);
			break;
		case CIRCULARSTRING:
			minimize((CircularString) geometry, maxX);
			break;
		case COMPOUNDCURVE:
			minimize((CompoundCurve) geometry, maxX);
			break;
		case CURVEPOLYGON:
			@SuppressWarnings("unchecked")
			CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
			minimize(curvePolygon, maxX);
			break;
		case POLYHEDRALSURFACE:
			minimize((PolyhedralSurface) geometry, maxX);
			break;
		case TIN:
			minimize((TIN) geometry, maxX);
			break;
		case TRIANGLE:
			minimize((Triangle) geometry, maxX);
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			for (Geometry subGeometry : geomCollection.getGeometries()) {
				minimize(subGeometry, maxX);
			}
			break;
		default:
			break;

		}
	}

	/**
	 * Minimize the line string
	 * 
	 * @param lineString
	 *            line string
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(LineString lineString, double maxX) {

		List<Point> points = lineString.getPoints();
		if (points.size() > 1) {
			Point point = points.get(0);
			for (int i = 1; i < points.size(); i++) {
				Point nextPoint = points.get(i);
				if (point.getX() < nextPoint.getX()) {
					if (nextPoint.getX() - point.getX() > point.getX()
							- nextPoint.getX() + (maxX * 2.0)) {
						nextPoint.setX(nextPoint.getX() - (maxX * 2.0));
					}
				} else if (point.getX() > nextPoint.getX()) {
					if (point.getX() - nextPoint.getX() > nextPoint.getX()
							- point.getX() + (maxX * 2.0)) {
						nextPoint.setX(nextPoint.getX() + (maxX * 2.0));
					}
				}
			}
		}
	}

	/**
	 * Minimize the multi line string
	 * 
	 * @param multiLineString
	 *            multi line string
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(MultiLineString multiLineString, double maxX) {

		List<LineString> lineStrings = multiLineString.getLineStrings();
		for (LineString lineString : lineStrings) {
			minimize(lineString, maxX);
		}
	}

	/**
	 * Minimize the polygon
	 * 
	 * @param polygon
	 *            polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(Polygon polygon, double maxX) {

		for (LineString ring : polygon.getRings()) {
			minimize(ring, maxX);
		}
	}

	/**
	 * Minimize the multi polygon
	 * 
	 * @param multiPolygon
	 *            multi polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(MultiPolygon multiPolygon, double maxX) {

		List<Polygon> polygons = multiPolygon.getPolygons();
		for (Polygon polygon : polygons) {
			minimize(polygon, maxX);
		}
	}

	/**
	 * Minimize the compound curve
	 * 
	 * @param compoundCurve
	 *            compound curve
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(CompoundCurve compoundCurve, double maxX) {

		for (LineString lineString : compoundCurve.getLineStrings()) {
			minimize(lineString, maxX);
		}
	}

	/**
	 * Minimize the curve polygon
	 * 
	 * @param curvePolygon
	 *            curve polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(CurvePolygon<Curve> curvePolygon,
			double maxX) {

		for (Curve ring : curvePolygon.getRings()) {
			minimize(ring, maxX);
		}
	}

	/**
	 * Minimize the polyhedral surface
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void minimize(PolyhedralSurface polyhedralSurface,
			double maxX) {

		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			minimize(polygon, maxX);
		}
	}

	/**
	 * Normalize the WGS84 geometry using the shortest x distance between each
	 * connected set of points. Resulting x values will be in the range: -180.0
	 * &lt;= x &lt;= 180.0
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void normalizeWGS84(Geometry geometry) {
		normalize(geometry, GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH);
	}

	/**
	 * Normalize the Web Mercator geometry using the shortest x distance between
	 * each connected set of points. Resulting x values will be in the range:
	 * -20037508.342789244 &lt;= x &lt;= 20037508.342789244
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void normalizeWebMercator(Geometry geometry) {
		normalize(geometry, GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH);
	}

	/**
	 * Normalize the geometry so all points outside of the min and max value
	 * range are adjusted to fall within the range.
	 *
	 * Example: For WGS84 provide a max x of
	 * {@link GeometryConstants#WGS84_HALF_WORLD_LON_WIDTH}. Resulting x values
	 * will be in the range: -180.0 &lt;= x &lt;= 180.0
	 *
	 * Example: For web mercator provide a world width of
	 * {@link GeometryConstants#WEB_MERCATOR_HALF_WORLD_WIDTH}. Resulting x
	 * values will be in the range: -20037508.342789244 &lt;= x &lt;=
	 * 20037508.342789244
	 *
	 * @param geometry
	 *            geometry
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	public static void normalizeGeometry(Geometry geometry, double maxX) {
		normalize(geometry, maxX);
	}

	/**
	 * Normalize the geometry so all points outside of the min and max value
	 * range are adjusted to fall within the range.
	 *
	 * Example: For WGS84 provide a max x of
	 * {@link GeometryConstants#WGS84_HALF_WORLD_LON_WIDTH}. Resulting x values
	 * will be in the range: -180.0 &lt;= x &lt;= 180.0
	 *
	 * Example: For web mercator provide a world width of
	 * {@link GeometryConstants#WEB_MERCATOR_HALF_WORLD_WIDTH}. Resulting x
	 * values will be in the range: -20037508.342789244 &lt;= x &lt;=
	 * 20037508.342789244
	 *
	 * @param geometry
	 *            geometry
	 * @param maxX
	 *            max positive x value in the geometry projection
	 * @since 2.2.0
	 */
	public static void normalize(Geometry geometry, double maxX) {

		GeometryType geometryType = geometry.getGeometryType();
		switch (geometryType) {
		case POINT:
			normalize((Point) geometry, maxX);
			break;
		case LINESTRING:
			normalize((LineString) geometry, maxX);
			break;
		case POLYGON:
			normalize((Polygon) geometry, maxX);
			break;
		case MULTIPOINT:
			normalize((MultiPoint) geometry, maxX);
			break;
		case MULTILINESTRING:
			normalize((MultiLineString) geometry, maxX);
			break;
		case MULTIPOLYGON:
			normalize((MultiPolygon) geometry, maxX);
			break;
		case CIRCULARSTRING:
			normalize((CircularString) geometry, maxX);
			break;
		case COMPOUNDCURVE:
			normalize((CompoundCurve) geometry, maxX);
			break;
		case CURVEPOLYGON:
			@SuppressWarnings("unchecked")
			CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
			normalize(curvePolygon, maxX);
			break;
		case POLYHEDRALSURFACE:
			normalize((PolyhedralSurface) geometry, maxX);
			break;
		case TIN:
			normalize((TIN) geometry, maxX);
			break;
		case TRIANGLE:
			normalize((Triangle) geometry, maxX);
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			for (Geometry subGeometry : geomCollection.getGeometries()) {
				normalize(subGeometry, maxX);
			}
			break;
		default:
			break;

		}

	}

	/**
	 * Normalize the point
	 * 
	 * @param point
	 *            point
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(Point point, double maxX) {
		point.setX(normalize(point.getX(), maxX));
	}

	/**
	 * Normalize the x value
	 * 
	 * @param x
	 *            x value
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static double normalize(double x, double maxX) {
		if (x < -maxX) {
			x = x + (maxX * 2.0);
		} else if (x > maxX) {
			x = x - (maxX * 2.0);
		}
		return x;
	}

	/**
	 * Normalize the multi point
	 * 
	 * @param multiPoint
	 *            multi point
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(MultiPoint multiPoint, double maxX) {

		List<Point> points = multiPoint.getPoints();
		for (Point point : points) {
			normalize(point, maxX);
		}
	}

	/**
	 * Normalize the line string
	 * 
	 * @param lineString
	 *            line string
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(LineString lineString, double maxX) {

		for (Point point : lineString.getPoints()) {
			normalize(point, maxX);
		}
	}

	/**
	 * Normalize the multi line string
	 * 
	 * @param multiLineString
	 *            multi line string
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(MultiLineString multiLineString,
			double maxX) {

		List<LineString> lineStrings = multiLineString.getLineStrings();
		for (LineString lineString : lineStrings) {
			normalize(lineString, maxX);
		}
	}

	/**
	 * Normalize the polygon
	 * 
	 * @param polygon
	 *            polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(Polygon polygon, double maxX) {

		for (LineString ring : polygon.getRings()) {
			normalize(ring, maxX);
		}
	}

	/**
	 * Normalize the multi polygon
	 * 
	 * @param multiPolygon
	 *            multi polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(MultiPolygon multiPolygon, double maxX) {

		List<Polygon> polygons = multiPolygon.getPolygons();
		for (Polygon polygon : polygons) {
			normalize(polygon, maxX);
		}
	}

	/**
	 * Normalize the compound curve
	 * 
	 * @param compoundCurve
	 *            compound curve
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(CompoundCurve compoundCurve, double maxX) {

		for (LineString lineString : compoundCurve.getLineStrings()) {
			normalize(lineString, maxX);
		}
	}

	/**
	 * Normalize the curve polygon
	 * 
	 * @param curvePolygon
	 *            curve polygon
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(CurvePolygon<Curve> curvePolygon,
			double maxX) {

		for (Curve ring : curvePolygon.getRings()) {
			normalize(ring, maxX);
		}
	}

	/**
	 * Normalize the polyhedral surface
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface
	 * @param maxX
	 *            max positive x value in the geometry projection
	 */
	private static void normalize(PolyhedralSurface polyhedralSurface,
			double maxX) {

		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			normalize(polygon, maxX);
		}
	}

	/**
	 * Simplify the ordered points (representing a line, polygon, etc) using the
	 * Douglas Peucker algorithm to create a similar curve with fewer points.
	 * Points should be in a meters unit type projection. The tolerance is the
	 * minimum tolerated distance between consecutive points.
	 *
	 * @param points
	 *            geometry points
	 * @param tolerance
	 *            minimum tolerance in meters for consecutive points
	 * @return simplified points
	 * @since 1.0.4
	 */
	public static List<Point> simplifyPoints(List<Point> points,
			double tolerance) {
		return simplifyPoints(points, tolerance, 0, points.size() - 1);
	}

	/**
	 * Simplify the ordered points (representing a line, polygon, etc) using the
	 * Douglas Peucker algorithm to create a similar curve with fewer points.
	 * Points should be in a meters unit type projection. The tolerance is the
	 * minimum tolerated distance between consecutive points.
	 *
	 * @param points
	 *            geometry points
	 * @param tolerance
	 *            minimum tolerance in meters for consecutive points
	 * @param startIndex
	 *            start index
	 * @param endIndex
	 *            end index
	 * @return simplified points
	 */
	private static List<Point> simplifyPoints(List<Point> points,
			double tolerance, int startIndex, int endIndex) {

		List<Point> result = null;

		double dmax = 0.0;
		int index = 0;

		Point startPoint = points.get(startIndex);
		Point endPoint = points.get(endIndex);

		for (int i = startIndex + 1; i < endIndex; i++) {
			Point point = points.get(i);

			double d = perpendicularDistance(point, startPoint, endPoint);

			if (d > dmax) {
				index = i;
				dmax = d;
			}
		}

		if (dmax > tolerance) {

			List<Point> recResults1 = simplifyPoints(points, tolerance,
					startIndex, index);
			List<Point> recResults2 = simplifyPoints(points, tolerance, index,
					endIndex);

			result = recResults1.subList(0, recResults1.size() - 1);
			result.addAll(recResults2);

		} else {
			result = new ArrayList<Point>();
			result.add(startPoint);
			result.add(endPoint);
		}

		return result;
	}

	/**
	 * Create a geodesic path of a line string in degrees with a max distance
	 * between any two path points
	 * 
	 * @param lineString
	 *            line string in degrees
	 * @param maxDistance
	 *            max distance allowed between path points
	 * @return geodesic path of points
	 * @since 2.2.2
	 */
	public static List<Point> geodesicPath(LineString lineString,
			double maxDistance) {
		return geodesicPath(lineString.getPoints(), maxDistance);
	}

	/**
	 * Create a geodesic path of a points in degrees with a max distance between
	 * any two path points
	 * 
	 * @param points
	 *            points in degrees
	 * @param maxDistance
	 *            max distance allowed between path points
	 * @return geodesic path of points
	 * @since 2.2.2
	 */
	public static List<Point> geodesicPath(List<Point> points,
			double maxDistance) {
		List<Point> path = new ArrayList<>();
		if (!points.isEmpty()) {
			path.add(points.get(0));
			for (int i = 0; i < points.size() - 1; i++) {
				List<Point> subPath = geodesicPath(points.get(i),
						points.get(i + 1), maxDistance);
				path.addAll(subPath.subList(1, subPath.size()));
			}
		}
		return path;
	}

	/**
	 * Create a geodesic path between the two points in degrees with a max
	 * distance between any two path points
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @param maxDistance
	 *            max distance allowed between path points
	 * @return geodesic path of points
	 * @since 2.2.2
	 */
	public static List<Point> geodesicPath(Point point1, Point point2,
			double maxDistance) {
		List<Point> path = new ArrayList<>();
		path.add(point1);
		geodesicPath(point1, point2, maxDistance, path);
		path.add(point2);
		return path;
	}

	/**
	 * Populate a geodesic path between the two points in degrees with a max
	 * distance between any two path points
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @param maxDistance
	 *            max distance allowed between path points
	 * @param path
	 *            geodesic path of points
	 */
	private static void geodesicPath(Point point1, Point point2,
			double maxDistance, List<Point> path) {
		double distance = distanceHaversine(point1, point2);
		if (distance > maxDistance) {
			Point midpoint = geodesicMidpoint(point1, point2);
			geodesicPath(point1, midpoint, maxDistance, path);
			path.add(midpoint);
			geodesicPath(midpoint, point2, maxDistance, path);
		}
	}

	/**
	 * Expand the vertical bounds of a geometry envelope in degrees by including
	 * geodesic bounds
	 * 
	 * @param envelope
	 *            geometry envelope in degrees
	 * @return geodesic expanded geometry envelope in degrees
	 * @since 2.2.2
	 */
	public static GeometryEnvelope geodesicEnvelope(GeometryEnvelope envelope) {
		GeometryEnvelope geodesic = envelope.copy();
		if (envelope.getMinY() < 0) {
			Point left = envelope.getBottomLeft();
			Point right = envelope.getBottomRight();
			Point midpoint = geodesicMidpoint(left, right);
			double y = midpoint.getY();
			if (y < geodesic.getMinY()) {
				geodesic.setMinY(y);
			}

		}
		if (envelope.getMaxY() > 0) {
			Point left = envelope.getTopLeft();
			Point right = envelope.getTopRight();
			Point midpoint = geodesicMidpoint(left, right);
			double y = midpoint.getY();
			if (y > geodesic.getMaxY()) {
				geodesic.setMaxY(y);
			}

		}
		return geodesic;
	}

	/**
	 * Calculate the perpendicular distance between the point and the line
	 * represented by the start and end points. Points should be in a meters
	 * unit type projection.
	 *
	 * @param point
	 *            point
	 * @param lineStart
	 *            point representing the line start
	 * @param lineEnd
	 *            point representing the line end
	 * @return distance in meters
	 * @since 1.0.4
	 */
	public static double perpendicularDistance(Point point, Point lineStart,
			Point lineEnd) {

		double x = point.getX();
		double y = point.getY();
		double startX = lineStart.getX();
		double startY = lineStart.getY();
		double endX = lineEnd.getX();
		double endY = lineEnd.getY();

		double vX = endX - startX;
		double vY = endY - startY;
		double wX = x - startX;
		double wY = y - startY;
		double c1 = wX * vX + wY * vY;
		double c2 = vX * vX + vY * vY;

		double x2;
		double y2;
		if (c1 <= 0) {
			x2 = startX;
			y2 = startY;
		} else if (c2 <= c1) {
			x2 = endX;
			y2 = endY;
		} else {
			double b = c1 / c2;
			x2 = startX + b * vX;
			y2 = startY + b * vY;
		}

		double distance = Math.sqrt(Math.pow(x2 - x, 2) + Math.pow(y2 - y, 2));

		return distance;
	}

	/**
	 * Check if the point is in the polygon
	 * 
	 * @param point
	 *            point
	 * @param polygon
	 *            polygon
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, Polygon polygon) {
		return pointInPolygon(point, polygon,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is in the polygon
	 * 
	 * @param point
	 *            point
	 * @param polygon
	 *            polygon
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, Polygon polygon,
			double epsilon) {

		boolean contains = false;
		List<LineString> rings = polygon.getRings();
		if (!rings.isEmpty()) {
			contains = pointInPolygon(point, rings.get(0), epsilon);
			if (contains) {
				// Check the holes
				for (int i = 1; i < rings.size(); i++) {
					if (pointInPolygon(point, rings.get(i), epsilon)) {
						contains = false;
						break;
					}
				}
			}
		}

		return contains;
	}

	/**
	 * Check if the point is in the polygon ring
	 * 
	 * @param point
	 *            point
	 * @param ring
	 *            polygon ring
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, LineString ring) {
		return pointInPolygon(point, ring,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is in the polygon ring
	 * 
	 * @param point
	 *            point
	 * @param ring
	 *            polygon ring
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, LineString ring,
			double epsilon) {
		return pointInPolygon(point, ring.getPoints(), epsilon);
	}

	/**
	 * Check if the point is in the polygon points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            polygon points
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, List<Point> points) {
		return pointInPolygon(point, points,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is in the polygon points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            polygon points
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if in the polygon
	 * @since 1.0.5
	 */
	public static boolean pointInPolygon(Point point, List<Point> points,
			double epsilon) {

		boolean contains = false;

		int i = 0;
		int j = points.size() - 1;
		if (closedPolygon(points)) {
			j = i++;
		}

		for (; i < points.size(); j = i++) {
			Point point1 = points.get(i);
			Point point2 = points.get(j);

			// Shortcut check if polygon contains the point within tolerance
			if (Math.abs(point1.getX() - point.getX()) <= epsilon
					&& Math.abs(point1.getY() - point.getY()) <= epsilon) {
				contains = true;
				break;
			}

			if (((point1.getY() > point.getY()) != (point2.getY() > point
					.getY()))
					&& (point.getX() < (point2.getX() - point1.getX())
							* (point.getY() - point1.getY())
							/ (point2.getY() - point1.getY())
							+ point1.getX())) {
				contains = !contains;
			}
		}

		if (!contains) {
			// Check the polygon edges
			contains = pointOnPolygonEdge(point, points);
		}

		return contains;
	}

	/**
	 * Check if the point is on the polygon edge
	 * 
	 * @param point
	 *            point
	 * @param polygon
	 *            polygon
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, Polygon polygon) {
		return pointOnPolygonEdge(point, polygon,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the polygon edge
	 * 
	 * @param point
	 *            point
	 * @param polygon
	 *            polygon
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, Polygon polygon,
			double epsilon) {
		return polygon.numRings() > 0 && pointOnPolygonEdge(point,
				polygon.getRings().get(0), epsilon);
	}

	/**
	 * Check if the point is on the polygon ring edge
	 * 
	 * @param point
	 *            point
	 * @param ring
	 *            polygon ring
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, LineString ring) {
		return pointOnPolygonEdge(point, ring,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the polygon ring edge
	 * 
	 * @param point
	 *            point
	 * @param ring
	 *            polygon ring
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, LineString ring,
			double epsilon) {
		return pointOnPolygonEdge(point, ring.getPoints(), epsilon);
	}

	/**
	 * Check if the point is on the polygon ring edge points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            polygon points
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, List<Point> points) {
		return pointOnPolygonEdge(point, points,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the polygon ring edge points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            polygon points
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the polygon edge
	 * @since 1.0.5
	 */
	public static boolean pointOnPolygonEdge(Point point, List<Point> points,
			double epsilon) {
		return pointOnPath(point, points, epsilon, !closedPolygon(points));
	}

	/**
	 * Check if the polygon outer ring is explicitly closed, where the first and
	 * last point are the same
	 * 
	 * @param polygon
	 *            polygon
	 * @return true if the first and last points are the same
	 * @since 1.0.5
	 */
	public static boolean closedPolygon(Polygon polygon) {
		return polygon.numRings() > 0
				&& closedPolygon(polygon.getRings().get(0));
	}

	/**
	 * Check if the polygon ring is explicitly closed, where the first and last
	 * point are the same
	 * 
	 * @param ring
	 *            polygon ring
	 * @return true if the first and last points are the same
	 * @since 1.0.5
	 */
	public static boolean closedPolygon(LineString ring) {
		return closedPolygon(ring.getPoints());
	}

	/**
	 * Check if the polygon ring points are explicitly closed, where the first
	 * and last point are the same
	 * 
	 * @param points
	 *            polygon ring points
	 * @return true if the first and last points are the same
	 * @since 1.0.5
	 */
	public static boolean closedPolygon(List<Point> points) {
		boolean closed = false;
		if (!points.isEmpty()) {
			Point first = points.get(0);
			Point last = points.get(points.size() - 1);
			closed = first.getX() == last.getX() && first.getY() == last.getY();
		}
		return closed;
	}

	/**
	 * Check if the point is on the line
	 * 
	 * @param point
	 *            point
	 * @param line
	 *            line
	 * @return true if on the line
	 * @since 1.0.5
	 */
	public static boolean pointOnLine(Point point, LineString line) {
		return pointOnLine(point, line, GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the line
	 * 
	 * @param point
	 *            point
	 * @param line
	 *            line
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the line
	 * @since 1.0.5
	 */
	public static boolean pointOnLine(Point point, LineString line,
			double epsilon) {
		return pointOnLine(point, line.getPoints(), epsilon);
	}

	/**
	 * Check if the point is on the line represented by the points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            line points
	 * @return true if on the line
	 * @since 1.0.5
	 */
	public static boolean pointOnLine(Point point, List<Point> points) {
		return pointOnLine(point, points,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the line represented by the points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            line points
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the line
	 * @since 1.0.5
	 */
	public static boolean pointOnLine(Point point, List<Point> points,
			double epsilon) {
		return pointOnPath(point, points, epsilon, false);
	}

	/**
	 * Check if the point is on the path between point 1 and point 2
	 * 
	 * @param point
	 *            point
	 * @param point1
	 *            path point 1
	 * @param point2
	 *            path point 2
	 * @return true if on the path
	 * @since 1.0.5
	 */
	public static boolean pointOnPath(Point point, Point point1, Point point2) {
		return pointOnPath(point, point1, point2,
				GeometryConstants.DEFAULT_LINE_EPSILON);
	}

	/**
	 * Check if the point is on the path between point 1 and point 2
	 * 
	 * @param point
	 *            point
	 * @param point1
	 *            path point 1
	 * @param point2
	 *            path point 2
	 * @param epsilon
	 *            epsilon line tolerance
	 * @return true if on the path
	 * @since 1.0.5
	 */
	public static boolean pointOnPath(Point point, Point point1, Point point2,
			double epsilon) {

		boolean contains = false;

		double x21 = point2.getX() - point1.getX();
		double y21 = point2.getY() - point1.getY();
		double xP1 = point.getX() - point1.getX();
		double yP1 = point.getY() - point1.getY();

		double dp = xP1 * x21 + yP1 * y21;
		if (dp >= 0.0) {

			double lengthP1 = xP1 * xP1 + yP1 * yP1;
			double length21 = x21 * x21 + y21 * y21;

			if (lengthP1 <= length21) {
				contains = Math.abs(dp * dp - lengthP1 * length21) <= epsilon;
			}
		}

		return contains;
	}

	/**
	 * Check if the point is on the path between the points
	 * 
	 * @param point
	 *            point
	 * @param points
	 *            path points
	 * @param epsilon
	 *            epsilon line tolerance
	 * @param circular
	 *            true if a path exists between the first and last point (a non
	 *            explicitly closed polygon)
	 * @return true if on the path
	 */
	private static boolean pointOnPath(Point point, List<Point> points,
			double epsilon, boolean circular) {

		boolean onPath = false;

		int i = 0;
		int j = points.size() - 1;
		if (!circular) {
			j = i++;
		}

		for (; i < points.size(); j = i++) {
			Point point1 = points.get(i);
			Point point2 = points.get(j);
			if (pointOnPath(point, point1, point2, epsilon)) {
				onPath = true;
				break;
			}
		}

		return onPath;
	}

	/**
	 * Get the point intersection between two lines
	 * 
	 * @param line1
	 *            first line
	 * @param line2
	 *            second line
	 * @return intersection point or null if no intersection
	 * @since 2.1.0
	 */
	public static Point intersection(Line line1, Line line2) {
		return intersection(line1.startPoint(), line1.endPoint(),
				line2.startPoint(), line2.endPoint());
	}

	/**
	 * Get the point intersection between end points of two lines
	 * 
	 * @param line1Point1
	 *            first point of the first line
	 * @param line1Point2
	 *            second point of the first line
	 * @param line2Point1
	 *            first point of the second line
	 * @param line2Point2
	 *            second point of the second line
	 * @return intersection point or null if no intersection
	 * @since 2.1.0
	 */
	public static Point intersection(Point line1Point1, Point line1Point2,
			Point line2Point1, Point line2Point2) {

		Point intersection = null;

		double a1 = line1Point2.getY() - line1Point1.getY();
		double b1 = line1Point1.getX() - line1Point2.getX();
		double c1 = a1 * (line1Point1.getX()) + b1 * (line1Point1.getY());

		double a2 = line2Point2.getY() - line2Point1.getY();
		double b2 = line2Point1.getX() - line2Point2.getX();
		double c2 = a2 * (line2Point1.getX()) + b2 * (line2Point1.getY());

		double determinant = a1 * b2 - a2 * b1;

		if (determinant != 0) {
			double x = (b2 * c1 - b1 * c2) / determinant;
			double y = (a1 * c2 - a2 * c1) / determinant;
			intersection = new Point(x, y);
		}

		return intersection;
	}

	/**
	 * Convert a geometry in degrees to a geometry in meters
	 * 
	 * @param geometry
	 *            geometry in degrees
	 * @return geometry in meters
	 * @since 2.2.0
	 */
	public static Geometry degreesToMeters(Geometry geometry) {

		Geometry meters = null;

		switch (geometry.getGeometryType()) {
		case POINT:
			meters = degreesToMeters((Point) geometry);
			break;
		case LINESTRING:
			meters = degreesToMeters((LineString) geometry);
			break;
		case POLYGON:
			meters = degreesToMeters((Polygon) geometry);
			break;
		case MULTIPOINT:
			meters = degreesToMeters((MultiPoint) geometry);
			break;
		case MULTILINESTRING:
			meters = degreesToMeters((MultiLineString) geometry);
			break;
		case MULTIPOLYGON:
			meters = degreesToMeters((MultiPolygon) geometry);
			break;
		case CIRCULARSTRING:
			meters = degreesToMeters((CircularString) geometry);
			break;
		case COMPOUNDCURVE:
			meters = degreesToMeters((CompoundCurve) geometry);
			break;
		case CURVEPOLYGON:
			@SuppressWarnings("unchecked")
			CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
			meters = degreesToMeters(curvePolygon);
			break;
		case POLYHEDRALSURFACE:
			meters = degreesToMeters((PolyhedralSurface) geometry);
			break;
		case TIN:
			meters = degreesToMeters((TIN) geometry);
			break;
		case TRIANGLE:
			meters = degreesToMeters((Triangle) geometry);
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			GeometryCollection<Geometry> metersCollection = new GeometryCollection<>();
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			for (Geometry subGeometry : geomCollection.getGeometries()) {
				metersCollection.addGeometry(degreesToMeters(subGeometry));
			}
			meters = metersCollection;
			break;
		default:
			break;

		}

		return meters;
	}

	/**
	 * Convert a point in degrees to a point in meters
	 * 
	 * @param point
	 *            point in degrees
	 * @return point in meters
	 * @since 2.1.0
	 */
	public static Point degreesToMeters(Point point) {
		Point value = degreesToMeters(point.getX(), point.getY());
		value.setZ(point.getZ());
		value.setM(point.getM());
		return value;
	}

	/**
	 * Convert a coordinate in degrees to a point in meters
	 * 
	 * @param x
	 *            x value in degrees
	 * @param y
	 *            y value in degrees
	 * @return point in meters
	 * @since 2.1.0
	 */
	public static Point degreesToMeters(double x, double y) {
		x = normalize(x, GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH);
		y = Math.min(y, GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT);
		y = Math.max(y, GeometryConstants.DEGREES_TO_METERS_MIN_LAT);
		double xValue = x * GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH
				/ GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH;
		double yValue = Math.log(Math.tan(
				(GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT + y) * Math.PI
						/ (2 * GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH)))
				/ (Math.PI / GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH);
		yValue = yValue * GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH
				/ GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH;
		return new Point(xValue, yValue);
	}

	/**
	 * Convert a multi point in degrees to a multi point in meters
	 * 
	 * @param multiPoint
	 *            multi point in degrees
	 * @return multi point in meters
	 * @since 2.2.0
	 */
	public static MultiPoint degreesToMeters(MultiPoint multiPoint) {
		MultiPoint meters = new MultiPoint(multiPoint.hasZ(),
				multiPoint.hasM());
		for (Point point : multiPoint.getPoints()) {
			meters.addPoint(degreesToMeters(point));
		}
		return meters;
	}

	/**
	 * Convert a line string in degrees to a line string in meters
	 * 
	 * @param lineString
	 *            line string in degrees
	 * @return line string in meters
	 * @since 2.2.0
	 */
	public static LineString degreesToMeters(LineString lineString) {
		LineString meters = new LineString(lineString.hasZ(),
				lineString.hasM());
		for (Point point : lineString.getPoints()) {
			meters.addPoint(degreesToMeters(point));
		}
		return meters;
	}

	/**
	 * Convert a line in degrees to a line in meters
	 * 
	 * @param line
	 *            line in degrees
	 * @return line in meters
	 * @since 2.2.0
	 */
	public static Line degreesToMeters(Line line) {
		Line meters = new Line(line.hasZ(), line.hasM());
		for (Point point : line.getPoints()) {
			meters.addPoint(degreesToMeters(point));
		}
		return meters;
	}

	/**
	 * Convert a multi line string in degrees to a multi line string in meters
	 * 
	 * @param multiLineString
	 *            multi line string in degrees
	 * @return multi line string in meters
	 * @since 2.2.0
	 */
	public static MultiLineString degreesToMeters(
			MultiLineString multiLineString) {
		MultiLineString meters = new MultiLineString(multiLineString.hasZ(),
				multiLineString.hasM());
		for (LineString lineString : multiLineString.getLineStrings()) {
			meters.addLineString(degreesToMeters(lineString));
		}
		return meters;
	}

	/**
	 * Convert a polygon in degrees to a polygon in meters
	 * 
	 * @param polygon
	 *            polygon in degrees
	 * @return polygon in meters
	 * @since 2.2.0
	 */
	public static Polygon degreesToMeters(Polygon polygon) {
		Polygon meters = new Polygon(polygon.hasZ(), polygon.hasM());
		for (LineString ring : polygon.getRings()) {
			meters.addRing(degreesToMeters(ring));
		}
		return meters;
	}

	/**
	 * Convert a multi polygon in degrees to a multi polygon in meters
	 * 
	 * @param multiPolygon
	 *            multi polygon in degrees
	 * @return multi polygon in meters
	 * @since 2.2.0
	 */
	public static MultiPolygon degreesToMeters(MultiPolygon multiPolygon) {
		MultiPolygon meters = new MultiPolygon(multiPolygon.hasZ(),
				multiPolygon.hasM());
		for (Polygon polygon : multiPolygon.getPolygons()) {
			meters.addPolygon(degreesToMeters(polygon));
		}
		return meters;
	}

	/**
	 * Convert a circular string in degrees to a circular string in meters
	 * 
	 * @param circularString
	 *            circular string in degrees
	 * @return circular string in meters
	 * @since 2.2.0
	 */
	public static CircularString degreesToMeters(
			CircularString circularString) {
		CircularString meters = new CircularString(circularString.hasZ(),
				circularString.hasM());
		for (Point point : circularString.getPoints()) {
			meters.addPoint(degreesToMeters(point));
		}
		return meters;
	}

	/**
	 * Convert a compound curve in degrees to a compound curve in meters
	 * 
	 * @param compoundCurve
	 *            compound curve in degrees
	 * @return compound curve in meters
	 * @since 2.2.0
	 */
	public static CompoundCurve degreesToMeters(CompoundCurve compoundCurve) {
		CompoundCurve meters = new CompoundCurve(compoundCurve.hasZ(),
				compoundCurve.hasM());
		for (LineString lineString : compoundCurve.getLineStrings()) {
			meters.addLineString(degreesToMeters(lineString));
		}
		return meters;
	}

	/**
	 * Convert a curve polygon in degrees to a curve polygon in meters
	 * 
	 * @param curvePolygon
	 *            curve polygon in degrees
	 * @return curve polygon in meters
	 * @since 2.2.0
	 */
	public static CurvePolygon<Curve> degreesToMeters(
			CurvePolygon<Curve> curvePolygon) {
		CurvePolygon<Curve> meters = new CurvePolygon<>(curvePolygon.hasZ(),
				curvePolygon.hasM());
		for (Curve ring : curvePolygon.getRings()) {
			meters.addRing((Curve) degreesToMeters(ring));
		}
		return meters;
	}

	/**
	 * Convert a polyhedral surface in degrees to a polyhedral surface in meters
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface in degrees
	 * @return polyhedral surface in meters
	 * @since 2.2.0
	 */
	public static PolyhedralSurface degreesToMeters(
			PolyhedralSurface polyhedralSurface) {
		PolyhedralSurface meters = new PolyhedralSurface(
				polyhedralSurface.hasZ(), polyhedralSurface.hasM());
		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			meters.addPolygon(degreesToMeters(polygon));
		}
		return meters;
	}

	/**
	 * Convert a TIN in degrees to a TIN in meters
	 * 
	 * @param tin
	 *            TIN in degrees
	 * @return TIN in meters
	 * @since 2.2.0
	 */
	public static TIN degreesToMeters(TIN tin) {
		TIN meters = new TIN(tin.hasZ(), tin.hasM());
		for (Polygon polygon : tin.getPolygons()) {
			meters.addPolygon(degreesToMeters(polygon));
		}
		return meters;
	}

	/**
	 * Convert a triangle in degrees to a triangle in meters
	 * 
	 * @param triangle
	 *            triangle in degrees
	 * @return triangle in meters
	 * @since 2.2.0
	 */
	public static Triangle degreesToMeters(Triangle triangle) {
		Triangle meters = new Triangle(triangle.hasZ(), triangle.hasM());
		for (LineString ring : triangle.getRings()) {
			meters.addRing(degreesToMeters(ring));
		}
		return meters;
	}

	/**
	 * Convert a geometry in meters to a geometry in degrees
	 * 
	 * @param geometry
	 *            geometry in meters
	 * @return geometry in degrees
	 * @since 2.2.0
	 */
	public static Geometry metersToDegrees(Geometry geometry) {

		Geometry degrees = null;

		switch (geometry.getGeometryType()) {
		case POINT:
			degrees = metersToDegrees((Point) geometry);
			break;
		case LINESTRING:
			degrees = metersToDegrees((LineString) geometry);
			break;
		case POLYGON:
			degrees = metersToDegrees((Polygon) geometry);
			break;
		case MULTIPOINT:
			degrees = metersToDegrees((MultiPoint) geometry);
			break;
		case MULTILINESTRING:
			degrees = metersToDegrees((MultiLineString) geometry);
			break;
		case MULTIPOLYGON:
			degrees = metersToDegrees((MultiPolygon) geometry);
			break;
		case CIRCULARSTRING:
			degrees = metersToDegrees((CircularString) geometry);
			break;
		case COMPOUNDCURVE:
			degrees = metersToDegrees((CompoundCurve) geometry);
			break;
		case CURVEPOLYGON:
			@SuppressWarnings("unchecked")
			CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
			degrees = metersToDegrees(curvePolygon);
			break;
		case POLYHEDRALSURFACE:
			degrees = metersToDegrees((PolyhedralSurface) geometry);
			break;
		case TIN:
			degrees = metersToDegrees((TIN) geometry);
			break;
		case TRIANGLE:
			degrees = metersToDegrees((Triangle) geometry);
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			GeometryCollection<Geometry> degreesCollection = new GeometryCollection<>();
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			for (Geometry subGeometry : geomCollection.getGeometries()) {
				degreesCollection.addGeometry(metersToDegrees(subGeometry));
			}
			degrees = degreesCollection;
			break;
		default:
			break;

		}

		return degrees;
	}

	/**
	 * Convert a point in meters to a point in degrees
	 * 
	 * @param point
	 *            point in meters
	 * @return point in degrees
	 * @since 2.1.0
	 */
	public static Point metersToDegrees(Point point) {
		Point value = metersToDegrees(point.getX(), point.getY());
		value.setZ(point.getZ());
		value.setM(point.getM());
		return value;
	}

	/**
	 * Convert a coordinate in meters to a point in degrees
	 * 
	 * @param x
	 *            x value in meters
	 * @param y
	 *            y value in meters
	 * @return point in degrees
	 * @since 2.1.0
	 */
	public static Point metersToDegrees(double x, double y) {
		double xValue = x * GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH
				/ GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH;
		double yValue = y * GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH
				/ GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH;
		yValue = Math.atan(Math.exp(yValue
				* (Math.PI / GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH)))
				/ Math.PI * (2 * GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH)
				- GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT;
		return new Point(xValue, yValue);
	}

	/**
	 * Convert a multi point in meters to a multi point in degrees
	 * 
	 * @param multiPoint
	 *            multi point in meters
	 * @return multi point in degrees
	 * @since 2.2.0
	 */
	public static MultiPoint metersToDegrees(MultiPoint multiPoint) {
		MultiPoint degrees = new MultiPoint(multiPoint.hasZ(),
				multiPoint.hasM());
		for (Point point : multiPoint.getPoints()) {
			degrees.addPoint(metersToDegrees(point));
		}
		return degrees;
	}

	/**
	 * Convert a line string in meters to a line string in degrees
	 * 
	 * @param lineString
	 *            line string in meters
	 * @return line string in degrees
	 * @since 2.2.0
	 */
	public static LineString metersToDegrees(LineString lineString) {
		LineString degrees = new LineString(lineString.hasZ(),
				lineString.hasM());
		for (Point point : lineString.getPoints()) {
			degrees.addPoint(metersToDegrees(point));
		}
		return degrees;
	}

	/**
	 * Convert a line in meters to a line in degrees
	 * 
	 * @param line
	 *            line in meters
	 * @return line in degrees
	 * @since 2.2.0
	 */
	public static Line metersToDegrees(Line line) {
		Line degrees = new Line(line.hasZ(), line.hasM());
		for (Point point : line.getPoints()) {
			degrees.addPoint(metersToDegrees(point));
		}
		return degrees;
	}

	/**
	 * Convert a multi line string in meters to a multi line string in degrees
	 * 
	 * @param multiLineString
	 *            multi line string in meters
	 * @return multi line string in degrees
	 * @since 2.2.0
	 */
	public static MultiLineString metersToDegrees(
			MultiLineString multiLineString) {
		MultiLineString degrees = new MultiLineString(multiLineString.hasZ(),
				multiLineString.hasM());
		for (LineString lineString : multiLineString.getLineStrings()) {
			degrees.addLineString(metersToDegrees(lineString));
		}
		return degrees;
	}

	/**
	 * Convert a polygon in meters to a polygon in degrees
	 * 
	 * @param polygon
	 *            polygon in meters
	 * @return polygon in degrees
	 * @since 2.2.0
	 */
	public static Polygon metersToDegrees(Polygon polygon) {
		Polygon degrees = new Polygon(polygon.hasZ(), polygon.hasM());
		for (LineString ring : polygon.getRings()) {
			degrees.addRing(metersToDegrees(ring));
		}
		return degrees;
	}

	/**
	 * Convert a multi polygon in meters to a multi polygon in degrees
	 * 
	 * @param multiPolygon
	 *            multi polygon in meters
	 * @return multi polygon in degrees
	 * @since 2.2.0
	 */
	public static MultiPolygon metersToDegrees(MultiPolygon multiPolygon) {
		MultiPolygon degrees = new MultiPolygon(multiPolygon.hasZ(),
				multiPolygon.hasM());
		for (Polygon polygon : multiPolygon.getPolygons()) {
			degrees.addPolygon(metersToDegrees(polygon));
		}
		return degrees;
	}

	/**
	 * Convert a circular string in meters to a circular string in degrees
	 * 
	 * @param circularString
	 *            circular string in meters
	 * @return circular string in degrees
	 * @since 2.2.0
	 */
	public static CircularString metersToDegrees(
			CircularString circularString) {
		CircularString degrees = new CircularString(circularString.hasZ(),
				circularString.hasM());
		for (Point point : circularString.getPoints()) {
			degrees.addPoint(metersToDegrees(point));
		}
		return degrees;
	}

	/**
	 * Convert a compound curve in meters to a compound curve in degrees
	 * 
	 * @param compoundCurve
	 *            compound curve in meters
	 * @return compound curve in degrees
	 * @since 2.2.0
	 */
	public static CompoundCurve metersToDegrees(CompoundCurve compoundCurve) {
		CompoundCurve degrees = new CompoundCurve(compoundCurve.hasZ(),
				compoundCurve.hasM());
		for (LineString lineString : compoundCurve.getLineStrings()) {
			degrees.addLineString(metersToDegrees(lineString));
		}
		return degrees;
	}

	/**
	 * Convert a curve polygon in meters to a curve polygon in degrees
	 * 
	 * @param curvePolygon
	 *            curve polygon in meters
	 * @return curve polygon in degrees
	 * @since 2.2.0
	 */
	public static CurvePolygon<Curve> metersToDegrees(
			CurvePolygon<Curve> curvePolygon) {
		CurvePolygon<Curve> degrees = new CurvePolygon<>(curvePolygon.hasZ(),
				curvePolygon.hasM());
		for (Curve ring : curvePolygon.getRings()) {
			degrees.addRing((Curve) metersToDegrees(ring));
		}
		return degrees;
	}

	/**
	 * Convert a polyhedral surface in meters to a polyhedral surface in degrees
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface in meters
	 * @return polyhedral surface in degrees
	 * @since 2.2.0
	 */
	public static PolyhedralSurface metersToDegrees(
			PolyhedralSurface polyhedralSurface) {
		PolyhedralSurface degrees = new PolyhedralSurface(
				polyhedralSurface.hasZ(), polyhedralSurface.hasM());
		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			degrees.addPolygon(metersToDegrees(polygon));
		}
		return degrees;
	}

	/**
	 * Convert a TIN in meters to a TIN in degrees
	 * 
	 * @param tin
	 *            TIN in meters
	 * @return TIN in degrees
	 * @since 2.2.0
	 */
	public static TIN metersToDegrees(TIN tin) {
		TIN degrees = new TIN(tin.hasZ(), tin.hasM());
		for (Polygon polygon : tin.getPolygons()) {
			degrees.addPolygon(metersToDegrees(polygon));
		}
		return degrees;
	}

	/**
	 * Convert a triangle in meters to a triangle in degrees
	 * 
	 * @param triangle
	 *            triangle in meters
	 * @return triangle in degrees
	 * @since 2.2.0
	 */
	public static Triangle metersToDegrees(Triangle triangle) {
		Triangle degrees = new Triangle(triangle.hasZ(), triangle.hasM());
		for (LineString ring : triangle.getRings()) {
			degrees.addRing(metersToDegrees(ring));
		}
		return degrees;
	}

	/**
	 * Get a WGS84 bounded geometry envelope
	 * 
	 * @return geometry envelope
	 * @since 2.2.0
	 */
	public static GeometryEnvelope wgs84Envelope() {
		return new GeometryEnvelope(
				-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				-GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT);
	}

	/**
	 * Get a WGS84 bounded geometry envelope used for projection transformations
	 * (degrees to meters)
	 * 
	 * @return geometry envelope
	 * @since 2.2.0
	 */
	public static GeometryEnvelope wgs84TransformableEnvelope() {
		return new GeometryEnvelope(
				-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				GeometryConstants.DEGREES_TO_METERS_MIN_LAT,
				GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT);
	}

	/**
	 * Get a Web Mercator bounded geometry envelope
	 * 
	 * @return geometry envelope
	 * @since 2.2.0
	 */
	public static GeometryEnvelope webMercatorEnvelope() {
		return new GeometryEnvelope(
				-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH);
	}

	/**
	 * Get a WGS84 geometry envelope with Web Mercator bounds
	 * 
	 * @return geometry envelope
	 * @since 2.2.0
	 */
	public static GeometryEnvelope wgs84EnvelopeWithWebMercator() {
		return new GeometryEnvelope(
				-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE);
	}

	/**
	 * Crop the geometry in meters by web mercator world bounds. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param geometry
	 *            geometry in meters
	 * @return cropped geometry in meters or null
	 * @since 2.2.0
	 */
	public static Geometry cropWebMercator(Geometry geometry) {
		return crop(geometry, webMercatorEnvelope());
	}

	/**
	 * Crop the geometry in meters by the envelope bounds in meters. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param geometry
	 *            geometry in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped geometry in meters or null
	 * @since 2.2.0
	 */
	public static Geometry crop(Geometry geometry, GeometryEnvelope envelope) {

		Geometry crop = null;

		if (contains(envelope, geometry.getEnvelope())) {
			crop = geometry;
		} else {

			switch (geometry.getGeometryType()) {
			case POINT:
				crop = crop((Point) geometry, envelope);
				break;
			case LINESTRING:
				crop = crop((LineString) geometry, envelope);
				break;
			case POLYGON:
				crop = crop((Polygon) geometry, envelope);
				break;
			case MULTIPOINT:
				crop = crop((MultiPoint) geometry, envelope);
				break;
			case MULTILINESTRING:
				crop = crop((MultiLineString) geometry, envelope);
				break;
			case MULTIPOLYGON:
				crop = crop((MultiPolygon) geometry, envelope);
				break;
			case CIRCULARSTRING:
				crop = crop((CircularString) geometry, envelope);
				break;
			case COMPOUNDCURVE:
				crop = crop((CompoundCurve) geometry, envelope);
				break;
			case CURVEPOLYGON:
				@SuppressWarnings("unchecked")
				CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
				crop = crop(curvePolygon, envelope);
				break;
			case POLYHEDRALSURFACE:
				crop = crop((PolyhedralSurface) geometry, envelope);
				break;
			case TIN:
				crop = crop((TIN) geometry, envelope);
				break;
			case TRIANGLE:
				crop = crop((Triangle) geometry, envelope);
				break;
			case GEOMETRYCOLLECTION:
			case MULTICURVE:
			case MULTISURFACE:
				GeometryCollection<Geometry> cropCollection = new GeometryCollection<>();
				@SuppressWarnings("unchecked")
				GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
				for (Geometry subGeometry : geomCollection.getGeometries()) {
					cropCollection.addGeometry(crop(subGeometry, envelope));
				}
				crop = cropCollection;
				break;
			default:
				break;

			}

		}

		return crop;
	}

	/**
	 * Crop the point by the envelope bounds.
	 * 
	 * @param point
	 *            point
	 * @param envelope
	 *            envelope
	 * @return cropped point or null
	 * @since 2.2.0
	 */
	public static Point crop(Point point, GeometryEnvelope envelope) {
		Point crop = null;
		if (contains(envelope, point)) {
			crop = new Point(point);
		}
		return crop;
	}

	/**
	 * Crop the list of consecutive points in meters by the envelope bounds in
	 * meters. Cropping removes points outside the envelope and creates new
	 * points on the line intersections with the envelope.
	 * 
	 * @param points
	 *            consecutive points
	 * @param envelope
	 *            envelope in meters
	 * @return cropped points in meters or null
	 * @since 2.2.0
	 */
	public static List<Point> crop(List<Point> points,
			GeometryEnvelope envelope) {

		List<Point> crop = new ArrayList<>();

		Line left = envelope.getLeft();
		Line bottom = envelope.getBottom();
		Line right = envelope.getRight();
		Line top = envelope.getTop();

		Point previousPoint = null;
		boolean previousContains = false;
		for (Point point : points) {
			boolean contains = contains(envelope, point);

			if (previousPoint != null && (!contains || !previousContains)) {

				Line line = new Line(previousPoint, point);
				double bearing = bearing(metersToDegrees(line));

				boolean westBearing = isWestBearing(bearing);
				boolean eastBearing = isEastBearing(bearing);
				boolean southBearing = isSouthBearing(bearing);
				boolean northBearing = isNorthBearing(bearing);

				Line vertLine = null;
				if (point.getX() > envelope.getMaxX()) {
					if (eastBearing) {
						vertLine = right;
					}
				} else if (point.getX() < envelope.getMinX()) {
					if (westBearing) {
						vertLine = left;
					}
				} else if (eastBearing) {
					vertLine = left;
				} else if (westBearing) {
					vertLine = right;
				}

				Line horizLine = null;
				if (point.getY() > envelope.getMaxY()) {
					if (northBearing) {
						horizLine = top;
					}
				} else if (point.getY() < envelope.getMinY()) {
					if (southBearing) {
						horizLine = bottom;
					}
				} else if (northBearing) {
					horizLine = bottom;
				} else if (southBearing) {
					horizLine = top;
				}

				Point vertIntersection = null;
				if (vertLine != null) {
					vertIntersection = intersection(line, vertLine);
					if (vertIntersection != null
							&& !contains(envelope, vertIntersection)) {
						vertIntersection = null;
					}
				}

				Point horizIntersection = null;
				if (horizLine != null) {
					horizIntersection = intersection(line, horizLine);
					if (horizIntersection != null
							&& !contains(envelope, horizIntersection)) {
						horizIntersection = null;
					}
				}

				Point intersection1 = null;
				Point intersection2 = null;
				if (vertIntersection != null && horizIntersection != null) {
					double vertDistance = distance(previousPoint,
							vertIntersection);
					double horizDistance = distance(previousPoint,
							horizIntersection);
					if (vertDistance <= horizDistance) {
						intersection1 = vertIntersection;
						intersection2 = horizIntersection;
					} else {
						intersection1 = horizIntersection;
						intersection2 = vertIntersection;
					}
				} else if (vertIntersection != null) {
					intersection1 = vertIntersection;
				} else {
					intersection1 = horizIntersection;
				}

				if (intersection1 != null && !isEqual(intersection1, point)
						&& !isEqual(intersection1, previousPoint)) {

					crop.add(intersection1);

					if (!contains && !previousContains && intersection2 != null
							&& !isEqual(intersection2, intersection1)) {
						crop.add(intersection2);
					}
				}

			}

			if (contains) {
				crop.add(point);
			}

			previousPoint = point;
			previousContains = contains;
		}

		if (crop.isEmpty()) {
			crop = null;
		} else if (crop.size() > 1) {

			if (points.get(0).equals(points.get(points.size() - 1))
					&& !crop.get(0).equals(crop.get(crop.size() - 1))) {
				crop.add(new Point(crop.get(0)));
			}

			if (crop.size() > 2) {

				List<Point> simplified = new ArrayList<>();
				simplified.add(crop.get(0));
				for (int i = 1; i < crop.size() - 1; i++) {
					Point previous = simplified.get(simplified.size() - 1);
					Point point = crop.get(i);
					Point next = crop.get(i + 1);
					if (!pointOnPath(point, previous, next)) {
						simplified.add(point);
					}
				}
				simplified.add(crop.get(crop.size() - 1));
				crop = simplified;

			}

		}

		return crop;
	}

	/**
	 * Crop the multi point by the envelope bounds.
	 * 
	 * @param multiPoint
	 *            multi point
	 * @param envelope
	 *            envelope
	 * @return cropped multi point or null
	 * @since 2.2.0
	 */
	public static MultiPoint crop(MultiPoint multiPoint,
			GeometryEnvelope envelope) {
		MultiPoint crop = null;
		List<Point> cropPoints = new ArrayList<>();
		for (Point point : multiPoint.getPoints()) {
			Point cropPoint = crop(point, envelope);
			if (cropPoint != null) {
				cropPoints.add(cropPoint);
			}
		}
		if (!cropPoints.isEmpty()) {
			crop = new MultiPoint(multiPoint.hasZ(), multiPoint.hasM());
			crop.setPoints(cropPoints);
		}
		return crop;
	}

	/**
	 * Crop the line string in meters by the envelope bounds in meters. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param lineString
	 *            line string in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped line string in meters or null
	 * @since 2.2.0
	 */
	public static LineString crop(LineString lineString,
			GeometryEnvelope envelope) {
		LineString crop = null;
		List<Point> cropPoints = crop(lineString.getPoints(), envelope);
		if (cropPoints != null) {
			crop = new LineString(lineString.hasZ(), lineString.hasM());
			crop.setPoints(cropPoints);
		}
		return crop;
	}

	/**
	 * Crop the line in meters by the envelope bounds in meters. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param line
	 *            line in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped line in meters or null
	 * @since 2.2.0
	 */
	public static Line crop(Line line, GeometryEnvelope envelope) {
		Line crop = null;
		List<Point> cropPoints = crop(line.getPoints(), envelope);
		if (cropPoints != null) {
			crop = new Line(line.hasZ(), line.hasM());
			crop.setPoints(cropPoints);
		}
		return crop;
	}

	/**
	 * Crop the multi line string in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param multiLineString
	 *            multi line string in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped multi line string in meters or null
	 * @since 2.2.0
	 */
	public static MultiLineString crop(MultiLineString multiLineString,
			GeometryEnvelope envelope) {
		MultiLineString crop = null;
		List<LineString> cropLineStrings = new ArrayList<>();
		for (LineString lineString : multiLineString.getLineStrings()) {
			LineString cropLineString = crop(lineString, envelope);
			if (cropLineString != null) {
				cropLineStrings.add(cropLineString);
			}
		}
		if (!cropLineStrings.isEmpty()) {
			crop = new MultiLineString(multiLineString.hasZ(),
					multiLineString.hasM());
			crop.setLineStrings(cropLineStrings);
		}
		return crop;
	}

	/**
	 * Crop the polygon in meters by the envelope bounds in meters. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param polygon
	 *            polygon in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped polygon in meters or null
	 * @since 2.2.0
	 */
	public static Polygon crop(Polygon polygon, GeometryEnvelope envelope) {
		Polygon crop = null;
		List<LineString> cropRings = new ArrayList<>();
		for (LineString ring : polygon.getRings()) {
			List<Point> points = ring.getPoints();
			if (!ring.isClosed()) {
				points.add(new Point(points.get(0)));
			}
			List<Point> cropPoints = crop(points, envelope);
			if (cropPoints != null) {
				LineString cropRing = new LineString(ring.hasZ(), ring.hasM());
				cropRing.setPoints(cropPoints);
				cropRings.add(cropRing);
			}
		}
		if (!cropRings.isEmpty()) {
			crop = new Polygon(polygon.hasZ(), polygon.hasM());
			crop.setRings(cropRings);
		}
		return crop;
	}

	/**
	 * Crop the multi polygon in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param multiPolygon
	 *            multi polygon in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped multi polygon in meters or null
	 * @since 2.2.0
	 */
	public static MultiPolygon crop(MultiPolygon multiPolygon,
			GeometryEnvelope envelope) {
		MultiPolygon crop = null;
		List<Polygon> cropPolygons = new ArrayList<>();
		for (Polygon polygon : multiPolygon.getPolygons()) {
			Polygon cropPolygon = crop(polygon, envelope);
			if (cropPolygon != null) {
				cropPolygons.add(cropPolygon);
			}
		}
		if (!cropPolygons.isEmpty()) {
			crop = new MultiPolygon(multiPolygon.hasZ(), multiPolygon.hasM());
			crop.setPolygons(cropPolygons);
		}
		return crop;
	}

	/**
	 * Crop the circular string in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param circularString
	 *            circular string in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped circular string in meters or null
	 * @since 2.2.0
	 */
	public static CircularString crop(CircularString circularString,
			GeometryEnvelope envelope) {
		CircularString crop = null;
		List<Point> cropPoints = crop(circularString.getPoints(), envelope);
		if (cropPoints != null) {
			crop = new CircularString(circularString.hasZ(),
					circularString.hasM());
			crop.setPoints(cropPoints);
		}
		return crop;
	}

	/**
	 * Crop the compound curve in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param compoundCurve
	 *            compound curve in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped compound curve in meters or null
	 * @since 2.2.0
	 */
	public static CompoundCurve crop(CompoundCurve compoundCurve,
			GeometryEnvelope envelope) {
		CompoundCurve crop = null;
		List<LineString> cropLineStrings = new ArrayList<>();
		for (LineString lineString : compoundCurve.getLineStrings()) {
			LineString cropLineString = crop(lineString, envelope);
			if (cropLineString != null) {
				cropLineStrings.add(cropLineString);
			}
		}
		if (!cropLineStrings.isEmpty()) {
			crop = new CompoundCurve(compoundCurve.hasZ(),
					compoundCurve.hasM());
			crop.setLineStrings(cropLineStrings);
		}
		return crop;
	}

	/**
	 * Crop the curve polygon in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param curvePolygon
	 *            curve polygon in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped curve polygon in meters or null
	 * @since 2.2.0
	 */
	public static CurvePolygon<Curve> crop(CurvePolygon<Curve> curvePolygon,
			GeometryEnvelope envelope) {
		CurvePolygon<Curve> crop = null;
		List<Curve> cropRings = new ArrayList<>();
		for (Curve ring : curvePolygon.getRings()) {
			Geometry cropRing = crop(ring, envelope);
			if (cropRing != null) {
				cropRings.add((Curve) cropRing);
			}
		}
		if (!cropRings.isEmpty()) {
			crop = new CurvePolygon<>(curvePolygon.hasZ(), curvePolygon.hasM());
			crop.setRings(cropRings);
		}
		return crop;
	}

	/**
	 * Crop the polyhedral surface in meters by the envelope bounds in meters.
	 * Cropping removes points outside the envelope and creates new points on
	 * the line intersections with the envelope.
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped polyhedral surface in meters or null
	 * @since 2.2.0
	 */
	public static PolyhedralSurface crop(PolyhedralSurface polyhedralSurface,
			GeometryEnvelope envelope) {
		PolyhedralSurface crop = null;
		List<Polygon> cropPolygons = new ArrayList<>();
		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			Polygon cropPolygon = crop(polygon, envelope);
			if (cropPolygon != null) {
				cropPolygons.add(cropPolygon);
			}
		}
		if (!cropPolygons.isEmpty()) {
			crop = new PolyhedralSurface(polyhedralSurface.hasZ(),
					polyhedralSurface.hasM());
			crop.setPolygons(cropPolygons);
		}
		return crop;
	}

	/**
	 * Crop the TIN in meters by the envelope bounds in meters. Cropping removes
	 * points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param tin
	 *            TIN in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped TIN in meters or null
	 * @since 2.2.0
	 */
	public static TIN crop(TIN tin, GeometryEnvelope envelope) {
		TIN crop = null;
		List<Polygon> cropPolygons = new ArrayList<>();
		for (Polygon polygon : tin.getPolygons()) {
			Polygon cropPolygon = crop(polygon, envelope);
			if (cropPolygon != null) {
				cropPolygons.add(cropPolygon);
			}
		}
		if (!cropPolygons.isEmpty()) {
			crop = new TIN(tin.hasZ(), tin.hasM());
			crop.setPolygons(cropPolygons);
		}
		return crop;
	}

	/**
	 * Crop the triangle in meters by the envelope bounds in meters. Cropping
	 * removes points outside the envelope and creates new points on the line
	 * intersections with the envelope.
	 * 
	 * @param triangle
	 *            triangle in meters
	 * @param envelope
	 *            envelope in meters
	 * @return cropped triangle in meters or null
	 * @since 2.2.0
	 */
	public static Triangle crop(Triangle triangle, GeometryEnvelope envelope) {
		Triangle crop = null;
		List<LineString> cropRings = new ArrayList<>();
		for (LineString ring : triangle.getRings()) {
			List<Point> points = ring.getPoints();
			if (!ring.isClosed()) {
				points.add(new Point(points.get(0)));
			}
			List<Point> cropPoints = crop(points, envelope);
			if (cropPoints != null) {
				LineString cropRing = new LineString(ring.hasZ(), ring.hasM());
				cropRing.setPoints(cropPoints);
				cropRings.add(cropRing);
			}
		}
		if (!cropRings.isEmpty()) {
			crop = new Triangle(triangle.hasZ(), triangle.hasM());
			crop.setRings(cropRings);
		}
		return crop;
	}

	/**
	 * Determine if the points are equal within the default tolerance of
	 * {@link GeometryConstants#DEFAULT_EQUAL_EPSILON}. For exact equality, use
	 * {@link Point#equals(Object)}.
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @return true if equal
	 * @since 2.2.0
	 */
	public static boolean isEqual(Point point1, Point point2) {
		return isEqual(point1, point2, GeometryConstants.DEFAULT_EQUAL_EPSILON);
	}

	/**
	 * Determine if the points are equal within the tolerance. For exact
	 * equality, use {@link Point#equals(Object)}.
	 * 
	 * @param point1
	 *            point 1
	 * @param point2
	 *            point 2
	 * @param epsilon
	 *            epsilon equality tolerance
	 * @return true if equal
	 * @since 2.2.0
	 */
	public static boolean isEqual(Point point1, Point point2, double epsilon) {
		boolean equal = Math.abs(point1.getX() - point2.getX()) <= epsilon
				&& Math.abs(point1.getY() - point2.getY()) <= epsilon
				&& point1.hasZ() == point2.hasZ()
				&& point1.hasM() == point2.hasM();
		if (equal) {
			if (point1.hasZ()) {
				equal = Math.abs(point1.getZ() - point2.getZ()) <= epsilon;
			}
			if (equal && point1.hasM()) {
				equal = Math.abs(point1.getM() - point2.getM()) <= epsilon;
			}
		}
		return equal;
	}

	/**
	 * Determine if the envelope contains the point within the default tolerance
	 * of {@link GeometryConstants#DEFAULT_EQUAL_EPSILON}. For exact equality,
	 * use {@link GeometryEnvelope#contains(Point)}.
	 * 
	 * @param envelope
	 *            envelope
	 * @param point
	 *            point
	 * @return true if contains
	 * @since 2.2.0
	 */
	public static boolean contains(GeometryEnvelope envelope, Point point) {
		return envelope.contains(point,
				GeometryConstants.DEFAULT_EQUAL_EPSILON);
	}

	/**
	 * Determine if the first envelope contains the second within the default
	 * tolerance of {@link GeometryConstants#DEFAULT_EQUAL_EPSILON}. For exact
	 * equality, use {@link GeometryEnvelope#contains(GeometryEnvelope)}.
	 * 
	 * @param envelope1
	 *            envelope 1
	 * @param envelope2
	 *            envelope 2
	 * @return true if contains
	 * @since 2.2.0
	 */
	public static boolean contains(GeometryEnvelope envelope1,
			GeometryEnvelope envelope2) {
		return envelope1.contains(envelope2,
				GeometryConstants.DEFAULT_EQUAL_EPSILON);
	}

	/**
	 * Bound all points in the geometry to be within WGS84 limits.
	 * 
	 * To perform a geometry crop using line intersections, see
	 * {@link #degreesToMeters(Geometry)} and
	 * {@link #crop(Geometry, GeometryEnvelope)}.
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void boundWGS84(Geometry geometry) {
		bound(geometry, wgs84Envelope());
	}

	/**
	 * Bound all points in the geometry to be within WGS84 projection
	 * transformable (degrees to meters) limits.
	 * 
	 * To perform a geometry crop using line intersections, see
	 * {@link #degreesToMeters(Geometry)} and
	 * {@link #crop(Geometry, GeometryEnvelope)}.
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void boundWGS84Transformable(Geometry geometry) {
		bound(geometry, wgs84TransformableEnvelope());
	}

	/**
	 * Bound all points in the geometry to be within Web Mercator limits.
	 * 
	 * To perform a geometry crop using line intersections, see
	 * {@link #cropWebMercator(Geometry)}.
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void boundWebMercator(Geometry geometry) {
		bound(geometry, webMercatorEnvelope());
	}

	/**
	 * Bound all points in the WGS84 geometry to be within degree Web Mercator
	 * limits.
	 * 
	 * To perform a geometry crop using line intersections, see
	 * {@link #degreesToMeters(Geometry)} and
	 * {@link #cropWebMercator(Geometry)}.
	 * 
	 * @param geometry
	 *            geometry
	 * @since 2.2.0
	 */
	public static void boundWGS84WithWebMercator(Geometry geometry) {
		bound(geometry, wgs84EnvelopeWithWebMercator());
	}

	/**
	 * Bound all points in the geometry to be within the geometry envelope.
	 * Point x and y values are bounded by the min and max envelope values.
	 * 
	 * To perform a geometry crop using line intersections, see
	 * {@link #crop(Geometry, GeometryEnvelope)} (requires geometry in meters).
	 * 
	 * @param geometry
	 *            geometry
	 * @param envelope
	 *            geometry envelope
	 * @since 2.2.0
	 */
	public static void bound(Geometry geometry, GeometryEnvelope envelope) {

		GeometryType geometryType = geometry.getGeometryType();
		switch (geometryType) {
		case POINT:
			bound((Point) geometry, envelope);
			break;
		case LINESTRING:
			bound((LineString) geometry, envelope);
			break;
		case POLYGON:
			bound((Polygon) geometry, envelope);
			break;
		case MULTIPOINT:
			bound((MultiPoint) geometry, envelope);
			break;
		case MULTILINESTRING:
			bound((MultiLineString) geometry, envelope);
			break;
		case MULTIPOLYGON:
			bound((MultiPolygon) geometry, envelope);
			break;
		case CIRCULARSTRING:
			bound((CircularString) geometry, envelope);
			break;
		case COMPOUNDCURVE:
			bound((CompoundCurve) geometry, envelope);
			break;
		case CURVEPOLYGON:
			@SuppressWarnings("unchecked")
			CurvePolygon<Curve> curvePolygon = (CurvePolygon<Curve>) geometry;
			bound(curvePolygon, envelope);
			break;
		case POLYHEDRALSURFACE:
			bound((PolyhedralSurface) geometry, envelope);
			break;
		case TIN:
			bound((TIN) geometry, envelope);
			break;
		case TRIANGLE:
			bound((Triangle) geometry, envelope);
			break;
		case GEOMETRYCOLLECTION:
		case MULTICURVE:
		case MULTISURFACE:
			@SuppressWarnings("unchecked")
			GeometryCollection<Geometry> geomCollection = (GeometryCollection<Geometry>) geometry;
			for (Geometry subGeometry : geomCollection.getGeometries()) {
				bound(subGeometry, envelope);
			}
			break;
		default:
			break;

		}

	}

	/**
	 * Bound the point by the geometry envelope
	 * 
	 * @param point
	 *            point
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(Point point, GeometryEnvelope envelope) {
		double x = point.getX();
		double y = point.getY();
		if (x < envelope.getMinX()) {
			point.setX(envelope.getMinX());
		} else if (x > envelope.getMaxX()) {
			point.setX(envelope.getMaxX());
		}
		if (y < envelope.getMinY()) {
			point.setY(envelope.getMinY());
		} else if (y > envelope.getMaxY()) {
			point.setY(envelope.getMaxY());
		}
	}

	/**
	 * Bound the multi point by the geometry envelope
	 * 
	 * @param multiPoint
	 *            multi point
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(MultiPoint multiPoint,
			GeometryEnvelope envelope) {
		for (Point point : multiPoint.getPoints()) {
			bound(point, envelope);
		}
	}

	/**
	 * Bound the line string by the geometry envelope
	 * 
	 * @param lineString
	 *            line string
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(LineString lineString,
			GeometryEnvelope envelope) {
		for (Point point : lineString.getPoints()) {
			bound(point, envelope);
		}
	}

	/**
	 * Bound the multi line string by the geometry envelope
	 * 
	 * @param multiLineString
	 *            multi line string
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(MultiLineString multiLineString,
			GeometryEnvelope envelope) {
		for (LineString lineString : multiLineString.getLineStrings()) {
			bound(lineString, envelope);
		}
	}

	/**
	 * Bound the polygon by the geometry envelope
	 * 
	 * @param polygon
	 *            polygon
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(Polygon polygon, GeometryEnvelope envelope) {
		for (LineString ring : polygon.getRings()) {
			bound(ring, envelope);
		}
	}

	/**
	 * Bound the multi polygon by the geometry envelope
	 * 
	 * @param multiPolygon
	 *            multi polygon
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(MultiPolygon multiPolygon,
			GeometryEnvelope envelope) {
		for (Polygon polygon : multiPolygon.getPolygons()) {
			bound(polygon, envelope);
		}
	}

	/**
	 * Bound the compound curve by the geometry envelope
	 * 
	 * @param compoundCurve
	 *            compound curve
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(CompoundCurve compoundCurve,
			GeometryEnvelope envelope) {
		for (LineString lineString : compoundCurve.getLineStrings()) {
			bound(lineString, envelope);
		}
	}

	/**
	 * Bound the curve polygon by the geometry envelope
	 * 
	 * @param curvePolygon
	 *            curve polygon
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(CurvePolygon<Curve> curvePolygon,
			GeometryEnvelope envelope) {
		for (Curve ring : curvePolygon.getRings()) {
			bound(ring, envelope);
		}
	}

	/**
	 * Bound the polyhedral surface by the geometry envelope
	 * 
	 * @param polyhedralSurface
	 *            polyhedral surface
	 * @param envelope
	 *            geometry envelope
	 */
	private static void bound(PolyhedralSurface polyhedralSurface,
			GeometryEnvelope envelope) {
		for (Polygon polygon : polyhedralSurface.getPolygons()) {
			bound(polygon, envelope);
		}
	}

	/**
	 * Determine if the geometries contain a Z value
	 * 
	 * @param geometries
	 *            list of geometries
	 * @param <T>
	 *            geometry type
	 * @return true if has z
	 */
	public static <T extends Geometry> boolean hasZ(List<T> geometries) {
		boolean hasZ = false;
		for (Geometry geometry : geometries) {
			if (geometry.hasZ()) {
				hasZ = true;
				break;
			}
		}
		return hasZ;
	}

	/**
	 * Determine if the geometries contain a M value
	 * 
	 * @param geometries
	 *            list of geometries
	 * @param <T>
	 *            geometry type
	 * @return true if has m
	 */
	public static <T extends Geometry> boolean hasM(List<T> geometries) {
		boolean hasM = false;
		for (Geometry geometry : geometries) {
			if (geometry.hasM()) {
				hasM = true;
				break;
			}
		}
		return hasM;
	}

	/**
	 * Get the parent type hierarchy of the provided geometry type starting with
	 * the immediate parent. If the argument is GEOMETRY, an empty list is
	 * returned, else the final type in the list will be GEOMETRY.
	 * 
	 * @param geometryType
	 *            geometry type
	 * @return list of increasing parent types
	 * @since 2.0.1
	 */
	public static List<GeometryType> parentHierarchy(
			GeometryType geometryType) {

		List<GeometryType> hierarchy = new ArrayList<>();

		GeometryType parentType = parentType(geometryType);
		while (parentType != null) {
			hierarchy.add(parentType);
			parentType = parentType(parentType);
		}

		return hierarchy;
	}

	/**
	 * Get the parent Geometry Type of the provided geometry type
	 * 
	 * @param geometryType
	 *            geometry type
	 * @return parent geometry type or null if argument is GEOMETRY (no parent
	 *         type)
	 * @since 2.0.1
	 */
	public static GeometryType parentType(GeometryType geometryType) {

		GeometryType parentType = null;

		switch (geometryType) {

		case GEOMETRY:
			break;
		case POINT:
			parentType = GeometryType.GEOMETRY;
			break;
		case LINESTRING:
			parentType = GeometryType.CURVE;
			break;
		case POLYGON:
			parentType = GeometryType.CURVEPOLYGON;
			break;
		case MULTIPOINT:
			parentType = GeometryType.GEOMETRYCOLLECTION;
			break;
		case MULTILINESTRING:
			parentType = GeometryType.MULTICURVE;
			break;
		case MULTIPOLYGON:
			parentType = GeometryType.MULTISURFACE;
			break;
		case GEOMETRYCOLLECTION:
			parentType = GeometryType.GEOMETRY;
			break;
		case CIRCULARSTRING:
			parentType = GeometryType.LINESTRING;
			break;
		case COMPOUNDCURVE:
			parentType = GeometryType.CURVE;
			break;
		case CURVEPOLYGON:
			parentType = GeometryType.SURFACE;
			break;
		case MULTICURVE:
			parentType = GeometryType.GEOMETRYCOLLECTION;
			break;
		case MULTISURFACE:
			parentType = GeometryType.GEOMETRYCOLLECTION;
			break;
		case CURVE:
			parentType = GeometryType.GEOMETRY;
			break;
		case SURFACE:
			parentType = GeometryType.GEOMETRY;
			break;
		case POLYHEDRALSURFACE:
			parentType = GeometryType.SURFACE;
			break;
		case TIN:
			parentType = GeometryType.POLYHEDRALSURFACE;
			break;
		case TRIANGLE:
			parentType = GeometryType.POLYGON;
			break;
		default:
			throw new SFException(
					"Geometry Type not supported: " + geometryType);
		}

		return parentType;
	}

	/**
	 * Get the child type hierarchy of the provided geometry type.
	 * 
	 * @param geometryType
	 *            geometry type
	 * @return child type hierarchy, null if no children
	 * @since 2.0.1
	 */
	public static Map<GeometryType, Map<GeometryType, ?>> childHierarchy(
			GeometryType geometryType) {

		Map<GeometryType, Map<GeometryType, ?>> hierarchy = null;

		List<GeometryType> childTypes = childTypes(geometryType);

		if (!childTypes.isEmpty()) {
			hierarchy = new HashMap<>();

			for (GeometryType childType : childTypes) {
				hierarchy.put(childType, childHierarchy(childType));
			}
		}

		return hierarchy;
	}

	/**
	 * Get the immediate child Geometry Types of the provided geometry type
	 * 
	 * @param geometryType
	 *            geometry type
	 * @return child geometry types, empty list if no child types
	 * @since 2.0.1
	 */
	public static List<GeometryType> childTypes(GeometryType geometryType) {

		List<GeometryType> childTypes = new ArrayList<>();

		switch (geometryType) {

		case GEOMETRY:
			childTypes.add(GeometryType.POINT);
			childTypes.add(GeometryType.GEOMETRYCOLLECTION);
			childTypes.add(GeometryType.CURVE);
			childTypes.add(GeometryType.SURFACE);
			break;
		case POINT:
			break;
		case LINESTRING:
			childTypes.add(GeometryType.CIRCULARSTRING);
			break;
		case POLYGON:
			childTypes.add(GeometryType.TRIANGLE);
			break;
		case MULTIPOINT:
			break;
		case MULTILINESTRING:
			break;
		case MULTIPOLYGON:
			break;
		case GEOMETRYCOLLECTION:
			childTypes.add(GeometryType.MULTIPOINT);
			childTypes.add(GeometryType.MULTICURVE);
			childTypes.add(GeometryType.MULTISURFACE);
			break;
		case CIRCULARSTRING:
			break;
		case COMPOUNDCURVE:
			break;
		case CURVEPOLYGON:
			childTypes.add(GeometryType.POLYGON);
			break;
		case MULTICURVE:
			childTypes.add(GeometryType.MULTILINESTRING);
			break;
		case MULTISURFACE:
			childTypes.add(GeometryType.MULTIPOLYGON);
			break;
		case CURVE:
			childTypes.add(GeometryType.LINESTRING);
			childTypes.add(GeometryType.COMPOUNDCURVE);
			break;
		case SURFACE:
			childTypes.add(GeometryType.CURVEPOLYGON);
			childTypes.add(GeometryType.POLYHEDRALSURFACE);
			break;
		case POLYHEDRALSURFACE:
			childTypes.add(GeometryType.TIN);
			break;
		case TIN:
			break;
		case TRIANGLE:
			break;
		default:
			throw new SFException(
					"Geometry Type not supported: " + geometryType);
		}

		return childTypes;
	}

	/**
	 * Serialize the geometry to bytes
	 * 
	 * @param geometry
	 *            geometry
	 * @return serialized bytes
	 * @since 2.0.1
	 */
	public static byte[] serialize(Geometry geometry) {

		byte[] bytes = null;

		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		ObjectOutput out = null;
		try {
			out = new ObjectOutputStream(bos);
			out.writeObject(geometry);
			out.flush();
			bytes = bos.toByteArray();
		} catch (IOException e) {
			throw new SFException("Failed to serialize geometry into bytes", e);
		} finally {
			try {
				bos.close();
			} catch (IOException e) {
				logger.log(Level.WARNING, "Failed to close stream", e);
			}
			if (out != null) {
				try {
					out.close();
				} catch (IOException e) {
					logger.log(Level.WARNING, "Failed to close stream", e);
				}
			}
		}

		return bytes;
	}

	/**
	 * Deserialize the bytes into a geometry
	 * 
	 * @param bytes
	 *            serialized bytes
	 * @return geometry
	 * @since 2.0.1
	 */
	public static Geometry deserialize(byte[] bytes) {

		Geometry geometry = null;

		ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
		ObjectInput in = null;
		try {
			in = new ObjectInputStream(bis);
			geometry = (Geometry) in.readObject();
		} catch (IOException | ClassNotFoundException e) {
			throw new SFException("Failed to deserialize geometry into bytes",
					e);
		} finally {
			try {
				bis.close();
			} catch (IOException e) {
				logger.log(Level.WARNING, "Failed to close stream", e);
			}
			if (in != null) {
				try {
					in.close();
				} catch (IOException e) {
					logger.log(Level.WARNING, "Failed to close stream", e);
				}
			}
		}

		return geometry;
	}

}
