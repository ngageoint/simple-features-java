package mil.nga.sf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import junit.framework.TestCase;
import mil.nga.sf.util.GeometryConstants;
import mil.nga.sf.util.GeometryEnvelopeBuilder;
import mil.nga.sf.util.GeometryUtils;

/**
 * Test Geometry Utilities
 * 
 * @author osbornb
 */
public class GeometryUtilsTest {

	/**
	 * Number of random geometries to create for each test
	 */
	private static final int GEOMETRIES_PER_TEST = 10;

	/**
	 * Constructor
	 */
	public GeometryUtilsTest() {

	}

	/**
	 * Test point centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testPointCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a point
			Point point = SFTestUtils.createPoint(SFTestUtils.coinFlip(),
					SFTestUtils.coinFlip());
			TestCase.assertEquals(0, GeometryUtils.getDimension(point));
			TestCase.assertEquals(0, point.getDimension());
			geometryCentroidTester(point);
		}

	}

	/**
	 * Test line string centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testLineStringCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a line string
			LineString lineString = SFTestUtils.createLineString(
					SFTestUtils.coinFlip(), SFTestUtils.coinFlip());
			TestCase.assertEquals(1, GeometryUtils.getDimension(lineString));
			TestCase.assertEquals(1, lineString.getDimension());
			geometryCentroidTester(lineString);
		}

	}

	/**
	 * Test polygon centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testPolygonCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a polygon
			Polygon polygon = createPolygon();
			TestCase.assertEquals(2, GeometryUtils.getDimension(polygon));
			TestCase.assertEquals(2, polygon.getDimension());
			geometryCentroidTester(polygon);
		}

	}

	/**
	 * Test multi point centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testMultiPointCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a multi point
			MultiPoint multiPoint = SFTestUtils.createMultiPoint(
					SFTestUtils.coinFlip(), SFTestUtils.coinFlip());
			TestCase.assertEquals(0, GeometryUtils.getDimension(multiPoint));
			TestCase.assertEquals(0, multiPoint.getDimension());
			geometryCentroidTester(multiPoint);
		}

	}

	/**
	 * Test multi line string centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testMultiLineStringCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a multi line string
			MultiLineString multiLineString = SFTestUtils.createMultiLineString(
					SFTestUtils.coinFlip(), SFTestUtils.coinFlip());
			TestCase.assertEquals(1,
					GeometryUtils.getDimension(multiLineString));
			TestCase.assertEquals(1, multiLineString.getDimension());
			geometryCentroidTester(multiLineString);
		}

	}

	/**
	 * Test multi polygon centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testMultiPolygonCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a multi polygon
			MultiPolygon multiPolygon = createMultiPolygon();
			TestCase.assertEquals(2, GeometryUtils.getDimension(multiPolygon));
			TestCase.assertEquals(2, multiPolygon.getDimension());
			geometryCentroidTester(multiPolygon);
		}

	}

	/**
	 * Test geometry collection centroid
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testGeometryCollectionCentroid() throws IOException {

		for (int i = 0; i < GEOMETRIES_PER_TEST; i++) {
			// Create and test a geometry collection
			GeometryCollection<Geometry> geometryCollection = createGeometryCollection(
					SFTestUtils.coinFlip(), SFTestUtils.coinFlip());
			TestCase.assertEquals(
					GeometryUtils.getDimension(geometryCollection),
					geometryCollection.getDimension());
			geometryCentroidTester(geometryCollection);
		}

	}

	/**
	 * Test polygon centroid with and without hole
	 * 
	 * @throws IOException
	 *             upon error
	 */
	@Test
	public void testPolygonCentroidWithAndWithoutHole() throws IOException {

		Polygon polygon = new Polygon();
		LineString lineString = new LineString();
		lineString.addPoint(new Point(-90, 45));
		lineString.addPoint(new Point(-90, -45));
		lineString.addPoint(new Point(90, -45));
		lineString.addPoint(new Point(90, 45));
		polygon.addRing(lineString);

		TestCase.assertEquals(2, GeometryUtils.getDimension(polygon));
		TestCase.assertEquals(2, polygon.getDimension());
		Point centroid = geometryCentroidTester(polygon);

		TestCase.assertEquals(0.0, centroid.getX());
		TestCase.assertEquals(0.0, centroid.getY());

		LineString holeLineString = new LineString();
		holeLineString.addPoint(new Point(0, 45));
		holeLineString.addPoint(new Point(0, 0));
		holeLineString.addPoint(new Point(90, 0));
		holeLineString.addPoint(new Point(90, 45));
		polygon.addRing(holeLineString);

		TestCase.assertEquals(2, GeometryUtils.getDimension(polygon));
		TestCase.assertEquals(2, polygon.getDimension());
		centroid = geometryCentroidTester(polygon);

		TestCase.assertEquals(-15.0, centroid.getX());
		TestCase.assertEquals(-7.5, centroid.getY());
	}

	/**
	 * Test the geometry centroid
	 * 
	 * @param geometry
	 * @throws IOException
	 *             upon error
	 */
	private Point geometryCentroidTester(Geometry geometry) throws IOException {

		Point point = GeometryUtils.getCentroid(geometry);
		TestCase.assertEquals(point, geometry.getCentroid());

		GeometryEnvelope envelope = geometry.getEnvelope();

		if (geometry.getGeometryType() == GeometryType.POINT) {
			TestCase.assertEquals(envelope.getMinX(), point.getX());
			TestCase.assertEquals(envelope.getMaxX(), point.getX());
			TestCase.assertEquals(envelope.getMinY(), point.getY());
			TestCase.assertEquals(envelope.getMaxY(), point.getY());
		}

		TestCase.assertTrue(point.getX() >= envelope.getMinX());
		TestCase.assertTrue(point.getX() <= envelope.getMaxX());
		TestCase.assertTrue(point.getY() >= envelope.getMinY());
		TestCase.assertTrue(point.getY() <= envelope.getMaxY());

		Point envelopeCentroid1 = envelope.buildGeometry().getCentroid();
		Point envelopeCentroid2 = envelope.getCentroid();
		TestCase.assertEquals(envelopeCentroid1.getX(),
				envelopeCentroid2.getX(), 0.0000000000001);
		TestCase.assertEquals(envelopeCentroid1.getY(),
				envelopeCentroid2.getY(), 0.0000000000001);

		return point;
	}

	/**
	 * Create a polygon
	 * 
	 * @return polygon
	 */
	private static Polygon createPolygon() {

		Polygon polygon = new Polygon();
		LineString lineString = new LineString();
		lineString.addPoint(createPoint(-180.0, 45.0, 90.0, 45.0));
		lineString.addPoint(createPoint(-180.0, -90.0, 90.0, 45.0));
		lineString.addPoint(createPoint(90.0, -90.0, 90.0, 45.0));
		lineString.addPoint(createPoint(90.0, 45.0, 90.0, 45.0));
		polygon.addRing(lineString);

		LineString holeLineString = new LineString();
		holeLineString.addPoint(createPoint(-90.0, 0.0, 90.0, 45.0));
		holeLineString.addPoint(createPoint(-90.0, -45.0, 90.0, 45.0));
		holeLineString.addPoint(createPoint(0.0, -45.0, 90.0, 45.0));
		holeLineString.addPoint(createPoint(0.0, 0.0, 90.0, 45.0));
		polygon.addRing(holeLineString);

		return polygon;
	}

	/**
	 * Create a point
	 * 
	 * @param minX
	 *            min x
	 * @param minY
	 *            min y
	 * @param xRange
	 *            x range
	 * @param yRange
	 *            y range
	 * @return point
	 */
	private static Point createPoint(double minX, double minY, double xRange,
			double yRange) {

		double x = minX + (Math.random() * xRange);
		double y = minY + (Math.random() * yRange);

		Point point = new Point(x, y);

		return point;
	}

	/**
	 * Create a multi polygon
	 * 
	 * @return multi polygon
	 */
	private static MultiPolygon createMultiPolygon() {

		MultiPolygon multiPolygon = new MultiPolygon();

		int num = 1 + ((int) (Math.random() * 5));

		for (int i = 0; i < num; i++) {
			multiPolygon.addPolygon(createPolygon());
		}

		return multiPolygon;
	}

	/**
	 * Create a geometry collection
	 * 
	 * @param hasZ
	 *            has z
	 * @param hasM
	 *            has m
	 * @return geometry collection
	 */
	private static GeometryCollection<Geometry> createGeometryCollection(
			boolean hasZ, boolean hasM) {

		GeometryCollection<Geometry> geometryCollection = new GeometryCollection<Geometry>(
				hasZ, hasM);

		int num = 1 + ((int) (Math.random() * 5));

		for (int i = 0; i < num; i++) {

			Geometry geometry = null;
			int randomGeometry = (int) (Math.random() * 6);

			switch (randomGeometry) {
			case 0:
				geometry = SFTestUtils.createPoint(hasZ, hasM);
				break;
			case 1:
				geometry = SFTestUtils.createLineString(hasZ, hasM);
				break;
			case 2:
				geometry = createPolygon();
				break;
			case 3:
				geometry = SFTestUtils.createMultiPoint(hasZ, hasM);
				break;
			case 4:
				geometry = SFTestUtils.createMultiLineString(hasZ, hasM);
				break;
			case 5:
				geometry = createMultiPolygon();
				break;
			}

			geometryCollection.addGeometry(geometry);
		}

		return geometryCollection;
	}

	/**
	 * Test copy minimize and normalize
	 */
	@Test
	public void testCopyMinimizeAndNormalize() {

		Polygon polygon = new Polygon();
		LineString ring = new LineString();
		double random = Math.random();
		if (random < .5) {
			ring.addPoint(createPoint(90.0, 0.0, 90.0, 90.0));
			ring.addPoint(createPoint(90.0, -90.0, 90.0, 90.0));
			ring.addPoint(createPoint(-180.0, -90.0, 89.0, 90.0));
			ring.addPoint(createPoint(-180.0, 0.0, 89.0, 90.0));
		} else {
			ring.addPoint(createPoint(-180.0, 0.0, 89.0, 90.0));
			ring.addPoint(createPoint(-180.0, -90.0, 89.0, 90.0));
			ring.addPoint(createPoint(90.0, -90.0, 90.0, 90.0));
			ring.addPoint(createPoint(90.0, 0.0, 90.0, 90.0));
		}
		polygon.addRing(ring);

		Polygon polygon2 = (Polygon) polygon.copy();
		GeometryUtils.minimizeWGS84(polygon2);

		Polygon polygon3 = (Polygon) polygon2.copy();
		GeometryUtils.normalizeWGS84(polygon3);

		List<Point> points = ring.getPoints();
		LineString ring2 = polygon2.getRings().get(0);
		List<Point> points2 = ring2.getPoints();
		LineString ring3 = polygon3.getRings().get(0);
		List<Point> points3 = ring3.getPoints();

		for (int i = 0; i < points.size(); i++) {

			Point point = points.get(i);
			Point point2 = points2.get(i);
			Point point3 = points3.get(i);

			TestCase.assertEquals(point.getY(), point2.getY(), .0000000001);
			TestCase.assertEquals(point.getY(), point3.getY(), .0000000001);
			TestCase.assertEquals(point.getX(), point3.getX(), .0000000001);
			if (i < 2) {
				TestCase.assertEquals(point.getX(), point2.getX(), .0000000001);
			} else {
				double point2Value = point2.getX();
				if (random < .5) {
					point2Value -= 360.0;
				} else {
					point2Value += 360.0;
				}
				TestCase.assertEquals(point.getX(), point2Value, .0000000001);
			}
		}

	}

	/**
	 * Test simplify points
	 */
	@Test
	public void testSimplifyPoints() {

		final double halfWorldWidth = GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH;

		List<Point> points = new ArrayList<>();
		List<Double> distances = new ArrayList<>();

		double x = (Math.random() * halfWorldWidth * 2) - halfWorldWidth;
		double y = (Math.random() * halfWorldWidth * 2) - halfWorldWidth;
		Point point = new Point(x, y);
		points.add(point);

		for (int i = 1; i < 100; i++) {

			double xChange = 100000.0 * Math.random()
					* (Math.random() < .5 ? 1 : -1);
			x += xChange;

			double yChange = 100000.0 * Math.random()
					* (Math.random() < .5 ? 1 : -1);
			y += yChange;
			if (y > halfWorldWidth || y < -halfWorldWidth) {
				y -= 2 * yChange;
			}

			Point previousPoint = point;
			point = new Point(x, y);
			points.add(point);

			double distance = GeometryUtils.distance(previousPoint, point);
			distances.add(distance);

		}

		List<Double> sortedDistances = new ArrayList<>(distances);
		Collections.sort(sortedDistances);
		double tolerance = sortedDistances.get(sortedDistances.size() / 2);

		List<Point> simplifiedPoints = GeometryUtils.simplifyPoints(points,
				tolerance);
		TestCase.assertTrue(simplifiedPoints.size() <= points.size());

		Point firstPoint = points.get(0);
		Point lastPoint = points.get(points.size() - 1);
		Point firstSimplifiedPoint = simplifiedPoints.get(0);
		Point lastSimplifiedPoint = simplifiedPoints
				.get(simplifiedPoints.size() - 1);

		TestCase.assertEquals(firstPoint.getX(), firstSimplifiedPoint.getX());
		TestCase.assertEquals(firstPoint.getY(), firstSimplifiedPoint.getY());
		TestCase.assertEquals(lastPoint.getX(), lastSimplifiedPoint.getX());
		TestCase.assertEquals(lastPoint.getY(), lastSimplifiedPoint.getY());

		int pointIndex = 0;
		for (int i = 1; i < simplifiedPoints.size(); i++) {
			Point simplifiedPoint = simplifiedPoints.get(i);
			double simplifiedDistance = GeometryUtils
					.distance(simplifiedPoints.get(i - 1), simplifiedPoint);
			TestCase.assertTrue(simplifiedDistance >= tolerance);

			for (pointIndex++; pointIndex < points.size(); pointIndex++) {
				Point newPoint = points.get(pointIndex);
				if (newPoint.getX() == simplifiedPoint.getX()
						&& newPoint.getY() == simplifiedPoint.getY()) {
					break;
				}
			}
			TestCase.assertTrue(pointIndex < points.size());
		}

	}

	/**
	 * Test point in polygon
	 */
	@Test
	public void testPointInPolygon() {

		List<Point> points = new ArrayList<>();
		points.add(new Point(0, 5));
		points.add(new Point(5, 0));
		points.add(new Point(10, 5));
		points.add(new Point(5, 10));

		TestCase.assertFalse(GeometryUtils.closedPolygon(points));

		double deviation = 0.000000000000001;

		for (Point point : points) {
			TestCase.assertTrue(GeometryUtils.pointInPolygon(point, points));
		}

		TestCase.assertTrue(GeometryUtils
				.pointInPolygon(new Point(0 + deviation, 5), points));
		TestCase.assertTrue(GeometryUtils
				.pointInPolygon(new Point(5, 0 + deviation), points));
		TestCase.assertTrue(GeometryUtils
				.pointInPolygon(new Point(10 - deviation, 5), points));
		TestCase.assertTrue(GeometryUtils
				.pointInPolygon(new Point(5, 10 - deviation), points));

		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(5, 5), points));

		TestCase.assertTrue(GeometryUtils.pointInPolygon(
				new Point(2.5 + deviation, 7.5 - deviation), points));
		TestCase.assertTrue(GeometryUtils.pointInPolygon(
				new Point(2.5 + deviation, 2.5 + deviation), points));
		TestCase.assertTrue(GeometryUtils.pointInPolygon(
				new Point(7.5 - deviation, 2.5 + deviation), points));
		TestCase.assertTrue(GeometryUtils.pointInPolygon(
				new Point(7.5 - deviation, 7.5 - deviation), points));

		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(2.5, 7.5), points));
		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(2.5, 2.5), points));
		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(7.5, 2.5), points));
		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(7.5, 7.5), points));

		deviation = .0000001;

		TestCase.assertFalse(
				GeometryUtils.pointInPolygon(new Point(0, 0), points));
		TestCase.assertFalse(GeometryUtils
				.pointInPolygon(new Point(0 - deviation, 5), points));
		TestCase.assertFalse(GeometryUtils
				.pointInPolygon(new Point(5, 0 - deviation), points));
		TestCase.assertFalse(GeometryUtils
				.pointInPolygon(new Point(10 + deviation, 5), points));
		TestCase.assertFalse(GeometryUtils
				.pointInPolygon(new Point(5, 10 + deviation), points));

		TestCase.assertFalse(GeometryUtils.pointInPolygon(
				new Point(2.5 - deviation, 7.5 + deviation), points));
		TestCase.assertFalse(GeometryUtils.pointInPolygon(
				new Point(2.5 - deviation, 2.5 - deviation), points));
		TestCase.assertFalse(GeometryUtils.pointInPolygon(
				new Point(7.5 + deviation, 2.5 - deviation), points));
		TestCase.assertFalse(GeometryUtils.pointInPolygon(
				new Point(7.5 + deviation, 7.5 + deviation), points));

		Point firstPoint = points.get(0);
		points.add(new Point(firstPoint.getX(), firstPoint.getY()));

		TestCase.assertTrue(GeometryUtils.closedPolygon(points));

		for (Point point : points) {
			TestCase.assertTrue(GeometryUtils.pointInPolygon(point, points));
		}
		TestCase.assertTrue(GeometryUtils.pointInPolygon(
				new Point(2.5 + deviation, 7.5 - deviation), points));
		TestCase.assertTrue(
				GeometryUtils.pointInPolygon(new Point(2.5, 7.5), points));
		TestCase.assertFalse(GeometryUtils.pointInPolygon(
				new Point(2.5 - deviation, 7.5 + deviation), points));

	}

	/**
	 * Test closed polygon
	 */
	@Test
	public void testClosedPolygon() {

		List<Point> points = new ArrayList<>();
		points.add(new Point(0.1, 0.2));
		points.add(new Point(5.3, 0.4));
		points.add(new Point(5.5, 5.6));

		TestCase.assertFalse(GeometryUtils.closedPolygon(points));

		Point firstPoint = points.get(0);
		points.add(new Point(firstPoint.getX(), firstPoint.getY()));

		TestCase.assertTrue(GeometryUtils.closedPolygon(points));
	}

	/**
	 * Test point on line
	 */
	@Test
	public void testPointOnLine() {

		List<Point> points = new ArrayList<>();
		points.add(new Point(0, 0));
		points.add(new Point(5, 0));
		points.add(new Point(5, 5));

		for (Point point : points) {
			TestCase.assertTrue(GeometryUtils.pointOnLine(point, points));
		}
		TestCase.assertTrue(
				GeometryUtils.pointOnLine(new Point(2.5, 0), points));
		TestCase.assertTrue(
				GeometryUtils.pointOnLine(new Point(5, 2.5), points));
		TestCase.assertTrue(
				GeometryUtils.pointOnLine(new Point(2.5, 0.00000001), points));
		TestCase.assertFalse(
				GeometryUtils.pointOnLine(new Point(2.5, 0.0000001), points));
		TestCase.assertTrue(
				GeometryUtils.pointOnLine(new Point(5, 2.500000001), points));
		TestCase.assertFalse(
				GeometryUtils.pointOnLine(new Point(5, 2.50000001), points));
		TestCase.assertTrue(GeometryUtils
				.pointOnLine(new Point(-0.0000000000000001, 0), points));
		TestCase.assertFalse(GeometryUtils
				.pointOnLine(new Point(-0.000000000000001, 0), points));
		TestCase.assertTrue(GeometryUtils
				.pointOnLine(new Point(5, 5.0000000000000001), points));
		TestCase.assertFalse(GeometryUtils
				.pointOnLine(new Point(5, 5.000000000000001), points));

	}

	/**
	 * Test line intersection
	 */
	@Test
	public void testIntersection() {

		Line line1 = new Line(new Point(-40.0, -40.0), new Point(40.0, 40.0));
		Line line2 = new Line(new Point(-40.0, 40.0), new Point(40.0, -40.0));

		Point point = GeometryUtils.intersection(line1, line2);
		assertEquals(0.0, point.getX(), 0.0);
		assertEquals(0.0, point.getY(), 0.0);

		line1 = new Line(GeometryUtils.degreesToMeters(line1).getPoints());
		line2 = new Line(GeometryUtils.degreesToMeters(line2).getPoints());

		point = GeometryUtils.intersection(line1, line2);
		assertEquals(0.0, point.getX(), 0.0);
		assertEquals(0.0, point.getY(), 0.0);

		line1 = new Line(new Point(-40.0, -10.0), new Point(20.0, 70.0));
		line2 = new Line(new Point(-40.0, 70.0), new Point(20.0, -10.0));

		point = GeometryUtils.intersection(line1, line2);
		assertEquals(-10.0, point.getX(), 0.0);
		assertEquals(30.0, point.getY(), 0.0);

		line1 = GeometryUtils.degreesToMeters(line1);
		line2 = GeometryUtils.degreesToMeters(line2);

		point = GeometryUtils.intersection(line1, line2);
		assertEquals(-1113194.9079327362, point.getX(), 0.0);
		assertEquals(4974912.842260765, point.getY(), 0.0);

		point = GeometryUtils.metersToDegrees(point);
		assertEquals(-10.0, point.getX(), 0.00000000000001);
		assertEquals(40.745756618323014, point.getY(), 0.0);

		line1 = new Line(
				new Point(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
						GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2),
				new Point(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
						GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH));
		line2 = new Line(
				new Point(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
						GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH),
				new Point(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
						GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2));

		point = GeometryUtils.intersection(line1, line2);
		assertEquals(0.75 * -GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				point.getX(), 0.00000001);
		assertEquals(0.75 * GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				point.getY(), 0.00000001);

		point = GeometryUtils.metersToDegrees(point);
		assertEquals(-135.0, point.getX(), 0.0);
		assertEquals(79.17133464081945, point.getY(), 0.0);

	}

	/**
	 * Test point conversion
	 */
	@Test
	public void testConversion() {

		Point point = new Point(-112.500003, 21.943049);

		Point point2 = GeometryUtils.degreesToMeters(point);
		assertEquals(-12523443.048201751, point2.getX(), 0.0);
		assertEquals(2504688.958883909, point2.getY(), 0.0);

		Point point3 = GeometryUtils.metersToDegrees(point2);
		assertEquals(-112.500003, point3.getX(), 0.0000000000001);
		assertEquals(21.943049, point3.getY(), 0.0000000000001);

	}

	/**
	 * Test crop
	 */
	@Test
	public void testCrop() {

		double side = 40;
		double sideC = Math.sqrt(2 * Math.pow(side, 2));

		double min = 10;
		double max = min + side;
		double mid = (min + max) / 2;

		double extraWidth = (sideC - side) / 2;
		double starLength = 11.7157287525381;

		// Test with two squares forming a star (Star of Lakshmi)
		// Polygon as the diamond, square as the crop bounds

		Polygon polygon = new Polygon();
		LineString ring = new LineString();
		ring.addPoint(new Point(min - extraWidth, mid));
		ring.addPoint(new Point(mid, min - extraWidth));
		ring.addPoint(new Point(max + extraWidth, mid));
		ring.addPoint(new Point(mid, max + extraWidth));
		polygon.addRing(ring);

		GeometryEnvelope envelope = new GeometryEnvelope(min, min, max, max);

		Polygon crop = GeometryUtils.crop(polygon, envelope);

		LineString cropRing = crop.getRing(0);
		assertEquals(9, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(min + starLength, cropRing.getPoint(0).getY(),
				0.00000000001);

		assertEquals(min + starLength, cropRing.getPoint(1).getX(),
				0.00000000001);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.000001);

		assertEquals(max - starLength, cropRing.getPoint(2).getX(),
				0.00000000001);
		assertEquals(min, cropRing.getPoint(2).getY(), 0.000001);

		assertEquals(max, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(min + starLength, cropRing.getPoint(3).getY(),
				0.00000000001);

		assertEquals(max, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(max - starLength, cropRing.getPoint(4).getY(),
				0.00000000001);

		assertEquals(max - starLength, cropRing.getPoint(5).getX(),
				0.00000000001);
		assertEquals(max, cropRing.getPoint(5).getY(), 0.0);

		assertEquals(min + starLength, cropRing.getPoint(6).getX(),
				0.00000000001);
		assertEquals(max, cropRing.getPoint(6).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(7).getX(), 0.0);
		assertEquals(max - starLength, cropRing.getPoint(7).getY(),
				0.00000000001);

		assertEquals(min, cropRing.getPoint(8).getX(), 0.0);
		assertEquals(min + starLength, cropRing.getPoint(8).getY(),
				0.00000000001);

		crop = GeometryUtils.crop(GeometryUtils.degreesToMeters(polygon),
				GeometryUtils.degreesToMeters(envelope.buildGeometry())
						.getEnvelope());
		crop = GeometryUtils.metersToDegrees(crop);

		cropRing = crop.getRing(0);
		assertEquals(9, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(22.181521688501903, cropRing.getPoint(0).getY(), 0.0);

		assertEquals(22.077332337134834, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.000001);

		assertEquals(37.922667662865166, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(2).getY(), 0.000001);

		assertEquals(max, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(22.181521688501903, cropRing.getPoint(3).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(39.74197667744292, cropRing.getPoint(4).getY(), 0.0);

		assertEquals(39.88507567296453, cropRing.getPoint(5).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(5).getY(), 0.0);

		assertEquals(20.114924327035485, cropRing.getPoint(6).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(6).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(7).getX(), 0.0);
		assertEquals(39.74197667744289, cropRing.getPoint(7).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(8).getX(), 0.0);
		assertEquals(22.181521688501903, cropRing.getPoint(8).getY(), 0.0);

		// Test with a diamond fully fitting within the crop bounds

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(min, mid));
		ring.addPoint(new Point(mid, min));
		ring.addPoint(new Point(max, mid));
		ring.addPoint(new Point(mid, max));
		polygon.addRing(ring);

		crop = GeometryUtils.crop(polygon, envelope);

		cropRing = crop.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(0).getY(), 0.0);

		assertEquals(mid, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(mid, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(3).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(4).getY(), 0.0);

		crop = GeometryUtils.crop(GeometryUtils.degreesToMeters(polygon),
				GeometryUtils.degreesToMeters(envelope.buildGeometry())
						.getEnvelope());
		crop = GeometryUtils.metersToDegrees(crop);

		cropRing = crop.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(0).getY(), 0.0);

		assertEquals(mid, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(mid, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(3).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(mid, cropRing.getPoint(4).getY(), 0.0);

		// Test with a star (Star of Lakshmi outer border) polygon and square as
		// the crop bounds

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(min - extraWidth, mid));
		ring.addPoint(new Point(min, min + extraWidth));
		ring.addPoint(new Point(min, min));
		ring.addPoint(new Point(min + extraWidth, min));
		ring.addPoint(new Point(mid, min - extraWidth));
		ring.addPoint(new Point(max - extraWidth, min));
		ring.addPoint(new Point(max, min));
		ring.addPoint(new Point(max, min + extraWidth));
		ring.addPoint(new Point(max + extraWidth, mid));
		ring.addPoint(new Point(max, max - extraWidth));
		ring.addPoint(new Point(max, max));
		ring.addPoint(new Point(max - extraWidth, max));
		ring.addPoint(new Point(mid, max + extraWidth));
		ring.addPoint(new Point(min + extraWidth, max));
		ring.addPoint(new Point(min, max));
		ring.addPoint(new Point(min, max - extraWidth));
		polygon.addRing(ring);

		crop = GeometryUtils.crop(polygon, envelope);

		cropRing = crop.getRing(0);
		assertEquals(6, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(min + extraWidth, cropRing.getPoint(0).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(3).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(4).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(5).getX(), 0.0);
		assertEquals(min + extraWidth, cropRing.getPoint(5).getY(), 0.0);

		crop = GeometryUtils.crop(GeometryUtils.degreesToMeters(polygon),
				GeometryUtils.degreesToMeters(envelope.buildGeometry())
						.getEnvelope());
		crop = GeometryUtils.metersToDegrees(crop);

		cropRing = crop.getRing(0);
		assertEquals(9, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(min, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(min + extraWidth, cropRing.getPoint(0).getY(),
				0.00000000001);

		assertEquals(min, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(max - extraWidth, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(min, cropRing.getPoint(3).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(max - extraWidth, cropRing.getPoint(4).getY(), 0.0);

		assertEquals(max, cropRing.getPoint(5).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(5).getY(), 0.0);

		assertEquals(min + extraWidth, cropRing.getPoint(6).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(6).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(7).getX(), 0.0);
		assertEquals(max, cropRing.getPoint(7).getY(), 0.0);

		assertEquals(min, cropRing.getPoint(8).getX(), 0.0);
		assertEquals(min + extraWidth, cropRing.getPoint(8).getY(),
				0.00000000001);

	}

	/**
	 * Test crop over the international date line
	 */
	@Test
	public void testCropIDL() {

		Polygon polygon = new Polygon();
		LineString ring = new LineString();
		ring.addPoint(new Point(-168.967, 67.0));
		ring.addPoint(new Point(-168.967, 90.0));
		ring.addPoint(new Point(-180.0000000001, 90.0000001148));
		ring.addPoint(new Point(-180.0000000001, 67.0));
		ring.addPoint(new Point(-168.967, 67.0));
		polygon.addRing(ring);

		Polygon meters = GeometryUtils.degreesToMeters(polygon);
		Polygon crop = (Polygon) GeometryUtils.cropWebMercator(meters);
		Polygon degrees = GeometryUtils.metersToDegrees(crop);
		GeometryUtils.minimizeWGS84(degrees);

		LineString cropRing = degrees.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(-168.967, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(67.0, cropRing.getPoint(0).getY(), 0.00000000001);

		assertEquals(-168.967, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE,
				cropRing.getPoint(1).getY(), 0.00000000001);

		assertEquals(-180.0000000001, cropRing.getPoint(2).getX(),
				0.00000000001);
		assertEquals(GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE,
				cropRing.getPoint(2).getY(), 0.00000000001);

		assertEquals(-180.0000000001, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(67.0, cropRing.getPoint(3).getY(), 0.00000000001);

		assertEquals(-168.967, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(67.0, cropRing.getPoint(4).getY(), 0.00000000001);

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(-18809320.400867056, 10156058.722522344));
		ring.addPoint(new Point(-18809320.400867056, 238107693.26496765));
		ring.addPoint(new Point(-20037508.342800375, 238107693.26496765));
		ring.addPoint(new Point(-20037508.342800375, 10156058.722522344));
		ring.addPoint(new Point(-18809320.400867056, 10156058.722522344));
		polygon.addRing(ring);

		GeometryEnvelope envelope = GeometryUtils.webMercatorEnvelope();
		envelope.setMinX(-20037508.342800375);

		crop = GeometryUtils.crop(polygon, envelope);

		cropRing = crop.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(-18809320.400867056, cropRing.getPoint(0).getX(), 0.0);
		assertEquals(10156058.722522344, cropRing.getPoint(0).getY(),
				0.00000000001);

		assertEquals(-18809320.400867056, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(1).getY(), 0.00000001);

		assertEquals(-20037508.342800375, cropRing.getPoint(2).getX(),
				0.00000000001);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(2).getY(), 0.00000001);

		assertEquals(-20037508.342800375, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(10156058.722522344, cropRing.getPoint(3).getY(),
				0.00000000001);

		assertEquals(-18809320.400867056, cropRing.getPoint(4).getX(), 0.0);
		assertEquals(10156058.722522344, cropRing.getPoint(4).getY(),
				0.00000000001);

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(-120.0, -90.0));
		ring.addPoint(new Point(-120.0, 0.0));
		ring.addPoint(new Point(-180.0, 0.0));
		ring.addPoint(new Point(-180.0, -90.0));
		ring.addPoint(new Point(-120.0, -90.0));
		polygon.addRing(ring);

		meters = GeometryUtils.degreesToMeters(polygon);
		crop = (Polygon) GeometryUtils.cropWebMercator(meters);
		degrees = GeometryUtils.metersToDegrees(crop);
		GeometryUtils.minimizeWGS84(degrees);

		cropRing = degrees.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(-120.0, cropRing.getPoint(0).getX(), 0.00000000001);
		assertEquals(GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				cropRing.getPoint(0).getY(), 0.0);

		assertEquals(-120.0, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(0.0, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(-180.0, cropRing.getPoint(2).getX(), 0.0);
		assertEquals(0.0, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(-180.0, cropRing.getPoint(3).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				cropRing.getPoint(3).getY(), 0.0);

		assertEquals(-120.0, cropRing.getPoint(4).getX(), 0.00000000001);
		assertEquals(GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				cropRing.getPoint(4).getY(), 0.0);

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(-13358338.89519283, -233606567.09255272));
		ring.addPoint(new Point(-13358338.89519283, 0.0));
		ring.addPoint(new Point(
				-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH, 0.0));
		ring.addPoint(
				new Point(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
						-233606567.09255272));
		ring.addPoint(new Point(-13358338.89519283, -233606567.09255272));
		polygon.addRing(ring);

		crop = (Polygon) GeometryUtils.cropWebMercator(polygon);

		cropRing = crop.getRing(0);
		assertEquals(5, cropRing.numPoints());
		assertTrue(cropRing.isClosed());

		assertEquals(-13358338.89519283, cropRing.getPoint(0).getX(),
				0.00000001);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(0).getY(), 0.00000001);

		assertEquals(-13358338.89519283, cropRing.getPoint(1).getX(), 0.0);
		assertEquals(0.0, cropRing.getPoint(1).getY(), 0.0);

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(2).getX(), 0.0);
		assertEquals(0.0, cropRing.getPoint(2).getY(), 0.0);

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(3).getX(), 0.0);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(3).getY(), 0.00000001);

		assertEquals(-13358338.89519283, cropRing.getPoint(4).getX(),
				0.00000001);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				cropRing.getPoint(4).getY(), 0.00000001);

	}

	/**
	 * Test bound
	 */
	@Test
	public void testBound() {

		Polygon polygon = new Polygon();
		LineString ring = new LineString();
		ring.addPoint(new Point(-180.01, 90.01));
		ring.addPoint(new Point(-90.0, 0.0));
		ring.addPoint(new Point(-181.0, -91.0));
		ring.addPoint(new Point(0, -45.0));
		ring.addPoint(new Point(180.00000000001, -90.00000000001));
		ring.addPoint(new Point(90.0, 0.0));
		ring.addPoint(new Point(180.0, 90.0));
		ring.addPoint(new Point(0, 45.0));
		ring.addPoint(new Point(-180.01, 90.01));
		polygon.addRing(ring);

		assertFalse(
				GeometryUtils.wgs84Envelope().contains(polygon.getEnvelope()));

		Polygon bounded = (Polygon) polygon.copy();

		GeometryUtils.boundWGS84(bounded);

		LineString boundedRing = bounded.getRing(0);
		assertEquals(ring.numPoints(), boundedRing.numPoints());
		assertTrue(boundedRing.isClosed());
		assertTrue(
				GeometryUtils.wgs84Envelope().contains(bounded.getEnvelope()));

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(0).getX(), 0.0);
		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				boundedRing.getPoint(0).getY(), 0.0);

		assertEquals(-90.0, boundedRing.getPoint(1).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(1).getY(), 0.0);

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(2).getX(), 0.0);
		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				boundedRing.getPoint(2).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(3).getX(), 0.0);
		assertEquals(-45.0, boundedRing.getPoint(3).getY(), 0.0);

		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(4).getX(), 0.0);
		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				boundedRing.getPoint(4).getY(), 0.0);

		assertEquals(90.0, boundedRing.getPoint(5).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(5).getY(), 0.0);

		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(6).getX(), 0.0);
		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				boundedRing.getPoint(6).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(7).getX(), 0.0);
		assertEquals(45.0, boundedRing.getPoint(7).getY(), 0.0);

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(8).getX(), 0.0);
		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LAT_HEIGHT,
				boundedRing.getPoint(8).getY(), 0.0);

		assertFalse(GeometryUtils.wgs84EnvelopeWithWebMercator()
				.contains(polygon.getEnvelope()));

		bounded = (Polygon) polygon.copy();

		GeometryUtils.boundWGS84WithWebMercator(bounded);

		boundedRing = bounded.getRing(0);
		assertEquals(ring.numPoints(), boundedRing.numPoints());
		assertTrue(boundedRing.isClosed());
		assertTrue(GeometryUtils.wgs84EnvelopeWithWebMercator()
				.contains(bounded.getEnvelope()));

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(0).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE,
				boundedRing.getPoint(0).getY(), 0.0);

		assertEquals(-90.0, boundedRing.getPoint(1).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(1).getY(), 0.0);

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(2).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				boundedRing.getPoint(2).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(3).getX(), 0.0);
		assertEquals(-45.0, boundedRing.getPoint(3).getY(), 0.0);

		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(4).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MIN_LAT_RANGE,
				boundedRing.getPoint(4).getY(), 0.0);

		assertEquals(90.0, boundedRing.getPoint(5).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(5).getY(), 0.0);

		assertEquals(GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(6).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE,
				boundedRing.getPoint(6).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(7).getX(), 0.0);
		assertEquals(45.0, boundedRing.getPoint(7).getY(), 0.0);

		assertEquals(-GeometryConstants.WGS84_HALF_WORLD_LON_WIDTH,
				boundedRing.getPoint(8).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_MAX_LAT_RANGE,
				boundedRing.getPoint(8).getY(), 0.0);

		polygon = new Polygon();
		ring = new LineString();
		ring.addPoint(new Point(-20037508.35, 20037508.35));
		ring.addPoint(new Point(-10018754.171394622, 0.0));
		ring.addPoint(new Point(-20037509, -20037509));
		ring.addPoint(new Point(0, -10018754.171394622));
		ring.addPoint(new Point(20037508.34278925, -20037508.34278925));
		ring.addPoint(new Point(10018754.171394622, 0.0));
		ring.addPoint(new Point(20037508.342789244, 20037508.342789244));
		ring.addPoint(new Point(0, 10018754.171394622));
		ring.addPoint(new Point(-20037508.35, 20037508.35));
		polygon.addRing(ring);

		assertFalse(GeometryUtils.webMercatorEnvelope()
				.contains(polygon.getEnvelope()));

		bounded = (Polygon) polygon.copy();

		GeometryUtils.boundWebMercator(bounded);

		boundedRing = bounded.getRing(0);
		assertEquals(ring.numPoints(), boundedRing.numPoints());
		assertTrue(boundedRing.isClosed());
		assertTrue(GeometryUtils.webMercatorEnvelope()
				.contains(bounded.getEnvelope()));

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(0).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(0).getY(), 0.0);

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
				boundedRing.getPoint(1).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(1).getY(), 0.00000001);

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(2).getX(), 0.0);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(2).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(3).getX(), 0.0);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
				boundedRing.getPoint(3).getY(), 0.0);

		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(4).getX(), 0.0);
		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(4).getY(), 0.0);

		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
				boundedRing.getPoint(5).getX(), 0.0);
		assertEquals(0.0, boundedRing.getPoint(5).getY(), 0.0);

		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(6).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(6).getY(), 0.0);

		assertEquals(0.0, boundedRing.getPoint(7).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH / 2,
				boundedRing.getPoint(7).getY(), 0.0);

		assertEquals(-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(8).getX(), 0.0);
		assertEquals(GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				boundedRing.getPoint(8).getY(), 0.0);

	}

	/**
	 * Test the geometry type parent and child hierarchy methods
	 */
	@Test
	public void testHierarchy() {

		for (GeometryType geometryType : GeometryType.values()) {

			GeometryType parentType = GeometryUtils.parentType(geometryType);
			List<GeometryType> parentHierarchy = GeometryUtils
					.parentHierarchy(geometryType);

			GeometryType previousParentType = null;

			while (parentType != null) {
				TestCase.assertEquals(parentType, parentHierarchy.get(0));

				if (previousParentType != null) {
					List<GeometryType> childTypes = GeometryUtils
							.childTypes(parentType);
					TestCase.assertTrue(
							childTypes.contains(previousParentType));
					Map<GeometryType, Map<GeometryType, ?>> childHierarchy = GeometryUtils
							.childHierarchy(parentType);
					TestCase.assertTrue(
							childHierarchy.containsKey(previousParentType));
				}

				previousParentType = parentType;
				parentType = GeometryUtils.parentType(previousParentType);
				parentHierarchy = GeometryUtils
						.parentHierarchy(previousParentType);

			}
			TestCase.assertTrue(parentHierarchy.isEmpty());

			Map<GeometryType, Map<GeometryType, ?>> childHierarchy = GeometryUtils
					.childHierarchy(geometryType);
			testChildHierarchy(geometryType, childHierarchy);

		}

	}

	/**
	 * Test the child hierarchy recursively
	 * 
	 * @param geometryType
	 *            geometry type
	 * @param childHierarchy
	 *            child hierarchy
	 */
	private void testChildHierarchy(GeometryType geometryType,
			Map<GeometryType, ?> childHierarchy) {
		List<GeometryType> childTypes = GeometryUtils.childTypes(geometryType);
		if (childTypes.isEmpty()) {
			TestCase.assertNull(childHierarchy);
		} else {
			TestCase.assertEquals(childTypes.size(), childHierarchy.size());
			for (GeometryType childType : childTypes) {
				TestCase.assertTrue(childHierarchy.containsKey(childType));

				TestCase.assertEquals(geometryType,
						GeometryUtils.parentType(childType));
				TestCase.assertEquals(geometryType,
						GeometryUtils.parentHierarchy(childType).get(0));

				testChildHierarchy(childType,
						GeometryUtils.childHierarchy(childType));
			}
		}
	}

	/**
	 * Test centroid and degrees centroid
	 */
	@Test
	public void testCentroid() {

		Point point = new Point(15, 35);

		Point centroid = point.getCentroid();

		TestCase.assertEquals(15.0, centroid.getX());
		TestCase.assertEquals(35.0, centroid.getY());

		Point degreesCentroid = point.getDegreesCentroid();

		TestCase.assertEquals(15.0, degreesCentroid.getX());
		TestCase.assertEquals(35.0, degreesCentroid.getY());

		LineString lineString = new LineString();
		lineString.addPoint(new Point(0, 5));
		lineString.addPoint(point);

		centroid = lineString.getCentroid();

		TestCase.assertEquals(7.5, centroid.getX());
		TestCase.assertEquals(20.0, centroid.getY());

		degreesCentroid = lineString.getDegreesCentroid();

		TestCase.assertEquals(6.764392425440724, degreesCentroid.getX());
		TestCase.assertEquals(20.157209770845522, degreesCentroid.getY());

		lineString.addPoint(new Point(2, 65));

		centroid = lineString.getCentroid();

		TestCase.assertEquals(7.993617921179541, centroid.getX());
		TestCase.assertEquals(34.808537635386266, centroid.getY());

		degreesCentroid = lineString.getDegreesCentroid();

		TestCase.assertEquals(5.85897989020252, degreesCentroid.getX());
		TestCase.assertEquals(35.20025371999032, degreesCentroid.getY());

		Polygon polygon = new Polygon(lineString);

		centroid = polygon.getCentroid();

		TestCase.assertEquals(5.666666666666667, centroid.getX());
		TestCase.assertEquals(35.0, centroid.getY());

		degreesCentroid = polygon.getDegreesCentroid();

		TestCase.assertEquals(5.85897989020252, degreesCentroid.getX());
		TestCase.assertEquals(35.20025371999032, degreesCentroid.getY());

		lineString.addPoint(new Point(-20, 40));
		lineString.addPoint(new Point(0, 5));

		centroid = polygon.getCentroid();

		TestCase.assertEquals(-1.3554502369668247, centroid.getX());
		TestCase.assertEquals(36.00315955766193, centroid.getY());

		degreesCentroid = polygon.getDegreesCentroid();

		TestCase.assertEquals(-0.6891904581641471, degreesCentroid.getX());
		TestCase.assertEquals(37.02524099014426, degreesCentroid.getY());

	}

	/**
	 * Test distance Haversine
	 */
	@Test
	public void testDistanceHaversine() {

		testDistanceHaversine(-73.779, 40.640, 103.989, 1.359,
				15356717.55865963);

		testDistanceHaversine(-61.207542, 15.526518, -18.124573, 27.697002,
				4633776.207109179);

		testDistanceHaversine(-115.49, 39.64, 52.98, -22.69,
				17858784.720537618);

	}

	/**
	 * Test distance Haversine
	 * 
	 * @param lon1
	 *            longitude 1
	 * @param lat1
	 *            latitude 1
	 * @param lon2
	 *            longitude 2
	 * @param lat2
	 *            latitude 2
	 * @param expectedDistance
	 *            expected distance
	 */
	private void testDistanceHaversine(double lon1, double lat1, double lon2,
			double lat2, double expectedDistance) {

		Point point1 = new Point(lon1, lat1);
		Point point2 = new Point(lon2, lat2);

		double distance = GeometryUtils.distanceHaversine(point1, point2);

		assertEquals(expectedDistance, distance, 0.0);

	}

	/**
	 * Test bearing
	 */
	@Test
	public void testBearing() {

		testBearing(-73.779, 40.640, 103.989, 1.359, 3.3326543286976857);

		testBearing(-61.207542, 15.526518, -18.124573, 27.697002,
				65.56992873258116);

		testBearing(-115.49, 39.64, 52.98, -22.69, 33.401404803852586);

	}

	/**
	 * Test bearing
	 * 
	 * @param lon1
	 *            longitude 1
	 * @param lat1
	 *            latitude 1
	 * @param lon2
	 *            longitude 2
	 * @param lat2
	 *            latitude 2
	 * @param expectedBearing
	 *            expected bearing
	 */
	private void testBearing(double lon1, double lat1, double lon2, double lat2,
			double expectedBearing) {

		Point point1 = new Point(lon1, lat1);
		Point point2 = new Point(lon2, lat2);

		double bearing = GeometryUtils.bearing(point1, point2);

		assertEquals(expectedBearing, bearing, 0.0);

	}

	/**
	 * Test midpoint
	 */
	@Test
	public void testMidpoint() {

		testMidpoint(-73.779, 40.640, 103.989, 1.359, 97.01165658499957,
				70.180706801119);

		testMidpoint(-61.207542, 15.526518, -18.124573, 27.697002,
				-40.62120498446578, 23.06700073523901);

		testMidpoint(-115.49, 39.64, 52.98, -22.69, 10.497130585764902,
				47.89844929382955);

	}

	/**
	 * Test midpoint
	 * 
	 * @param lon1
	 *            longitude 1
	 * @param lat1
	 *            latitude 1
	 * @param lon2
	 *            longitude 2
	 * @param lat2
	 *            latitude 2
	 * @param expectedLon
	 *            expected longitude
	 * @param expectedLat
	 *            expected latitude
	 */
	private void testMidpoint(double lon1, double lat1, double lon2,
			double lat2, double expectedLon, double expectedLat) {

		Point point1 = new Point(lon1, lat1);
		Point point2 = new Point(lon2, lat2);

		Point midpoint = GeometryUtils.geodesicMidpoint(point1, point2);

		assertEquals(expectedLon, midpoint.getX(), 0.0);
		assertEquals(expectedLat, midpoint.getY(), 0.0);

		Point point1Radians = GeometryUtils.degreesToRadians(point1);
		Point point2Radians = GeometryUtils.degreesToRadians(point2);

		Point midpointRadians = GeometryUtils
				.geodesicMidpointRadians(point1Radians, point2Radians);

		assertEquals(GeometryUtils.degreesToRadians(midpoint), midpointRadians);

	}

	/**
	 * Test geodesic path
	 */
	@Test
	public void testGeodesicPath() {

		final double METERS_TEST = 2500000;

		Point point1 = new Point(-73.779, 40.640);
		Point point2 = new Point(103.989, 1.359);

		List<Point> path = GeometryUtils.geodesicPath(point1, point2,
				100000000);

		assertEquals(2, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));

		path = GeometryUtils.geodesicPath(point1, point2, 10000000);

		assertEquals(3, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(97.01165658499957, path.get(1).getX(), 0.0);
		assertEquals(70.180706801119, path.get(1).getY(), 0.0);

		path = GeometryUtils.geodesicPath(point1, point2, METERS_TEST);

		final int PATH_COUNT_1 = 9;
		assertEquals(PATH_COUNT_1, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(-71.92354211648598, path.get(1).getX(), 0.0);
		assertEquals(57.84299258409111, path.get(1).getY(), 0.0);
		assertEquals(-66.48824612217787, path.get(2).getX(), 0.0);
		assertEquals(74.96658102555037, path.get(2).getY(), 0.0);
		assertEquals(57.819970305247146, path.get(3).getX(), 0.0);
		assertEquals(86.50086370806171, path.get(3).getY(), 0.0);
		assertEquals(97.01165658499957, path.get(4).getX(), 0.0);
		assertEquals(70.180706801119, path.get(4).getY(), 0.0);
		assertEquals(100.68758599469604, path.get(5).getX(), 0.0);
		assertEquals(53.01802965780041, path.get(5).getY(), 0.0);
		assertEquals(102.22354481789006, path.get(6).getX(), 0.0);
		assertEquals(35.80797373447713, path.get(6).getY(), 0.0);
		assertEquals(103.19828968828496, path.get(7).getX(), 0.0);
		assertEquals(18.585518760662953, path.get(7).getY(), 0.0);

		point1 = new Point(-61.207542, 15.526518);
		point2 = new Point(-18.124573, 27.697002);

		path = GeometryUtils.geodesicPath(point1, point2, METERS_TEST);

		final int PATH_COUNT_2 = 3;
		assertEquals(PATH_COUNT_2, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(-40.62120498446578, path.get(1).getX(), 0.0);
		assertEquals(23.06700073523901, path.get(1).getY(), 0.0);

		path = GeometryUtils.geodesicPath(point1, point2, 1500000);

		assertEquals(5, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(-51.154455380800286, path.get(1).getX(), 0.0);
		assertEquals(19.58837936536169, path.get(1).getY(), 0.0);
		assertEquals(-40.62120498446578, path.get(2).getX(), 0.0);
		assertEquals(23.06700073523901, path.get(2).getY(), 0.0);
		assertEquals(-29.591449129559173, path.get(3).getX(), 0.0);
		assertEquals(25.814850212565023, path.get(3).getY(), 0.0);

		point1 = new Point(-115.49, 39.64);
		point2 = new Point(52.98, -22.69);

		path = GeometryUtils.geodesicPath(point1, point2, 10000000);

		assertEquals(3, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(10.497130585764902, path.get(1).getX(), 0.0);
		assertEquals(47.89844929382955, path.get(1).getY(), 0.0);

		path = GeometryUtils.geodesicPath(point1, point2, METERS_TEST);

		final int PATH_COUNT_3 = 9;
		assertEquals(PATH_COUNT_3, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));
		assertEquals(-96.24706669925085, path.get(1).getX(), 0.0);
		assertEquals(55.05735117318652, path.get(1).getY(), 0.0);
		assertEquals(-60.22371275346116, path.get(2).getX(), 0.0);
		assertEquals(64.43472917406754, path.get(2).getY(), 0.0);
		assertEquals(-16.117214263945066, path.get(3).getX(), 0.0);
		assertEquals(61.05440885911931, path.get(3).getY(), 0.0);
		assertEquals(10.497130585764902, path.get(4).getX(), 0.0);
		assertEquals(47.89844929382955, path.get(4).getY(), 0.0);
		assertEquals(25.189402921803275, path.get(5).getX(), 0.0);
		assertEquals(31.25647851927572, path.get(5).getY(), 0.0);
		assertEquals(35.259366393531955, path.get(6).getX(), 0.0);
		assertEquals(13.465910465448717, path.get(6).getY(), 0.0);
		assertEquals(43.88449185422331, path.get(7).getX(), 0.0);
		assertEquals(-4.667462203083873, path.get(7).getY(), 0.0);

		LineString lineString = new LineString();

		int expectedPoints = 0;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);

		point1 = new Point(-73.779, 40.640);
		lineString.addPoint(point1);
		expectedPoints++;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));

		point2 = new Point(103.989, 1.359);
		lineString.addPoint(point2);
		expectedPoints += PATH_COUNT_1 - 1;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(path.size() - 1));

		Point point3 = new Point(-61.207542, 15.526518);
		lineString.addPoint(point3);
		final int PATH_COUNT_1_2 = GeometryUtils
				.geodesicPath(point2, point3, METERS_TEST).size();
		expectedPoints += PATH_COUNT_1_2 - 1;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(PATH_COUNT_1 - 1));
		assertEquals(point3, path.get(path.size() - 1));

		Point point4 = new Point(-18.124573, 27.697002);
		lineString.addPoint(point4);
		expectedPoints += PATH_COUNT_2 - 1;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(PATH_COUNT_1 - 1));
		assertEquals(point3, path.get(PATH_COUNT_1 + PATH_COUNT_1_2 - 2));
		assertEquals(point4, path.get(path.size() - 1));

		Point point5 = new Point(-115.49, 39.64);
		lineString.addPoint(point5);
		final int PATH_COUNT_2_3 = GeometryUtils
				.geodesicPath(point4, point5, METERS_TEST).size();
		expectedPoints += PATH_COUNT_2_3 - 1;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(PATH_COUNT_1 - 1));
		assertEquals(point3, path.get(PATH_COUNT_1 + PATH_COUNT_1_2 - 2));
		assertEquals(point4,
				path.get(PATH_COUNT_1 + PATH_COUNT_1_2 + PATH_COUNT_2 - 3));
		assertEquals(point5, path.get(path.size() - 1));

		Point point6 = new Point(52.98, -22.69);
		lineString.addPoint(point6);
		expectedPoints += PATH_COUNT_3 - 1;

		path = GeometryUtils.geodesicPath(lineString, METERS_TEST);
		assertEquals(expectedPoints, path.size());
		pathDuplicateCheck(path);
		assertEquals(point1, path.get(0));
		assertEquals(point2, path.get(PATH_COUNT_1 - 1));
		assertEquals(point3, path.get(PATH_COUNT_1 + PATH_COUNT_1_2 - 2));
		assertEquals(point4,
				path.get(PATH_COUNT_1 + PATH_COUNT_1_2 + PATH_COUNT_2 - 3));
		assertEquals(point5, path.get(PATH_COUNT_1 + PATH_COUNT_1_2
				+ PATH_COUNT_2 + PATH_COUNT_2_3 - 4));
		assertEquals(point6, path.get(path.size() - 1));

	}

	/**
	 * Test geodesic envelope
	 */
	@Test
	public void testGeodesicEnvelope() {

		testGeodesicEnvelope(new GeometryEnvelope());

		testGeodesicEnvelope(new GeometryEnvelope(-85, 0, 85, 0));

		testGeodesicEnvelope(new GeometryEnvelope(0, -45, 0, 45));

		testGeodesicEnvelope(new GeometryEnvelope(-85, -45, 85, 45),
				new GeometryEnvelope(-85, -85.0189306062998, 85,
						85.0189306062998));

		testGeodesicEnvelope(new GeometryEnvelope(0, 40, 60, 60),
				new GeometryEnvelope(0, 40, 60, 63.43494882292201));

		testGeodesicEnvelope(
				new GeometryEnvelope(-116.564009, 52.257876, 21.002792,
						55.548544),
				new GeometryEnvelope(-116.564009, 52.257876, 21.002792,
						76.05697069912907));

		testGeodesicEnvelope(
				new GeometryEnvelope(-0.118092, 1.290270, 103.851959,
						51.509865),
				new GeometryEnvelope(-0.118092, 1.290270, 103.851959,
						63.908548725225884));

		testGeodesicEnvelope(
				new GeometryEnvelope(-71.038887, -33.92584, 18.42322,
						42.364506),
				new GeometryEnvelope(-71.038887, -43.43480368259327, 18.42322,
						52.08227546634191));

		testGeodesicEnvelope(
				new GeometryEnvelope(-65.116900, -54.656860, 13.008587,
						-9.120679),
				new GeometryEnvelope(-65.116900, -61.16106506177795, 13.008587,
						-9.120679));

		testGeodesicEnvelope(
				new GeometryEnvelope(-69.001773, -51.614743, 120.316646,
						22.794475),
				new GeometryEnvelope(-69.001773, -86.31825003835286, 120.316646,
						79.0603064734963));

	}

	/**
	 * Test geodesic geometry envelope expansion
	 * 
	 * @param envelope
	 *            geometry envelope and expected geometry envelope
	 */
	private void testGeodesicEnvelope(GeometryEnvelope envelope) {
		testGeodesicEnvelope(envelope, envelope);
	}

	/**
	 * Test geodesic geometry envelope expansion
	 * 
	 * @param envelope
	 *            geometry envelope
	 * @param expected
	 *            expected geometry envelope
	 */
	private void testGeodesicEnvelope(GeometryEnvelope envelope,
			GeometryEnvelope expected) {

		final double distancePixels = 512.0;
		final double delta = 0.0000000000001;

		GeometryEnvelope geodesic = GeometryUtils.geodesicEnvelope(envelope);
		compareEnvelopes(expected, geodesic, delta);

		if (envelope.getMaxX() - envelope.getMinX() <= 180.0) {

			// Top line
			double distance = GeometryUtils.distanceHaversine(
					envelope.getTopLeft(), envelope.getTopRight());
			double maxDistance = distance / distancePixels;
			List<Point> path = GeometryUtils.geodesicPath(envelope.getTopLeft(),
					envelope.getTopRight(), maxDistance);
			GeometryEnvelope pathEnvelope = GeometryEnvelopeBuilder
					.buildEnvelope(new LineString(path));

			// Bottom line
			distance = GeometryUtils.distanceHaversine(envelope.getBottomLeft(),
					envelope.getBottomRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getBottomLeft(),
					envelope.getBottomRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

			// The rest of the line tests below are not really needed, but
			// included as extra sanity checks

			// Left line
			distance = GeometryUtils.distanceHaversine(envelope.getTopLeft(),
					envelope.getBottomLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopLeft(),
					envelope.getBottomLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			// Right line
			distance = GeometryUtils.distanceHaversine(envelope.getTopRight(),
					envelope.getBottomRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopRight(),
					envelope.getBottomRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

			// Diagonal line
			distance = GeometryUtils.distanceHaversine(envelope.getTopLeft(),
					envelope.getBottomRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopLeft(),
					envelope.getBottomRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			// Other diagonal line
			distance = GeometryUtils.distanceHaversine(envelope.getTopRight(),
					envelope.getBottomLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopRight(),
					envelope.getBottomLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

			// Mid horizontal line
			distance = GeometryUtils.distanceHaversine(envelope.getLeftMid(),
					envelope.getRightMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getLeftMid(),
					envelope.getRightMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			// Mid vertical line
			distance = GeometryUtils.distanceHaversine(envelope.getTopMid(),
					envelope.getBottomMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopMid(),
					envelope.getBottomMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

			// Mid to neighbor mid lines

			distance = GeometryUtils.distanceHaversine(envelope.getLeftMid(),
					envelope.getTopMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getLeftMid(),
					envelope.getTopMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getTopMid(),
					envelope.getRightMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopMid(),
					envelope.getRightMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getRightMid(),
					envelope.getBottomMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getRightMid(),
					envelope.getBottomMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getBottomMid(),
					envelope.getLeftMid());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getBottomMid(),
					envelope.getLeftMid(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

			// Mid to corner lines

			distance = GeometryUtils.distanceHaversine(envelope.getLeftMid(),
					envelope.getTopRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getLeftMid(),
					envelope.getTopRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getLeftMid(),
					envelope.getBottomRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getLeftMid(),
					envelope.getBottomRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getTopMid(),
					envelope.getBottomLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopMid(),
					envelope.getBottomLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getTopMid(),
					envelope.getBottomRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getTopMid(),
					envelope.getBottomRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getRightMid(),
					envelope.getTopLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getRightMid(),
					envelope.getTopLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getRightMid(),
					envelope.getBottomLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getRightMid(),
					envelope.getBottomLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getBottomMid(),
					envelope.getTopLeft());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getBottomMid(),
					envelope.getTopLeft(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			distance = GeometryUtils.distanceHaversine(envelope.getBottomMid(),
					envelope.getTopRight());
			maxDistance = distance / distancePixels;
			path = GeometryUtils.geodesicPath(envelope.getBottomMid(),
					envelope.getTopRight(), maxDistance);
			GeometryEnvelopeBuilder.buildEnvelope(new LineString(path),
					pathEnvelope);

			compareEnvelopes(pathEnvelope, geodesic, delta);

		}

	}

	/**
	 * Compare the envelopes for equality
	 * 
	 * @param expected
	 *            expected geometry envelope
	 * @param envelope
	 *            geometry envelope
	 * @param delta
	 *            comparison delta
	 */
	private void compareEnvelopes(GeometryEnvelope expected,
			GeometryEnvelope envelope, double delta) {
		if (!envelope.equals(expected)) {
			assertEquals(expected.getMinX(), envelope.getMinX(), delta);
			assertEquals(expected.getMinY(), envelope.getMinY(), delta);
			assertEquals(expected.getMaxX(), envelope.getMaxX(), delta);
			assertEquals(expected.getMaxY(), envelope.getMaxY(), delta);
		}
	}

	/**
	 * Check that there are no back to back duplicate points
	 * 
	 * @param path
	 *            point path
	 */
	private void pathDuplicateCheck(List<Point> path) {

		for (int i = 0; i < path.size() - 1; i++) {
			assertNotEquals(path.get(i), path.get(i + 1));
		}

	}

}
