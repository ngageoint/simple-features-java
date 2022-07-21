package mil.nga.sf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import junit.framework.TestCase;
import mil.nga.sf.util.GeometryConstants;
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
		GeometryUtils.minimizeWGS84Geometry(polygon2);

		Polygon polygon3 = (Polygon) polygon2.copy();
		GeometryUtils.normalizeWGS84Geometry(polygon3);

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

		GeometryEnvelope envelope = new GeometryEnvelope(
				-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				-GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH,
				GeometryConstants.WEB_MERCATOR_HALF_WORLD_WIDTH);

		Polygon meters = GeometryUtils.degreesToMeters(polygon);
		Polygon crop = GeometryUtils.crop(meters, envelope);
		Polygon degrees = GeometryUtils.metersToDegrees(crop);
		GeometryUtils.minimizeWGS84Geometry(degrees);

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

}
