package mil.nga.sf.util.sweep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import junit.framework.TestCase;
import mil.nga.sf.LineString;
import mil.nga.sf.Point;
import mil.nga.sf.Polygon;

/**
 * Test Shamos-Hoey simple polygon methods
 * 
 * @author osbornb
 */
public class ShamosHoeyTest {

	@Test
	public void testSimple() throws IOException {

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 1, 0);
		addPoint(points, .5, 1);

		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		simple(new Polygon(new LineString(points)));
		TestCase.assertEquals(3, points.size());

		addPoint(points, 0, 0);

		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		simple(new Polygon(new LineString(points)));
		TestCase.assertEquals(4, points.size());

		points.clear();

		addPoint(points, 0, 100);
		addPoint(points, 100, 0);
		addPoint(points, 200, 100);
		addPoint(points, 100, 200);
		addPoint(points, 0, 100);

		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		simple(new Polygon(new LineString(points)));
		TestCase.assertEquals(5, points.size());

		points.clear();

		addPoint(points, -104.8384094, 39.753657);
		addPoint(points, -104.8377228, 39.7354422);
		addPoint(points, -104.7930908, 39.7364983);
		addPoint(points, -104.8233891, 39.7440222);
		addPoint(points, -104.7930908, 39.7369603);
		addPoint(points, -104.808197, 39.7541849);
		addPoint(points, -104.8383236, 39.753723);

		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		simple(new Polygon(new LineString(points)));
		TestCase.assertEquals(7, points.size());

		points.clear();

		addPoint(points, -106.3256836, 40.2962865);
		addPoint(points, -105.6445313, 38.5911138);
		addPoint(points, -105.0842285, 40.3046654);
		addPoint(points, -105.6445313, 38.5911139);

		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		simple(new Polygon(new LineString(points)));
		TestCase.assertEquals(4, points.size());

	}

	@Test
	public void testNonSimple() throws IOException {

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(1, points.size());

		addPoint(points, 1, 0);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(2, points.size());

		addPoint(points, 0, 0);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(3, points.size());

		points.clear();

		addPoint(points, 0, 100);
		addPoint(points, 100, 0);
		addPoint(points, 200, 100);
		addPoint(points, 100, 200);
		addPoint(points, 100.01, 200);
		addPoint(points, 0, 100);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(6, points.size());

		points.clear();

		addPoint(points, -104.8384094, 39.753657);
		addPoint(points, -104.8377228, 39.7354422);
		addPoint(points, -104.7930908, 39.7364983);
		addPoint(points, -104.8233891, 39.7440222);
		addPoint(points, -104.8034763, 39.7387424);
		addPoint(points, -104.7930908, 39.7369603);
		addPoint(points, -104.808197, 39.7541849);
		addPoint(points, -104.8383236, 39.753723);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(8, points.size());

		points.clear();

		addPoint(points, -106.3256836, 40.2962865);
		addPoint(points, -105.6445313, 38.5911138);
		addPoint(points, -105.0842285, 40.3046654);
		addPoint(points, -105.6445313, 38.5911138);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(4, points.size());

		points.clear();

		addPoint(points, 1, 0);
		addPoint(points, 0, 1);
		addPoint(points, 1, 0);
		addPoint(points, 2, 2);

		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		complex(new Polygon(new LineString(points)));
		TestCase.assertEquals(4, points.size());

	}

	@Test
	public void testSimpleHole() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints = new ArrayList<>();

		addPoint(holePoints, 1, 1);
		addPoint(holePoints, 9, 1);
		addPoint(holePoints, 5, 9);

		LineString hole = new LineString();
		hole.setPoints(holePoints);

		polygon.addRing(hole);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());
	}

	@Test
	public void testNonSimpleHole() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints = new ArrayList<>();

		addPoint(holePoints, 1, 1);
		addPoint(holePoints, 9, 1);
		addPoint(holePoints, 5, 9);
		addPoint(holePoints, 5.000001, 9);

		LineString hole = new LineString();
		hole.setPoints(holePoints);

		polygon.addRing(hole);

		TestCase.assertFalse(ShamosHoey.simplePolygon(polygon));
		complex(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(4, polygon.getRings().get(1).numPoints());
	}

	@Test
	public void testIntersectingHole() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints = new ArrayList<>();

		addPoint(holePoints, 1, 1);
		addPoint(holePoints, 9, 1);
		addPoint(holePoints, 5, 10);

		LineString hole = new LineString();
		hole.setPoints(holePoints);

		polygon.addRing(hole);

		TestCase.assertFalse(ShamosHoey.simplePolygon(polygon));
		complex(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());
	}

	@Test
	public void testIntersectingHoles() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints1 = new ArrayList<>();

		addPoint(holePoints1, 1, 1);
		addPoint(holePoints1, 9, 1);
		addPoint(holePoints1, 5, 9);

		LineString hole1 = new LineString();
		hole1.setPoints(holePoints1);

		polygon.addRing(hole1);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());

		List<Point> holePoints2 = new ArrayList<>();

		addPoint(holePoints2, 5.0, 0.1);
		addPoint(holePoints2, 6.0, 0.1);
		addPoint(holePoints2, 5.5, 1.00001);

		LineString hole2 = new LineString();
		hole2.setPoints(holePoints2);

		polygon.addRing(hole2);

		TestCase.assertFalse(ShamosHoey.simplePolygon(polygon));
		complex(polygon);
		TestCase.assertEquals(3, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(2).numPoints());
	}

	@Test
	public void testHoleInsideHole() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints1 = new ArrayList<>();

		addPoint(holePoints1, 1, 1);
		addPoint(holePoints1, 9, 1);
		addPoint(holePoints1, 5, 9);

		LineString hole1 = new LineString();
		hole1.setPoints(holePoints1);

		polygon.addRing(hole1);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());

		List<Point> holePoints2 = new ArrayList<>();

		addPoint(holePoints2, 2, 2);
		addPoint(holePoints2, 8, 2);
		addPoint(holePoints2, 5, 8);

		LineString hole2 = new LineString();
		hole2.setPoints(holePoints2);

		polygon.addRing(hole2);

		TestCase.assertFalse(ShamosHoey.simplePolygon(polygon));
		complex(polygon);
		TestCase.assertEquals(3, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(2).numPoints());
	}

	@Test
	public void testExternalHole() throws IOException {

		Polygon polygon = new Polygon();

		List<Point> points = new ArrayList<>();

		addPoint(points, 0, 0);
		addPoint(points, 10, 0);
		addPoint(points, 5, 10);

		LineString ring = new LineString();
		ring.setPoints(points);

		polygon.addRing(ring);

		TestCase.assertTrue(ShamosHoey.simplePolygon(polygon));
		simple(polygon);
		TestCase.assertEquals(1, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());

		List<Point> holePoints = new ArrayList<>();

		addPoint(holePoints, -1, 1);
		addPoint(holePoints, -1, 3);
		addPoint(holePoints, -2, 1);

		LineString hole = new LineString();
		hole.setPoints(holePoints);

		polygon.addRing(hole);

		TestCase.assertFalse(ShamosHoey.simplePolygon(polygon));
		complex(polygon);
		TestCase.assertEquals(2, polygon.numRings());
		TestCase.assertEquals(3, polygon.getRings().get(0).numPoints());
		TestCase.assertEquals(3, polygon.getRings().get(1).numPoints());
	}

	@Test
	public void testLargeSimple() throws IOException {

		double increment = .001;
		double radius = 1250;
		double x = -radius + increment;
		double y = 0;

		List<Point> points = new ArrayList<>();

		while (x <= radius) {
			if (x <= 0) {
				y -= increment;
			} else {
				y += increment;
			}
			addPoint(points, x, y);
			x += increment;
		}

		x = radius - increment;
		while (x >= -radius) {
			if (x >= 0) {
				y += increment;
			} else {
				y -= increment;
			}
			addPoint(points, x, y);
			x -= increment;
		}

		// Date before = new Date();
		TestCase.assertTrue(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertTrue(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertTrue(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertTrue((new LineString(points)).isSimple());
		TestCase.assertTrue((new Polygon(new LineString(points))).isSimple());
		// Date after = new Date();
		// long time = after.getTime() - before.getTime();
		// System.out.println("Points: " + points.size() + ", Time: " + time);
		TestCase.assertEquals((int) (radius / increment * 4), points.size());

	}

	@Test
	public void testLargeNonSimple() throws IOException {

		double increment = .001;
		double radius = 1250;
		double x = -radius + increment;
		double y = 0;

		List<Point> points = new ArrayList<>();

		while (x <= radius) {
			if (x <= 0) {
				y -= increment;
			} else {
				y += increment;
			}
			addPoint(points, x, y);
			x += increment;
		}

		Point previousPoint = points.get(points.size() - 2);
		int invalidIndex = points.size();
		addPoint(points, previousPoint.getX(),
				previousPoint.getY() - .000000000000001);

		x = radius - increment;
		while (x >= -radius) {
			if (x >= 0) {
				y += increment;
			} else {
				y -= increment;
			}
			addPoint(points, x, y);
			x -= increment;
		}

		// Date before = new Date();
		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		TestCase.assertFalse((new Polygon(new LineString(points))).isSimple());
		// Date after = new Date();
		// long time = after.getTime() - before.getTime();
		// System.out.println("Points: " + points.size() + ", Time: " + time);
		TestCase.assertEquals(1 + (int) (radius / increment * 4),
				points.size());

		points.remove(invalidIndex);
		previousPoint = points.get(points.size() - 3);
		addPoint(points, previousPoint.getX(),
				previousPoint.getY() + .000000000000001);

		// Date before2 = new Date();
		TestCase.assertFalse(ShamosHoey.simplePolygonPoints(points));
		TestCase.assertFalse(ShamosHoey.simplePolygon(new LineString(points)));
		TestCase.assertFalse(
				ShamosHoey.simplePolygon(new Polygon(new LineString(points))));
		TestCase.assertFalse((new LineString(points)).isSimple());
		TestCase.assertFalse((new Polygon(new LineString(points))).isSimple());
		// Date after2 = new Date();
		// long time2 = after2.getTime() - before2.getTime();
		// System.out.println("Points: " + points.size() + ", Time: " + time2);
		TestCase.assertEquals(1 + (int) (radius / increment * 4),
				points.size());
	}

	@Test
	public void testPolygons1() {

		LineString ring = new LineString();
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(0, 4));
		ring.addPoint(new Point(4, 0));
		Polygon polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(4, 0));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(0, 4));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(0, 4));
		ring.addPoint(new Point(4, 0));
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(4, 4));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(0, 4));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(4, 0));
		ring.addPoint(new Point(0, 0));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(1, 3));
		ring.addPoint(new Point(3, 1));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(3, 1));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(1, 3));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(1, 3));
		ring.addPoint(new Point(3, 1));
		ring.addPoint(new Point(0, 0));
		ring.addPoint(new Point(4, 4));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(1, 3));
		ring.addPoint(new Point(4, 4));
		ring.addPoint(new Point(3, 1));
		ring.addPoint(new Point(0, 0));
		polygon = new Polygon(ring);

		simple(polygon);

	}

	@Test
	public void testPolygons2() {

		LineString ring = new LineString();
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		ring.addPoint(new Point(120.94208471603682, 241.47274889005215));
		ring.addPoint(new Point(120.57389187028015, 240.42380619065557));
		ring.addPoint(new Point(120.40553233952696, 239.3249423106921));
		ring.addPoint(new Point(120.44278575100797, 238.2138802324909));
		ring.addPoint(new Point(120.68437322950616, 237.12876169126298));
		ring.addPoint(new Point(121.12200129996195, 236.1068378045508));
		ring.addPoint(new Point(121.74064659481407, 235.18319027576013));
		ring.addPoint(new Point(122.51907159233552, 234.38952707167576));
		Polygon polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		ring.addPoint(new Point(120.94208471603682, 241.47274889005215));
		// ring.addPoint(new Point(120.57389187028015, 240.42380619065557));
		// ring.addPoint(new Point(120.40553233952696, 239.3249423106921));
		ring.addPoint(new Point(120.44278575100797, 238.2138802324909));
		ring.addPoint(new Point(120.68437322950616, 237.12876169126298));
		ring.addPoint(new Point(121.12200129996195, 236.1068378045508));
		ring.addPoint(new Point(121.74064659481407, 235.18319027576013));
		ring.addPoint(new Point(122.51907159233552, 234.38952707167576));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		ring.addPoint(new Point(120.94208471603682, 241.47274889005215));
		// ring.addPoint(new Point(120.57389187028015, 240.42380619065557));
		// ring.addPoint(new Point(120.40553233952696, 239.3249423106921));
		ring.addPoint(new Point(120.44278575100797, 238.2138802324909));
		ring.addPoint(new Point(120.68437322950616, 237.12876169126298));
		ring.addPoint(new Point(121.12200129996195, 236.1068378045508));
		ring.addPoint(new Point(121.74064659481407, 235.18319027576013));
		ring.addPoint(new Point(122.51907159233552, 234.38952707167576));
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		ring.addPoint(new Point(120.94208471603682, 241.47274889005215));
		ring.addPoint(new Point(120.57389187028015, 240.42380619065557));
		// ring.addPoint(new Point(120.40553233952696, 239.3249423106921));
		ring.addPoint(new Point(120.44278575100797, 238.2138802324909));
		ring.addPoint(new Point(120.68437322950616, 237.12876169126298));
		ring.addPoint(new Point(121.12200129996195, 236.1068378045508));
		ring.addPoint(new Point(121.74064659481407, 235.18319027576013));
		ring.addPoint(new Point(122.51907159233552, 234.38952707167576));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(119.65450502825215, 234.97190110269844));
		ring.addPoint(new Point(120.94208471603682, 241.47274889005215));
		// ring.addPoint(new Point(120.57389187028015, 240.42380619065557));
		ring.addPoint(new Point(120.40553233952696, 239.3249423106921));
		ring.addPoint(new Point(120.44278575100797, 238.2138802324909));
		ring.addPoint(new Point(120.68437322950616, 237.12876169126298));
		ring.addPoint(new Point(121.12200129996195, 236.1068378045508));
		ring.addPoint(new Point(121.74064659481407, 235.18319027576013));
		ring.addPoint(new Point(122.51907159233552, 234.38952707167576));
		polygon = new Polygon(ring);

		complex(polygon);

	}

	@Test
	public void testPolygons3() {

		LineString ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(2.642, 4.22));
		ring.addPoint(new Point(4.41, 2.91));
		ring.addPoint(new Point(4.69, 4.42));
		ring.addPoint(new Point(2.68, 4.90));
		Polygon polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(2.641, 4.22));
		ring.addPoint(new Point(4.41, 2.91));
		ring.addPoint(new Point(4.69, 4.42));
		ring.addPoint(new Point(2.68, 4.90));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(2.642, 4.22));
		ring.addPoint(new Point(3.60, 4.68));
		ring.addPoint(new Point(4.69, 4.42));
		ring.addPoint(new Point(2.68, 4.90));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(2.642, 4.22));
		ring.addPoint(new Point(3.60, 4.681));
		ring.addPoint(new Point(4.69, 4.42));
		ring.addPoint(new Point(2.68, 4.90));
		polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(6.34, 2.15));
		ring.addPoint(new Point(5.85, 3.34));
		ring.addPoint(new Point(5.98, 3.07));
		ring.addPoint(new Point(5.09, 4.98));
		polygon = new Polygon(ring);

		simple(polygon);

		ring = new LineString();
		ring.addPoint(new Point(2.48, 1.38));
		ring.addPoint(new Point(6.34, 2.15));
		ring.addPoint(new Point(5.855, 3.34));
		ring.addPoint(new Point(5.98, 3.07));
		ring.addPoint(new Point(5.09, 4.98));
		polygon = new Polygon(ring);

		complex(polygon);

	}

	@Test
	public void testPolygons4() {

		LineString ring = new LineString();
		ring.addPoint(new Point(4.96, 4.18));
		ring.addPoint(new Point(2.68, 4.90));
		ring.addPoint(new Point(3.50, 2.64));
		ring.addPoint(new Point(5.20, 1.86));
		ring.addPoint(new Point(8.00, 2.83));
		ring.addPoint(new Point(3.50, 2.64));
		Polygon polygon = new Polygon(ring);

		complex(polygon);

		ring = new LineString();
		ring.addPoint(new Point(4.96, 4.18));
		ring.addPoint(new Point(2.68, 4.90));
		ring.addPoint(new Point(3.50, 2.64));
		ring.addPoint(new Point(5.20, 1.86));
		ring.addPoint(new Point(8.00, 2.83));
		ring.addPoint(new Point(3.500000000000001, 2.64));
		polygon = new Polygon(ring);

		simple(polygon);

	}

	private void simple(Polygon polygon) {
		validate(polygon, true);
	}

	private void complex(Polygon polygon) {
		validate(polygon, false);
	}

	private void validate(Polygon polygon, boolean simple) {

		TestCase.assertEquals(simple, polygon.isSimple());

		Polygon copy = new Polygon(polygon);
		List<Point> points = copy.getRing(0).getPoints();

		Point first = points.get(0);
		Point last = points.get(points.size() - 1);
		if (first.equalsXY(last)) {
			points.remove(points.size() - 1);
		}

		for (int i = 1; i < points.size(); i++) {

			points.add(points.remove(0));

			TestCase.assertEquals(simple, copy.isSimple());
		}

	}

	private void addPoint(List<Point> points, double x, double y) {
		points.add(new Point(x, y));
	}

}
