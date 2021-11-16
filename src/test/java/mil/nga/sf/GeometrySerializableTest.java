package mil.nga.sf;

import java.io.IOException;

import mil.nga.sf.Geometry;
import mil.nga.sf.util.GeometryUtils;

import org.junit.Test;

/**
 * Geometry Collection tests
 * 
 * @author osbornb
 */
public class GeometrySerializableTest {

	@Test
	public void testPolygon() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createPolygon(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testLineString() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createLineString(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testPoint() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createPoint(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testGeometryCollection() throws ClassNotFoundException,
			IOException {

		testSerializable(SFTestUtils.createGeometryCollection(
				Math.random() < .5, Math.random() < .5));

	}

	@Test
	public void testMultiPolygon() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createMultiPolygon(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testMultiLineString() throws ClassNotFoundException,
			IOException {

		testSerializable(SFTestUtils.createMultiLineString(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testMultiPoint() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createMultiPoint(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testCurvePolygon() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createCurvePolygon(Math.random() < .5,
				Math.random() < .5));

	}

	@Test
	public void testCompoundCurve() throws ClassNotFoundException, IOException {

		testSerializable(SFTestUtils.createCompoundCurve(Math.random() < .5,
				Math.random() < .5));

	}

	private void testSerializable(Geometry geometry) throws IOException,
			ClassNotFoundException {

		byte[] bytes = GeometryUtils.serialize(geometry);
		Geometry deserializedGeometry = GeometryUtils.deserialize(bytes);

		SFTestUtils.compareGeometries(geometry, deserializedGeometry);

	}

}
