package mil.nga.sf.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Write a byte array
 * 
 * @author osbornb
 */
public class ByteWriter {

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(ByteWriter.class.getName());

	/**
	 * Default write byte order
	 * 
	 * @since 2.0.3
	 */
	public static final ByteOrder DEFAULT_BYTE_ORDER = ByteOrder.BIG_ENDIAN;

	/**
	 * Output stream to write bytes to
	 */
	private final OutputStream outputStream;

	/**
	 * Byte order
	 */
	private ByteOrder byteOrder = DEFAULT_BYTE_ORDER;

	/**
	 * Constructor
	 */
	public ByteWriter() {
		outputStream = new ByteArrayOutputStream();
	}

	/**
	 * Constructor
	 * 
	 * @param outputStream
	 *            output stream
	 * @since 2.0.3
	 */
	public ByteWriter(OutputStream outputStream) {
		this.outputStream = outputStream;
	}

	/**
	 * Constructor
	 * 
	 * @param byteOrder
	 *            byte order
	 * @since 2.0.3
	 */
	public ByteWriter(ByteOrder byteOrder) {
		this();
		this.byteOrder = byteOrder;
	}

	/**
	 * Constructor
	 * 
	 * @param outputStream
	 *            output stream
	 * @param byteOrder
	 *            byte order
	 * @since 2.0.3
	 */
	public ByteWriter(OutputStream outputStream, ByteOrder byteOrder) {
		this(outputStream);
		this.byteOrder = byteOrder;
	}

	/**
	 * Get the output stream
	 * 
	 * @return output stream
	 * @since 2.0.3
	 */
	public OutputStream getOutputStream() {
		return outputStream;
	}

	/**
	 * Get the output stream
	 * 
	 * @return output stream
	 * @since 2.0.3
	 */
	public ByteArrayOutputStream getByteArrayOutputStream() {
		if (!(outputStream instanceof ByteArrayOutputStream)) {
			throw new SFException(
					"Output Stream is not a ByteArrayOutputStream: "
							+ outputStream.getClass().getName());
		}
		return (ByteArrayOutputStream) outputStream;
	}

	/**
	 * Close the byte writer
	 */
	public void close() {
		try {
			outputStream.close();
		} catch (IOException e) {
			logger.log(Level.WARNING,
					"Failed to close byte writer output stream", e);
		}
	}

	/**
	 * Get the byte order
	 * 
	 * @return byte order
	 */
	public ByteOrder getByteOrder() {
		return byteOrder;
	}

	/**
	 * Set the byte order
	 * 
	 * @param byteOrder
	 *            byte order
	 */
	public void setByteOrder(ByteOrder byteOrder) {
		this.byteOrder = byteOrder;
	}

	/**
	 * Get the written bytes
	 * 
	 * @return written bytes
	 */
	public byte[] getBytes() {
		return getByteArrayOutputStream().toByteArray();
	}

	/**
	 * Get the current size in bytes written
	 * 
	 * @return bytes written
	 */
	public int size() {
		return getByteArrayOutputStream().size();
	}

	/**
	 * Write a String
	 * 
	 * @param value
	 *            string value
	 * @throws IOException
	 *             upon error
	 */
	public void writeString(String value) throws IOException {
		byte[] valueBytes = value.getBytes();
		outputStream.write(valueBytes);
	}

	/**
	 * Write a byte
	 * 
	 * @param value
	 *            byte
	 * @throws IOException
	 *             upon error
	 */
	public void writeByte(byte value) throws IOException {
		outputStream.write(value);
	}

	/**
	 * Write an integer
	 * 
	 * @param value
	 *            int value
	 * @throws IOException
	 *             upon error
	 */
	public void writeInt(int value) throws IOException {
		byte[] valueBytes = new byte[4];
		ByteBuffer byteBuffer = ByteBuffer.allocate(4).order(byteOrder)
				.putInt(value);
		byteBuffer.flip();
		byteBuffer.get(valueBytes);
		outputStream.write(valueBytes);
	}

	/**
	 * Write a double
	 * 
	 * @param value
	 *            double
	 * @throws IOException
	 *             upon error
	 */
	public void writeDouble(double value) throws IOException {
		byte[] valueBytes = new byte[8];
		ByteBuffer byteBuffer = ByteBuffer.allocate(8).order(byteOrder)
				.putDouble(value);
		byteBuffer.flip();
		byteBuffer.get(valueBytes);
		outputStream.write(valueBytes);
	}

}
