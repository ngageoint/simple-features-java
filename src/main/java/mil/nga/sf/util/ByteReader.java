package mil.nga.sf.util;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Read through a byte array
 * 
 * @author osbornb
 */
public class ByteReader {

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(ByteReader.class.getName());

	/**
	 * Default read byte order
	 * 
	 * @since 2.0.3
	 */
	public static final ByteOrder DEFAULT_BYTE_ORDER = ByteOrder.BIG_ENDIAN;

	/**
	 * Character set
	 */
	private static final String CHAR_SET = "UTF-8";

	/**
	 * Bytes to read from
	 */
	private final byte[] bytes;

	/**
	 * Input stream to read bytes from
	 */
	private final InputStream inputStream;

	/**
	 * Next byte index to read
	 */
	private int nextByte = 0;

	/**
	 * Byte order
	 */
	private ByteOrder byteOrder = DEFAULT_BYTE_ORDER;

	/**
	 * Constructor
	 * 
	 * @param bytes
	 *            bytes
	 */
	public ByteReader(byte[] bytes) {
		this.bytes = bytes;
		inputStream = new ByteArrayInputStream(bytes);
	}

	/**
	 * Constructor
	 * 
	 * @param inputStream
	 *            input stream
	 * @since 2.0.3
	 */
	public ByteReader(InputStream inputStream) {
		this.bytes = null;
		this.inputStream = inputStream;
	}

	/**
	 * Constructor
	 * 
	 * @param bytes
	 *            bytes
	 * @param byteOrder
	 *            byte order
	 * @since 2.0.3
	 */
	public ByteReader(byte[] bytes, ByteOrder byteOrder) {
		this(bytes);
		this.byteOrder = byteOrder;
	}

	/**
	 * Constructor
	 * 
	 * @param inputStream
	 *            input stream
	 * @param byteOrder
	 *            byte order
	 * @since 2.0.3
	 */
	public ByteReader(InputStream inputStream, ByteOrder byteOrder) {
		this(inputStream);
		this.byteOrder = byteOrder;
	}

	/**
	 * Get the bytes
	 * 
	 * @return bytes
	 * 
	 * @since 2.0.3
	 */
	public byte[] getBytes() {
		return bytes;
	}

	/**
	 * Get the output stream
	 * 
	 * @return output stream
	 * 
	 * @since 2.0.3
	 */
	public InputStream getInputStream() {
		return inputStream;
	}

	/**
	 * Close the byte reader
	 * 
	 * @since 2.0.3
	 */
	public void close() {
		try {
			inputStream.close();
		} catch (IOException e) {
			logger.log(Level.WARNING,
					"Failed to close byte reader intput stream", e);
		}
	}

	/**
	 * Get the next byte to be read
	 * 
	 * @return next byte to be read
	 */
	public int getNextByte() {
		return nextByte;
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
	 * Read a String from the provided number of bytes
	 * 
	 * @param num
	 *            number of bytes
	 * @return String
	 * @throws IOException
	 *             upon error
	 */
	public String readString(int num) throws IOException {
		return new String(readBytes(num), CHAR_SET);
	}

	/**
	 * Read the number of bytes
	 * 
	 * @param num
	 *            number of bytes
	 * @return bytes
	 * @throws IOException
	 *             upon error
	 */
	public byte[] readBytes(int num) throws IOException {
		verifyRemainingBytes(num);
		byte[] bytes = new byte[num];
		inputStream.read(bytes);
		nextByte += num;
		return bytes;
	}

	/**
	 * Read a byte
	 * 
	 * @return byte
	 * @throws IOException
	 *             upon error
	 */
	public byte readByte() throws IOException {
		return readBytes(1)[0];
	}

	/**
	 * Read an integer
	 * 
	 * @return integer
	 * @throws IOException
	 *             upon error
	 */
	public int readInt() throws IOException {
		return ByteBuffer.wrap(readBytes(4)).order(byteOrder).getInt();
	}

	/**
	 * Read an unsigned integer
	 * 
	 * @return unsigned integer
	 * @throws IOException
	 *             upon error
	 */
	public long readUnsignedInt() throws IOException {
		int intValue = readInt();
		long value = intValue & 0xffffffffL;
		return value;
	}

	/**
	 * Read a double
	 * 
	 * @return double
	 * @throws IOException
	 *             upon error
	 */
	public double readDouble() throws IOException {
		return ByteBuffer.wrap(readBytes(8)).order(byteOrder).getDouble();
	}

	/**
	 * Verify with the remaining bytes that there are enough remaining to read
	 * the provided amount
	 * 
	 * @param bytesToRead
	 *            number of bytes to read
	 */
	private void verifyRemainingBytes(int bytesToRead) {
		if (bytes != null && nextByte + bytesToRead > bytes.length) {
			throw new SFException(
					"No more remaining bytes to read. Total Bytes: "
							+ bytes.length + ", Bytes already read: " + nextByte
							+ ", Attempted to read: " + bytesToRead);
		}
	}

}
