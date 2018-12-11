package mil.nga.sf.util;

/**
 * Simple Features exception
 * 
 * @author osbornb
 */
public class SFException extends RuntimeException {

	/**
	 * Serial version id
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Constructor
	 */
	public SFException() {
		super();
	}

	/**
	 * Constructor
	 * 
	 * @param message
	 *            error message
	 */
	public SFException(String message) {
		super(message);
	}

	/**
	 * Constructor
	 * 
	 * @param message
	 *            error message
	 * @param throwable
	 *            throwable
	 */
	public SFException(String message, Throwable throwable) {
		super(message, throwable);
	}

	/**
	 * Constructor
	 * 
	 * @param throwable
	 *            throwable
	 */
	public SFException(Throwable throwable) {
		super(throwable);
	}

}
