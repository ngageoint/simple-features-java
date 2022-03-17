package mil.nga.sf.util;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Read through text string
 * 
 * @author osbornb
 * @since 2.0.3
 */
public class TextReader {

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(TextReader.class.getName());

	/**
	 * Reader
	 */
	private final Reader reader;

	/**
	 * Next token cache for peeks
	 */
	private String nextToken;

	/**
	 * Next character number cache for between token caching
	 */
	private Integer nextCharacterNum;

	/**
	 * Constructor
	 * 
	 * @param text
	 *            text
	 */
	public TextReader(String text) {
		this(new StringReader(text));
	}

	/**
	 * Constructor
	 * 
	 * @param reader
	 *            reader
	 */
	public TextReader(Reader reader) {
		this.reader = reader;
	}

	/**
	 * Get the reader
	 * 
	 * @return reader
	 */
	public Reader getReader() {
		return reader;
	}

	/**
	 * Close the text reader
	 */
	public void close() {
		try {
			reader.close();
		} catch (IOException e) {
			logger.log(Level.WARNING, "Failed to close text reader", e);
		}
	}

	/**
	 * Read the next token. Ignores whitespace until a non whitespace character
	 * is encountered. Returns a contiguous block of token characters ( [a-z] |
	 * [A-Z] | [0-9] | - | . | + ) or a non whitespace single character.
	 * 
	 * @return token
	 * @throws IOException
	 *             upon read error
	 */
	public String readToken() throws IOException {

		String token = null;

		// Get the next token, cached or read
		if (nextToken != null) {
			token = nextToken;
			nextToken = null;
		} else {

			StringBuilder builder = null;

			// Get the next character, cached or read
			int characterNum;
			if (nextCharacterNum != null) {
				characterNum = nextCharacterNum;
				nextCharacterNum = null;
			} else {
				characterNum = reader.read();
			}

			// Continue while characters are left
			while (characterNum != -1) {

				char character = (char) characterNum;

				// Check if not the first character in the token
				if (builder != null) {

					// Append token characters
					if (isTokenCharacter(character)) {
						builder.append(character);
					} else {
						// Complete the token before this character and cache
						// the character
						if (!isWhitespace(character)) {
							nextCharacterNum = characterNum;
						}
						break;
					}

				} else if (!isWhitespace(character)) {

					// First non whitespace character in the token
					builder = new StringBuilder();
					builder.append(character);

					// Complete token if a single character token
					if (!isTokenCharacter(character)) {
						break;
					}

				}

				// Read the next character
				characterNum = reader.read();
			}

			if (builder != null) {
				token = builder.toString();
			}

		}
		return token;
	}

	/**
	 * Peek at the next token without reading past it
	 * 
	 * @return next token
	 * @throws IOException
	 *             upon read error
	 */
	public String peekToken() throws IOException {
		if (nextToken == null) {
			nextToken = readToken();
		}
		return nextToken;
	}

	/**
	 * Read a double
	 * 
	 * @return double
	 * @throws IOException
	 *             upon read error
	 */
	public double readDouble() throws IOException {
		String token = readToken();
		if (token == null) {
			throw new SFException("Failed to read expected double value");
		}
		return Double.parseDouble(token);
	}

	/**
	 * Check if the character is a contiguous block token character: ( [a-z] |
	 * [A-Z] | [0-9] | - | . | + )
	 * 
	 * @param c
	 * @return
	 */
	private static boolean isTokenCharacter(char c) {
		return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
				|| (c >= '0' && c <= '9') || c == '-' || c == '.' || c == '+';
	}

	/**
	 * Check if the character is whitespace or a space character
	 * 
	 * @param c
	 *            character
	 * @return true if whitespace
	 */
	private static boolean isWhitespace(char c) {
		return Character.isWhitespace(c) || Character.isSpaceChar(c);
	}

}
