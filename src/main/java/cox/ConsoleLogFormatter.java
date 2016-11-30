package cox;

import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.stream.Collectors;

public class ConsoleLogFormatter extends Formatter {

	@Override
	public String format(LogRecord record) {
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss", Locale.getDefault());
		Date date = new Date(record.getMillis());
		String d = sdf.format(date);
		String level = record.getLevel().getName();
		String msg = formatMessage(record);
		String method = record.getSourceMethodName();
		String clas = record.getSourceClassName();

		String[] details = { d, level, clas, method, msg };
		return Arrays.stream(details).collect(Collectors.joining(" | ")) + System.lineSeparator();
	}

    public static String center(String s, int size, char pad) {
        if (s == null || size <= s.length())
            return s;

        StringBuilder sb = new StringBuilder(size);
        for (int i = 0; i < (size - s.length()) / 2; i++) {
            sb.append(pad);
        }
        sb.append(s);
        while (sb.length() < size) {
            sb.append(pad);
        }
        return sb.toString();
    }
}
