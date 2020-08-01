package net.finmath.time.businessdaycalendar;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

/**
 * A business day calendar, where every day is a business day, except for weekends and London holidays
 *
 * @author Niklas Rodi
 * @version 1.0
 */
public class BusinessdayCalendarExcludingLONHolidays extends BusinessdayCalendarExcludingGivenHolidays {

	/**
	 *
	 */
	private static final long serialVersionUID = 7451923974528392081L;
	/*
	 * Details of this calendar.
	 * If you like to create a similar calendar, just duplicate this class and
	 * modify the following two lines.
	 */
	private static final String NAME = "London";

	private static final String[] HOLIDAYLISTASSTRINGS = new String[] { "03/01/2000", "21/04/2000", "24/04/2000",
			"01/05/2000", "29/05/2000", "28/08/2000", "25/12/2000", "26/12/2000", "01/01/2001", "13/04/2001",
			"16/04/2001", "07/05/2001", "28/05/2001", "27/08/2001", "25/12/2001", "26/12/2001", "01/01/2002",
			"29/03/2002", "01/04/2002", "06/05/2002", "03/06/2002", "04/06/2002", "26/08/2002", "25/12/2002",
			"26/12/2002", "01/01/2003", "18/04/2003", "21/04/2003", "05/05/2003", "26/05/2003", "25/08/2003",
			"25/12/2003", "26/12/2003", "01/01/2004", "09/04/2004", "12/04/2004", "03/05/2004", "31/05/2004",
			"30/08/2004", "27/12/2004", "28/12/2004", "03/01/2005", "25/03/2005", "28/03/2005", "02/05/2005",
			"30/05/2005", "29/08/2005", "26/12/2005", "27/12/2005", "02/01/2006", "14/04/2006", "17/04/2006",
			"01/05/2006", "29/05/2006", "28/08/2006", "25/12/2006", "26/12/2006", "01/01/2007", "06/04/2007",
			"09/04/2007", "07/05/2007", "28/05/2007", "27/08/2007", "25/12/2007", "26/12/2007", "01/01/2008",
			"21/03/2008", "24/03/2008", "05/05/2008", "26/05/2008", "25/08/2008", "25/12/2008", "26/12/2008",
			"01/01/2009", "10/04/2009", "13/04/2009", "04/05/2009", "25/05/2009", "31/08/2009", "25/12/2009",
			"28/12/2009", "01/01/2010", "02/04/2010", "05/04/2010", "03/05/2010", "31/05/2010", "30/08/2010",
			"27/12/2010", "28/12/2010", "03/01/2011", "22/04/2011", "25/04/2011", "29/04/2011", "02/05/2011",
			"30/05/2011", "29/08/2011", "26/12/2011", "27/12/2011", "02/01/2012", "06/04/2012", "09/04/2012",
			"07/05/2012", "04/06/2012", "05/06/2012", "27/08/2012", "25/12/2012", "26/12/2012", "01/01/2013",
			"29/03/2013", "01/04/2013", "06/05/2013", "27/05/2013", "26/08/2013", "25/12/2013", "26/12/2013",
			"01/01/2014", "18/04/2014", "21/04/2014", "05/05/2014", "26/05/2014", "25/08/2014", "25/12/2014",
			"26/12/2014", "01/01/2015", "03/04/2015", "06/04/2015", "04/05/2015", "25/05/2015", "31/08/2015",
			"25/12/2015", "28/12/2015", "01/01/2016", "25/03/2016", "28/03/2016", "02/05/2016", "30/05/2016",
			"29/08/2016", "26/12/2016", "27/12/2016", "02/01/2017", "14/04/2017", "17/04/2017", "01/05/2017",
			"29/05/2017", "28/08/2017", "25/12/2017", "26/12/2017", "01/01/2018", "30/03/2018", "02/04/2018",
			"07/05/2018", "28/05/2018", "27/08/2018", "25/12/2018", "26/12/2018", "01/01/2019", "19/04/2019",
			"22/04/2019", "06/05/2019", "27/05/2019", "26/08/2019", "25/12/2019", "26/12/2019", "01/01/2020",
			"10/04/2020", "13/04/2020", "04/05/2020", "25/05/2020", "31/08/2020", "25/12/2020", "28/12/2020",
			"01/01/2021", "02/04/2021", "05/04/2021", "03/05/2021", "31/05/2021", "30/08/2021", "27/12/2021",
			"28/12/2021", "03/01/2022", "15/04/2022", "18/04/2022", "02/05/2022", "30/05/2022", "29/08/2022",
			"26/12/2022", "27/12/2022", "02/01/2023", "07/04/2023", "10/04/2023", "01/05/2023", "29/05/2023",
			"28/08/2023", "25/12/2023", "26/12/2023", "01/01/2024", "29/03/2024", "01/04/2024", "06/05/2024",
			"27/05/2024", "26/08/2024", "25/12/2024", "26/12/2024", "01/01/2025", "18/04/2025", "21/04/2025",
			"05/05/2025", "26/05/2025", "25/08/2025", "25/12/2025", "26/12/2025", "01/01/2026", "03/04/2026",
			"06/04/2026", "04/05/2026", "25/05/2026", "31/08/2026", "25/12/2026", "28/12/2026", "01/01/2027",
			"26/03/2027", "29/03/2027", "03/05/2027", "31/05/2027", "30/08/2027", "27/12/2027", "28/12/2027",
			"03/01/2028", "14/04/2028", "17/04/2028", "01/05/2028", "29/05/2028", "28/08/2028", "25/12/2028",
			"26/12/2028", "01/01/2029", "30/03/2029", "02/04/2029", "07/05/2029", "28/05/2029", "27/08/2029",
			"25/12/2029", "26/12/2029", "01/01/2030", "19/04/2030", "22/04/2030", "06/05/2030", "27/05/2030",
			"26/08/2030", "25/12/2030", "26/12/2030", "01/01/2031", "11/04/2031", "14/04/2031", "05/05/2031",
			"26/05/2031", "25/08/2031", "25/12/2031", "26/12/2031", "01/01/2032", "26/03/2032", "29/03/2032",
			"03/05/2032", "31/05/2032", "30/08/2032", "27/12/2032", "28/12/2032", "03/01/2033", "15/04/2033",
			"18/04/2033", "02/05/2033", "30/05/2033", "29/08/2033", "26/12/2033", "27/12/2033", "02/01/2034",
			"07/04/2034", "10/04/2034", "01/05/2034", "29/05/2034", "28/08/2034", "25/12/2034", "26/12/2034",
			"01/01/2035", "23/03/2035", "26/03/2035", "07/05/2035", "28/05/2035", "27/08/2035", "25/12/2035",
			"26/12/2035", "01/01/2036", "11/04/2036", "14/04/2036", "05/05/2036", "26/05/2036", "25/08/2036",
			"25/12/2036", "26/12/2036", "01/01/2037", "03/04/2037", "06/04/2037", "04/05/2037", "25/05/2037",
			"31/08/2037", "25/12/2037", "28/12/2037", "01/01/2038", "23/04/2038", "26/04/2038", "03/05/2038",
			"31/05/2038", "30/08/2038", "27/12/2038", "28/12/2038", "03/01/2039", "08/04/2039", "11/04/2039",
			"02/05/2039", "30/05/2039", "29/08/2039", "26/12/2039", "27/12/2039", "02/01/2040", "30/03/2040",
			"02/04/2040", "07/05/2040", "28/05/2040", "27/08/2040", "25/12/2040", "26/12/2040", "01/01/2041",
			"19/04/2041", "22/04/2041", "06/05/2041", "27/05/2041", "26/08/2041", "25/12/2041", "26/12/2041",
			"01/01/2042", "04/04/2042", "07/04/2042", "05/05/2042", "26/05/2042", "25/08/2042", "25/12/2042",
			"26/12/2042", "01/01/2043", "27/03/2043", "30/03/2043", "04/05/2043", "25/05/2043", "31/08/2043",
			"25/12/2043", "28/12/2043", "01/01/2044", "15/04/2044", "18/04/2044", "02/05/2044", "30/05/2044",
			"29/08/2044", "26/12/2044", "27/12/2044", "02/01/2045", "07/04/2045", "10/04/2045", "01/05/2045",
			"29/05/2045", "28/08/2045", "25/12/2045", "26/12/2045", "01/01/2046", "23/03/2046", "26/03/2046",
			"07/05/2046", "28/05/2046", "27/08/2046", "25/12/2046", "26/12/2046", "01/01/2047", "12/04/2047",
			"15/04/2047", "06/05/2047", "27/05/2047", "26/08/2047", "25/12/2047", "26/12/2047", "01/01/2048",
			"03/04/2048", "06/04/2048", "04/05/2048", "25/05/2048", "31/08/2048", "25/12/2048", "28/12/2048",
			"01/01/2049", "16/04/2049", "19/04/2049", "03/05/2049", "31/05/2049", "30/08/2049", "27/12/2049",
			"28/12/2049", "03/01/2050", "08/04/2050", "11/04/2050", "02/05/2050", "30/05/2050", "29/08/2050",
			"26/12/2050", "27/12/2050", "02/01/2051", "31/03/2051", "03/04/2051", "01/05/2051", "29/05/2051",
			"28/08/2051", "25/12/2051", "26/12/2051", "01/01/2052", "19/04/2052", "22/04/2052", "06/05/2052",
			"27/05/2052", "26/08/2052", "25/12/2052", "26/12/2052", "01/01/2053", "04/04/2053", "07/04/2053",
			"05/05/2053", "26/05/2053", "25/08/2053", "25/12/2053", "26/12/2053", "01/01/2054", "27/03/2054",
			"30/03/2054", "04/05/2054", "25/05/2054", "31/08/2054", "25/12/2054", "28/12/2054", "01/01/2055",
			"16/04/2055", "19/04/2055", "03/05/2055", "31/05/2055", "30/08/2055", "27/12/2055", "28/12/2055",
			"03/01/2056", "31/03/2056", "03/04/2056", "01/05/2056", "29/05/2056", "28/08/2056", "25/12/2056",
			"26/12/2056", "01/01/2057", "20/04/2057", "23/04/2057", "07/05/2057", "28/05/2057", "27/08/2057",
			"25/12/2057", "26/12/2057", "01/01/2058", "12/04/2058", "15/04/2058", "06/05/2058", "27/05/2058",
			"26/08/2058", "25/12/2058", "26/12/2058", "01/01/2059", "28/03/2059", "31/03/2059", "05/05/2059",
			"26/05/2059", "25/08/2059", "25/12/2059", "26/12/2059", "01/01/2060", "16/04/2060", "19/04/2060",
			"03/05/2060", "31/05/2060", "30/08/2060", "27/12/2060", "28/12/2060" };

	/*
	 * Static initializer for holidays set from dd/MM/yyyy spec
	 */
	private static final Set<LocalDate> HOLIDAYS;
	static {
		final DateTimeFormatter formatter = DateTimeFormatter.ofPattern("dd/MM/yyyy");
		final Set<LocalDate> holidaysSet = new TreeSet<>();
		for(final String holidayAsString : HOLIDAYLISTASSTRINGS) {
			holidaysSet.add(LocalDate.parse(holidayAsString, formatter));
		}
		HOLIDAYS = Collections.unmodifiableSet(holidaysSet);
	}

	/**
	 * Create LONDON business day calendar.
	 */
	public BusinessdayCalendarExcludingLONHolidays() {
		this(null);
	}

	/**
	 * Create LONDON business day calendar using a given business day calendar as basis.
	 *
	 * @param baseCalendar Calendar of business days.
	 */
	public BusinessdayCalendarExcludingLONHolidays(final BusinessdayCalendar baseCalendar) {
		super(NAME, baseCalendar, true);
	}

	@Override
	public Set<LocalDate> getHolidays() { return HOLIDAYS; }
}
