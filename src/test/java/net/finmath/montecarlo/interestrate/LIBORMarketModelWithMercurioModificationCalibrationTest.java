/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 16.01.2015
 */
package net.finmath.montecarlo.interestrate;

import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.Function;

import org.junit.Assert;
import org.junit.Test;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.calibration.ParameterObject;
import net.finmath.marketdata.calibration.Solver;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelFromCurvesAndVols;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveFromDiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionView;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModelWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.BlendedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.VolatilityReductionMercurioModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecayWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelPiecewiseConstantWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.Caplet.ValueUnit;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.interestrate.products.CapletOnBackwardLookingRate;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.OptimizerFactory;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.Schedule;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 */
public class LIBORMarketModelWithMercurioModificationCalibrationTest {

	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterParam		= new DecimalFormat(" #0.000;-#0.000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	private final RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	private CalibrationProduct createCalibrationItem(double weight, double exerciseDate, double swapPeriodLength, int numberOfPeriods, double moneyness, double targetVolatility, String targetVolatilityType, ForwardCurve forwardCurve, DiscountCurve discountCurve) throws CalculationException { 

		final double[]	fixingDates			= new double[numberOfPeriods];
		final double[]	paymentDates		= new double[numberOfPeriods];
		final double[]	swapTenor			= new double[numberOfPeriods + 1];

		for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
			fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
			paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex + 1) * swapPeriodLength;
			swapTenor[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
		}
		swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;

		// Swaptions swap rate
		final double swaprate = moneyness + getParSwaprate(forwardCurve, discountCurve, swapTenor);

		// Set swap rates for each period
		final double[] swaprates = new double[numberOfPeriods];
		Arrays.fill(swaprates, swaprate);

		/*
		 * We use Monte-Carlo calibration on implied volatility.
		 * Alternatively you may change here to Monte-Carlo valuation on price or
		 * use an analytic approximation formula, etc.
		 */
		final SwaptionSimple swaptionMonteCarlo = new SwaptionSimple(swaprate, swapTenor, SwaptionSimple.ValueUnit.valueOf(targetVolatilityType));
		//		double targetValuePrice = AnalyticFormulas.blackModelSwaptionValue(swaprate, targetVolatility, fixingDates[0], swaprate, getSwapAnnuity(discountCurve, swapTenor));
		return new CalibrationProduct(swaptionMonteCarlo, targetVolatility, weight);
	}

	
//	testSwaptionSmileCalibration()  <-- REMOVED

	// 						PART 2 
	/**
	 * Brute force Monte-Carlo calibration of swaptions.
	 *
	 * The test also performs a test on the serialization of the LMM. It serialized the calibrated model into a byte array,
	 * reads the model back and compares a simulation using the serialized model with the original one.
	 *
	 * @throws CalculationException Thrown if the model fails to calibrate.
	 * @throws SolverException Thrown if the solver fails to find a solution.
	 */
	@Test
	public void testATMSwaptionCalibration() throws CalculationException, SolverException {

		final int numberOfPaths		= 1000; // ideale è 10'000 paths!
		final int numberOfFactors	= 1;

		final long millisCurvesStart = System.currentTimeMillis();

		/*
		 * Calibration test
		 */
		System.out.println("Calibration to Swaptions.\n");

		/*
		 * Calibration of rate curves
		 */
		System.out.println("Calibration of rate curves:");

		final AnalyticModel curveModel = getCalibratedCurve();

		// Create the forward curve (initial value of the LIBOR market model)
		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,1D)");
		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");
		//		curveModel.addCurve(discountCurve.getName(), discountCurve);

		final long millisCurvesEnd = System.currentTimeMillis();
		System.out.println();

		/*
		 * Calibration of model volatilities
		 */
		System.out.println("Brute force Monte-Carlo calibration of model volatilities:");

		/*
		 * Create a set of calibration products.
		 */
		final ArrayList<String>				calibrationItemNames	= new ArrayList<>();
		final ArrayList<CalibrationProduct>	calibrationProducts		= new ArrayList<>();

		final double	swapPeriodLength	= 0.5;

		final String[] atmExpiries = {
				"1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M", "1M",
//				"2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", "2M", 
				"3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M", "3M",
				"6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M", "6M",
//				"9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M", "9M",
				"1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", "1Y", 
//				"18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M", "18M",
				"2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y", "2Y",
				"3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "3Y", "4Y",
				"4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "4Y", "5Y", "5Y", "5Y", "5Y",
				"5Y", "5Y", "5Y", "5Y", "5Y", "5Y", "5Y", "5Y", "5Y", "5Y", "7Y", "7Y", "7Y", "7Y", "7Y", "7Y", "7Y",
				"7Y", "7Y", "7Y", "7Y", "7Y", "7Y", "7Y", "10Y", "10Y", "10Y", "10Y", "10Y", "10Y", "10Y", "10Y", "10Y",
				"10Y", "10Y", "10Y", "10Y", "10Y", "15Y", "15Y", "15Y", "15Y", "15Y", "15Y", "15Y", "15Y", "15Y", "15Y",
				"15Y", "15Y", "15Y", "15Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y", "20Y",
				"20Y", "20Y", "20Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y", "25Y",
				"25Y", "25Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y", "30Y" };

		final String[] atmTenors = {
				"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", 
				
//				"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y",
//				"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y",
//				"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y",

				"1Y", "2Y","3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y",
				"1Y", "2Y", "3Y", "4Y","5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y",
				"7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y",
				"9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y",
				"15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y",
				"25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y",
				"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y",
				"3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y",
				"5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y",
				"7Y", "8Y", "9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y",
				"9Y", "10Y", "15Y", "20Y", "25Y", "30Y", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y",
				"15Y", "20Y", "25Y", "30Y" };

		final double[] atmNormalVolatilities = {
				0.0015335, 0.0015179, 0.0019499, 0.0024161, 0.0027817, 0.0031067, 0.0033722, 0.0035158, 0.0036656, 0.0037844, 0.00452, 0.0050913, 0.0054071, 0.0056496,
				//next is 2M
//				0.0016709, 0.0016287, 0.0020182, 0.0024951, 0.002827, 0.0031023, 0.0034348, 0.0036183, 0.0038008, 0.0039155, 0.0046602, 0.0051981, 0.0055116, 0.0057249,
				
				0.0015543, 0.0016509, 0.0020863, 0.002587, 0.002949, 0.0032105, 0.0035338, 0.0037133, 0.0038475, 0.0040674, 0.0047458, 0.005276, 0.005476, 0.005793,
				0.0016777, 0.001937, 0.0023423, 0.0027823, 0.0031476, 0.0034569, 0.0037466, 0.0039852, 0.0041802, 0.0043221, 0.0049649, 0.0054206, 0.0057009, 0.0059071,
				//next is 9M
//				0.0017809, 0.0020951, 0.0024978, 0.0029226, 0.0032379, 0.0035522, 0.0038397, 0.0040864, 0.0043122, 0.0044836, 0.0050939, 0.0054761, 0.0057374, 0.0059448,
				
				0.0020129, 0.0022865, 0.0027082, 0.0030921, 0.0033849, 0.0037107, 0.0039782, 0.0042058, 0.0044272, 0.0046082, 0.0051564, 0.0055307, 0.0057924, 0.0059811,
				//next is 18M
//				0.0022824, 0.0025971, 0.0029895, 0.0033299, 0.0036346, 0.0039337, 0.0042153, 0.0044347, 0.0046686, 0.0048244, 0.0052739, 0.005604, 0.0058311, 0.0060011,
			
				0.0026477, 0.0029709, 0.0033639, 0.0036507, 0.0039096, 0.0041553, 0.0044241, 0.00462, 0.0048265, 0.004989, 0.005361, 0.0056565, 0.0058529, 0.0060102,
				0.003382, 0.0036593, 0.0039353, 0.0041484, 0.0043526, 0.0045677, 0.004775, 0.0049506, 0.0051159, 0.0052722, 0.0055185, 0.0057089, 0.0058555, 0.0059432,
				0.0040679, 0.0042363, 0.0044602, 0.0046206, 0.0047527, 0.0048998, 0.0050513, 0.0051928, 0.0053439, 0.0054657, 0.0056016, 0.0057244, 0.0058153, 0.0058793,
				0.0045508, 0.0046174, 0.0047712, 0.0048999, 0.0050364, 0.0051504, 0.0052623, 0.0053821, 0.0054941, 0.0055918, 0.0056569, 0.0057283, 0.0057752, 0.0058109,
				0.0051385, 0.0051373, 0.0052236, 0.005312, 0.0053793, 0.0054396, 0.0055037, 0.0055537, 0.0056213, 0.0056943, 0.005671, 0.0056707, 0.0056468, 0.0056423,
				0.0055069, 0.0054836, 0.0055329, 0.0055696, 0.005605, 0.0056229, 0.0056562, 0.005655, 0.0056679, 0.0057382, 0.0056494, 0.0055831, 0.0055096, 0.0054526,
				0.0054486, 0.0054057, 0.0054439, 0.005462, 0.0054915, 0.0054993, 0.0055134, 0.0054985, 0.0055318, 0.0055596, 0.005369, 0.0052513, 0.0051405, 0.0050416,
				0.005317, 0.005268, 0.005312, 0.0053112, 0.0053417, 0.0053556, 0.0053323, 0.0053251, 0.0053233, 0.0053126, 0.0050827, 0.004922, 0.0047924, 0.0046666,
				0.0051198, 0.0051013, 0.0051421, 0.0051418, 0.0051538, 0.005133, 0.0051081, 0.0050552, 0.005055, 0.0050473, 0.0048161, 0.0045965, 0.0044512, 0.0043099,
				0.0049482, 0.004947, 0.0049805, 0.0049951, 0.0050215, 0.0049849, 0.0049111, 0.0048498, 0.0047879, 0.0047688, 0.0044943, 0.0042786, 0.0041191, 0.0039756};

		final LocalDate referenceDate = LocalDate.of(2020, Month.JULY, 31);
		final BusinessdayCalendarExcludingTARGETHolidays cal = new BusinessdayCalendarExcludingTARGETHolidays();
		final DayCountConvention_ACT_365 modelDC = new DayCountConvention_ACT_365();
		for(int i=0; i<atmNormalVolatilities.length; i++ ) {

			final LocalDate exerciseDate = cal.getDateFromDateAndOffsetCode(referenceDate, atmExpiries[i]);
			final LocalDate tenorEndDate = cal.getDateFromDateAndOffsetCode(exerciseDate, atmTenors[i]);
			double	exercise		= modelDC.getDaycountFraction(referenceDate, exerciseDate);
			double	tenor			= modelDC.getDaycountFraction(exerciseDate, tenorEndDate);

			// We consider an idealized tenor grid (alternative: adapt the model grid)
			exercise	= Math.round(exercise/0.25)*0.25;
			tenor		= Math.round(tenor/0.25)*0.25;

			if(exercise < 1.0) {
				continue;
			}

			final int numberOfPeriods = (int)Math.round(tenor / swapPeriodLength);

			final double	moneyness			= 0.0;
			final double	targetVolatility	= atmNormalVolatilities[i];

			final String	targetVolatilityType = "VOLATILITYNORMAL";

			final double	weight = 1.0;

			calibrationProducts.add(createCalibrationItem(weight, exercise, swapPeriodLength, numberOfPeriods, moneyness, targetVolatility, targetVolatilityType, forwardCurve, discountCurve));
			calibrationItemNames.add(atmExpiries[i]+"\t"+atmTenors[i]);
		}

		/*
		 * Create a simulation time discretization
		 */
		// If simulation time is below libor time, exceptions will be hard to track.
		final double lastTime	= 21.0;
		final double dtLibor	= 0.5;
		final double dt	= 0.25;
		final TimeDiscretization timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dt), dt);
		final TimeDiscretization liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dtLibor), dtLibor);
		
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, 31415 /* seed */);

		// needed if you use LIBORVolatilityModelPiecewiseConstantWithMercurioModification: TimeDiscretization optionMaturityDiscretization = new TimeDiscretizationFromArray(0.0, 0.25, 0.50, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0);
		TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9,00, 10.0, 12.5, 15.0, 21.0);
		double[] arrayValues = new double [timeToMaturityDiscretization.getNumberOfTimes()];
		for (int i=0; i<timeToMaturityDiscretization.getNumberOfTimes(); i++) {arrayValues[i]= 0.2/100;}

		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, timeToMaturityDiscretization, arrayValues);
		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, 0.20, 0.05, 0.10, 0.20, true);
		//final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization,optionMaturityDiscretization,timeToMaturityDiscretization, 0.50 / 100);
		
		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecayWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);
		// Create a covariance model
		//AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, new double[] { 0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 0.10} );
		final AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter (0=lognormal, > 1 = almost a normal model).
		final AbstractLIBORCovarianceModelParametric covarianceModelDisplaced = new DisplacedLocalVolatilityModel(covarianceModelParametric, 1.0/0.25, false /* isCalibrateable */);
		final AbstractLIBORCovarianceModelParametric covarianceModelReducedVolatility = new VolatilityReductionMercurioModel(covarianceModelDisplaced);
		
		// Set model properties
		final Map<String, Object> properties = new HashMap<>();
	
		System.out.println("Number of volatility parameters: " + volatilityModel.getParameter().length);

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModelFromCovarianceModelWithMercurioModification.Measure.SPOT.name());

		// Choose normal state space for the Euler scheme (the covariance model above carries a linear local volatility model, such that the resulting model is log-normal).
		properties.put("stateSpace", LIBORMarketModelFromCovarianceModelWithMercurioModification.StateSpace.NORMAL.name());

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).
		final Double accuracy = new Double(1E-8);	// Lower accuracy to reduce runtime of the unit test
		final int maxIterations = 200;
		final int numberOfThreads = 1;
		final OptimizerFactory optimizerFactory = new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		//penso che getParameterAsDouble().length mi ritorni il numero di tutti i paratemtri da calibrare cioè tutta la piecewise volatiliy e anche correlazion
		final double[] parameterStandardDeviation = new double[covarianceModelParametric.getParameterAsDouble().length];
		final double[] parameterLowerBound = new double[covarianceModelParametric.getParameterAsDouble().length];
		final double[] parameterUpperBound = new double[covarianceModelParametric.getParameterAsDouble().length];
		Arrays.fill(parameterStandardDeviation, 0.20/100.0);
		Arrays.fill(parameterLowerBound, 0.0);
		Arrays.fill(parameterUpperBound, Double.POSITIVE_INFINITY);

		//		optimizerFactory = new OptimizerFactoryCMAES(accuracy, maxIterations, parameterLowerBound, parameterUpperBound, parameterStandardDeviation);

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).
		final Map<String, Object> calibrationParameters = new HashMap<>();
		calibrationParameters.put("accuracy", accuracy);
		calibrationParameters.put("brownianMotion", brownianMotion);
		calibrationParameters.put("optimizerFactory", optimizerFactory);
		calibrationParameters.put("parameterStep", new Double(1E-4));
		properties.put("calibrationParameters", calibrationParameters);

		final long millisCalibrationStart = System.currentTimeMillis();

		/*
		 * Create corresponding LIBOR Market Model
		 */
		final CalibrationProduct[] calibrationItemsLMM = new CalibrationProduct[calibrationItemNames.size()];
		for(int i=0; i<calibrationItemNames.size(); i++) {
			calibrationItemsLMM[i] = new CalibrationProduct(calibrationProducts.get(i).getProduct(),calibrationProducts.get(i).getTargetValue(),calibrationProducts.get(i).getWeight());
		}
		final LIBORMarketModel mercurioModelCalibrated = LIBORMarketModelFromCovarianceModelWithMercurioModification.of(
				liborPeriodDiscretization,
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				randomVariableFactory,
				covarianceModelReducedVolatility,
				calibrationItemsLMM, properties);

		final long millisCalibrationEnd = System.currentTimeMillis();
		
//-------------------------------------------------------------------------------- fine calibrazione volatility------------------------------
		System.out.println("\nCalibrated parameters are:");
		final double[] param = ((AbstractLIBORCovarianceModelParametric)((LIBORMarketModelFromCovarianceModelWithMercurioModification) mercurioModelCalibrated).getCovarianceModel()).getParameterAsDouble();
		for (final double p : param) {
			System.out.println(p);
		}

		final EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(brownianMotion);
		final LIBORModelMonteCarloSimulationModel simulationCalibrated = new LIBORMonteCarloSimulationFromLIBORModel(mercurioModelCalibrated, process);

		System.out.println("\nValuation on calibrated model:");
		double deviationSum			= 0.0;
		double deviationSquaredSum	= 0.0;
		for (int i = 0; i < calibrationProducts.size(); i++) {
			final AbstractLIBORMonteCarloProduct calibrationProduct = calibrationProducts.get(i).getProduct();
			try {
				final double valueModel = calibrationProduct.getValue(simulationCalibrated);
				final double valueTarget = calibrationProducts.get(i).getTargetValue().getAverage();
				final double error = valueModel-valueTarget;
				deviationSum += error;
				deviationSquaredSum += error*error;
				System.out.println(calibrationItemNames.get(i) + "\t" + "Model: " + formatterValue.format(valueModel) + "\t Target: " + formatterValue.format(valueTarget) + "\t Deviation: " + formatterDeviation.format(valueModel-valueTarget));// + "\t" + calibrationProduct.toString());
			}
			catch(final Exception e) {
			}
		}


		/*
		 * Checking serilization
		 */
		byte[] lmmSerialized = null;
		try {
			final ByteArrayOutputStream baos = new ByteArrayOutputStream();
			final ObjectOutputStream oos = new ObjectOutputStream( baos );
			oos.writeObject(mercurioModelCalibrated.getCloneWithModifiedData(null));
			lmmSerialized = baos.toByteArray();
		} catch (final IOException e) {
			fail("Serialization failed with exception " + e.getMessage());
		}

		LIBORMarketModelFromCovarianceModelWithMercurioModification liborMarketModelFromSerialization = null;
		try {
			final ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(lmmSerialized) );
			liborMarketModelFromSerialization = (LIBORMarketModelFromCovarianceModelWithMercurioModification)ois.readObject();
		} catch (IOException | ClassNotFoundException e) {
			fail("Deserialization failed with exception " + e.getMessage());
		}

		/*
		 * Check if the deserialized model and the original calibrated model give the same valuations
		 */
		if(liborMarketModelFromSerialization != null) {
			final LIBORModelMonteCarloSimulationModel simulationFromSerialization = new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModelFromSerialization, new EulerSchemeFromProcessModel(brownianMotion));

			System.out.println("\nValuation on calibrated model:");
			for (int i = 0; i < calibrationProducts.size(); i++) {
				final AbstractLIBORMonteCarloProduct calibrationProduct = calibrationProducts.get(i).getProduct();
				try {
					final double valueFromCalibratedModel = calibrationProduct.getValue(simulationCalibrated);
					final double valueFromSerializedModel = calibrationProduct.getValue(simulationFromSerialization);
					final double error = valueFromSerializedModel-valueFromCalibratedModel;
					Assert.assertEquals("Valuation using deserilized model.", valueFromCalibratedModel, valueFromSerializedModel, 1E-12);
					System.out.println(calibrationItemNames.get(i) + "\t" + formatterDeviation.format(error));
				}
				catch(final Exception e) {
				}
			}
		}
		
		// CAPLET ON BACKWARD LOOKING RATE SEMESTRALI
		DecimalFormat formatterTimeValue = new DecimalFormat("##0.00;");
		DecimalFormat formatterVolValue = new DecimalFormat("##0.00000;");
		DecimalFormat formatterAnalytic = new DecimalFormat("##0.000;");
		DecimalFormat formatterPercentage = new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));

		double strike = 0.004783;
	
		int liborIndex=0;
		
		//Results with CALIBRATED model
		System.out.println("\n results on CALIBRATED model \n");
		while(liborIndex < liborPeriodDiscretization.getNumberOfTimes()) {
			if (liborIndex<5) {		//da i valori del caplet per maturity 1.5Y, 2Y, 2.5Y,.
				double maturityMinusLengthLibor =liborPeriodDiscretization.getTime(liborIndex);
				double fixBackwardTime=liborPeriodDiscretization.getTime(liborIndex+1);
				Caplet capletCassical = new Caplet(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
				//set to default ValueUnit.NORMALVOLATILITY for CapletOnBackwardLookingRate
				CapletOnBackwardLookingRate capletBackward = new CapletOnBackwardLookingRate(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
				double impliedVolClassical = capletCassical.getValue(simulationCalibrated);
				double impliedVolBackward = capletBackward.getValue(simulationCalibrated);
				double analyticFormulaPaper = Math.sqrt(1+0.5/(fixBackwardTime*3));
				double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
				System.out.println("Caplet on B(" + formatterTimeValue.format(maturityMinusLengthLibor) + ", "
						+ formatterTimeValue.format(fixBackwardTime) + "; "
						+ "" + formatterTimeValue.format(fixBackwardTime) + "). Implied vols: Backward" 
						+ formatterPercentage.format(impliedVolBackward) + "  Classic"
						+ formatterPercentage.format(impliedVolClassical)+ "  Analytic result: " 
						+ formatterAnalytic.format(analyticFormulaPaper)+  "  Ratio vol: "
						+ formatterAnalytic.format(ratioImpliedVol) + "  Error:"
						+ formatterPercentage.format((ratioImpliedVol-analyticFormulaPaper)/ratioImpliedVol) );
				liborIndex+=1;
			}
			else {	//secondo loop da i valori del caplet per maturity 4Y,5Y,...,21
				double maturityMinusLengthLibor =liborPeriodDiscretization.getTime(liborIndex);
				double fixBackwardTime=liborPeriodDiscretization.getTime(liborIndex+1);
				Caplet capletCassical = new Caplet(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
				CapletOnBackwardLookingRate capletBackward = new CapletOnBackwardLookingRate(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
				double impliedVolClassical = capletCassical.getValue(simulationCalibrated);
				double impliedVolBackward = capletBackward.getValue(simulationCalibrated);
				double analyticFormulaPaper = Math.sqrt(1+0.5/(fixBackwardTime*3));
				double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
				System.out.println("Caplet on B(" + formatterTimeValue.format(maturityMinusLengthLibor) + ", "
						+ formatterTimeValue.format(fixBackwardTime) + "; "
						+ "" + formatterTimeValue.format(fixBackwardTime) + "). Implied vols: Backward" 
						+ formatterPercentage.format(impliedVolBackward) + "  Classic"
						+ formatterPercentage.format(impliedVolClassical)+ "  Analytic result: " 
						+ formatterAnalytic.format(analyticFormulaPaper)+  "  Ratio vol: "
						+ formatterAnalytic.format(ratioImpliedVol) + "  Error:"
						+ formatterPercentage.format((ratioImpliedVol-analyticFormulaPaper)/ratioImpliedVol) );
				liborIndex+=2;
			}
		}
		
		// Set model properties
		Map<String, String> properties2 = new HashMap<String, String >();
		properties2.put("measure", LIBORMarketModelFromCovarianceModelWithMercurioModification.Measure.SPOT.name());
		properties2.put("stateSpace", LIBORMarketModelFromCovarianceModelWithMercurioModification.StateSpace.NORMAL.name());

		final LIBORMarketModel MercurioModelNONcalibrated=  new LIBORMarketModelFromCovarianceModelWithMercurioModification(
				liborPeriodDiscretization,
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				randomVariableFactory,
				covarianceModelReducedVolatility,
				properties2
				);

		final EulerSchemeFromProcessModel process2 = new EulerSchemeFromProcessModel(brownianMotion);
		final LIBORModelMonteCarloSimulationModel simulationMercurioModelNONcalibrated = new LIBORMonteCarloSimulationFromLIBORModel(MercurioModelNONcalibrated, process2);
	
		
		//Results with NON calibrated model
		liborIndex=0;
		System.out.println("\n results on NON-CALIBRATED model \n");
		while(liborIndex < liborPeriodDiscretization.getNumberOfTimes()) {
			if (liborIndex<5) {		//da i valori del caplet per maturity 1.5Y, 2Y, 2.5Y,.
				double maturityMinusLengthLibor =liborPeriodDiscretization.getTime(liborIndex);
				double fixBackwardTime=liborPeriodDiscretization.getTime(liborIndex+1);
				Caplet capletCassical = new Caplet(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
				//set to default ValueUnit.NORMALVOLATILITY for CapletOnBackwardLookingRate
				CapletOnBackwardLookingRate capletBackward = new CapletOnBackwardLookingRate(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
				double impliedVolClassical = capletCassical.getValue(simulationMercurioModelNONcalibrated);
				double impliedVolBackward = capletBackward.getValue(simulationMercurioModelNONcalibrated);
				double analyticFormulaPaper = Math.sqrt(1+0.5/(fixBackwardTime*3));
				double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
				System.out.println("Caplet on B(" + formatterTimeValue.format(maturityMinusLengthLibor) + ", "
						+ formatterTimeValue.format(fixBackwardTime) + "; "
						+ "" + formatterTimeValue.format(fixBackwardTime) + "). Implied vols: Backward" 
						+ formatterPercentage.format(impliedVolBackward) + "  Classic"
						+ formatterPercentage.format(impliedVolClassical)+ "  Analytic result: " 
						+ formatterAnalytic.format(analyticFormulaPaper)+  "  Ratio vol: "
						+ formatterAnalytic.format(ratioImpliedVol) + "  Error:"
						+ formatterPercentage.format((ratioImpliedVol-analyticFormulaPaper)/ratioImpliedVol) );
				liborIndex+=1;
			}
			else {	//secondo loop da i valori del caplet per maturity 4Y,5Y,...,21
				double maturityMinusLengthLibor =liborPeriodDiscretization.getTime(liborIndex);
				double fixBackwardTime=liborPeriodDiscretization.getTime(liborIndex+1);
				Caplet capletCassical = new Caplet(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
				CapletOnBackwardLookingRate capletBackward = new CapletOnBackwardLookingRate(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
				double impliedVolClassical = capletCassical.getValue(simulationMercurioModelNONcalibrated);
				double impliedVolBackward = capletBackward.getValue(simulationMercurioModelNONcalibrated);
				double analyticFormulaPaper = Math.sqrt(1+0.5/(fixBackwardTime*3));
				double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
				System.out.println("Caplet on B(" + formatterTimeValue.format(maturityMinusLengthLibor) + ", "
						+ formatterTimeValue.format(fixBackwardTime) + "; "
						+ "" + formatterTimeValue.format(fixBackwardTime) + "). Implied vols: Backward" 
						+ formatterPercentage.format(impliedVolBackward) + "  Classic"
						+ formatterPercentage.format(impliedVolClassical)+ "  Analytic result: " 
						+ formatterAnalytic.format(analyticFormulaPaper)+  "  Ratio vol: "
						+ formatterAnalytic.format(ratioImpliedVol) + "  Error:"
						+ formatterPercentage.format((ratioImpliedVol-analyticFormulaPaper)/ratioImpliedVol) );
				liborIndex+=2;
			}
		}
		

		int timeIndex =  2;
		System.out.println(" \n Backward looking rate evaluated at the END of the accrual period:");
		for( liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
			double evaluationTime=timeDiscretizationFromArray.getTime(timeIndex);
			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
			RandomVariable backwardLookingRate =  simulationCalibrated.getLIBOR(timeIndex, liborIndex);
			double avgBackwardLookingRate =backwardLookingRate.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(evaluationTime) + ", avg. value " + avgBackwardLookingRate);
			timeIndex +=  2;
		}	
		
		int timeIndex2 =  0;
		System.out.println("\n Backward looking rate evaluated at the BEGIN of the accrual period:");
		for( liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
			double evaluationTime=timeDiscretizationFromArray.getTime(timeIndex2);
			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
			RandomVariable backwardLookingRate =  simulationCalibrated.getLIBOR(timeIndex2, liborIndex);
			double avgBackwardLookingRate =backwardLookingRate.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(evaluationTime) + ", avg. value " + avgBackwardLookingRate);
			timeIndex2 +=  2;
		}
		
		
	//CORRETTO, CALIBRA LA IR CURVE E BASTA, L(0) è CORRETTO!	
//		timeIndex2 =  0;
//		System.out.println("\n CHECK OF NON-CALIBRATED MODEL :");
//		for( liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
//			double evaluationTime=timeDiscretizationFromArray.getTime(timeIndex2);
//			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
//			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
//			RandomVariable backwardLookingRate =  simulationMercurioModelNONcalibrated.getLIBOR(timeIndex2, liborIndex);
//			double avgBackwardLookingRate =backwardLookingRate.getAverage();
//			System.out.println("Backward B(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(evaluationTime) + ", avg. value " + avgBackwardLookingRate);
//			timeIndex2 +=  2;
//		}
		
		
		System.out.println("Time required for calibration of curves.........: " + (millisCurvesEnd-millisCurvesStart)/1000.0 + " s.");
		System.out.println("Time required for calibration of volatilities...: " + (millisCalibrationEnd-millisCalibrationStart)/1000.0 + " s.");

		final double averageDeviation = deviationSum/calibrationProducts.size();
		System.out.println("Mean Deviation:" + formatterValue.format(averageDeviation));
		System.out.println("RMS Error.....:" + formatterValue.format(Math.sqrt(deviationSquaredSum/calibrationProducts.size())));
		System.out.println("__________________________________________________________________________________________\n");

		Assert.assertTrue(Math.abs(averageDeviation) < 2E-4);
	}

	
	
	
//------------------------------- getCalibratedCurve() -------------------------------
	public AnalyticModel getCalibratedCurve() throws SolverException {
		final String[] maturity					= { "7D", "14D", "21D", "1M", "2M", "3M", "4M", "5M", "6M", "7M", "8M", "9M", "12M", "15M", "18M", "21M", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "11Y", "12Y", "15Y", "20Y", "25Y", "30Y", "40Y", "50Y" };
		final String[] frequency				= { "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
		final String[] frequencyFloat			= { "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
		final String[] daycountConventions	    = { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
		final String[] daycountConventionsFloat	= { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
		final double[] rates					= { -0.00553, -0.00553, -0.00553, -0.00553, -0.00555, -0.00556, -0.00559, -0.00564, -0.00568, -0.00572, -0.00577, -0.00581, -0.00592, -0.00601, -0.00608, -0.00613, -0.00619, -0.00627, -0.00622, -0.00606, -0.00582, -0.00553, -0.00519, -0.00482, -0.00442, -0.00402, -0.00362, -0.00261, -0.00189, -0.00197, -0.0023, -0.00286, -0.00333};
		final HashMap<String, Object> parameters = new HashMap<>();

		parameters.put("referenceDate", LocalDate.of(2020, Month.JULY, 31));
		parameters.put("currency", "EUR");
		parameters.put("forwardCurveTenor", "1D");
		parameters.put("maturities", maturity);
		parameters.put("fixLegFrequencies", frequency);
		parameters.put("floatLegFrequencies", frequencyFloat);
		parameters.put("fixLegDaycountConventions", daycountConventions);
		parameters.put("floatLegDaycountConventions", daycountConventionsFloat);
		parameters.put("rates", rates);

		return getCalibratedCurve(null, parameters);
	}

	private static AnalyticModel getCalibratedCurve(final AnalyticModel model2, final Map<String, Object> parameters) throws SolverException {

		final LocalDate	referenceDate		= (LocalDate) parameters.get("referenceDate");
		final String	currency			= (String) parameters.get("currency");
		final String	forwardCurveTenor	= (String) parameters.get("forwardCurveTenor");
		final String[]	maturities			= (String[]) parameters.get("maturities");
		final String[]	frequency			= (String[]) parameters.get("fixLegFrequencies");
		final String[]	frequencyFloat		= (String[]) parameters.get("floatLegFrequencies");
		final String[]	daycountConventions	= (String[]) parameters.get("fixLegDaycountConventions");
		final String[]	daycountConventionsFloat	= (String[]) parameters.get("floatLegDaycountConventions");
		final double[]	rates						= (double[]) parameters.get("rates");

		Assert.assertEquals(maturities.length, frequency.length);
		Assert.assertEquals(maturities.length, daycountConventions.length);
		Assert.assertEquals(maturities.length, rates.length);

		Assert.assertEquals(frequency.length, frequencyFloat.length);
		Assert.assertEquals(daycountConventions.length, daycountConventionsFloat.length);

		final int		spotOffsetDays = 2;
		final String	forwardStartPeriod = "0D";

		final String curveNameDiscount = "discountCurve-" + currency;

		/*
		 * We create a forward curve by referencing the same discount curve, since
		 * this is a single curve setup.
		 *
		 * Note that using an independent NSS forward curve with its own NSS parameters
		 * would result in a problem where both, the forward curve and the discount curve
		 * have free parameters.
		 */
		final ForwardCurve forwardCurve		= new ForwardCurveFromDiscountCurve(curveNameDiscount, referenceDate, forwardCurveTenor);

		// Create a collection of objective functions (calibration products)
		final Vector<AnalyticProduct> calibrationProducts = new Vector<>();
		final double[] curveMaturities	= new double[rates.length+1];
		final double[] curveValue			= new double[rates.length+1];
		final boolean[] curveIsParameter	= new boolean[rates.length+1];
		curveMaturities[0] = 0.0;
		curveValue[0] = 1.0;
		curveIsParameter[0] = false;
		for(int i=0; i<rates.length; i++) {

			final Schedule schedulePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequency[i], daycountConventions[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
			final Schedule scheduleRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequencyFloat[i], daycountConventionsFloat[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);

			curveMaturities[i+1] = Math.max(schedulePay.getPayment(schedulePay.getNumberOfPeriods()-1),scheduleRec.getPayment(scheduleRec.getNumberOfPeriods()-1));
			curveValue[i+1] = 1.0;
			curveIsParameter[i+1] = true;
			calibrationProducts.add(new Swap(schedulePay, null, rates[i], curveNameDiscount, scheduleRec, forwardCurve.getName(), 0.0, curveNameDiscount));
		}

		final InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;

		// Create a discount curve
		final DiscountCurveInterpolation discountCurveInterpolation = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				curveNameDiscount								/* name */,
				referenceDate	/* referenceDate */,
				curveMaturities	/* maturities */,
				curveValue		/* discount factors */,
				curveIsParameter,
				interpolationMethod ,
				ExtrapolationMethod.CONSTANT,
				InterpolationEntity.LOG_OF_VALUE
				);

		/*
		 * Model consists of the two curves, but only one of them provides free parameters.
		 */
		AnalyticModel model = new AnalyticModelFromCurvesAndVols(new Curve[] { discountCurveInterpolation, forwardCurve });

		/*
		 * Create a collection of curves to calibrate
		 */
		final Set<ParameterObject> curvesToCalibrate = new HashSet<>();
		curvesToCalibrate.add(discountCurveInterpolation);

		/*
		 * Calibrate the curve
		 */
		final Solver solver = new Solver(model, calibrationProducts, 0.0, 1E-4 /* target accuracy */);
		final AnalyticModel calibratedModel = solver.getCalibratedModel(curvesToCalibrate);
		System.out.println("Solver reported acccurary....: " + solver.getAccuracy());

		Assert.assertEquals("Calibration accurarcy", 0.0, solver.getAccuracy(), 1E-3);

		// Get best parameters
		final double[] parametersBest = calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getParameter();

		// Test calibration
		model			= calibratedModel;

		double squaredErrorSum = 0.0;
		for(final AnalyticProduct c : calibrationProducts) {
			final double value = c.getValue(0.0, model);
			final double valueTaget = 0.0;
			final double error = value - valueTaget;
			squaredErrorSum += error*error;
		}
		final double rms = Math.sqrt(squaredErrorSum/calibrationProducts.size());

		System.out.println("Independent checked acccurary: " + rms);

		System.out.println("Calibrated discount curve: ");
		for(int i=0; i<curveMaturities.length; i++) {
			final double maturity = curveMaturities[i];
			System.out.println(maturity + "\t" + calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getDiscountFactor(maturity));
		}
		return model;
	}

	private static double getParSwaprate(final ForwardCurve forwardCurve, final DiscountCurve discountCurve, final double[] swapTenor) {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretizationFromArray(swapTenor), new TimeDiscretizationFromArray(swapTenor), forwardCurve, discountCurve);
	}
}
