package net.finmath.montecarlo.interestrate;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.Month;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.junit.Assert;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.calibration.ParameterObject;
import net.finmath.marketdata.calibration.Solver;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelFromCurvesAndVols;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveFromDiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.marketdata.model.volatilities.SwaptionMarketData;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionLazyInit;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModelWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecayWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelPiecewiseConstantWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.VolatilityReductionMercurioModel;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.Caplet.ValueUnit;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.Schedule;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.TimeDiscretizationFromArray.ShortPeriodLocation;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;

/**
 * This class creates a LIBOR market model, basing on the classes of the Finmath library
 * 
 * @author Andrea Mazzon
 *
 */
public class LIBORMarketModelTestingMercurio {
	
	public static void main(final String[] args) throws CalculationException, SolverException {
		RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

		final int numberOfPaths		= 500;
		final int numberOfFactors	= 1;
		final AnalyticModel curveModel = getCalibratedCurve();
		// Create the forward curve (initial value of the LIBOR market model)
		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,3M)");
		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");
		final double lastTime	= 40.0;
		final double dt		= 0.1;	
		final double liborDt		= 0.3;	
//		double[] fixingForForwards = { 0.5, 1.0, 3.0, 4.0, 20.0 }; // fixing for the forwards
//		double[] forwardsForCurve = { 0.05, 0.05, 0.05, 0.05, 0.05 };
//		ForwardCurve forwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards(
//				"forwardCurve", // name of the curve
//				fixingForForwards, // fixings of the forward 
//				forwardsForCurve, // forwards 
//				liborDt 	);

			
		final TimeDiscretization timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, 22, 0.33, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		final TimeDiscretization liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0, 22, 0.43, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		
		final TimeDiscretization optionMaturityDiscretization = new TimeDiscretizationFromArray(0.0, 0.25, 0.50, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0);
		TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0);
		
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, 31415 /* seed */);

		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization,optionMaturityDiscretization,timeToMaturityDiscretization, 0.50 / 100);
		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecayWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);
		final AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);
		final AbstractLIBORCovarianceModelParametric covarianceModelDisplaced = new DisplacedLocalVolatilityModel(covarianceModelParametric, 1.0/0.25, false /* isCalibrateable */);
		final AbstractLIBORCovarianceModelParametric covarianceModelReducedVolatility = new VolatilityReductionMercurioModel(covarianceModelDisplaced);
		
		
		// Set model properties
		Map<String, String> properties = new HashMap<String, String >();
		properties.put("measure", LIBORMarketModelFromCovarianceModelWithMercurioModification.Measure.SPOT.name());
		properties.put("stateSpace", LIBORMarketModelFromCovarianceModelWithMercurioModification.StateSpace.NORMAL.name());
		SwaptionMarketData swaptionMarketData = null;
	

	
	final LIBORMarketModel liborMarketModel=  new LIBORMarketModelFromCovarianceModelWithMercurioModification(
			liborPeriodDiscretization,
			curveModel,
			forwardCurve,
			new DiscountCurveFromForwardCurve(forwardCurve),
			randomVariableFactory,
			covarianceModelReducedVolatility,
			properties);


		EulerSchemeFromProcessModel process = new 
				EulerSchemeFromProcessModel(brownianMotion,
						EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);
		EulerSchemeFromProcessModel process2 = new 
				EulerSchemeFromProcessModel(brownianMotion,
						EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR);

		LIBORMonteCarloSimulationFromLIBORModel simulationLiborModel = new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModel, process);
		
		final LIBORMarketModel liborMarketModelNoVolReduction =  new LIBORMarketModelFromCovarianceModelWithMercurioModification(
				liborPeriodDiscretization,
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				randomVariableFactory,
				covarianceModelDisplaced,
				properties);

		LIBORMonteCarloSimulationFromLIBORModel simulationLiborModelNoVolReduction = new LIBORMonteCarloSimulationFromLIBORModel(liborMarketModelNoVolReduction, process2);

		RandomVariable backwardLookingRate =  simulationLiborModelNoVolReduction.getLIBOR(3, 10);
		System.out.println(backwardLookingRate);
		RandomVariable backwardLookingRate2 =  simulationLiborModelNoVolReduction.getLIBOR(4, 10);
		System.out.println(backwardLookingRate2);
		RandomVariable backwardLookingRate3 =  simulationLiborModelNoVolReduction.getLIBOR(5, 10);
		System.out.println(backwardLookingRate3);
//		DecimalFormat formatterTimeValue = new DecimalFormat("##0.00;");
//
//		for(int liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
//			int fixingBackwardIndex =  (liborIndex+1)*3;
//			double fixingTime=timeDiscretizationFromArray.getTime(fixingBackwardIndex);
//			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
//			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
//			RandomVariable backwardLookingRate =  simulationLiborModel.getLIBOR(fixingBackwardIndex, liborIndex);
//			double volatilityBackwardLookingRate = backwardLookingRate.getVariance();
//	
//			System.out.println("liborIndex " + liborIndex +" i.e. libor L(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(fixingTime) + ", vol. value " + volatilityBackwardLookingRate);
//	    }
//		System.out.println("\n WITHOUT VOLATILITY REDUCTION: \n ");
//		
//		for(int liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
//			int fixingBackwardIndex =  (liborIndex+1)*3;
//			double fixingTime=timeDiscretizationFromArray.getTime(fixingBackwardIndex);
//			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
//			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
//			RandomVariable backwardLookingRate =  simulationLiborModelNoVolReduction.getLIBOR(fixingBackwardIndex, liborIndex);
//			double volatilityBackwardLookingRate = backwardLookingRate.getVariance();
//	
//			System.out.println("liborIndex " + liborIndex +" i.e. libor L(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(fixingTime) + ", vol. value " + volatilityBackwardLookingRate);
//	    }
//		System.out.println("\n check correct fixing, value should be the same as before: \n ");
//		
//		for(int liborIndex=0; liborIndex<liborPeriodDiscretization.getNumberOfTimes()-1; liborIndex++) {
//			int fixingBackwardIndexClassical=  (liborIndex+1)*3;
//			int fixingBackwardIndex =  (liborIndex+1)*4;
//			double fixingTime=timeDiscretizationFromArray.getTime(fixingBackwardIndex);
//			double liborStartingTime=liborPeriodDiscretization.getTime(liborIndex);
//			double liborEndingTime=liborPeriodDiscretization.getTime(liborIndex+1);
//			RandomVariable backwardLookingRate =  simulationLiborModelNoVolReduction.getLIBOR(fixingBackwardIndexClassical, liborIndex);
//			RandomVariable backwardLookingRateClassical =  simulationLiborModelNoVolReduction.getLIBOR(fixingBackwardIndex, liborIndex);
//			double volatilityBackwardLookingRateClassical = backwardLookingRateClassical.getVariance();
//
//			double volatilityBackwardLookingRate = backwardLookingRate.getVariance();
//	
//			System.out.println("liborIndex " + liborIndex +" i.e. libor L(" + formatterTimeValue.format(liborStartingTime) +  ", " + formatterTimeValue.format(liborEndingTime) + ") evaluated in t= " + formatterTimeValue.format(fixingTime)+ ", back. vol. value " + volatilityBackwardLookingRate + "is =? to LMM vol " + volatilityBackwardLookingRateClassical );
//
//			}
		
		}
	

   
	
	public static AnalyticModel getCalibratedCurve() throws SolverException {
		final String[] maturity					= { "12M", "15M", "18M", "21M", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "12Y", "15Y", "20Y", "25Y", "30Y", "40Y", "50Y" };
		final String[] frequency				= { "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
		final String[] frequencyFloat			= { "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly", "quarterly"};
		final String[] daycountConventions		= { "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360", "E30/360"};
		final String[] daycountConventionsFloat	= { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
		final double[] rates					= { -0.00478, -0.00484, -0.00489, -0.00489, -0.00492, -0.00492, -0.0048, -0.00459, -0.00431, -0.00398, -0.00363, -0.00323, -0.00281, -0.00199, -0.00098, -0.00025, -0.00033, -0.00068, -0.00129, -0.00182 };

		final HashMap<String, Object> parameters = new HashMap<>();

		parameters.put("referenceDate", LocalDate.of(2020, Month.JULY, 31));
		parameters.put("currency", "EUR");
		parameters.put("forwardCurveTenor", "3M");
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
}
 

	