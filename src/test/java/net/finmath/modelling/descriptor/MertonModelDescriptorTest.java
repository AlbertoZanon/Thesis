package net.finmath.modelling.descriptor;

import java.time.LocalDate;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.modelling.DescribedModel;
import net.finmath.modelling.Product;
import net.finmath.modelling.ProductDescriptor;
import net.finmath.modelling.modelfactory.AssetModelFourierMethodFactory;
import net.finmath.modelling.modelfactory.AssetModelMonteCarloFactory;
import net.finmath.montecarlo.MertonJumpProcess;
import net.finmath.time.FloatingpointDate;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

/**
 * Unit test creating a Merton model and a European option from corresponding model descriptors and product descriptors
 * using two different factories: Fourier versus Monte-Carlo.
 *
 * @author Alessandro Gnoatto
 */
public class MertonModelDescriptorTest {

	// Model properties
	private static final LocalDate referenceDate = LocalDate.of(2017,8,15);

	private final double	initialValue   = 1.0;
	private final double	riskFreeRate   = 0.05;
	private final double	volatility     = 0.30;

	private final double m = 1.0;
	private final double nu = 0.15;

	private final double lambda = 0.4;
	private final double jumpSizeStdDev = nu;
	private final double jumpSizeMean = Math.log(m);

	// Process discretization properties
	private final int		numberOfPaths		= 1000000;
	private final int		numberOfTimeSteps	= 100;
	private final double	deltaT				= 0.02;

	private final int		seed				= 3141;

	// Product properties
	private static final double maturity 			= 1.0;
	private static final LocalDate maturityDate	= FloatingpointDate.getDateFromFloatingPointDate(referenceDate, maturity);
	private static final double strike				= 0.95;

	@Test
	public void test() {
		MertonModelDescriptor mertonModelDescriptor = new MertonModelDescriptor(referenceDate, initialValue,
				getDiscountCurve("forward curve", referenceDate, riskFreeRate),
				getDiscountCurve("discount curve", referenceDate, riskFreeRate),
				volatility, lambda, jumpSizeMean, jumpSizeStdDev);
		/*
		 * Create European option descriptor
		 */
		String underlyingName = "eurostoxx";
		ProductDescriptor europeanOptionDescriptor = (new SingleAssetEuropeanOptionProductDescriptor(underlyingName, maturityDate, strike));

		/*
		 * Create Fourier implementation of model and product
		 */

		// Create Fourier implementation of Heston model
		DescribedModel<?> mertonModelFourier = (new AssetModelFourierMethodFactory()).getModelFromDescriptor(mertonModelDescriptor);

		// Create product implementation compatible with Heston model
		Product europeanOptionFourier = mertonModelFourier.getProductFromDescriptor(europeanOptionDescriptor)
				;
		// Evaluate product
		double evaluationTime = 0.0;
		Map<String, Object> valueFourier = europeanOptionFourier.getValues(evaluationTime, mertonModelFourier);

		System.out.println(valueFourier);

		/*
		 * Create Monte Carlo implementation of model and product
		 */
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0 /* initial */, numberOfTimeSteps, deltaT);
		
		MertonJumpProcess mertonProcess = new MertonJumpProcess(lambda, jumpSizeMean, jumpSizeStdDev, timeDiscretization, numberOfPaths, seed);

		// Create Fourier implementation of Heston model
		DescribedModel<?> mertonModelMonteCarlo = (new AssetModelMonteCarloFactory(mertonProcess)).getModelFromDescriptor(mertonModelDescriptor);
		
		// Create product implementation compatible with Variance Gamma model
		Product europeanOptionMonteCarlo = mertonModelMonteCarlo.getProductFromDescriptor(europeanOptionDescriptor);
		
		Map<String, Object> valueMonteCarlo = europeanOptionMonteCarlo.getValues(evaluationTime, mertonModelMonteCarlo);

		System.out.println(valueMonteCarlo);

		double deviation = (Double)valueMonteCarlo.get("value") - (Double)valueFourier.get("value");
		Assert.assertEquals("Difference of Fourier and Monte-Carlo valuation", 0.0, deviation, 1E-2);
	}

	/**
	 * Get the discount curve using the riskFreeRate.
	 *
	 * @param name Name of the curve
	 * @param referenceDate Date corresponding to t=0.
	 * @param riskFreeRate Constant continuously compounded rate
	 *
	 * @return the discount curve using the riskFreeRate.
	 */
	private static DiscountCurve getDiscountCurve(String name, LocalDate referenceDate, double riskFreeRate) {
		double[] times = new double[] { 1.0 };
		double[] givenAnnualizedZeroRates = new double[] { riskFreeRate };
		InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;
		InterpolationEntity interpolationEntity = InterpolationEntity.LOG_OF_VALUE_PER_TIME;
		ExtrapolationMethod extrapolationMethod = ExtrapolationMethod.CONSTANT;
		DiscountCurve discountCurve = DiscountCurveInterpolation.createDiscountCurveFromAnnualizedZeroRates(name, referenceDate, times, givenAnnualizedZeroRates, interpolationMethod, extrapolationMethod, interpolationEntity);
		return discountCurve;
	}
}
