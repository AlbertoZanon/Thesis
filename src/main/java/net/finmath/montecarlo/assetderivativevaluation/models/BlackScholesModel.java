/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 20.01.2004
 */
package net.finmath.montecarlo.assetderivativevaluation.models;

import java.util.Arrays;
import java.util.Map;

import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.stochastic.RandomVariable;

/**
 * This class implements a Black Scholes Model, that is, it provides the drift and volatility specification
 * and performs the calculation of the numeraire (consistent with the dynamics, i.e. the drift).
 *
 * The model is
 * \[
 * 	dS = r S dt + \sigma S dW, \quad S(0) = S_{0},
 * \]
 * \[
 * 	dN = r N dt, \quad N(0) = N_{0},
 * \]
 *
 * The class provides the model of S to an <code>{@link net.finmath.montecarlo.process.MonteCarloProcess}</code> via the specification of
 * \( f = exp \), \( \mu = r - \frac{1}{2} \sigma^2 \), \( \lambda_{1,1} = \sigma \), i.e.,
 * of the SDE
 * \[
 * 	dX = \mu dt + \lambda_{1,1} dW, \quad X(0) = \log(S_{0}),
 * \]
 * with \( S = f(X) \). See {@link net.finmath.montecarlo.process.MonteCarloProcess} for the notation.
 *
 * @author Christian Fries
 * @see net.finmath.montecarlo.process.MonteCarloProcess The interface for numerical schemes.
 * @see net.finmath.montecarlo.model.ProcessModel The interface for models provinding parameters to numerical schemes.
 * @version 1.0
 */
public class BlackScholesModel extends AbstractProcessModel {

	private final RandomVariable initialValue;
	private final RandomVariable riskFreeRate;
	private final RandomVariable volatility;

	private final RandomVariableFactory abstractRandomVariableFactory;

	// Cache for arrays provided though AbstractProcessModel
	private final RandomVariable[]	initialState;
	private final RandomVariable[]	drift;
	private final RandomVariable[]	factorLoadings;

	/**
	 * Create a Black-Scholes specification implementing AbstractProcessModel.
	 *
	 * @param initialValue Spot value.
	 * @param riskFreeRate The risk free rate.
	 * @param volatility The log volatility.
	 * @param abstractRandomVariableFactory The random variable factory used to create random variables from constants.
	 */
	public BlackScholesModel(
			final RandomVariable initialValue,
			final RandomVariable riskFreeRate,
			final RandomVariable volatility,
			final RandomVariableFactory abstractRandomVariableFactory) {
		super();

		this.initialValue = initialValue;
		this.volatility = volatility;
		this.riskFreeRate	= riskFreeRate;
		this.abstractRandomVariableFactory = abstractRandomVariableFactory;

		// Cache
		initialState = new RandomVariable[] { initialValue.log() };
		drift = new RandomVariable[] { riskFreeRate.sub(volatility.squared().div(2)) };
		factorLoadings = new RandomVariable[] { volatility };
	}

	/**
	 * Create a Monte-Carlo simulation using given time discretization.
	 *
	 * @param initialValue Spot value.
	 * @param riskFreeRate The risk free rate.
	 * @param volatility The log volatility.
	 * @param abstractRandomVariableFactory The random variable factory used to create random variables from constants.
	 */
	public BlackScholesModel(
			final double initialValue,
			final double riskFreeRate,
			final double volatility,
			final RandomVariableFactory abstractRandomVariableFactory) {
		this(abstractRandomVariableFactory.createRandomVariable(initialValue), abstractRandomVariableFactory.createRandomVariable(riskFreeRate), abstractRandomVariableFactory.createRandomVariable(volatility), abstractRandomVariableFactory);
	}

	/**
	 * Create a Black-Scholes model from given parameters.
	 *
	 * @param initialValue Spot value.
	 * @param riskFreeRate The risk free rate.
	 * @param volatility The log volatility.
	 */
	public BlackScholesModel(
			final double initialValue,
			final double riskFreeRate,
			final double volatility) {
		this(initialValue, riskFreeRate, volatility, new RandomVariableFromArrayFactory());
	}

	@Override
	public RandomVariable[] getInitialState() {
		return initialState;
	}

	@Override
	public RandomVariable[] getDrift(final int timeIndex, final RandomVariable[] realizationAtTimeIndex, final RandomVariable[] realizationPredictor) {
		return drift;
	}

	@Override
	public RandomVariable[] getFactorLoading(final int timeIndex, final int component, final RandomVariable[] realizationAtTimeIndex) {
		return factorLoadings;
	}

	@Override
	public RandomVariable applyStateSpaceTransform(final int componentIndex, final RandomVariable randomVariable) {
		return randomVariable.exp();
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(final int componentIndex, final RandomVariable randomVariable) {
		return randomVariable.log();
	}

	@Override
	public RandomVariable getNumeraire(final double time) {
		return riskFreeRate.mult(time).exp();
	}

	@Override
	public int getNumberOfComponents() {
		return 1;
	}

	@Override
	public RandomVariable getRandomVariableForConstant(final double value) {
		return abstractRandomVariableFactory.createRandomVariable(value);
	}

	@Override
	public BlackScholesModel getCloneWithModifiedData(final Map<String, Object> dataModified) {
		/*
		 * Determine the new model parameters from the provided parameter map.
		 */
		final double	newInitialValue	= dataModified.get("initialValue") != null	? ((Number)dataModified.get("initialValue")).doubleValue() 	: initialValue.getAverage();
		final double	newRiskFreeRate	= dataModified.get("riskFreeRate") != null	? ((Number)dataModified.get("riskFreeRate")).doubleValue()	: getRiskFreeRate().getAverage();
		final double	newVolatility	= dataModified.get("volatility") != null	? ((Number)dataModified.get("volatility")).doubleValue()	: getVolatility().getAverage();

		return new BlackScholesModel(newInitialValue, newRiskFreeRate, newVolatility, abstractRandomVariableFactory);
	}

	/**
	 * Return the initial value of this model.
	 *
	 * @return the initial value of this model.
	 */
	@Override
	public RandomVariable[] getInitialValue() {
		return new RandomVariable[] { initialValue };
	}

	/**
	 * Returns the risk free rate parameter of this model.
	 *
	 * @return Returns the riskFreeRate.
	 */
	public RandomVariable getRiskFreeRate() {
		return riskFreeRate;
	}

	/**
	 * Returns the volatility parameter of this model.
	 *
	 * @return Returns the volatility.
	 */
	public RandomVariable getVolatility() {
		return factorLoadings[0];
	}

	@Override
	public String toString() {
		return "BlackScholesModel [initialValue=" + initialValue + ", riskFreeRate=" + riskFreeRate + ", volatility="
				+ volatility + ", abstractRandomVariableFactory=" + abstractRandomVariableFactory + ", initialState="
				+ Arrays.toString(initialState) + ", drift=" + Arrays.toString(drift) + ", factorLoadings="
				+ Arrays.toString(factorLoadings) + ", getProcess()=" + getProcess() + "]";
	}
}
