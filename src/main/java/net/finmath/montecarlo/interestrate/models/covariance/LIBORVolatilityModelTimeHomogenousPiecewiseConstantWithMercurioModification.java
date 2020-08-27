/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 08.08.2005
 */
package net.finmath.montecarlo.interestrate.models.covariance;

import java.util.Map;

import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

/**
 * Implements a piecewise constant volatility model, where
 * \( \sigma(t,T) = sigma_{i} \) where \( i = \max \{ j : \tau_{j} \leq T-t \} \) and
 * \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \) is a given time discretization.
 *
 * @author Christian Fries
 * @version 1.0
 */
public class LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification extends LIBORVolatilityModel {

	private static final long serialVersionUID = -1942151065049237807L;

	private final RandomVariableFactory	abstractRandomVariableFactory;

	private final TimeDiscretization timeToMaturityDiscretization;
	private final RandomVariable[] volatility;

	/**
	 * Create a piecewise constant volatility model, where
	 * \( \sigma(t,T) = sigma_{i} \) where \( i = \max \{ j : \tau_{j} \leq T-t \} \) and
	 * \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \) is a given time discretization.
	 *
	 * @param abstractRandomVariableFactory The random variable factor used to construct random variables from the parameters.
	 * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
	 * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
	 * @param timeToMaturityDiscretization The discretization \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \)  of the piecewise constant volatility function.
	 * @param volatility The values \( \sigma_{0}, \sigma_{1}, \ldots, \sigma_{n-1} \) of the piecewise constant volatility function.
	 */
	public LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(final RandomVariableFactory abstractRandomVariableFactory, final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization timeToMaturityDiscretization, final RandomVariable[] volatility) {
		super(timeDiscretization, liborPeriodDiscretization);

		if(timeToMaturityDiscretization.getTime(0) != 0) {
			throw new IllegalArgumentException("timeToMaturityDiscretization should start with 0 as first time point.");
		}
		if(timeToMaturityDiscretization.getNumberOfTimes() != volatility.length) {
			throw new IllegalArgumentException("volatility.length should equal timeToMaturityDiscretization.getNumberOfTimes() .");
		}

		this.abstractRandomVariableFactory = abstractRandomVariableFactory;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.volatility = volatility;
	}

	/**
	 * Create a piecewise constant volatility model, where
	 * \( \sigma(t,T) = sigma_{i} \) where \( i = \max \{ j : \tau_{j} \leq T-t \} \) and
	 * \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \) is a given time discretization.
	 *
	 * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
	 * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
	 * @param timeToMaturityDiscretization The discretization \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \)  of the piecewise constant volatility function.
	 * @param volatility The values \( \sigma_{0}, \sigma_{1}, \ldots, \sigma_{n-1} \) of the piecewise constant volatility function.
	 */
	public LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization timeToMaturityDiscretization, final RandomVariable[] volatility) {
		this(null, timeDiscretization, liborPeriodDiscretization, timeToMaturityDiscretization, volatility);
	}

	/**
	 * Create a piecewise constant volatility model, where
	 * \( \sigma(t,T) = sigma_{i} \) where \( i = \max \{ j : \tau_{j} \leq T-t \} \) and
	 * \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \) is a given time discretization.
	 *
	 * @param abstractRandomVariableFactory The random variable factor used to construct random variables from the parameters.
	 * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
	 * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
	 * @param timeToMaturityDiscretization The discretization \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \)  of the piecewise constant volatility function.
	 * @param volatility The values \( \sigma_{0}, \sigma_{1}, \ldots, \sigma_{n-1} \) of the piecewise constant volatility function.
	 */
	public LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(final RandomVariableFactory abstractRandomVariableFactory, final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[] volatility) {
		super(timeDiscretization, liborPeriodDiscretization);

		if(timeToMaturityDiscretization.getTime(0) != 0) {
			throw new IllegalArgumentException("timeToMaturityDiscretization should start with 0 as first time point.");
		}
		if(timeToMaturityDiscretization.getNumberOfTimes() != volatility.length) {
			throw new IllegalArgumentException("volatility.length should equal timeToMaturityDiscretization.getNumberOfTimes() .");
		}

		this.abstractRandomVariableFactory = abstractRandomVariableFactory;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.volatility = abstractRandomVariableFactory.createRandomVariableArray(volatility);
	}

	/**
	 * Create a piecewise constant volatility model, where
	 * \( \sigma(t,T) = sigma_{i} \) where \( i = \max \{ j : \tau_{j} \leq T-t \} \) and
	 * \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \) is a given time discretization.
	 *
	 * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
	 * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
	 * @param timeToMaturityDiscretization The discretization \( \tau_{0}, \tau_{1}, \ldots, \tau_{n-1} \)  of the piecewise constant volatility function.
	 * @param volatility The values \( \sigma_{0}, \sigma_{1}, \ldots, \sigma_{n-1} \) of the piecewise constant volatility function.
	 */
	
	
// ----------> IMPORTANTE
	// praticamente timeToMaturityDiscretization ti dice a quale LIBOR ci stiamo riferendeo quindi se prendi una discretizzazione annuale significa che, se i libor sono trimestrali, partendo dAL TEMPO 0, i primi 4 libor avranno lo stesso valore della volatiltà poi i secondi 4, avranno un'altra volatilità è cosi via.. rappresenta l'asse vertical della tua idea delle volatilità del LIBOR, mentre la simulationDiscretization è l'asse orizzontale, cioè ti dice ogni quanto i valori della volatilità cambiano, se ad esempio trimestrale, significa che ogni volatilità resta la stesa per un trimestre di simulazione.
	// ma per avere un piecewise constant trimestrale non basterebbe specificare quindi un timeToMaturityDiscretization? e optionMaturityDiscretization dovrebbe essere tipo (0,40).. NO! questo praticamente ti da una time-homogeneous piecewise constatn!!
	// NB:timeToMaturityDiscretization è proprio time to maturity quindi assumendo sia trimestrale un L(0,1;0) prenderà il valore della volatilità con Index=4.
	public LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[] volatility) {
		this(new RandomVariableFromArrayFactory(), timeDiscretization, liborPeriodDiscretization, timeToMaturityDiscretization, volatility);
	}

	@Override
	public RandomVariable[] getParameter() {
		return volatility;
	}

	@Override
	public LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification getCloneWithModifiedParameter(final RandomVariable[] parameter) {
		return new LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(
				abstractRandomVariableFactory,
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				timeToMaturityDiscretization,
				parameter
				);
	}

	@Override
	public RandomVariable getVolatility(final int timeIndex, final int liborIndex) {
		// Create a very simple volatility model here
		final double time             = getTimeDiscretization().getTime(timeIndex);
		
//---> MERCURIO:
		//final double maturity         = getLiborPeriodDiscretization().getTime(liborIndex);
		final double maturity         = getLiborPeriodDiscretization().getTime(liborIndex+1);
		final double timeToMaturity   = maturity-time;

		RandomVariable volatilityInstanteaneous;
		if(timeToMaturity <= 0)
		{
			volatilityInstanteaneous = new Scalar(0.0);   // This forward rate is already fixed, no volatility
		}
		else
		{
			int timeIndexTimeToMaturity = timeToMaturityDiscretization.getTimeIndex(timeToMaturity);
			if(timeIndexTimeToMaturity < 0) {
				timeIndexTimeToMaturity = -timeIndexTimeToMaturity-1-1;
			}
			if(timeIndexTimeToMaturity < 0) {
				timeIndexTimeToMaturity = 0;
			}
			if(timeIndexTimeToMaturity >= timeToMaturityDiscretization.getNumberOfTimes()) {
				timeIndexTimeToMaturity--;
			}
			volatilityInstanteaneous = volatility[timeIndexTimeToMaturity];
		}

		return volatilityInstanteaneous;
	}

	@Override
	public Object clone() {
		return new LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				timeToMaturityDiscretization,
				volatility.clone()
				);
	}

	@Override
	public LIBORVolatilityModel getCloneWithModifiedData(final Map<String, Object> dataModified) {
		// TODO Auto-generated method stub
		return null;
	}
}
