"""
This module implements an Option class which the calculates the premium,
implied volatility, and Greeks of European- and American-style options using
the Black-Scholes and binomial tree models, respectively.
"""

from dataclasses import dataclass, field
import scipy.stats
from numpy import log, sqrt, exp, maximum, zeros, e
from enum import Enum

class OptionCategory(Enum):
    EQUITY = 'E'
    FUTURE = 'F'
    CURRENCY = 'C'

class OptionType(Enum):
    CALL = 'C'
    PUT = 'P'

class ExerciseStyle(Enum):
    EUROPEAN = 'European'
    AMERICAN = 'American'

@dataclass
class Option:
    category: OptionCategory
    option_type: OptionType
    stock_price: float
    strike_price: float
    time_to_expiration: float
    risk_free_rate: float
    volatility: float
    foreign_rate: float = field(default=None)
    dividend_rate: float = field(default=0.0)
    exercise_style: ExerciseStyle = field(default=ExerciseStyle.EUROPEAN)

    def __post_init__(self) -> None:
        self.normal_pdf = scipy.stats.norm.pdf
        self.normal_cdf = scipy.stats.norm.cdf
        self.calculate_d1_and_d2()
        self.calculate_cost_of_carry()

        if self.exercise_style == ExerciseStyle.EUROPEAN:
            self.calculate_premium_european()
        elif self.exercise_style == ExerciseStyle.AMERICAN:
            self.calculate_premium_american()

        self.calculate_greeks()

    def calculate_d1_and_d2(self) -> None:
        """
        Calculate the d1 and d2 parameters used in the Black-Scholes option pricing model.

        These parameters are used in the calculation of option premia and Greeks.
        """
        self.d1 = (log(self.stock_price / self.strike_price) \
                   + (self.risk_free_rate + (self.volatility ** 2) / 2) \
                   * self.time_to_expiration) \
                   / (self.volatility * sqrt(self.time_to_expiration))
        self.d2 = self.d1 - self.volatility * sqrt(self.time_to_expiration)

    def calculate_cost_of_carry(self) -> None:
        """
        Calcualtes the cost of carry based on the category of option.

        This method sets the cost of carry (commonly expressed as the 'b' term in the Black-Scholes model),
        which is determined by the category of topion (equity, future, or currency).

        This method raises a TypeError if the currency category is chosen without providing a foreign interest rate.
        """

        try:
            if self.category == OptionCategory.EQUITY:
                self.cost_of_carry = self.risk_free_rate
            elif self.category == OptionCategory.FUTURE:
                self.cost_of_carry = 0
            elif self.category == OptionCategory.CURRENCY:
                self.cost_of_carry = self.risk_free_rate - self.foreign_rate
            elif self.dividend_rate is not None:
                self.cost_of_carry = self.risk_free_rate - self.dividend_rate
        except TypeError:
            print('ERROR: FX Options Must Include a Foreign Interest-Rate')

    def calculate_premium_european(self) -> None:
        """
        Calculates the premium of a European-style option.

        This method uses the Black-Scholes formula to calculate the premium of
        a European-style call or put option. The calculated premium is then stored in the
        attribute `self.premium`.
        """
        if self.option_type == OptionType.CALL:
            self.premium = self.normal_cdf(self.d1) \
                * self.stock_price - self.normal_cdf(self.d2) \
                * self.strike_price \
                * exp(-self.risk_free_rate * self.time_to_expiration)
        elif self.option_type == OptionType.PUT:
            self.premium = self.normal_cdf(-self.d2) \
                * self.strike_price \
                * exp(-self.risk_free_rate * self.time_to_expiration) \
                - self.normal_cdf(-self.d1) \
                * self.stock_price

    def calculate_greeks(self) -> None:
        self.calculate_delta()
        self.calculate_gamma()
        self.calculate_theta()
        self.calculate_vega()
        self.calculate_rho()

    def calculate_delta(self) -> None:
        if self.option_type == OptionType.CALL:
            self.delta = exp(self.cost_of_carry - self.risk_free_rate) * self.time_to_expiration

    def calculate_gamma(self) -> None:
        self.gamma = exp((self.cost_of_carry - self.risk_free_rate) * self.time_to_expiration) \
            * self.normal_pdf(self.d1) / self.stock_price \
            * self.volatility \
            * sqrt(self.time_to_expiration)

    def calculate_theta(self) -> None:
        if self.option_type == OptionType.CALL:
            self.theta = -self.stock_price \
                * self.normal_pdf(self.d1) \
                * self.volatility / (2 * sqrt(self.time_to_expiration)) \
                - self.cost_of_carry * self.stock_price \
                * exp(-self.cost_of_carry * self.time_to_expiration) \
                * self.normal_cdf(self.d1) - self.risk_free_rate \
                * self.strike_price * exp(-self.risk_free_rate * self.time_to_expiration) \
                * self.normal_cdf(self.d2)
        elif self.option_type == OptionType.PUT:
            self.theta = -self.stock_price \
                * self.normal_pdf(self.d1) \
                * self.volatility / (2 * sqrt(self.time_to_expiration)) \
                + self.cost_of_carry * self.stock_price \
                * exp(-self.cost_of_carry * self.time_to_expiration) \
                * self.normal_cdf(-self.d1) + self.risk_free_rate \
                * self.strike_price * exp(-self.risk_free_rate * self.time_to_expiration) \
                * self.normal_cdf(-self.d2)
            
    def calculate_vega(self) -> None:
        self.vega = self.stock_price * self.normal_cdf(self.d1) * sqrt(self.time_to_expiration)

    def calculate_rho(self) -> None:
        if self.option_type == OptionType.CALL and self.cost_of_carry != 0:
            self.rho = self.time_to_expiration \
                * self.strike_price \
                * exp(-self.risk_free_rate * self.time_to_expiration) \
                * self.normal_cdf(self.d2)
        elif self.option_type == OptionType.PUT and self.cost_of_carry != 0:
            self.rho = -self.time_to_expiration \
                * self.strike_price \
                * e(-self.risk_free_rate * self.time_to_expiration) \
                * self.normal_cdf(-self.d2)
        else:
            self.rho = -self.time_to_expiration * self.premium()

    def calculate_premium_american(self) -> None:
        self.premium = self.binomial_tree()

    def binomial_tree(self, steps=100) -> float:
        """
        Calculates the premium of an American-style option.

        This method uses the binomial pricing model to calculate the premium of
        an American-style call or put option. The calculated premium is then stored in the
        attribute `self.premium`.
        """

        time_step = self.time_to_expiration / steps
        up_factor = exp(self.volatility * sqrt(time_step))
        down_factor = 1 / up_factor
        probability = (exp((self.cost_of_carry - self.dividend_rate) * time_step) - down_factor) / (up_factor - down_factor)
        discount_factor = exp(-self.risk_free_rate * time_step)

        stock_tree = zeros((steps + 1, steps + 1))
        option_tree = zeros((steps + 1, steps + 1))

        for i in range(steps + 1):
            for j in range(i + 1):
                stock_tree[j, i] = self.stock_price * (up_factor ** (i - j)) * (down_factor ** j)

        if self.option_type == OptionType.CALL:
            option_tree[:, steps] = maximum(zeros(steps + 1), stock_tree[:, steps] - self.strike_price)
        elif self.option_type == OptionType.PUT:
            option_tree[:, steps] = maximum(zeros(steps + 1), self.strike_price - stock_tree[:, steps])

        for i in range(steps - 1, -1, -1):
            for j in range(i + 1):
                option_tree[j, i] = discount_factor * (probability * option_tree[j, i + 1] + (1 - probability) * option_tree[j + 1, i + 1])
                if self.option_type == OptionType.CALL:
                    option_tree[j, i] = max(option_tree[j, i], stock_tree[j, i] - self.strike_price)
                elif self.option_type == OptionType.PUT:
                    option_tree[j, i] = max(option_tree[j, i], self.strike_price - stock_tree[j, i])

        return option_tree[0, 0]

    def implied_volatility(self, target_premium, initial_guess=0.2, max_iterations=100, tolerance=1e-6) -> None:
        """
        Calculates the implied volatility of an option, given a target option premium.

        This method uses the Newton-Raphson method; an iterative method for finding the root of a real-valued function.
        In seeking the minimize the difference between the calculated option premium and the target option premium,
        this method updates the current estimate of implied volatility by subtracting the ratio of the function value (`premium_difference`)
        to its derivative (vega) at the current estimate. The method then iterates until either `premium_difference`
        is within the specified tolerance or the maximum number of iterations has been reached.
        
        This method raises a RuntimeError if it does not converge within the specified number of iterations.
        """
        sigma = initial_guess
        for i in range(max_iterations):
            self.volatility = sigma
            self.calculate_d1_and_d2()
            self.calculate_premium_european()
            premium_difference = self.premium - target_premium
            if abs(premium_difference) < tolerance:
                return sigma

            self.calculate_vega()
            
            # Check if vega is very close to zero
            if abs(self.vega) < 1e-6:
                # Update sigma by a small fixed amount
                sigma += 1e-6 if premium_difference > 0 else -1e-6
            else:
                sigma -= premium_difference / self.vega

        raise RuntimeError("Implied volatility calculation did not converge")
