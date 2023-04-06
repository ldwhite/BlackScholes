from dataclasses import dataclass, field
import scipy.stats
from numpy import *
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
    dividend_rate: float = field(default=None)
    exercise_style: ExerciseStyle = field(default=ExerciseStyle.EUROPEAN)

    def __post_init__(self):
        self.normal_pdf = scipy.stats.norm.pdf
        self.normal_cdf = scipy.stats.norm.cdf
        self.calculate_d1_and_d2()
        self.calculate_cost_of_carry()

        if self.exercise_style == ExerciseStyle.EUROPEAN:
            self.calculate_premium_european()
        elif self.exercise_style == ExerciseStyle.AMERICAN:
            self.calculate_premium_american()

        self.calculate_greeks()

    def calculate_d1_and_d2(self):
        self.d1 = (log(self.stock_price / self.strike_price) + \
                   (self.risk_free_rate + (self.volatility ** 2) / 2) * \
                    self.time_to_expiration) / \
                    (self.volatility * sqrt(self.time_to_expiration))
        self.d2 = self.d1 - self.volatility * sqrt(self.time_to_expiration)

    def calculate_cost_of_carry(self):
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

    def calculate_premium_european(self):
        if self.option_type == OptionType.CALL:
            self.premium = self.normal_cdf(self.d1) * \
                self.stock_price - self.normal_cdf(self.d2) * \
                self.strike_price * \
                exp(-self.risk_free_rate * self.time_to_expiration)
        elif self.option_type == OptionType.PUT:
            self.premium = self.normal_cdf(-self.d2) * \
                self.strike_price * \
                exp(-self.risk_free_rate * self.time_to_expiration) - \
                self.normal_cdf(-self.d1) * \
                self.stock_price

    def calculate_greeks(self):
        self.calculate_delta()
        self.calculate_gamma()
        self.calculate_theta()
        self.calculate_vega()
        self.calculate_rho()

    def calculate_delta(self):
        if self.option_type == OptionType.CALL:
            self.delta = exp(self.cost_of_carry - self.risk_free_rate) * self.time_to_expiration

    def calculate_gamma(self):
        self.gamma = exp((self.cost_of_carry - self.risk_free_rate) * self.time_to_expiration) * \
            self.normal_pdf(self.d1) / self.stock_price * \
            self.volatility * \
            sqrt(self.time_to_expiration)

    def calculate_theta(self):
        if self.option_type == OptionType.CALL:
            self.theta = -self.stock_price * \
                self.normal_pdf(self.d1) * \
                self.volatility / (2 * sqrt(self.time_to_expiration)) - \
                self.cost_of_carry * self.stock_price * \
                exp(-self.cost_of_carry * self.time_to_expiration) * \
                self.normal_cdf(self.d1) - self.risk_free_rate * \
                self.strike_price * exp(-self.risk_free_rate * self.time_to_expiration) * \
                self.normal_cdf(self.d2)
        elif self.option_type == OptionType.PUT:
            self.theta = -self.stock_price * \
                self.normal_pdf(self.d1) * \
                self.volatility / (2 * sqrt(self.time_to_expiration)) + \
                self.cost_of_carry * self.stock_price * \
                exp(-self.cost_of_carry * self.time_to_expiration) * \
                self.normal_cdf(-self.d1) + self.risk_free_rate * \
                self.strike_price * exp(-self.risk_free_rate * self.time_to_expiration) * \
                self.normal_cdf(-self.d2)
            
    def calculate_vega(self):
        self.vega = self.stock_price * \
            self.normal_cdf(self.d1) * sqrt(self.time_to_expiration)

    def calculate_rho(self):
        if self.option_type == OptionType.CALL and self.cost_of_carry != 0:
            self.rho = self.time_to_expiration * \
                self.strike_price * \
                exp(-self.risk_free_rate * self.time_to_expiration) * \
                self.normal_cdf(self.d2)
        elif self.option_type == OptionType.PUT and self.cost_of_carry != 0:
            self.rho = -self.time_to_expiration * \
                self.strike_price * \
                e(-self.risk_free_rate * self.time_to_expiration) * \
                self.normal_cdf(-self.d2)
        else:
            self.rho = -self.time_to_expiration * self.premium()

    def calculate_premium_american(self):
        self.premium = self.binomial_tree()

    def binomial_tree(self, steps=100):
        dt = self.time_to_expiration / steps
        u = exp(self.volatility * sqrt(dt))
        d = 1 / u
        p = (exp((self.cost_of_carry - self.dividend_yield) * dt) - d) / (u - d)
        discount_factor = exp(-self.risk_free_rate * dt)

        stock_tree = zeros((steps + 1, steps + 1))
        option_tree = zeros((steps + 1, steps + 1))

        for i in range(steps + 1):
            for j in range(i + 1):
                stock_tree[j, i] = self.stock_price * (u ** (i - j)) * (d ** j)

        if self.option_type == OptionType.CALL:
            option_tree[:, steps] = maximum(zeros(steps + 1), stock_tree[:, steps] - self.strike_price)
        elif self.option_type == OptionType.PUT:
            option_tree[:, steps] = maximum(zeros(steps + 1), self.strike_price - stock_tree[:, steps])

        for i in range(steps - 1, -1, -1):
            for j in range(i + 1):
                option_tree[j, i] = discount_factor * (p * option_tree[j, i + 1] + (1 - p) * option_tree[j + 1, i + 1])
                if self.option_type == OptionType.CALL:
                    option_tree[j, i] = max(option_tree[j, i], stock_tree[j, i] - self.strike_price)
                elif self.option_type == OptionType.PUT:
                    option_tree[j, i] = max(option_tree[j, i], self.strike_price - stock_tree[j, i])

        return option_tree[0, 0]

    def implied_volatility(self, target_premium, initial_guess=0.2, max_iterations=100, tolerance=1e-6):
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
