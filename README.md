# Option Calculator

This Python script constructs a DataClass which, when provided with the Option category (Equity, FX, etc), type (call or put), strike price, price of the underlying, days to expiration (which the user must divide by total days in the year), interest rate, foreign interest rate (if an FX option), and annualized volatility. This script will then generate an Option object with the option's premium, as well as all of its greeks, as attributes to that object using either the Black-Scholes model (for a European option) or the binomial tree method (for an American option). Additionally, if one provide the model with a target premium for the option, the model will calculate its implied volatility using the Newton-Raphson method.

The option premium of a European call, C, is given by
$$C = S e^{(b - r)t} N(d_{1}) - X e^{-rt} N(d_{2})$$

and the premium of a European put, P, is given by

$$P = X e^{-rt} N(-d_{2}) - S e^{(b-r)t} N(-d_{1})$$

Where

$S =$ Spot Price  
$X =$ Exercise Price  
$t =$ Time to expiration, in years  
$r =$ domestic interest rate  
$\sigma =$ Annualized volatility (standard deviation) in percent  

Further,

$$d_{1} = \frac{\ln(\frac{S}{X}) + (b + \frac{\sigma^2}{2})t}{\sigma \sqrt{t}}$$

and $$d_{2} = d_{1} - \sigma \sqrt{t}$$

Finally, options are differentiated by the type of underlying asset with

$b = r:$ Equity options  
$b = r = 0:$ Options on futures, subject to futures-type settlement  
$b = 0:$ Options on futures, subject to equity-type settlement  
$b = r - r_{f}:$ Garman-Kohlhagen model for foreign currency options, where $r_f$ is the foreing interest rate  

Sources: Options Volatility and Pricing by Sheldon Natenburg (2e), The Complete Guide to Option Pricing Formulas by Espen Gaarder Haug (2e)

### Example

The following is an example using an equity call option with a strike price of 3870, underlying price of 3858.01, time to expiration of 24 days, interest rate of 0.0325, and annualized volatility of 0.2288 (note that these numbers are for example purposes only).

<p align="center">
<img src = "https://raw.githubusercontent.com/ldwhite/BlackScholes/main/options.png" style = "width:80%" />
</p>

Then, by running `option.premium` you will find that the Black-Scholes model determines a fair value of this option's premium to be `88.48`

Additionally, if you find that the market has priced that same option at 87, then you can calculate the implied probability of the option by setting `market_premium = 87` and running `option.implied_volatility(market_premium)`.
