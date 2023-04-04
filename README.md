# Black-Scholes Model Option Calculator

This Python script constructs a DataClass which, when provided with the Option category (Equity, FX, etc), type (call or put), strike price, price of the underlying, days to expiration (which the user must divide by total days in the year), interest rate, foreign interest rate (if an FX option), and annualized volatility.

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

<p align="center">
<img src = "https://raw.githubusercontent.com/ldwhite/BlackScholes/main/options.png" style = "width:80%" />
</p>
