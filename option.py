from dataclasses import dataclass, field
import scipy.stats
from numpy import *

@dataclass
class BSMOption:
    cat: str
    c_p: str
    S: float
    X: float
    t: float
    r: float
    sigma: float
    rf: float = field(default=None)
    q: float = field(default=None)

    def __post_init__(self):
        self.n = scipy.stats.norm.pdf
        self.N = scipy.stats.norm.cdf
        self.d1 = (log(self.S/self.X) + (self.r+self.sigma**2/2)*self.t) / (self.sigma*sqrt(self.t))
        self.d2 = self.d1 - self.sigma * sqrt(self.t)
        try:
            if self.cat == 'E':
                self.b = self.r
            elif self.cat == 'F':
                self.b = 0
            elif self.cat == 'C':
                self.b = self.r - self.rf
            elif self.div != None:
                self.b = self.r - self.q
        except TypeError:
            print('ERROR: FX Options Must Include a Foreign Interest-Rate')

        if self.c_p == 'C':
            self.price = self.N(self.d1) * self.S - self.N(self.d2) * self.X * exp(-self.r*self.t)
        elif self.c_p == 'P':
            self.price = self.N(-self.d2) * self.X * exp(-self.r*self.t) - self.N(-self.d1) * self.S

        if self.c_p == 'C':
            self.delta =  exp((self.b-self.r)*self.t)*(self.N(self.d1))
        elif self.c_p == 'P':
            self.delta = exp((self.b-self.r)*self.t)*(self.N(self.d1) - 1)

        self.gamma = exp((self.b-self.r)*self.t)*self.n(self.d1)/self.S*self.sigma*sqrt(self.t)

        if self.c_p == 'C':
            self.theta = -self.S*exp((self.b-self.r)*self.t)*self.n(self.d1)*self.sigma/2*sqrt(self.t) - (self.b-self.r)*self.S*exp((self.b-self.r)*self.t)*self.N(self.d1)-self.r*self.X*exp(-self.r*self.t)*self.N(self.d2)
        elif self.c_p == 'P':
            self.theta = -self.S*exp((self.b-self.r)*self.t)*self.n(self.d1)*self.sigma/2*sqrt(self.t) + (self.b-self.r)*self.S*exp((self.b-self.r)*self.t)*self.N(self.d1)+self.r*self.X*exp(-self.r*self.t)*self.N(-self.d2)

        self.vega = self.S*exp((self.b-self.r)*self.t)*self.b*(self.d1)*sqrt(self.t)

        if self.c_p == 'C' and self.b != 0:
            self.rho = self.t*self.X*exp(-self.r*self.t)*self.N(self.d2)
        elif self.c_p == 'P' and self.b != 0:
            self.rho -self.t*self.X*e(-self.r*self.t)*self.N(-self.d2)
        else:
            self.rho = -self.t*self.price()