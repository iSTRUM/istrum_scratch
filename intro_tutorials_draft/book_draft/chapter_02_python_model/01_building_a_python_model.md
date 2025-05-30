# Building a Python class inheritance model for the phenomenological models 

The form of the base creep function lends itself well to an object-oriented design. 
In this chapter, we'll identify the common characteristics of the various 
phenomenological models and construct an object-based data model in Python. 

## Building the base model

In the previous chapters, we described the following phenomelogical models:

* maxwell 
* andrade 
* burgers
* zener 
* extended burgers 

You'll notice that for all of these, there are common quantities that we set or 
calculated and plotted repeatedly. Looking just at the general time dependent creep 
function, 

$$
J(t) = J_u (1 + \Gamma(t) / J_u + \frac{t}{\tau_m}
$$

We can list a number of properties that **every** model will have. 
First, the unrelaxed modulus, $J_u$, and the maxwell time $\tau_m$ (recognizing 
that some models have multiple maxwell times). And every model calculates 
it $\Gamma(t)$ slightly differently. Furthermore, for all the models we 
want to be able to calculate the $J(t)$ but also the complex compliance, $J*(\omega)$ 
along with all the associated frequency-dependent variables ($J_1(\omega)$, $J_2(\omega)$).

So while we could write functions to do all these calculations for each
model, the problem lends itself very well to an object-oriented design.

What we want is a general idea of a model, the **abstract** base, on top 
of which we superimpose the specific behavior of each mechanical model.

In Python, you can achieve this behavior by creating a class that inherits from
the Abstract Base Class:

```python
import abc 

class MaterialModel(abc.ABC):
    pass
```

Mechanical properties (whether inputs or outputs) that are common to all 
our mechanical models can now be defined in a single location,the base `MaterialModel`,
after which we will add on our model-specific operations. 

To start, we want every model to at least accept two constants: $J_u$ and $\tau_m$:

```python 
import abc 

class MaterialModel(abc.ABC):
    def __init__(self, Ju, tau_m):
        self.Ju = Ju
        self.tau_m = tau_m
```

Next, because each model will implement the time- and frequency-dependent 
calculation differently, we want to define that functionality as an "abstract method".
By using the `@abc.abstractmethod` decorator, it will **require** that our 
mechanical models override and implement the function. We'll call our functions 
`J_t` for time-dependent compliance, `J1_w` and `J2_w` for the frequency-dependent 
storage and loss moduli:

```python 
import abc 
import numpy as np

class MaterialModel(abc.ABC):
    def __init__(self, Ju, tau_m):
        self.Ju = Ju
        self.tau_m = tau_m
        
    @abc.abstractmethod
    def J_t(self, t):
        pass 
    
    @abc.abstractmethod
    def J1_w(self, w):
        pass 
    
    @abc.abstractmethod    
    def J2_w(self, w):
        pass
```

Now we can add on all the definitions for all the calculations that all the mechanical models share.
For example -- regardless of how a method defines the `J_t` function, it will always 
calculate `M_t` in the same way (the inverse of `J_t`). 

```python 
import abc 
import numpy as np

class MaterialModel(abc.ABC):
    def __init__(self, Ju, tau_m):
        self.Ju = Ju
        self.tau_m = tau_m
        
    @abc.abstractmethod
    def J_t(self, t):
        pass 
    
    @abc.abstractmethod
    def J1_w(self, w):
        pass 
    
    @abc.abstractmethod    
    def J2_w(self, w):
        pass 
    
    def M_t(self, t):
        return 1 / self.J_t(t)
    
    def J_w(self, w):
        """The full complex compliance"""
        return np.complex(self.J1_w(w), self.J2_w(w))
        
    def M_w(self, w):
        """The full complex modulus"""
        return 1 / self.J_w(w)
        
    def M1_w(self, w):
        """Real part of the complex modulus"""
        return np.real(self.M_w(w))
    
    def M2_w(self,w):
        """Imaginary part of the complex modulus"""
        return np.imag(self.M_w(w))
    
    def Q_w_approx(self,w):
        """Q with the small Q approximation"""
        return self.J1_w(w) / self.J2_w(w)
    
    def Q_w(self, w):        
        Qfac = 1.0
        return self.Q_w_approx(w) * Qfac
```

## Inheriting from the base model

### Maxwell
So now that we have our base class, we can write versions of it for every
mechanical model. The maxwell model

```python 
class MaxwellModel(MaterialModel):
     
    def J_t(self, t):        
        return self.Ju * (1 + t / self.tau_m)     
    
    def J1_w(self, w):
        return self.Ju

    def J2_w(self, w):
        return self.Ju / (w * self.tau_m)
```

And that's all! 

For andrade, we also need to modify the initialization routine, `__init__` to 
accept the additional $\beta$ and $\alpha$ parameters, but in order to avoid copy/pasting
code, we will call the "parent" or "base" class's `__init__` method by using the `super()`
class call, which will identify the parent class for you:

### Andrade

```python 
class AndradeModel(MaterialModel):
     
    def __init__(self, Ju, tau_m, beta, alpha):
        super().__init__(Ju, tau_m)
        self.beta = beta 
        self.alpha = alpha 
     
```
Now let's add on our `J_t`, `J1_w` and `J2_w` definitions for `AndradeModel`:
```python 
from scipy.special import gamma
import numpy as np

class AndradeModel(MaterialModel):
     
    def __init__(self, Ju, tau_m, beta=1e-5, alpha=1./3):
        super().__init__(Ju, tau_m)
        self.beta = beta 
        self.alpha = alpha 
        
    def J_t(self, t):        
        return self.Ju + self.beta * t**self.alpha + self.Ju * t / self.tau_m
    
    def J1_w(self, w):
        alf = self.alpha
        J_fac = 1 + self.beta * gamma(1+alf) * np.cos(alf * np.pi /2) / (w**alf)
        return self.Ju * J_fac

    def J2_w(self, w):
        alf = self.alpha
        J_fac = 1. / (w * self.tau_m) + self.beta * gamma(1+alf)*np.sin(alf*np.pi/2)/(w**alf)
        return self.Ju * J_fac
```

### Side note on dependencies

In the above `AndradeModel` implementation, we make use of the gamma function implementation 
from `scipy` as well as trig functions from `numpy`. To do so, the import statements are placed above
the class definition -- if we were writing our classes in a module on disk, those imports 
would go up top at the start of the file following standard Python style guidelines, which generally recommend
placing all imports at the start of a file. While it is generally good practice to place your imports at the 
top of your file, there are some cases where it benefits performance or user experience to place imports within 
a function. 

For example, if you expected the above `AndradeModel` to be used infrequently, you could move the scipy import to 
within the functions that need it:

```python 
    def J1_w(self, w):
        from scipy.special import gamma
        alf = self.alpha
        J_fac = 1 + self.beta * gamma(1+alf) * np.cos(alf * np.pi /2) / (w**alf)
        return self.Ju * J_fac

    def J2_w(self, w):
        from scipy.special import gamma
        alf = self.alpha
        J_fac = 1. / (w * self.tau_m) + self.beta * gamma(1+alf)*np.sin(alf*np.pi/2)/(w**alf)
        return self.Ju * J_fac
```
Having the import within the function will delay the import of those libraries until the functionality is actually needed. 
And uses who never use the `AndradeModel` will not incur any of the cost of loading in another library. Additionally, when one 
of these functions is called for the first time, the import will load in the libraries and they will remain imported in the 
user's Python session, so subsequent calls to the function will not have to re-import. 

In addition to potentially saving on load times, this design would allow you to utilize `scipy` as an 
optional dependency, so that users only need to install scipy if they plan on using these functions. To do that you 
would need to adjust your package's metadata and also help the user out a bit -- the simplest approach is to add some 
error handling: 

```python 
    def J1_w(self, w):
        try: 
            from scipy.special import gamma
        except ImportError:
            raise ImportError("This functionality requires scipy, install it with: pip install scipy")
        alf = self.alpha
        J_fac = 1 + self.beta * gamma(1+alf) * np.cos(alf * np.pi /2) / (w**alf)
        return self.Ju * J_fac

    def J2_w(self, w):
        try: 
            from scipy.special import gamma
        except ImportError:
            raise ImportError("This functionality requires scipy, install it with: pip install scipy")
        alf = self.alpha
        J_fac = 1. / (w * self.tau_m) + self.beta * gamma(1+alf)*np.sin(alf*np.pi/2)/(w**alf)
        return self.Ju * J_fac
```
For the remaining examples, we'll keep it simple and keep the imports up top. 

### Using the models 
So let's actually put our models to use!

In this example, we'll assume that you've written the final models in a Pytho module (explain this). 

```python
import material_models as mm 

Ju = 65 * 1e9; 
tau_m = 1e3 * 3600 * 24 * 365
maxwell = mm.MaxwellModel(Ju, tau_m)
andrade = mm.AndradeModel(Ju, tau_m)  # using default beta, alpha

tau_m_2 = tau_m * 0.5
zener = mm.ZenerModel(Ju, tau_m, tau_m_2)
```

### visualizing the Material model inheritance structure 

While the above case is not **too** complicated, some codebases have very 
expansive class structures that are hard to get a handle on if you're new to the 
code. 

```python 
import inheritance_explorer 
```
