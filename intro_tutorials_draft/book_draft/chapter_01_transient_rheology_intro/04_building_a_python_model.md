## Building a Python class inheritance model for the phenomelogical models 

The form of the base creep function lends itself well to an object-oriented design. 
In this chapter, we'll identify the common characteristics of the various 
phenomenological models and construct an object-based data model in Python. 

### Building the base model

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

### Inheriting from the base model

#### Maxwell
So now that we have our base class, we can write versions of it for every
mechanical model. The maxwell model

```python 
class MaxwellModel(MaterialModel):
     
    def J_t(self, t):        
        return self.Ju * (1 + t / self.tau_m)     
    
    def J1_w(self, w):
        pass     
        
    def J2_w(self, w):
        pass 
```

And that's all! 

For andrade, we also need to modify the initialization routine, `__init__` to 
accept the additional $\beta$ and $\alpha$ parameters, but in order to avoid copy/pasting
code, we will call the "parent" or "base" class's `__init__` method by using the `super()`
class call, which will identify the parent class for you:

#### Andrade

```python 
class AndradeModel(MaterialModel):
     
    def __init__(self, Ju, tau_m, beta, alpha):
        super().__init__(Ju, tau_m)
        self.beta = beta 
        self.alpha = alpha 
     
```
Now let's add on our `J_t`, `J1_w` and `J2_w` defitions for `AndradeModel`:
```python 
class AndradeModel(MaterialModel):
     
    def __init__(self, Ju, tau_m, beta=1e-5, alpha=1./3):
        super().__init__(Ju, tau_m)
        self.beta = beta 
        self.alpha = alpha 
        
    def J_t(self, t):        
        return self.Ju + self.beta * t**self.alpha + self.Ju * t / self.tau_m
    
    def J1_w(self, w):
        pass     
        
    def J2_w(self, w):
        pass 
```

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
