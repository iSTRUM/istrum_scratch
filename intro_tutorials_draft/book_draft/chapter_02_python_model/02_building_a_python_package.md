# Building a Python class inheritance model for the phenomenological models

In this chapter, we'll walk through the steps to take our set of phenomenological Python classes and create a Python **package**. 

For more on the details of Python packaging, check out the [Scientific Python Library Development Guide ](https://learn.scientific-python.org/development/).

In this example we'll utilize a number of tools:

* a Python environment manager 
* a Cookiecutter template


In this, using pyenv virtualenv. Could use conda. 

``` 
$ pyenv virtualenv 3.10.12 phenom_dev
$ pyenv activate phenom_dev
```

## cookiecutters

```commandline
$ pip install cookiecutter
$ cookiecutter gh:scientific-python/cookie
```

![](assets/chapter_01/packaging-01-cookies.png)


After initializing with the cookiecutter, you'll have a skeleton of a Python package 
laid out for you, along with a bunch of helpful developer tools initialized.

```commandline
$ cd deformation_models
$ ls 
LICENSE		README.md	docs		noxfile.py	pyproject.toml	src		tests
```

## add dependencies 

in `pyproject.toml`:

``` 
dependencies = ["numpy",
                "scipy",]
```

## adding code 

Open a new file, `src/deformation_models/main.py` and copy in the Python classes for the 
phenomenological models. 

Open up `src/deformation_models/__init__.py` and add 

``` 
from .main import MaxwellModel, AndradeModel, SLS, Burgers
```

Adding this to `__init__.py` makes it so that all the mechanical models will be 
available at the top level of the package API. i.e., you'll be able to do:

```python
import deformation_models 
andrade = deformation_models.Andrade_Model(...) 
```

rather than 

```python
from deformation_models.main import Andrade_Model 
andrade = Andrade_Model(...)  
```

## installing a development version

We're now ready to install a development version! Normally, you could jump right into 
an editable isntall with `pip install -e .`, but remember how we selected Y for that 
"use version control for versioning"? Well when you 
do that, it sets up your `pyproject.toml` to automatically detect the version of your 
package when you run `pip install -e .` by looking at git tags. So since we're 
starting from a fresh directory, let's first initialize our git repository and add 
everything we've got in an initial commit:

```
git init 
git add . 
git commit -m "initial commit" 
```

And now we need to add a git tag. A [git tag](https://git-scm.com/book/en/v2/Git-Basics-Tagging) is 
a way of bookmarking a particular snapshot of code in the git history. And there are a number of tools 
for Python packages that look for tags, and if they find one that looks like a version number, they will 
run certain actions. So let's create our first version tag:

``` 
git tag v0.0.1
```

And now, when we create an editable installation (`pip install -e .`), the version 
will be pulled from the tag. Additionally, let's install the optional dependencies for
running tests:

``` 
pip install -e ".[test]"
```

And we can real quickly run the tests that the cookiecutter came with by running `pytest` (one of the 
optional dependencies for the `test` option, defined in the `pyproject.toml`): 

```commandline 
$ pytest .
============================================================================= test session starts =============================================================================
platform darwin -- Python 3.10.12, pytest-8.3.3, pluggy-1.5.0
rootdir: /Users/chavlin/src/vbr_/istrum/deformation_models
configfile: pyproject.toml
plugins: cov-5.0.0
collected 1 item                                                                                                                                                              

tests/test_package.py .                                                                                                                                                 [100%]

============================================================================== 1 passed in 2.39s ==============================================================================
```

## push up initial code 

Create the github repo, don't add readme or license file, make it public. 

``` 
git remote add origin git@github.com:iSTRUM/deformation_models.git
git branch -M main
git push -u origin main
```

Now ready to push up!

``` 
$ git remote add origin git@github.com:iSTRUM/deformation_models.git
$ git branch -M main
$  git push -u origin main
Enumerating objects: 29, done.
Counting objects: 100% (29/29), done.
Delta compression using up to 12 threads
Compressing objects: 100% (25/25), done.
Writing objects: 100% (29/29), 12.18 KiB | 4.06 MiB/s, done.
Total 29 (delta 0), reused 0 (delta 0), pack-reused 0
To github.com:iSTRUM/deformation_models.git
 * [new branch]      main -> main
branch 'main' set up to track 'origin/main'.
```

After the initial push, you can go to Actions page on the repo on github and you should see 
some tests running.

If you're going to be the only developer, you can stop here but if you hope that your 
project will receive contributions from other developers, you may also want to reconfigure 
your remote repositories. For the above example, I initially created the repository in a github 
organization (iSTRUM), but I don't want to do my own development there directly. So after
pushing the code, I went to github and created a fork at my own personal github account, from which 
I plan to develop and submit pull requests to the main iSTRUM repository. To set that up, 
I could delete my local copy of the repository and clone a fresh copy off my fork, or, you can 
simply adjust where your local git repoistory points to. To do that, first check where 
your local repository points to up on github:

```
$ git remote -v
origin	git@github.com:iSTRUM/deformation_models.git (fetch)
origin	git@github.com:iSTRUM/deformation_models.git (push)
```  

We'll copy that address (`git@github.com:iSTRUM/deformation_models.git`) and use it to 
create a new upstream remote repository:

```
$ git remote add upstream git@github.com:iSTRUM/deformation_models.git 
```

And now, adjust the url for the `origin` to point to your fork:
```
$ git remote set-url origin git@github.com:chrishavlin/deformation_models.git
``` 

Now, we can check our urls: 

``` 
$ git remote -v
origin	git@github.com:chrishavlin/deformation_models.git (fetch)
origin	git@github.com:chrishavlin/deformation_models.git (push)
upstream	git@github.com:iSTRUM/deformation_models.git (fetch)
upstream	git@github.com:iSTRUM/deformation_models.git (push)
```

and we can fetch all the changes from both repositories with 

```
$ git fetch --all
Fetching origin
Fetching upstream
From github.com:iSTRUM/deformation_models
 * [new branch]      main       -> upstream/main
 ```

# publishing to pypi 

The final step of sharing your package is to publish it to a package repository.  