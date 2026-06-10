# wrap_xspec

NOT READY FOR USE YET -- STILL IN DEVELOPMENT

Wrap_xspec is a wrapper module to facilitate work with pyXspec for users that are less experienced with regular-xspec and/or command-line tools in general. 

It supports an analysis workflow that emphasis clarity and reproducability through the use of python scripts or python notebooks. 

The main idea is that any modification of the analysis should be done in the python code and the whole notebook can be re-run again. Or at the very least, a notebook cell would be rerun again, instead of using pyXspec in a command-line style, interactive mode. 

Some of the wrapper functions also restrict xspec from second-guessing the user (for example loading a spectrum in a different group than the one specified -- see the design notes below), or having silent failures when e.g. a file does not exist. 

Some functions are also provided to facilitate making graphs with matplotlib. 

Finally, we provide some tutorials on the specific use of the wrap_xspec functionalities, and also more general pyXspec use (to complement the tutorials available at https://github.com/HEASARC/PyXspec-Jupyter-notebooks)



## Installation

Option 1:

You can clone the repository locally, and make sure that the wrap_xspec/wrap_xspec is in your PYTHON_PATH environment variable. 

The package is still in active development: please pull the repro often. 

Option 2:

You can use `pip install "git+https://github.com/veropetit/wrap_xspec`. 

The package is still in active development: please pull the repo often. 

## Requirements

This package assumes that you have xspec and pyXspec installed and acessible in your local environment. 

To verify:
* `xspec` in the command line should load xspec. 
* `python` and then `import xspec` in the command line should be working.
* `import xspec` in a jupyter notebook should be working. 

## Tutorials

Tutorials and their associated data are located in the Tutorial folder. If you are cloning the repository, I suggest that you make a copy of the Tutorial folder somewhere else to try them out (as to not clutter the repository and have to discard you changes when pulling the new version of the package). If you are using the pip method, simply download the Tutorial folder from the repo -- the folder contains all the necessary notebooks and example data. 

## Some notes about the design

This project's goal is to facilitate the use of pyXspec for novice users, who might be familiar with python and python notebooks, but not with regular Xspec and/or command-line tools in general. Therefore, some choices were made here in some of the wrapper functions, that are different than the default pyXspec (or general Xspec)

* In regular Xspec, the response and background can be loaded automatically from the keywords in the header of the spectrum file. In pyXspec, xspec try to load the files listed in the header keywords, even if other files are passed as keywords in the Spectrum() initialization. 

    Then, if for example the data files are in a different directory than the current working directory (ie where the python notebook is), pyXspec will display some error messages when trying to load the files listed in the fits keywords, then will succesffuly load the file passed as argument to the function. Because the output is rather long, and editors like VScode will truncate the output (unless one click on "scroll"), it might look to students that there is a real problem with their codes, when in fact the error displayed is not a real error. 
    
    Furthermore, if a file path/name that does not exist is passed to pyXspec, the python code *does not crash*. This will potentially hide some problems (for example a missing arf or missing background), unless a student really understands the xspec output -- which is difficult for a novice user. 

    Finally, loading multiple spectra in some "groups" is a bit cumbersome (because one has to pass a rather complicated string, and then assign the reponses and background manually -- emulating the command-line way). 

    Therefore, the implementation in this package is to provide a `SpectrumLoader` class that improves readability and reproducability by enforcing strict rules and file checks and by cleaning up the output. 

1. The default pyXspec behavior of the Spectrum() function is to have the `respFile`, `arfFile`, and `backFile` function-keywords set to `'USE_DEFAULT'` (which means, use the spectrum header-keywords). When passed path/files instead in the function-keywords, it loads the header-keywords first (and potentially display an error message if they are not found), then overwrite with the supplied path/files in the function-keywords. If a function-keywords is set to `None`, then again, it loads the header-keywords first, then removes them.

2. In SpectrumLoader, there are only two choices: either the `respFile` function-keyword is a path/file, or it is set to `None` (which means that we don't want a response file to be loaded). If a keyword is *not given*, the default will be `None`. (For the Xspec experts -- if you'd like the keywords-headers to be used, you can simply directly use the xpsec.Spectrum() instead of wrap_xspec.SpectrumLoader())

3. During the initialization of a `SpectrumLoader` object, the File/paths given will be checked, and the python code **will crash** if the files given do not exists (with an appropriate error message, of course). This way, the students *will know* there is a problem.

4. Once the SpectrumLoader object is initialized, then we can run SpectrumLoader.xspec_load() to actually get Xspec to load up the files. The intermediate Xspec output (when it first (try to) load the files in the header-keyword, etc) is captured but not displayed (but it can be by setting the `verbose` keyword to `True`, for debugging). At the end of the loading, SpectrumLoader.xspec_load() will run `xspec.AllData.show()` to display the current status of the all the loaded data. The idea is that the student then see a clean output with the final status, and not all of the garbled output from the intermediaty steps. 

5. It is always good practice to keep research files organized. Which means that increasingly, the spectrum files from various observing runs might be organized in various folders, and the jupyter notebook might be in another. Therefore, it would require the students to type in the full path for each of the axiliary files (response, arf, background) which are usually in the same folder as a given spectrum. Therefore SpectrumLoader has a `path` keyword that can be set. This path will apply to all the file names passes to SpectrumLoader. The paths can be absolute of relative. For example, these would all be valid:
    * `path='Data/obsID123', pi='MyStar.pi', rmf='MyStar.rmf'`
    * `path='$HOME/MyProject/Data/obsID123', pi='MyStar.pi', rmf='MyStar.rmf'`
    * `pi='Data/obsID123/MyStar.pi', rmf='Data/obsID123/MyStar.rmf'`
    * `pi='$HOME/MyProject/Data/obsID123/MyStar.pi', rmf='$HOME/MyProject/Data/obsID123/MyStar.rmf'`

6. Running SpectrumLoader.xspec_load() returns the associated xspec.Spectrum object. 

7. Loading multiple spectra can be done one at the time (each call to a SpectrumLoader.xspec_load()) will add a new spectrum to xspec.AllData(), and return the associated xspec.Spectrum object). The spectra will all be in the same `fitting group` by default. 

    However, it is a bit more complicated to assign spectra to different groups. This has to be done by calling the xspec.AllData() (i.e. it is not possible at the moment to assign a group with a simple call to xspec.Spectrum() -- to have a group, we need to pass a string to the xspec.AllData(). 
    
    The behavior of AllData() may be hard to predict for a novice user. For example:
    ```
    xspec.AllData("1:1 Spectrum1.pi")
    ```
    will load the first spectrum in group 1. Then
    ```
    xspec.AllData("2:2 Spectrum2.pi")
    ```
    will load spectrum 2 in group 2. Then
    ```
    xspec.AllData("1:3 Spectrum3.pi")
    ```
    will add spectrum 3 in group 1. But then
    ```
    xspec.AllData("1:1 New_spectrum1.pi")  
    ```
    will replace spectrum 1 but unload all of the other spectra. 

    The goal of this package is to enable the research done by novice users to be easily reproducable and clear. Therefore, we certainly want to avoid doing too much handling of the loaded spectra inside of a notebook. Afterall, if we want to make a change in our analysis, our workflow would usually be to modify our code and then run the notebook again. 

    Thus, we opt to have three options in SpectrumLoader for handling multiple spectra. 
    1. Multiple separate calls to SpectrumLoader.xspec_load(). This can be used when there is only one fitting group. 
        ```
        s1 = SpectrumLoader("spectrum1.pi").xspec_load()
        s2 = SpectrumLoader("spectrum2.pi").xspec_load()
        s3 = SpectrumLoader("spectrum3.pi").xspec_load()
        ```

    2. We can pass lists to SpectrumLoader instead. It will return a list of spectrum objects, that can be unpacked
        ```
        s1, s2, s3 = SpectrumLoader(["spectrum1", "spectrum2", spectrum3"]).xspec_load()
        ```
        or left unpacked
        ```
        list_s = SpectrumLoader(["spectrum1", "spectrum2", spectrum3"]).xspec_load()
        ```
        the `path`, `rmf`, etc can also be lists, and can include some `None` if e.g. some responses or backgrounds are not to be loaded. 

        Note: with this option, each spectrum is loaded with the "xspec.Spectrum()" command. Therefore each new spectrum in these lists will be *added* to the current xspec.AllData. In other words, the xspec.AllData will not be cleared. 

    3. [NOT IMPLEMENTED YET] The function-keyword `groups` accepts a list of integers that will list the group each spectrum belongs to. SpectrumLoader will then construct the appropriate string and pass it to xspec.AllData. 

        First note: When using the groups keyword, SpectrumLoader will first clear AllData, otherwise things can get a bit messy because AllData might (or might not) remove some of the already loaded spectra. 

        Second note: xspec is doing a lot of second guessing. For example, if one do this call to AllData (remember that the numbers are GroupNumber:SpectrumNumber):
        ```
        xspec.AllData('1:1 Spectrum1Group1.pi 3:2 Spectrum2Group3.pi 2:3 Spectrum3Group2.pi')
        ```

        one would think that there would be 3 groups. But no. Xspec will load Spectrum2Group2.pi into group 2 (with Spectrum3Group2.pi), instead of group 3 as we asked, because group 2 didn't exist. A novice user might not notice or be able to interpret the xspec output and might not know there is a problem because the code didn't crash. 

        Therefore, SpectrumLoader will make a check in the ordering of the Group keyword. If a higher group is asked before its predecessors are asked, then the code will crash and throw an error (and ask the user to reorganize their groups.)

        Yes -- SpectrumLoader could sort out the lists internally under the hood, however then the spectrum index and the group values will not match what the user think they are. It wouldn't be a problem if a user always access the Spectrum characteristics and functions via its SpectrumObject, but sometimes the pyXspec tutorials will show how to access via the integer index instead. So I decided that it would be best to have to throw the error and get the user to correct their input instead. 