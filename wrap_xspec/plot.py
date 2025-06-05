# Module plot
# This module was inspired by the pyXspec tutorials
# 

import copy
import numpy as np
import xspec

class PlotData:
    '''
    Helper class to hold the data arrays necessary to reproduce an XSPEC Plot Data graph in matplotlib.

    The plot is in the form of data points with error bars in the x and y axis. 
    The format of the plot (log, units, etc) depends on what was setup when loading XSPEC. 
    See the details of the 'get_plot_data()' function
    '''
    def __init__(self, x, y, dx, dy) -> None:
        '''
        Initialization of the plot_data object.

        :param x: The array of energy/wavelength of the datapoints
        :param y: The array of rates/fluxes
        :param dx: The interval between energy/wavelength datapoints
        :param dy: The error bar on the rates/fluxes 
        '''
        self.x=x
        self.y=y
        self.dx=dx
        self.dy=dy
        #self.labels=labels

    def __getitem__(self, key):
        """
        Returns a PlotData object with only the values at the specified index(s)

        :param key: the index or slice being checked
        :rtype: PlotData
        """
        return PlotData(self.x[key], 
                              self.y[key],
                              self.dx[key],
                              self.dy[key])
    
    def __setitem__(self, key, newval):
        """
        Sets all values of the PlotData at the specified location equal
        to the input PlotData's values.

        :param key: the index or slice being overwritten
        :param newval: PlotData whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, PlotData)):
            raise TypeError()
        else:
            self.x[key] = newval.x
            self.y[key] = newval.y
            self.dx[key] = newval.dx
            self.dy[key] = newval.dy

    def __len__(self):
        return len(self.x)

    def save(self, filename):
        '''
        Save a plot_data object to an ascii file. 
        Can be read with the read_plot_data() function. 

        :param filename: the path/filename of the output file. 
        '''
        with open(filename, mode='w') as f:
            f.write('x y dx dy\n')
            for i in range(0,len(self.x)):
                f.write('{} {} {} {}\n'.format(self.x[i], self.y[i], self.dx[i], self.dy[i]))

    def plot_err(self, ax, fmt='.', **kwargs):
        '''
        Function to add a plot_data object to an existing matplotlib `ax` object.

        It will be a datapoint graph (`ax.errorbar`) with error bars in x and y. 

        :param ax: a matplotlib `ax` object
        '''
        ax.errorbar(self.x,self.y,xerr=self.dx,yerr=self.dy,fmt=fmt,**kwargs)
        return ax
    
    def plot_fill_err(self, ax, color='0.5', label='', **kwargs):
        '''
        Function to add a plot_data object to an existing matplotlib `ax` object. Create a fill between value+error bar and value-errorbar

        :param ax: a matplotlib `ax` object
        '''
        for i in range(len(self)):
            if i == 0:
                the_label=label
            else:
                the_label=''
            item = self[i]
            ax.fill_between([item.x-item.dx, item.x+item.dx], [item.y+item.dy]*2,
                            y2=item.y-item.dy,
                            color=color, label=the_label, **kwargs)
        return ax
    
    def plot_fill(self, ax, color='0.5', label='', **kwargs):
        '''
        Function to add a plot_data object to an existing matplotlib `ax` object. Create a fill between value+error bar and value-errorbar

        :param ax: a matplotlib `ax` object
        '''
        for i in range(len(self)):
            if i == 0:
                the_label=label
            else:
                the_label=''
            item = self[i]
            ax.fill_between([item.x-item.dx, item.x+item.dx], [item.y]*2,
                            y2=0,
                            color=color, label=the_label, **kwargs)
        return ax
    
    def plot_step(self, ax, c='k', label='', **kwargs):
        '''
        Function to add a plot_data object to an existing matplotlib `ax` object. Create a fill between value+error bar and value-errorbar

        :param ax: a matplotlib `ax` object
        '''
        for i in range(len(self)):
            if i == 0:
                the_label=label
            else:
                the_label=''
            item = self[i]
            ax.plot([item.x-item.dx, item.x+item.dx], [item.y]*2,
                            c=c, label=the_label, **kwargs)
        return ax
 
    def to_step(self):
        '''
        Function to create a step function out of a plot_data object. Returns a plot_model object
        '''

        # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
        # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
        nE = len(self.energies)
        #stepenergies = list()
        stepenergies = np.array([])
        for i in range(nE):
            stepenergies = np.append(stepenergies, self.energies[i] - self.edeltas[i])
            #stepenergies.append(self.energies[i] - self.edeltas[i])
        stepenergies = np.append(stepenergies, self.energies[-1]+self.edeltas[-1])
        #stepenergies.append(self.energies[-1]+self.edeltas[-1])
        rates = np.append(self.rates, self.rates[-1])

        return plot_model(stepenergies, rates, [])

def get_SpectrumData(s=1, w=1):
    '''
    Helper function to gather the info necessary to make a Plot Data in matplotlib.
    
    This function assumes that a spectrum is loaded in xspec already. 
    The xaxis is set to keV. 
    NOTE: need to check whether I can re-initialize the status of the xspec.Plot 
    to its default values, because this is basically using global variables...
    The yaxis default is count/s/keV, but I need to make sure that it cannot be 
    set in a previous cell... 
    
    (I think that the log status of the xspec plot is irrelevant when accessing the xspec.Plot.x()...)

    :param s: (s=1) the index of the spectrum???
    :param w: (w=1) ???
    :rtype: a `plot_data` object. 
    '''
    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    y = np.array(xspec.Plot.y(s,w))
    dy = np.array(xspec.Plot.yErr(s,w))
    #labels = xspec.Plot.labels()
    return PlotData(x, y, dx, dy)

def get_ModelData(s=1, w=1):
    '''
    Helper function to gather the info necessary to make a Plot Data in matplotlib.
    
    This function assumes that a spectrum is loaded in xspec already. 
    The xaxis is set to keV. 
    NOTE: need to check whether I can re-initialize the status of the xspec.Plot 
    to its default values, because this is basically using global variables...
    The yaxis default is count/s/keV, but I need to make sure that it cannot be 
    set in a previous cell... 
    
    (I think that the log status of the xspec plot is irrelevant when accessing the xspec.Plot.x()...)

    :param s: (s=1) the index of the spectrum???
    :param w: (w=1) ???
    :rtype: a `plot_data` object. 
    '''
    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    y = np.array(xspec.Plot.model(s,w))
    dy = np.zeros(len(x))
    #labels = xspec.Plot.labels()

    return PlotData(x, y, dx, dy)

def get_ModelCompData(s=1, w=1):
    '''
    Helper function to gather the info necessary to make a Plot Data in matplotlib.
    
    This function assumes that a spectrum is loaded in xspec already. 
    The xaxis is set to keV. 
    NOTE: need to check whether I can re-initialize the status of the xspec.Plot 
    to its default values, because this is basically using global variables...
    The yaxis default is count/s/keV, but I need to make sure that it cannot be 
    set in a previous cell... 
    
    (I think that the log status of the xspec plot is irrelevant when accessing the xspec.Plot.x()...)

    :param s: (s=1) the index of the spectrum???
    :param w: (w=1) ???
    :rtype: a `plot_data` object. 
    '''
    ### Should check is Plot.add = True....

    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    dy = np.zeros(len(x))
    #labels = xspec.Plot.labels()
    ncomp = xspec.Plot.nAddComps(s,w)
    comps = []

    if ncomp > 0:
        for i in range(0,ncomp):
            compmodel = xspec.Plot.addComp(i+1,s,w)
            compmodel.append(compmodel[-1])
            comps.append(np.array(compmodel))

    else:
        return comps


def read_plot_data(filename):
    '''
    Function to read in a saved plot_data object. 

    :param filename: the path/filename of the saved plot_data object. 
    :rtype plot_data:
    '''
    data = np.genfromtxt(filename, names=True )
    return PlotData(data['x'], data['y'], data['dx'], data['dy'])

class plot_model:
    def __init__(self, stepenergies, foldedmodel, foldedcomps):
        '''Initialization of a plot_model object
        
        :param stepenergies: The energy/wavelength at the edges of the bins
        :param foldedmodel: The value of the model at the edges of the bins
        :param foldedcomps: The values of the model componenets at the edges of the bins
        '''
        self.stepenergies=stepenergies
        self.foldedmodel=foldedmodel
        self.foldedcomps=foldedcomps

    def save(self, filename):
        '''
        Save a plot_model object to an ascii file. 
        Can be read with the read_plot_model() function. 

        :param filename: the path/filename of the output file. 
        '''
        ncomps = len(self.foldedcomps)
        if ncomps > 0:
            labels = ' '.join( ['comp{}'.format(i+1) for i in range(0,ncomps)] )
        else:
            labels = ''

        with open(filename, mode='w') as f:
            f.write('stepenergies foldedmodel ' + labels + '\n')
            for i in range(0,len(self.stepenergies)):
                row = ' '.join( [ {}.format( self.foldedcomps[x][i]) for x in range(0,ncomps) ] )

                f.write('{} {} {}\n'.format(self.stepenergies[i], self.foldedmodel[i], row))

    def plot(self, ax, comp=None, **kwargs):
        '''
        Function to add a model to an existing matplotlib `ax` object. If `comp` is None (default),
        the total model is plotted as a step function. If `comp` is an integer, it will plot 
        the model component of that index (NOTE: compoenent indexes start at 1). 
        '''

        if comp == None:
            ax.step(self.stepenergies,self.foldedmodel, where='post', **kwargs)
        else:
            if comp <= len(self.foldedcomps):
                ax.step(self.stepenergies,self.foldedcomps[comp-1], where='post', **kwargs)
            else:
                print('Component index {} out of range for model'.format(comp))

        return ax
    
def get_plot_model(s=1,w=1):
    '''
    Helper function to gather the info necessary to make a step Plot Model in matplotlib.
    This function assumes that a spectrum and model are loaded in xspec already. 

    :param s: (=1) ??
    :param w: (=1) ??
    :rtype: a `plot_model` object. 
    '''
    xspec.Plot.device = "/null"
    xspec.Plot.add = True
    xspec.Plot("data")
    energies = xspec.Plot.x(s,w)
    edeltas = xspec.Plot.xErr(s,w)
    foldedmodel = xspec.Plot.model(s,w)
    # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
    # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
    nE = len(energies)
    stepenergies = np.array([])
    for i in range(nE):
        stepenergies = np.append(stepenergies, energies[i] - edeltas[i])
    stepenergies = np.append(stepenergies, energies[-1]+edeltas[-1])
    foldedmodel = np.append(foldedmodel, foldedmodel[-1])
    ncomp = xspec.Plot.nAddComps(s,w)
    comps = []

    if ncomp > 0:
        for i in range(0,ncomp):
            compmodel = xspec.Plot.addComp(i+1,s,w)
            compmodel.append(compmodel[-1])
            comps.append(np.array(compmodel))

    return plot_model(stepenergies, foldedmodel, comps)

def read_plot_model(filename):
    '''
    Function to read in a saved plot_model object. 

    :param filename: the path/filename of the saved plot_model object. 
    :rtype plot_model:
    '''
    data = np.genfromtxt(filename, names=True )
    ncomps = len(data.dtype.names) - 2
    if ncomps > 0:
        comps = [ data['comp{}'.format(i+1)] for i in range(0,ncomps) ]
    else:
        comps = []

    return plot_model(data['stepenergies'], data['foldedmodel'], foldedcomps=comps)
