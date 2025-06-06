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

    def split(self):
        '''
        Checking for ignored chanels and returns a list PlotData
        '''
        right = self[:-1].x+self[:-1].dx
        left = self[1:].x-self[1:].dx
        diff = np.abs(right-left) > 0
        split_indices = np.where(diff)[0] + 1  # Add 1 to split *after* the True
        # Use np.split to split the array at those indices
        
        nseg = len(split_indices)+1
        #print('{} segments'.format(nseg))
        split_indices = np.hstack([[0],split_indices,[len(diff)]])
        seg = []
        for i in range(0,nseg):
            seg.append(self[split_indices[i]:split_indices[i+1]])
        return(seg)

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
    
    def plot_step_h(self, ax, c='k', label='', **kwargs):
        '''
        Function to add a plot_data object to an existing matplotlib `ax` object. 
        Create a fill between value+error bar and value-errorbar

        :param ax: a matplotlib `ax` object
        '''
        split = self.split()
        for j, s in enumerate(split):
            if j == 0:
                the_label=label
            else:
                the_label=''

            nx = len(s.x)
            step_x = np.array([])
            step_y = np.array([])
            for i in range(nx):
                step_x = np.append(step_x, s.x[i] - s.dx[i])
                step_y = np.append(step_y, s.y[i])
                step_x = np.append(step_x, s.x[i] + s.dx[i])
                step_y = np.append(step_y, s.y[i])

            ax.plot(step_x, step_y, c=c, label=the_label, **kwargs)
        return ax

def get_SpectrumData(spec, w=1):
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
    s = spec.index

    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    y = np.array(xspec.Plot.y(s,w))
    dy = np.array(xspec.Plot.yErr(s,w))
    #labels = xspec.Plot.labels()
    return PlotData(x, y, dx, dy)

def get_ModelData(spec, w=1):
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
    s = spec.index
    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    y = np.array(xspec.Plot.model(s,w))
    dy = np.zeros(len(x))
    #labels = xspec.Plot.labels()

    return PlotData(x, y, dx, dy)

def get_ModelCompData(spec, w=1):
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
    #print(xspec.Plot.add)
    if xspec.Plot.add == False:
        print('Plot.add must be set to True')
        return([])

    s = spec.index

    # get the info for the plotting
    x = np.array(xspec.Plot.x(s,w))
    dx = np.array(xspec.Plot.xErr(s,w))
    dy = np.zeros(len(x))
    ncomp = xspec.Plot.nAddComps(s,w)
    #print(ncomp)
    comps = []

    if ncomp > 0:
        for i in range(0,ncomp):
            compmodel = xspec.Plot.addComp(i+1,s,w)
            comps.append(PlotData(x,np.array(compmodel),dx,dy))

    return comps

def read_PlotData(filename):
    '''
    Function to read in a saved PlotData object. 

    :param filename: the path/filename of the saved PlotData object. 
    :rtype PlotData:
    '''
    data = np.genfromtxt(filename, names=True )
    return PlotData(data['x'], data['y'], data['dx'], data['dy'])

