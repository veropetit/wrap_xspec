# Module model


def model_to_dict(m, custom_row_label=None):
    '''
    Function to create a dictionary from the current model values.
    A custom label can be added in the first column named "Model". 
    '''
    # create a dictinary object to store the values
    if custom_row_label is None:
        d = {'Model':m.expression}
    else:
        d = {'Model':custom_row_label}
    CompNames = m.componentNames  
    for n in CompNames:
        comp = getattr(m, n)
        ParamNames = comp.parameterNames
        for p in ParamNames:
            val = getattr(comp,p)
            d['{}:{} {}'.format(n,p,val.unit)]=val.values[0]
            d['{}:{} {} err'.format(n,p,val.unit)]=val.sigma if val.frozen == False else 'frozen'
    return d