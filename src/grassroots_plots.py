import json

import numpy as np
from functools import reduce

import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64

import os

import plotly.express as px

def dict_phenotypes(pheno, plots):
    """Extract traits of phenotypes 

    Args:
        pheno: list of phenotypes of a particular study 
        plots     : plots data of a particular study

    Returns:
        dictionary: keys: phenotypes names, values: traits
    """

    names = []
    traits = []
    for key in pheno:
        #print("-->", key)
         

        names.append(key)
        traits.append(pheno[key]['definition']['trait']['so:name'])

    phenoDict = dict(zip(names, traits))

    for j in range(len(plots)):
        if ( 'discard' in plots[j]['rows'][0] ):
            pass
        if ('observations' in plots[j]['rows'][0]):
            for k in range(len(plots[j]['rows'][0]['observations'])):
                if ('raw_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['raw_value']
                if ('corrected_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['corrected_value']
                if ( type(rawValue) == str):
                    name = plots[j]['rows'][0]['observations'][k]['phenotype']['variable']
                    if( name in phenoDict.keys() ):
                        #print("Remove pheonotypes with non-numeric value:", phenoDict[name])
                        del phenoDict[name]
    
    return phenoDict    


###################################################################
def lookup_keys(dictionary, keys, default=None):
     return reduce(lambda d, key: d.get(key, default) if isinstance(d, dict) else default, keys.split("."), dictionary)

###################################################################
def searchPhenotypeTrait(listPheno, value):

    name = listPheno[value]['definition']['trait']['so:name']

    return name

###################################################################
def searchPhenotypeUnit(listPheno, value):

    name = listPheno[value]['definition']['unit']['so:name']

    return name


###################################################################
def search_phenotype(list_observations, value):

    found = False
    for i in range(len(list_observations)):

        dic            = list_observations[i]
        phenotype_name = lookup_keys(dic, 'phenotype.variable')
        if  (phenotype_name == value ):
              return True
              break

    return found

############################################################################
def search_phenotype_index(list_observations, value):

    for i in range(len(list_observations)):

        dic            = list_observations[i]
        phenotype_name = lookup_keys(dic, 'phenotype.variable')
        if  (phenotype_name == value ):
              return i


####################################################################
def dict_phenotypes(pheno, plots):
    """Extract traits of phenotypes 

    Args:
        pheno: list of phenotypes of a particular study 
        plots     : plots data of a particular study

    Returns:
        dictionary: keys: phenotypes names, values: traits
    """

    names = []
    traits = []
    for key in pheno:
        #print("-->", key)
         

        names.append(key)
        traits.append(pheno[key]['definition']['trait']['so:name'])

    phenoDict = dict(zip(names, traits))

    for j in range(len(plots)):
        if ( 'discard' in plots[j]['rows'][0] ):
            pass
        if ('observations' in plots[j]['rows'][0]):
            for k in range(len(plots[j]['rows'][0]['observations'])):
                if ('raw_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['raw_value']
                if ('corrected_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['corrected_value']
                if ( type(rawValue) == str):
                    name = plots[j]['rows'][0]['observations'][k]['phenotype']['variable']
                    if( name in phenoDict.keys() ):
                        # print("Remove:", phenoDict[name])
                        del phenoDict[name]
    
    return phenoDict    




##############--------------------------------##########################
#### new plotly function. Reduce lines of code for jupyter notebook###
def plotly_heatmap(json_study, colormap, index):
    single_study = json.loads(json_study) # "Deserialising" data 

    phenotypes         = single_study['results'][0]['results'][0]['data']['phenotypes']
    phenotype_selected = list(phenotypes.keys())[index]
    arrays             = matrices(single_study, phenotype_selected)

    rows       = arrays[0] 
    columns    = arrays[1]
    raw_values = arrays[2]
    title      = arrays[3]
    units      = arrays[4]
    accession  = arrays[5]

    accession   = accession.reshape(rows,columns)
    plotly_plot(raw_values, accession, title, units, colormap)

##############--------------------------------##########################
#### new seaborn function. Reduce lines of code for jupyter notebook###
def seaborn_heatmap(json_study, colormap, index):

    single_study = json.loads(json_study) # "Deserialising" data 

    phenotypes         = single_study['results'][0]['results'][0]['data']['phenotypes']
    phenotype_selected = list(phenotypes.keys())[index]
    arrays             = matrices(single_study, phenotype_selected)

    rows       = arrays[0] 
    columns    = arrays[1]
    raw_values = arrays[2]
    title      = arrays[3]
    units      = arrays[4]

    matrix   = raw_values.reshape(rows,columns)
    seaborn_plot(matrix, title, units, phenotype_selected, colormap)

##############--------------------------------##########################
########### reduce lines of code for Jupyer notebook  ###########
def matrices(single_study, selected):
    plots      = single_study['results'][0]['results'][0]['data']['plots']
    phenotypes = single_study['results'][0]['results'][0]['data']['phenotypes']
    total_rows = single_study['results'][0]['results'][0]['data']['num_rows']
    total_cols = single_study['results'][0]['results'][0]['data']['num_columns']
    traits     = dict_phenotypes(phenotypes, plots)

    matrices  = create_matrices(plots, phenotypes, selected, total_rows, total_cols)
    return matrices

####################################################################
def create_matrices(json, pheno, current_name, total_rows, total_columns):
    """create numpy matrices for plotting

    Args:
        json     : Plots data of a particular study
        pheno    : Phenotypes of particular study
        name     : Name of current study

    Returns:
        matrices: matrix with numpy matrices...
    """


    traitName = searchPhenotypeTrait(pheno, current_name)
    unit      = searchPhenotypeUnit( pheno, current_name)

    dtID= np.dtype(('U', 4))

    row_raw   = np.array([])
    matrix    = np.array([])
    row_acc   = np.array([])
    accession = np.array([])
    plotsIds  = np.array([], dtype=dtID)  #format of strings

    matrices = []

    num_columns = 1
    row    = 1
    column = 1
    #loop throght observations in the same fashion as in old JS code. 
    for j in range(len(json)):
        if ( int( json[j]['row_index'] ) == row ):
            if  (int( json[j]['column_index'] ) == column):
               if column > num_columns:
                   num_columns = column

               if   ( 'discard' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  # use NaN for discarded plots
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
               elif ( 'blank' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  # use NaN for discarded plots
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
      
               elif ( 'observations' in json[j]['rows'][0] ):
                    if( search_phenotype(json[j]['rows'][0]['observations'], current_name) ):
                        indexCurrentPhenotype = search_phenotype_index (json[j]['rows'][0]['observations'], current_name)
                        if ('raw_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                        if ('corrected_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):    
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                        row_raw  = np.append(row_raw, rawValue) 
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession']) 
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
                    else:
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
               else:
                    if ( 'rows' in json[j] ):
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
         
  
               column+=1
               columns = json[j]['column_index']#

        elif ( int( json[j]['row_index'] ) > row  ):
            if column > num_columns:
                   num_columns = column

            if   ( 'discard' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
            elif   ( 'blank' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
        
            elif ( 'observations' in json[j]['rows'][0] ):
                    if( search_phenotype(json[j]['rows'][0]['observations'], current_name) ):
                        indexCurrentPhenotype = search_phenotype_index (json[j]['rows'][0]['observations'], current_name)
                        if ('raw_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                        if ('corrected_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):    
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                        row_raw  = np.append(row_raw, rawValue) 
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession']) 
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
                    else:
                        row_raw  = np.append(row_raw, np.inf )
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
            else:
                    if ( 'rows' in json[j] ):
                        ##print("rows with no observations------",json[j])
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
             

            row+=1
            column=2
            columns = json[j]['column_index']


    #column = columns # use actual number of columns instead of counter
    column = num_columns-1

    if column<columns:
        column=columns
    
    #######print("number of plots and shape check", len(json), row, column, row*(column) )
    if (len(json) != row*column):
        #print("NOT rectangular")
        if(total_columns!=None):
          if(column<total_columns):
             column=total_columns

        # fit odd shape plot into bigger rectangular plot.
        row_raw  = oddShapeValues(   json, row, column, current_name)
        row_acc  = oddShapeAccession(json, row, column, current_name)
        plotsIds = oddShapePlotID(   json, row, column, current_name)

    matrices.append(row)
    matrices.append(column)
    matrices.append(row_raw)
    #matrices.append(row_acc)
    matrices.append(traitName)
    matrices.append(unit)
    matrices.append(row_acc)
    matrices.append(plotsIds)
    
    return matrices

####################################################################################
def oddShapeValues(arraysJson, rows, columns, phenotype):

    matrix = np.zeros((rows,columns))
    matrix[:] = np.nan

    for r in range(len(arraysJson)):
        if  ( 'discard' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = np.nan
        elif  ( 'blank' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = np.nan

        elif ( 'observations' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            if( search_phenotype(arraysJson[r]['rows'][0]['observations'], phenotype) ):
                indexCurrentPhenotype = search_phenotype_index (arraysJson[r]['rows'][0]['observations'], phenotype)
                if ('raw_value' in arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]):
                    rawValue = arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                if ('corrected_value' in arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]):
                    rawValue = arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                matrix[i][j] = rawValue
            else:
                matrix[i][j] = np.inf

        else:
            if('rows' in arraysJson[r]):        #rows field exists but it has no observations!
               i = int( arraysJson[r]['row_index']    )
               j = int( arraysJson[r]['column_index'] )
               i=i-1
               j=j-1
               matrix[i][j] = np.inf  # consider it N/A instead as default discarded (nan)
    

    #matrix = np.flipud(matrix)
    #print(matrix)
    matrix  = matrix.flatten()

    return matrix
#######################################################################
def oddShapeAccession(arraysJson, rows, columns, phenotype):

    dt= np.dtype(('U', 50)) # define string type for of strings (accession names)
    matrix = np.empty((rows,columns), dtype=dt)
    matrix[:] = 'Discarded' #  hovering text in empty plots

    for r in range(len(arraysJson)):
        if  ( 'discard' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = np.nan         #discarded plot
        elif  ( 'blank' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = np.nan         #discarded plot

        elif ( 'observations' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
         #   if( search_phenotype(arraysJson[r]['rows'][0]['observations'], phenotype) ):
            matrix[i][j] = arraysJson[r]['rows'][0]['material']['accession']
        elif('rows' in arraysJson[r]):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = arraysJson[r]['rows'][0]['material']['accession']


    matrix  = matrix.flatten()

    return matrix


#######################################################################
def oddShapePlotID(arraysJson, rows, columns, phenotype):

    dt= np.dtype(('U', 40))
    matrix = np.empty((rows,columns), dtype=dt)
    matrix[:] = 'N/A'

    for r in range(len(arraysJson)):
        if  ( 'discard' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = arraysJson[r]['rows'][0]['study_index']

        elif ( 'observations' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            if( search_phenotype(arraysJson[r]['rows'][0]['observations'], phenotype) ):
                matrix[i][j] = arraysJson[r]['rows'][0]['study_index']
            else:
            #    matrix[i][j] = np.nan    # No values for that phenotype
                matrix[i][j] = arraysJson[r]['rows'][0]['study_index']



    matrix  = matrix.flatten()

    return matrix

#####################################################################################################
'''
test rendering seaborn image
'''
#def seaborn_plot(numpy_matrix, title, unit, uuid, name):
#def seaborn_plot(numpy_matrix, title, unit, name):
def seaborn_plot(numpy_matrix, title, unit, name, color_map):

    sns.set(rc={'figure.figsize':(15.5,6.5)})

    numpy_matrix = np.flipud(numpy_matrix)      # To Match order shown originally in JS code
    notAvailable = np.zeros_like(numpy_matrix)
    discarded    =  np.zeros_like(numpy_matrix)
    indexInf     =  np.where(np.isinf(numpy_matrix))
    indexDiscard =  np.where(np.isnan(numpy_matrix))
    notAvailable[indexInf]   = 1
    discarded[indexDiscard]  = 1
    NA        = np.where(notAvailable < 1, np.nan, notAvailable)
    discarded = np.where(   discarded < 1, np.nan, discarded)
    units = 'Units: '+ unit

    numpy_matrix[indexInf] = np.nan # Replace Inf by NaN

    # Reverse Y ticks and start them from 1
    size  = numpy_matrix.shape
    Y     = size[0]
    Yvals = np.arange(0.5, Y+0.5, 1.0)
    Yaxis = np.arange(1,Y+1)
    Yaxis = np.flip(Yaxis)

    X     = size[1]
    Xvals = np.arange(0, X)
    Xaxis = np.arange(1,X+1)


    maxVal = np.nanmax(numpy_matrix)
    minVal = np.nanmin(numpy_matrix)
    #print(minVal)
    colormap  = sns.light_palette(color_map, as_cmap=True)
    dark      = sns.dark_palette((260, 75, 60), input="husl")
    sns.heatmap(NA, linewidth=0.5,cmap=dark, cbar=False )

    g = sns.heatmap(numpy_matrix,  vmax=maxVal, vmin=minVal,linewidth=0.5,cmap=colormap, cbar_kws={'label': units}) 
    ##g.set_facecolor('xkcd:black')

    g.patch.set_facecolor('white')
    g.patch.set_edgecolor('black')
    g.patch.set_hatch('xx')


    g.set_xlabel("Columns", fontsize = 14)
    g.set_ylabel("Rows", fontsize = 14)
    g.set_title(title, fontsize = 20)
    g.set_yticks( Yvals )
    g.set_yticklabels(Yaxis, size=10)
    g.tick_params(    axis='y', rotation=0)
    g.set_xticks(Xvals)
    g.set_xticklabels(Xaxis, size=10)
   
##############################################################################################
'''
test rendering plotly interactive heatmap
'''
#def plotly_plot(numpy_matrix, accession, title, unit, IDs, treatments):
def plotly_plot(numpy_matrix, accession, title, unit, colormap):
    #colormap = "Hot"
    ##numpy_matrix = np.flipud(numpy_matrix)      # To Match order shown originally in JS code
    #plotID      = np.flipud(IDs)        
    size = accession.shape
    Y    = size[0]
    X    = size[1]

    indexInf     =  np.where(np.isinf(numpy_matrix))
    indexDiscard =  np.where(np.isnan(numpy_matrix))
    
    numpy_matrix[indexInf] = np.nan # Replace Inf by NaN

    strings     = np.array(["%s" % x for x in numpy_matrix])  #matrix has to be flattened for conversion to strings

    for i in range(len(strings)):   #remove decimal place when floats are integers
        string1 = strings[i]
        string_split = string1.split('.')
        if( len(string_split)==2):
            if(string_split[1]=='0'):
                integer = strings[i]
                integer = integer[:-2]
                strings[i] = integer
    

    strings[indexInf]     = 'N/A'    # use array of string for custumising hovering text
    strings[indexDiscard] = 'N/A'

    accession = accession.flatten()
    accession[indexDiscard] = 'Discarded'
    accession = accession.reshape(Y,X)                       

    s_matrix = strings.reshape(Y,X)                  
    s_matrix = np.flipud(s_matrix)      

    numpy_matrix[indexInf] = np.inf # but back inf (N/A values)
    numpy_matrix           = numpy_matrix.reshape(Y,X) 


    numpy_matrix = np.flipud(numpy_matrix)  # For matching order of JS table
    accession    = np.flipud(accession)        
    #plotID       = np.flipud(IDs)        

    # Reverse Y ticks and start them from 1
    Yvals = np.arange(0,Y)
    Yaxis = np.arange(1,Y+1)
    Yaxis = np.flip(Yaxis)
    
    Xvals = np.arange(0, X)
    Xaxis = np.arange(1,X+1)

    units = 'Units: '+unit
    color = "px.colors.sequential."+colormap
    #print(color)
    #CM=eval("px.colors.sequential.Greens")
    CM=eval(color)
    #print(CM)
    fig = px.imshow(numpy_matrix, aspect="auto",
            labels=dict(x="columns", y="rows", color=units),
            #color_continuous_scale=px.colors.sequential.Hot, height=800 )
            color_continuous_scale=CM, height=600 )
    #else:
    fig.update_traces(
    #customdata = np.moveaxis([accession, s_matrix, plotID], 0,-1),
    customdata  = np.moveaxis([accession, s_matrix], 0,-1),
    hovertemplate="Accession: %{customdata[0]}<br>Raw value: %{customdata[1]}<br> (column: %{x}, row:%{y})<extra></extra>")
    fig.update_layout(font=dict(family="Courier New, monospace",size=12,color="Black"),title={
    'text': title,
    'y':0.98,'x':0.5,
    'xanchor': 'center','yanchor': 'top'})

    fig.update_layout( yaxis = dict(tickmode = 'array', tickvals = Yvals, ticktext = Yaxis ) )
    fig.update_layout( xaxis = dict(tickmode = 'array', tickvals = Xvals, ticktext = Xaxis ) )


    fig.update_xaxes(showgrid=True, gridwidth=7, gridcolor='Black', zeroline=False)
    fig.update_yaxes(showgrid=True, gridwidth=7, gridcolor='Black', zeroline=False)
    fig['layout'].update(plot_bgcolor='black')

    #plot_div = plotlyOffline(fig, output_type='div')
    fig.show()   ## ADDED ONLY FOR DASH TEST
    
    #return fig

