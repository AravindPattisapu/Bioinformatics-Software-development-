####### Dependencies #######
# Install flask using command line before: pip install flask
# Install flask_wtf using command line: pip install flask-wtf
# Install plotly using command line: pip install plotly

# Import flask
# Import neccessary library required from flask
from flask import Flask, render_template, request, redirect, url_for, flash, session, send_file

# Pandas for handling and data/dataframe
import pandas as pd

# sqlite3 for accessing pandas
import sqlite3

# Import forms from form.py file
from forms import SNPForm

# IO for input/output operations
import io

# import functions to calculate genotype frequency and allele frequency function's
from GF_AF_functions import gtfreq, allefreq

# importing the functions for stats
from stats import *

# importing search and query functions
from Validate_Search_Functions import *

# imports json
import json

# importing graph libraries
import plotly
import plotly.express as px

# importing functions to produce multiple plotly graphs
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# importing combinations to allow all pairwise combination calculations of populations for FST 
from itertools import combinations

################################################################################################################

# Create a flask application object
app = Flask(__name__, template_folder='templates')

# A secret key is needed to create and process forms(inputs)
app.config['SECRET_KEY'] = 'f466c2ef7b41e5c65e9f'


# Link to the main page - base.html.
# Methods GET and POST allow the user to input information from forms.
@app.route("/", methods=['GET', 'POST'])
def basepage():
    return render_template('base.html')


# This will render our SNP browser page.
@app.route('/home', methods=['GET', 'POST'])
def home():
    # reset session
    session.clear()
    # get the form
    form = SNPForm()
    # return the snp browser template
    return render_template('home_template.html',
                           text="JIDA SNP browser",
                           form=form)  # formatted using a separate html template file



# Takes the input from forms.py and home_template and redirects the user

@app.route("/search", methods=['GET', 'POST'])

########################        Search Functions       ##########################

def search():

    # do the calculation only if method is POST
    if request.method == 'POST':

        # Connects to the database
        conn = sqlite3.connect('SNPB (1).db')

        # creates a cursor, something that lets us interact with the db in sql language
        cur = conn.cursor()

        # Assigns the user input to variables
        gene = str(request.form['gene'])
        rs = str(request.form['rs'])
        start = request.form['start']
        end = request.form['end']

        # Runs the inputs through a function that filters if the input is valid
        search = search_term(gene, rs, start, end)

        # If the input is valid, the inputs are assigned to the an appropriate query(There are 3 separate ones for gene, rs and positions)
        if search != None:
            query = query_select(search[0], search[1], gene, rs, start, end)

            # Creates a pandas dataframe form the sql query
            sql_query = pd.read_sql_query(query, conn)
 
            df = pd.DataFrame(sql_query, columns=['chromosome', 'position', 'rs_value', 'population',
                                                  'reference',
                                                  'alternate',
                                                  'sample_count', '0|0', '0|1', '1|0', '1|1'])
        
        #Checking/ coverting position values to integer 
        df["position"] = df["position"].astype(int)

        # creating empty columns for calculations
        df['af-0'] = ""
        df['af-1'] = ""
        df['gt-00'] = ""
        df['gt-01'] = ""
        df['gt-10'] = ""
        df['gt-11'] = ""

        # iter over each row and calculate af and gf values
        for i, row in df.iterrows():

            gt = gtfreq(row['0|0'], row['0|1'], row['1|0'], row['1|1'])
            af = allefreq(row['0|0'], row['0|1'], row['1|0'], row['1|1'])

            # converting the values from dict into a list
            ls = list(gt.values())
            ls2 = list(af.values())

            # storing the values in specified location
            df['gt-00'].iloc[i] = ls[0]
            df['gt-01'].iloc[i] = ls[1]
            df['gt-10'].iloc[i] = ls[2]
            df['gt-11'].iloc[i] = ls[3]

            df['af-0'].iloc[i] = ls2[0]
            df['af-1'].iloc[i] = ls2[1]

        # getting the user input data from html file home_template for population checkboxex
        population_list = []
        population_names = []
        if 'GBR' in request.form:
            population_list.append('British')
            population_names.append('British')
        if 'COL' in request.form:
            population_list.append('Colombian')
            population_names.append('Colombian')
        if 'FIN' in request.form:
            population_list.append('Finnish')
            population_names.append('Finnish')
        if 'PUN' in request.form:
            population_list.append('Punjabi')
            population_names.append('Punjabi')
        if 'TEL' in request.form:
            population_list.append('Telugu')
            population_names.append('Telugu')

        # getting the statistics in a list acording to user's input
        # i.e if user clicks the check box on home_template
        stats_list = []
        if 'nud' in request.form:
            stats_list.append('nud')
        if 'hapd' in request.form:
            stats_list.append('hapd')
        if 'td' in request.form:
            stats_list.append('td')
        if 'fst' in request.form:
            stats_list.append('fst')
        
        #dropping the columns and storing in a new variable that can be later used to while rendering output
        DF = df
        DF = DF.drop(['sample_count', '0|0', '0|1', '1|0', '1|1'], axis=1)

        # selecting only specific rows acording to input
        ind = pd.DataFrame([], columns=DF.columns.values)

        # indexing the list Populaion_list ; i - (indexes (1, 2 , 3)) ;j - population list
        for i, j in enumerate(population_list):

            # to show the populations one by one
            d = DF[DF["population"] == j]

            # concating the data into ind from d
            ind = pd.concat([ind, d], axis=0)

        # Show the summary statistics only when the user enters genome and gene name
        if (request.form['start'] != "" and request.form['end'] != "") or (request.form['gene'] != ""):

            # p_names will contain -> all the data for selected population
            p_names = {}

            # zipping 2 lists together 
            for pl, pn in zip(population_list, population_names):
                
                p_names[pn] = df[df['population'] == pl]

            #Booleans- later made true acording to the user input 
            check_nud = False
            check_hapd = False
            check_td = False
            check_fst = False


            #################### Statistics functions #####################

            # key -> E.g. : df_british (name)
            # val -> all the rows for british pop only


            # Taking the input for window size (default = 3)
            try:
                Window = request.form['Window']
                Window = int(Window)
            except:
                Window = 3

            # at this point the df has one row per population
            number_of_pops = len(population_list)
            number_of_variants = len(df) / number_of_pops

            if 'nud' in request.form:
                if number_of_variants >= 3:
                    nud = {}
                    nud_w = {}
                    for key, val in p_names.items():
                        nud[key] = SQLtoNucDiv(val, val.iloc[0]["position"], val.iloc[-1]["position"])
                        nud_w[key] = nuc_div_sliding(val, Window) 
                else:
                    nud = {}
                    for key, val in p_names.items():
                        nud[key] = SQLtoNucDiv(val, val.iloc[0]["position"], val.iloc[-1]["position"])
                check_nud = True

            if 'hapd' in request.form:
                if number_of_variants >= 3:
                    hapd = {}
                    hapd_w = {}
                    for key, val in p_names.items():
                        hapd[key] = SQLtoHapDiv(val)  
                        hapd_w[key] = SQLtoHapDiv_window(val, Window).tolist()
                else:
                    hapd = {}
                    for key, val in p_names.items():
                        hapd[key] = SQLtoHapDiv(val)  
                check_hapd = True

            if 'td' in request.form:
                if number_of_variants >= 3:
                    td = {}
                    td_w = {}
                    for key, val in p_names.items():
                        td[key] = SQLtoTD(val) 
                        td_w[key] = SQLtoTD_window(val, Window).tolist()  
                else:
                    td = {}
                    for key, val in p_names.items():
                        td[key] = SQLtoTD(val)  
                check_td = True

            if 'fst' in request.form:
                combos = list(combinations(population_list, 2))
                fst = {}
                for element in combos:
                    pop1 = element[0]
                    pop2 = element[1]
                    key = str("FST - " + pop1 + " - " + pop2)
                    fst[key] = SQLtoFST(df[df["population"] == pop1], df[df["population"] == pop2])
                    check_fst = True
               
                if number_of_variants >= 3:
                    fst_w = {}
                    for element in combos:
                        pop1 = element[0]
                        pop2 = element[1]
                        key = str("FST - " + pop1 + " - " + pop2)
                        fst_w[key] = SQLtoFST_window(df[df["population"] == pop1], df[df["population"] == pop2],
                                                     Window).tolist()
                check_fst = True

            # final_stat contains all the values calculated  by stat functions
            final_stat = {}
            if check_nud:
                final_stat['Nucleotide diversity'] = nud
                if number_of_variants >= 3:
                    final_stat['Nucleotide Diversity Window'] = nud_w
            if check_hapd:
                final_stat['Haplotype Diversity'] = hapd
                if number_of_variants >= 3:
                    final_stat['Haplotype Diversity Window'] = hapd_w
            if check_td:
                final_stat['Tajimas D'] = td
                if number_of_variants >= 3:
                    final_stat['Tajimas D Window'] = td_w
            if check_fst:
                final_stat['FST'] = fst
                if number_of_variants >= 3:
                    final_stat['FST Window'] = fst_w
            


            #converting final_stat into a global variable 
            #global variables can be accessed throughout the program body by all functions
            global final_stat_global
            final_stat_global = final_stat 

        else: 
            final_stat_global = {}

           
        # if data is empty than return simple empty arrray (for safety)
        if df.empty:
            return render_template('output.html', d=ind, zip=zip)

        # final return.
        return render_template('output.html', column_names= DF.columns.values, d=list(ind.values.tolist()),
                                final_stat=final_stat_global, col_name = population_names ,zip=zip)


########################################    App route for graphs    ###############################################


############################################### APP ROUTE FOR GRAPHS ###################################################

all_pop = ["British", "Colombian", "Finnish", "Punjabi", "Telugu"]
# muted blue, safety orange, cooked asparagus green, brick red, muted purple
colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
pop_colour = dict(zip(all_pop, colours))

############################################### NUCLEOTIDE DIVERSITY ###################################################

@app.route("/Nucleotide diversity")
def Nucleotide_graph():

    name = []
    val = []
    for i, j in final_stat_global["Nucleotide diversity"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["name"] = name
    df["val"] = val

    final_colours = []
    for n in range(len(name)):
        if name[n] in pop_colour.keys():
            final_colours.append(pop_colour[name[n]])
    print(final_colours)

    fig = px.scatter(df, x="name", y="val", color=final_colours)

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header = "Nucleotide diversity graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


@app.route("/Nucleotide Diversity Window")
def Nucleotide_w_graph():

    # making empty lists to store the ND window values and the population name
    ND_window_list = []
    name = []

    # from the stats dictionary extract the population and the ND window values for ND and append these to the lists
    for population, ND_window_values in final_stat_global["Nucleotide Diversity Window"].items():
        name.append(population)
        ND_window_list.append(ND_window_values)

    # creates the list of colours for the populations that the user has selected
    final_colours = []
    for n in range(len(name)):
        if name[n] in pop_colour.keys():
            final_colours.append(pop_colour[name[n]])
    print(final_colours)

    # if the user has only selected multiple populations then output multiple graphs and an overall graph
    if len(ND_window_list) > 1:
        fig = make_subplots(rows=1, cols=1)

        # all the graphs
        n = 0
        for ND_values in ND_window_list:
            x_window = list(range(1, len(ND_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=ND_values, name=name[n], mode='markers', marker=dict(color=final_colours[n])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Nucleotide Diversity", row=1, col=1)
            n += 1

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "ND graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

    # if the user has only selected one population only that specific graph will be present
    else:
        fig = make_subplots(rows=(len(ND_window_list)), cols=1, subplot_titles=tuple(name))

        for ND_values in ND_window_list:
            x_window = list(range(1, len(ND_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=ND_values, name=name[0], mode='markers', marker=dict(color=final_colours[0])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Nucleotide Diversity", row=1, col=1)

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "Nucleotide Diversity graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

################################################# HAPLOTYPE DIVERSITY #################################################

@app.route("/Haplotype Diversity")
def Haplotype_graph():

    name = []
    val = []
    for i, j in final_stat_global["Haplotype Diversity"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["Population"] = name
    df["Value"] = val
    fig = px.scatter(df, x="Population", y="Value")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header = "Haplotype Diversity graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


@app.route("/Haplotype Diversity Window")
def Haplotype_w_graph():

    # making empty lists to store the ND window values and the population name
    HD_window_list = []
    name = []

    # from the stats dictionary extract the population and the ND window values for ND and append these to the lists
    for population, HD_window_values in final_stat_global["Haplotype Diversity Window"].items():
        name.append(population)
        HD_window_list.append(HD_window_values)

    # creates the list of colours for the populations that the user has selected
    final_colours = []
    for n in range(len(name)):
        if name[n] in pop_colour.keys():
            final_colours.append(pop_colour[name[n]])
    print(final_colours)

    # if the user has only selected multiple populations then output multiple graphs and an overall graph
    if len(HD_window_list) > 1:
        fig = make_subplots(rows=1, cols=1)

        # all the graphs
        n = 0
        for HD_values in HD_window_list:
            x_window = list(range(1, len(HD_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=HD_values, name=name[n], mode='markers', marker=dict(color=final_colours[n])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Haplotype Diversity", row=1, col=1)
            n += 1

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "HD graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

    # if the user has only selected one population only that specific graph will be present
    else:
        fig = make_subplots(rows=(len(HD_window_list)), cols=1, subplot_titles=tuple(name))

        for HD_values in HD_window_list:
            x_window = list(range(1, len(HD_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=HD_values, name=name[0], mode='markers', marker=dict(color=final_colours[0])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Haplotype Diversity", row=1, col=1)

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "Haplotype Diversity graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

###################################################### TAJIMAS D ######################################################

@app.route("/Tajimas D")
def Tajimas_graph():

    name = []
    val = []
    for i, j in final_stat_global["Tajimas D"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["Population"] = name
    df["Value"] = val
    fig = px.scatter(df, x="Population", y="Value")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header = "Tajima's D graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


@app.route("/Tajimas D Window")
def Tajimas_w_graph():

    # making empty lists to store the ND window values and the population name
    TD_window_list = []
    name = []

    # from the stats dictionary extract the population and the ND window values for ND and append these to the lists
    for population, TD_window_values in final_stat_global["Tajimas D Window"].items():
        name.append(population)
        TD_window_list.append(TD_window_values)

    # creates the list of colours for the populations that the user has selected
    final_colours = []
    for n in range(len(name)):
        if name[n] in pop_colour.keys():
            final_colours.append(pop_colour[name[n]])
    print(final_colours)

    # if the user has only selected multiple populations then output multiple graphs and an overall graph
    if len(TD_window_list) > 1:
        fig = make_subplots(rows=1, cols=1)

        # all the graphs
        n = 0
        for TD_values in TD_window_list:
            x_window = list(range(1, len(TD_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=TD_values, name=name[n], mode='markers', marker=dict(color=final_colours[n])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Tajima's D", row=1, col=1)
            n += 1

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "TD graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

    # if the user has only selected one population only that specific graph will be present
    else:
        fig = make_subplots(rows=(len(TD_window_list)), cols=1, subplot_titles=tuple(name))

        for TD_values in TD_window_list:
            x_window = list(range(1, len(TD_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=TD_values, name=name[0], mode='markers', marker=dict(color=final_colours[0])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="Tajima's D", row=1, col=1)

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "Tajima's D graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


###################################################### FST #############################################################

@app.route("/FST")
def FST_graph():

    name = []
    val = []
    for i, j in final_stat_global["FST"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["Population"] = name
    df["Value"] = val
    fig = px.scatter(df, x="Population", y="Value")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header = "FST graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


@app.route("/FST Window")
def FST_w_graph():

    # making empty lists to store the ND window values and the population name
    FST_window_list = []
    name = []

    # from the stats dictionary extract the population and the ND window values for ND and append these to the lists
    for population, FST_window_values in final_stat_global["FST Window"].items():
        name.append(population)
        FST_window_list.append(FST_window_values)

    # if the user has only selected multiple populations then output multiple graphs and an overall graph
    if len(FST_window_list) > 1:
        fig = make_subplots(rows=1, cols=1)

        # all the graphs
        n = 0
        for FST_values in FST_window_list:
            x_window = list(range(1, len(FST_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=FST_values, name=name[n], mode='markers'), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="FST", row=1, col=1)
            n += 1

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "FST graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)

    # if the user has only selected one population only that specific graph will be present
    else:
        fig = make_subplots(rows=(len(FST_window_list)), cols=1, subplot_titles=tuple(name))

        for FST_values in FST_window_list:
            x_window = list(range(1, len(FST_values) + 1))
            fig.add_trace(go.Scatter(x=x_window, y=FST_values, name=name[0], mode='markers', marker=dict(color=final_colours[0])), row=1, col=1)
            fig.update_layout(height=600, width=1500)
            fig.update_xaxes(title_text="Window Number", row=1, col=1)
            fig.update_yaxes(title_text="FST", row=1, col=1)

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        header = "FST graph"
        description = """..."""
        return render_template('notdash2.html', graphJSON=graphJSON, header=header, description=description)


# Download button 
@app.route("/Download", methods=['GET', 'POST'])

#this function will download the summary statistics output as a text file
def export_txt():
    data = {}
    for i, j in final_stat_global.items():
        name = []  
        val = []
        for k, l in j.items():
            name.append(k)
            val.append(l)
        df = pd.DataFrame()
        df["name"] = name
        df["val"] = val
        data[i] = df

    # convert the dataframe into text file
    # buffer is used to writing a file
    buffer = io.BytesIO()

    for i, df in data.items():
        buffer.write(bytes("\n" + i + "\n", 'utf-8'))
        # converting into csv
        df.to_csv(buffer, index=False)

    buffer.seek(0)
    return send_file(buffer,
                     attachment_filename="test.txt",
                     mimetype='text/csv')


################################################# Start web server #####################################

if __name__ == '__main__':
    app.run(debug=True)
