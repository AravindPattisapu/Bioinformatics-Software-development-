####### Dependencies #######
# Install flask using command line before: pip install flask
# Install flask_wtf from the command line: pip install flask-wtf
# Install plotly using command line: pip install plotly

# Import flask
# Import neccessary library required from flask
from flask import Flask, render_template, request, redirect, url_for, flash, session, send_file

from msilib.schema import tables

# Pandas for handling and data/dataframe
import pandas as pd

# sqlite3 for accessing pandas
import sqlite3

# Import forms from form.py file 
from forms import SNPForm

# IO for input/output operations
import io

# importing graph libraries
import matplotlib.pyplot as plt

# import the returned values for genotype frequency and allele frequency function 
from GF_AF_functions import gtfreq, allefreq

# importing the functions for stats  
from UPDATED_STATS import *

# importing search and query functions
from Validate_Search_Functions import *

# imports json 
import json

# importing graph libraries
import plotly
import plotly.express as px

# importing package to allow all pairwise combinations of populations for FST calculation  
from itertools import combinations

# Create a flask application object
app = Flask(__name__, template_folder='templates')

# A secret key is needed to create and process forms(inputs)
app.config['SECRET_KEY']='f466c2ef7b41e5c65e9f'

####### we put all our code for the functions here ########

# Link to the main/home page. Takes info from home_template and displays a message.
#Methods GET and POST allow the user to input information from forms. 
@app.route("/", methods = ['GET','POST'])
def basepage():
    return render_template('base.html')

# This will render our SNP browser page.
@app.route('/home', methods = ['GET','POST'])
def home():
    # get the form
    form=SNPForm()
    # return the template
    return render_template('home_template.html',
     text = "Welcome to the QMUL SNP browser!",form=form) # formatted using a separate html template file

# this will return graph on button click
@app.route("/visualize", methods = ['GET','POST'])
def visualize():
    return render_template('graph.html')

# Function for the forms. Takes the input from forms.py and redirects the user

#####        Search Functions       #####

@app.route("/search", methods = ['GET','POST'])
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

        df["position"] = df["position"].astype(int)
        
         # creating empty columns for calculations
        df['af_0'] = ""
        df['af_1'] = ""
        df['gt_00'] = ""
        df['gt_01'] = ""
        df['gt_10'] = ""
        df['gt_11'] = ""
        # iter over each row adn calculate af and gf values
        for i,row in df.iterrows():
            print(i)
            gt = gtfreq(row['0|0'], row['0|1'], row['1|0'], row['1|1'])
            af = allefreq(row['0|0'], row['0|1'], row['1|0'], row['1|1'])

            # converting the values from dict into a list
            ls = list(gt.values())
            ls2 = list(af.values())

            # storing the values in specified location
            df['gt_00'].iloc[i] = ls[0]
            df['gt_01'].iloc[i] = ls[1]
            df['gt_10'].iloc[i] = ls[2]
            df['gt_11'].iloc[i] = ls[3]

            df['af_0'].iloc[i] = ls2[0]
            df['af_1'].iloc[i] = ls2[1]

        print(df)

        # getting all the statistics function list
        stats_list = []
        if 'nud' in request.form:
            stats_list.append('nud')
        if 'hapd' in request.form:
            stats_list.append('hapd')
        if 'td' in request.form:
            stats_list.append('td')
        if 'fst' in request.form:
            stats_list.append('fst')

        # getting the user input data from html file for population checkbox
        population_list = []
        population_names = []
        if 'GBR' in request.form:
            population_list.append('British')
            population_names.append('British:')
        if 'COL' in request.form:
            population_list.append('Colombian')
            population_names.append('Colombian:')
        if 'FIN' in request.form:
            population_list.append('Finnish')
            population_names.append('Finnish:')
        if 'PUN' in request.form:
            population_list.append('Punjabi')
            population_names.append('Punjabi:')
        if 'TEL' in request.form:
            population_list.append('Telugu')
            population_names.append('Telugu:')

    
        #print(df)
        df_temp = df
        df_temp = df_temp.drop(['sample_count','0|0', '0|1', '1|0', '1|1'], axis=1)
        #print(df_temp)
        
        # selecting only specific rows acording to input (with user given population)
        # will contain data for selected population 
        ind=pd.DataFrame([], columns=df_temp.columns.values)
        
        #indexing the list Populaion_list ; i is int (indexes (1, 2 , 3)) ;j is population list
        for i,j in enumerate(population_list):

            #to show the populations one by one 

            d = df_temp[df_temp["population"] == j]

            # concating - adding the data into ind fropm d

            ind = pd.concat([ind,d], axis=0)

        # Show the summary statistics only when the user enters rs value and gene name 
        if (request.form['start'] != "" and request.form['end'] != "") or (request.form['gene'] != ""):


            # p_names will contain -> all the data for selected population
            p_names = {}
            
            #zipping 2 lists together

            for pl, pn in zip(population_list, population_names):

                # pn - column name ; pl - 
                p_names[pn] = df[df['population'] == pl]
            

            # for statistics selected these will make it as true - booleans 

            check_nud = False
            check_hapd = False
            check_td = False
            check_fst = False


            # key -> E.g. : df_british (name) 
            # val -> all the rows for british pop only  

            
            #################### Statistics functions #####################

            if 'nud' in request.form: 
                nud = {}
                for key, val in p_names.items():
                    nud[key] = SQLtoNucDiv(val, val.iloc[0]["position"],val.iloc[-1]["position"])
                check_nud = True
        
            if 'hapd' in request.form:
                hapd = {}
                for key, val in p_names.items():
                    hapd[key] = SQLtoHapDiv(val)
                check_hapd = True
            
            if 'td' in request.form:
                td = {}
                td_w = {}
                for key, val in p_names.items():
                    td[key] = SQLtoTD(val)
                    td_w[key] = SQLtoTD_window(val, 5).tolist()
                print("heyyyy")
                print(td)
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
                    
            # final_stat contains all the values calculated  by stat functions

            final_stat = {}

            # if the check booleans are true it will feed the respective data into ther final stat dic 
            if check_nud:
                final_stat['Nucleotide diversity'] = nud
            if check_hapd:
                final_stat['Haplotype Diversity'] = hapd
            if check_td:
                final_stat["Tajima's D"] = td
                final_stat['Tajimas D Window'] = td_w
            if check_fst:
                final_stat['FST'] = fst

            print(final_stat)
        else:
            final_stat = {}

        #session - global variable that store the data , can be used through out the python file

        # if data is empty than return simple empty arrray (for safety)
        if df.empty:
            session["data"] = df.to_dict()
            return render_template('output.html', d=ind, zip=zip)
        
        #final return.
        df = df.drop(['sample_count','0|0', '0|1', '1|0', '1|1'], axis=1)
        session["data"] = df.to_dict()
        session["final_stat"] = final_stat
        return render_template('output.html', column_names= df_temp.columns.values, d=list(ind.values.tolist()), final_stat = final_stat, col_name = population_names ,zip=zip)
        # df.temp contains colum names 

        ############### app route for graphs #################

@app.route("/Nucleotide diversity")
def Nucleotide_graph():
    final_stat = session["final_stat"]

    #empty lists 
    name = []
    val = []

    #

    for i,j in final_stat["Nucleotide diversity"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["name"] = name
    df["val"] = val
    fig = px.scatter(df, x="name", y="val")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header="Nucleotide diversity graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header,description=description)

@app.route("/Haplotype Diversity")
def Haplotype_graph():
    final_stat = session["final_stat"]
    name = []
    val = []
    for i,j in final_stat["Haplotype Diversity"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["Population"] = name
    df["Value"] = val
    fig = px.scatter(df, x="Population", y="Value")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header="Haplotype Diversity graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header,description=description)

from plotly.subplots import make_subplots
import plotly.graph_objects as go


@app.route("/Tajimas D")
def Tajimas_graph():

    final_stat = session["final_stat"]

    name = []
    TD_window_list = []


    TD_windows_list = final_stat["Tajimas D Window"]["British:"]
    x_values = list(range(1, len(TD_windows_list)+1))
    my_dict = {"Window":x_values, "TD":TD_windows_list}
    df = pd.DataFrame(my_dict)
    fig = px.scatter(df, x="Window", y="TD")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header="Tajima's D graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header,description=description)


@app.route("/FST")
def FST_graph():
    final_stat = session["final_stat"]
    name = []
    val = []
    for i,j in final_stat["FST"].items():
        name.append(i)
        val.append(j)
    df = pd.DataFrame()
    df["Population"] = name
    df["Value"] = val
    fig = px.scatter(df, x="Population", y="Value")

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    header="FST graph"
    description = """..."""
    return render_template('notdash2.html', graphJSON=graphJSON, header=header,description=description)


#this function will download text file when click on download button
@app.route("/Download", methods=['GET', 'POST'])
def export_csv():
    print(session["data"])

    #converting session ndata into data frame
    dataframe = pd.DataFrame(session["data"])
    
    data = {}
    for i,j in session["final_stat"].items():
        name = []
        val = []

        for k,l in j.items():
            name.append(k)
            val.append(l)
        df = pd.DataFrame()
        df["name"] = name
        df["val"] = val
        data[i] = df


#dropping index column 

    dataframe = dataframe.reset_index(drop=True)
    print(dataframe)
    #convert the dataframe into text file

    # buffer is used to writing a file
    buffer = io.BytesIO()

    # converting it into csv 
    dataframe.to_csv(buffer,index=False)

    for i,df in data.items():
        buffer.write(bytes("\n" + i + "\n", 'utf-8'))
        df.to_csv(buffer,index=False)
        
    buffer.seek(0)
    return send_file(buffer,
                 attachment_filename="test.csv",
                 mimetype='text/csv')





####### start web server #######

if __name__ == '__main__':
    app.run(debug=True)