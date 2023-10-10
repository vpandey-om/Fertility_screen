import pandas as pd
import dash
from dash.dependencies import Input, Output
import dash_bio as dashbio
from dash import html, dcc
from sklearn.impute import KNNImputer
import numpy as np
import flask

# app = dash.Dash(__name__)

df1 = pd.read_csv('genes_phenotype_clutsering.txt',sep='\t').set_index('Gene ID')
imputer = KNNImputer(n_neighbors=5, weights="uniform")
X=imputer.fit_transform(df1.values)

#### gene discription 
geneanot = pd.read_csv('gene_des.txt',sep='\t')
gene_des=[]
for index in df1.index:
    tmp=geneanot[geneanot['new_id']==index]
    if tmp.empty:
        gene_des.append(index+'|'+'NA')
    else:
        gene_des.append(index+'|'+tmp['description'].to_list()[0])

df1.index=gene_des

###

df=pd.DataFrame(index=df1.index,columns=df1.columns,data=X)
df['Blood RGR' ]=np.log2( df['Blood RGR' ])
columns = list(df.columns.values)
rows = list(df.index)


server = flask.Flask(__name__)

app = dash.Dash(
  __name__,
  server=server,
  routes_pathname_prefix='/phenocluster/')


app.layout = html.Div([
    "Genes",
    dcc.Dropdown(
        id='my-default-clustergram-input',
        options=[
            {'label': row, 'value': row} for row in list(df.index)
        ],
        value=rows[0],
        multi=False
    ),

    html.Div(id='my-default-clustergram')
])

@app.callback(
    Output('my-default-clustergram', 'children'),
    Input('my-default-clustergram-input', 'value')
)
def update_clustergram(gene):
    # if len(rows) < 2:
    #     return "Please select at least two rows to display."
    # Add a horizontal rectangle
    # print('gene',gene)
    # if gene==None:
    #     gene=df.index.to_list()[0]
        
    figure=dashbio.Clustergram(
        data=df.values,
        column_labels=columns,
        cluster='row',
        standardize='none',
        row_labels=df.index.to_list(),
        color_threshold={
            'row': 250,
            'col': 700
        },
        hidden_labels='row',
        color_map= [
        [0.0, '#ff0000'],
        [0.25, '#ff00bf'],
        [0.5, '#FFFFFF'],
        [0.75, '#00ff80'],
        [1.0, '#80ff00']
        ],
        # row_group_marker=[
        # {'group': 2, 'annotation': 'cluster 2', 'color': '#AB63FA'},
        # {'group': 1, 'annotation': '', 'color': '#19D3F3'},
        # ],
        height=800,
        width=700
    )
    xx=figure.to_dict()
   
    genes=xx['layout']['yaxis11']['ticktext']

    tickvals=xx['layout']['yaxis11']['tickvals']
    index=genes.index(gene)
    print(gene)
    tickval=tickvals[index]
    figure.add_hline(y=tickval,line_width=1,opacity=0.5)
    g=dcc.Graph(figure=figure)
    
    return g

if __name__ == '__main__':
    app.run_server(debug=True,host = '127.0.0.1',port=8055)
