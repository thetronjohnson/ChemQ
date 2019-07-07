from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
import pandas as pd
import os
from flask import Flask, render_template, url_for, send_from_directory,request, redirect

BASE_DIR = os.getcwd()

MEDIA_FOLDER = os.path.join(BASE_DIR,'media')
ALLOWED_EXTENSIONS = set(['svg','png'])

app = Flask(__name__)
app.config['MEDIA_FOLDER'] = MEDIA_FOLDER

df = pd.read_csv(os.path.join(BASE_DIR,'dataset/sample.csv'))

DrawingOptions.atomLabelFontSize = 55
DrawingOptions.dotsPerAngstrom = 100
DrawingOptions.bondLineWidth = 3.0
width = 500


def imagegen(value):
    mol = Chem.MolFromSmiles(value)
    Draw.MolToFile( mol, f"media/{value}.svg" )

names = df['name']
struct = df['structure']
pairs = dict(zip(names, struct))
	

@app.route('/')
def index():
	return render_template('index.html',pairs=pairs)

@app.route('/media/<filename>')	
def media(filename):
    return send_from_directory(app.config['MEDIA_FOLDER'],filename)

@app.route("/selected",methods=['GET', 'POST'])
def selected():
	if request.method == 'POST':
		select = request.form
		imagegen(select['compound'])
		return render_template('index.html',pairs=pairs,img=f"{select['compound']}.svg")
	return render_template('index.html',pairs=pairs)




if __name__ == "__main__":
	app.jinja_env.auto_reload = True
	app.config['TEMPLATES_AUTO_RELOAD'] = True
	app.run(debug=True)
