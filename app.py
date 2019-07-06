from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
import pandas as pd
import os
from werkzeug.utils import secure_filename
from flask import Flask, render_template, url_for, send_from_directory

BASE_DIR = os.getcwd()

MEDIA_FOLDER = os.path.join(BASE_DIR,'media')
ALLOWED_EXTENSIONS = set(['svg'])

app = Flask(__name__)
app.config['MEDIA_FOLDER'] = MEDIA_FOLDER

df = pd.read_csv(os.path.join(BASE_DIR,'dataset/sample.csv'))

DrawingOptions.atomLabelFontSize = 55
DrawingOptions.dotsPerAngstrom = 100
DrawingOptions.bondLineWidth = 3.0
width = 500


def imagegen(value):
    mol = Chem.MolFromSmiles(value)
    Draw.MolToFile( mol, f"media/compound.svg" )



names = df['name']
struct = df['structure']
pairs = dict(zip(names, struct))
	
print(pairs)
imagegen(pairs['Ethane'])


@app.route('/')
def index():
	return render_template('index.html')

@app.route('/media/<filename>')	
def media(filename):
    return send_from_directory(app.config['MEDIA_FOLDER'],
                               filename)


if __name__ == "__main__":
	app.run(debug=True)