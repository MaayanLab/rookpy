from flask import Flask, request, jsonify, session, send_file, redirect
from flask_cors import CORS
from flask_session import Session 

import os
import h5py as h5
import urllib.request
import pandas as pd
import numpy as np


app = Flask(__name__)
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)

print("server started")

base_name = os.environ.get('BASE_NAME', 'rookpy')
data_url = os.environ.get('DATA_URL', 'https://mssm-seq-matrix.s3.amazonaws.com/rooky_data.pkl')
token = os.environ.get('TOKEN', 'NA')

urllib.request.urlretrieve(data_url, 'data.h5')

@app.route('/'+base_name+'/genes', methods=["POST", "GET"])
def listgenes():
    if request.method == 'GET':
        data = {"library": request.args.get("library")}
    else:
        data = request.get_json()

    f = h5.File("data.h5")
    if data["library"] in list(f["meta"]["genes"].keys()):
        genes = [x.decode() for x in list(f["meta"]["genes"][data["library"]])]
        response = { 'genes': genes}
    else:  
        response = { 'error': 'library not found'}
    f.close()
    return jsonify(response)

@app.route('/'+base_name+'/libraries', methods=["POST", "GET"])
def listlibraries():
    f = h5.File("matrix.h5")
    ids = list(f["meta"]["colid"].keys())
    response = { 'libraries': sorted(ids)}
    return jsonify(response)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)