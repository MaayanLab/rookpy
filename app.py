from flask import Flask, request, jsonify, session, send_file, redirect
from flask_cors import CORS
from flask_session import Session 

import os
import urllib.request
import pandas as pd
import numpy as np
import pickle
import qnorm

import traceback
import sys

app = Flask(__name__)
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)

print("server started")

base_name = os.environ.get('BASE_NAME', 'rookpy')
data_url = os.environ.get('DATA_URL', 'https://mssm-seq-matrix.s3.amazonaws.com/rooky_data.pkl')
token = os.environ.get('TOKEN', 'NA')

urllib.request.urlretrieve(data_url, 'data.pkl')
print("data loaded")

data = 0
with open('data.pkl', 'rb') as f:
    data = pickle.load(f)

@app.route('/'+base_name+'/signature', methods=["POST", "GET"])
def signature_search():
    print("search")
    query_data = request.get_json()
    print(query_data)

    signature_name = query_data.get("signatureName", "signature")
    upgenes = [x.upper() for x in query_data.get("upgenes", [])]
    downgenes = [x.upper() for x in query_data.get("downgenes",[])]
    #direction = query_data["direction"]
    species = query_data.get("species",[])
    signature_type = query_data.get("type", "")
    signature_values = query_data.get("signature", "")
    signature_genes = query_data.get("siggenes", "")

    if signature_type == "full_signature":
        genes = [x.upper() for x in list(data[species]["normalization"].index)]
        input_sig = pd.DataFrame(signature_values, index=signature_genes)
        signature = pd.DataFrame([0] * len(genes), index=genes)
        intergene = list(set(signature_genes).intersection(set(genes)))
        signature.loc[intergene, :] = input_sig.loc[intergene,:]

        qexp = pd.DataFrame(qnorm.quantile_normalize(np.log2(1+signature), target=data[species]["normalization"].loc[:,"target_quantiles"]))
        qexp = qexp.sub(np.array(data[species]["normalization"].loc[:,"mean"]), axis=0).div(np.array(data[species]["normalization"].loc[:,"std"]), axis=0)
        qexp = pd.DataFrame(np.dot(data[species]["transform"], qexp), dtype=np.float16)
        
        sigs = data[species]["signatures"].astype(np.float32)
        bs = pd.DataFrame(np.sqrt((sigs*sigs).sum(axis=0)))
        qs = np.sqrt((signature*signature).sum())
        ok = data[species]["signatures"].T.dot(signature)
        cosine_dist = (ok/bs/qs).sort_values(0, ascending=True)
        cosine_dist = (cosine_dist-cosine_dist.mean())/np.std(cosine_dist)
        significant = cosine_dist[cosine_dist[0] > 3]
        samples = [int(x.replace("GSM", "")) for x in significant.index]

        response = { 'name': signature_name}
        response["samples"] = samples
        response["similarity"] = list(significant[0])

    elif signature_type == "geneset":
        genes = [x.upper() for x in list(data[species]["normalization"].index)]
        signature = pd.DataFrame([0] * len(genes), index=genes)
        inter = list(set(genes).intersection(set(upgenes)))
        
        signature.loc[inter,:] = 1
        inter = list(set(genes).intersection(set(downgenes)))
        signature.loc[inter,:] = -1
        signature = pd.DataFrame(np.dot(data[species]["transform"], signature), dtype=np.float32)
        
        sigs = data[species]["signatures"].astype(np.float32)
        bs = pd.DataFrame(np.sqrt((sigs*sigs).sum(axis=0)))
        qs = np.sqrt((signature*signature).sum())
        ok = sigs.T.dot(signature)

        cosine_dist = (ok/bs/qs).sort_values(0, ascending=True)
        cosine_dist = (cosine_dist-cosine_dist.mean())/np.std(cosine_dist)
        significant = cosine_dist[cosine_dist[0] > 2.5]
        samples = [int(x.replace("GSM", "")) for x in significant.index]
        
        response = { 'name': signature_name}
        response["samples"] = samples

    return jsonify(response)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)